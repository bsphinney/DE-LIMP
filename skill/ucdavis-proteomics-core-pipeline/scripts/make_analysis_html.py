#!/usr/bin/env python3
"""
make_analysis_html.py -- one self-contained HTML file with the QC panels, the
figures and the analysis text, openable on any machine with no network.

WHY THIS IS THE DEFAULT DELIVERABLE
-----------------------------------
The .docx needs Word and drops you into a document you have to scroll; a folder of
PNGs needs the reader to already know which panel matters. Most people receiving
one of these runs open it on a laptop that has never had a proteomics tool
installed, often after copying the results off a cluster. A single HTML file
double-clicks open in whatever browser is already there, on Windows, macOS and
Linux alike -- and because every image is inlined as a data URI there is nothing
to lose in transit and nothing to fetch at read time.

That last property is the one that matters operationally: ONE file to copy. A
report that references ./figures/pca.png silently loses every figure the moment
someone copies just the report, which is exactly what people do.

QC PANELS COME FIRST, ABOVE THE RESULTS
---------------------------------------
Deliberate. A volcano plot is persuasive whether or not the run was any good, so a
reader who meets the biology first has already formed a conclusion by the time
they reach the evidence about whether to trust it. Detected-vs-inferred and the
per-sample protein counts decide how much weight the DE table can carry, so they
are placed where they get read.

Usage
  python3 make_analysis_html.py --session <session dir> --out report.html
  python3 make_analysis_html.py --report AI_Analysis_Report.md --figures ./figures \\
      --tables ./tables --out report.html [--title "..."] [--quality SAMPLE_QUALITY.md]
"""
import argparse, base64, csv, html, json, mimetypes, os, re, sys

# QC first, then overview, then per-contrast results. Anything not listed still
# gets rendered -- an unknown figure is shown rather than silently dropped.
FIGURE_ORDER = [
    ("qc_detected_vs_inferred", "Quality control"),
    ("qc_protein_counts", "Quality control"),
    ("pca", "Overview"),
    ("heatmap_top", "Overview"),
    ("volcano", "Differential expression"),
    ("pvalue", "Differential expression"),
]
SECTION_ORDER = ["Quality control", "Overview", "Differential expression", "Other figures"]


def data_uri(path):
    mime = mimetypes.guess_type(path)[0] or "image/png"
    with open(path, "rb") as fh:
        return f"data:{mime};base64,{base64.b64encode(fh.read()).decode('ascii')}"


def classify(name):
    stem = os.path.splitext(os.path.basename(name))[0]
    for key, section in FIGURE_ORDER:
        if stem.startswith(key):
            return section, FIGURE_ORDER.index((key, section))
    return "Other figures", len(FIGURE_ORDER)


def md_inline(t):
    t = html.escape(t)
    t = re.sub(r"`([^`]+)`", r"<code>\1</code>", t)
    t = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", t)
    t = re.sub(r"(?<![*\w])\*([^*]+)\*(?!\w)", r"<em>\1</em>", t)
    return t


def md_to_html(md):
    """Enough Markdown for the report the model writes: headings, lists, tables,
    fenced code, blockquotes. Not a general converter -- deliberately small so it
    has no dependencies to install on a cluster."""
    out, lines, i, n = [], md.splitlines(), 0, len(md.splitlines())
    while i < n:
        ln = lines[i]
        if ln.startswith("```"):
            block = []
            i += 1
            while i < n and not lines[i].startswith("```"):
                block.append(lines[i]); i += 1
            i += 1
            out.append("<pre><code>" + html.escape("\n".join(block)) + "</code></pre>")
            continue
        m = re.match(r"^(#{1,6})\s+(.*)", ln)
        if m:
            lvl = len(m.group(1))
            out.append(f"<h{lvl}>{md_inline(m.group(2))}</h{lvl}>")
            i += 1
            continue
        # table: header row, separator, body
        if ln.strip().startswith("|") and i + 1 < n and re.match(r"^\s*\|[\s:|-]+\|\s*$", lines[i + 1]):
            def cells(r):
                return [c.strip() for c in r.strip().strip("|").split("|")]
            head = cells(ln)
            i += 2
            body = []
            while i < n and lines[i].strip().startswith("|"):
                body.append(cells(lines[i])); i += 1
            t = ['<div class="tablewrap"><table><thead><tr>']
            t += [f"<th>{md_inline(c)}</th>" for c in head]
            t.append("</tr></thead><tbody>")
            for r in body:
                t.append("<tr>" + "".join(f"<td>{md_inline(c)}</td>" for c in r) + "</tr>")
            t.append("</tbody></table></div>")
            out.append("".join(t))
            continue
        if re.match(r"^\s*[-*+]\s+", ln):
            items = []
            while i < n and re.match(r"^\s*[-*+]\s+", lines[i]):
                items.append(re.sub(r"^\s*[-*+]\s+", "", lines[i])); i += 1
            out.append("<ul>" + "".join(f"<li>{md_inline(x)}</li>" for x in items) + "</ul>")
            continue
        if re.match(r"^\s*\d+[.)]\s+", ln):
            items = []
            while i < n and re.match(r"^\s*\d+[.)]\s+", lines[i]):
                items.append(re.sub(r"^\s*\d+[.)]\s+", "", lines[i])); i += 1
            out.append("<ol>" + "".join(f"<li>{md_inline(x)}</li>" for x in items) + "</ol>")
            continue
        if ln.strip().startswith(">"):
            q = []
            while i < n and lines[i].strip().startswith(">"):
                q.append(re.sub(r"^\s*>\s?", "", lines[i])); i += 1
            out.append(f"<blockquote>{md_inline(' '.join(q))}</blockquote>")
            continue
        if not ln.strip():
            i += 1
            continue
        para = []
        while i < n and lines[i].strip() and not re.match(r"^(#{1,6}\s|```|\s*[-*+]\s|\s*\d+[.)]\s|\s*>)", lines[i]) \
                and not lines[i].strip().startswith("|"):
            para.append(lines[i]); i += 1
        out.append(f"<p>{md_inline(' '.join(para))}</p>")
    return "\n".join(out)


def de_summary(tables_dir, adjp=0.05, logfc=1.0):
    """Count significant proteins per contrast. Read from the DE CSVs rather than
    re-stating whatever the prose claimed -- if the two disagree, the reader can
    see it."""
    rows = []
    if not tables_dir or not os.path.isdir(tables_dir):
        return rows
    for fn in sorted(os.listdir(tables_dir)):
        m = re.match(r"^DE_(\w+)_(.+)\.csv$", fn)
        if not m:
            continue
        up = dn = tot = 0
        try:
            with open(os.path.join(tables_dir, fn), newline="") as fh:
                for r in csv.DictReader(fh):
                    tot += 1
                    p = r.get("adj.P.Val") or r.get("padj") or r.get("FDR")
                    lf = r.get("logFC") or r.get("log2FoldChange")
                    try:
                        p, lf = float(p), float(lf)
                    except (TypeError, ValueError):
                        continue
                    if p < adjp and abs(lf) >= logfc:
                        up += lf > 0
                        dn += lf < 0
        except Exception:
            continue
        rows.append({"contrast": m.group(2).replace(".", " vs "), "method": m.group(1),
                     "tested": tot, "up": up, "down": dn, "file": fn})
    return rows


CSS = """
:root{--bg:#fff;--fg:#16191d;--mut:#5b6470;--line:#e2e6eb;--card:#f7f9fb;--accent:#2a78d6;--warn:#b4531a}
@media (prefers-color-scheme:dark){:root{--bg:#14171a;--fg:#e8eaed;--mut:#9aa4b0;--line:#2b3138;--card:#1b1f24;--accent:#5fa3f0;--warn:#e08c4e}}
:root[data-theme=dark]{--bg:#14171a;--fg:#e8eaed;--mut:#9aa4b0;--line:#2b3138;--card:#1b1f24;--accent:#5fa3f0;--warn:#e08c4e}
:root[data-theme=light]{--bg:#fff;--fg:#16191d;--mut:#5b6470;--line:#e2e6eb;--card:#f7f9fb;--accent:#2a78d6;--warn:#b4531a}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:16px/1.65 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;overflow-x:hidden}
.wrap{max-width:60rem;margin:0 auto;padding:2.5rem 1.25rem 5rem}
h1{font-size:1.9rem;line-height:1.25;margin:0 0 .3rem}
h2{font-size:1.35rem;margin:2.6rem 0 .8rem;padding-bottom:.35rem;border-bottom:2px solid var(--line)}
h3{font-size:1.08rem;margin:1.8rem 0 .5rem}
p,li{color:var(--fg)}
code{background:var(--card);border:1px solid var(--line);border-radius:4px;padding:.1em .35em;font-size:.88em}
pre{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:.9rem 1rem;overflow-x:auto}
pre code{background:none;border:0;padding:0}
blockquote{margin:1rem 0;padding:.6rem 1rem;border-left:3px solid var(--accent);background:var(--card);color:var(--mut)}
.sub{color:var(--mut);margin:.2rem 0 0}
.tablewrap{overflow-x:auto;margin:1rem 0}
table{border-collapse:collapse;width:100%;font-size:.93rem}
th,td{text-align:left;padding:.5rem .7rem;border-bottom:1px solid var(--line);vertical-align:top}
th{font-weight:600;color:var(--mut);font-size:.82rem;text-transform:uppercase;letter-spacing:.03em}
figure{margin:1.6rem 0;padding:1rem;background:var(--card);border:1px solid var(--line);border-radius:10px}
figure img{width:100%;height:auto;display:block;border-radius:6px;background:#fff}
figcaption{color:var(--mut);font-size:.9rem;margin-top:.7rem}
.qc{border-left:3px solid var(--accent)}
.banner{background:var(--card);border:1px solid var(--line);border-left:3px solid var(--warn);border-radius:8px;padding:.9rem 1.1rem;margin:1.4rem 0;font-size:.94rem}
.toc{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:1rem 1.2rem;margin:1.6rem 0}
.toc ul{margin:.4rem 0 0;padding-left:1.1rem}
.toc a{color:var(--accent);text-decoration:none}
.toc a:hover{text-decoration:underline}
.toggle{position:fixed;top:1rem;right:1rem;background:var(--card);color:var(--fg);border:1px solid var(--line);border-radius:6px;padding:.4rem .7rem;font-size:.85rem;cursor:pointer;z-index:9}
@media print{.toggle{display:none}figure{break-inside:avoid}}
"""

JS = """
(function(){
 var b=document.getElementById('tt');
 function cur(){var a=document.documentElement.getAttribute('data-theme');
  if(a)return a;return matchMedia('(prefers-color-scheme:dark)').matches?'dark':'light';}
 function set(t){document.documentElement.setAttribute('data-theme',t);
  b.textContent=t==='dark'?'\\u2600 Light':'\\u263e Dark';}
 set(cur()); b.addEventListener('click',function(){set(cur()==='dark'?'light':'dark');});
})();
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--session", help="session dir; infers report/figures/tables under output/")
    ap.add_argument("--report", help="analysis report markdown")
    ap.add_argument("--figures", help="figures dir")
    ap.add_argument("--tables", help="DE tables dir")
    ap.add_argument("--quality", help="SAMPLE_QUALITY.md")
    ap.add_argument("--audit", help="AUDIT.md")
    ap.add_argument("--title", default="Proteomics Analysis Report")
    ap.add_argument("--adjp", type=float, default=0.05)
    ap.add_argument("--logfc", type=float, default=1.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    if a.session:
        o = os.path.join(a.session, "output")
        a.report = a.report or os.path.join(o, "AI_Analysis_Report.md")
        a.figures = a.figures or os.path.join(o, "figures")
        a.tables = a.tables or os.path.join(o, "tables")
        for attr, fn in (("quality", "SAMPLE_QUALITY.md"), ("audit", "AUDIT.md")):
            p = os.path.join(o, fn)
            if not getattr(a, attr) and os.path.exists(p):
                setattr(a, attr, p)

    caps = {}
    if a.figures and os.path.exists(os.path.join(a.figures, "figures.json")):
        try:
            fj = json.load(open(os.path.join(a.figures, "figures.json")))
            entries = fj if isinstance(fj, list) else fj.get("figures", [])
            for e in entries:
                if isinstance(e, dict):
                    k = e.get("file") or e.get("filename") or e.get("name")
                    if k:
                        caps[os.path.basename(k)] = e.get("caption") or e.get("title") or ""
        except Exception:
            pass

    figs = {}
    if a.figures and os.path.isdir(a.figures):
        for fn in sorted(os.listdir(a.figures)):
            if not fn.lower().endswith((".png", ".jpg", ".jpeg", ".svg", ".webp")):
                continue
            sec, rank = classify(fn)
            figs.setdefault(sec, []).append((rank, fn, os.path.join(a.figures, fn)))

    body, toc = [], []

    def sect(anchor, title):
        toc.append(f'<li><a href="#{anchor}">{html.escape(title)}</a></li>')
        body.append(f'<h2 id="{anchor}">{html.escape(title)}</h2>')

    de = de_summary(a.tables, a.adjp, a.logfc)
    if de:
        sect("summary", "Results at a glance")
        body.append(
            f"<p class='sub'>Counted directly from the DE tables at adjusted "
            f"p &lt; {a.adjp} and |log2 fold change| &ge; {a.logfc}.</p>")
        t = ["<div class='tablewrap'><table><thead><tr><th>Contrast</th><th>Proteins tested</th>"
             "<th>Higher</th><th>Lower</th><th>Total changed</th></tr></thead><tbody>"]
        for r in de:
            t.append(f"<tr><td>{html.escape(r['contrast'])}</td><td>{r['tested']:,}</td>"
                     f"<td>{r['up']:,}</td><td>{r['down']:,}</td><td>{r['up']+r['down']:,}</td></tr>")
        t.append("</tbody></table></div>")
        body.append("".join(t))

    # QC before results: see module docstring.
    for sec in SECTION_ORDER:
        if sec not in figs:
            continue
        anchor = re.sub(r"[^a-z0-9]+", "-", sec.lower()).strip("-")
        sect(anchor, sec)
        if sec == "Quality control":
            body.append("<div class='banner'>Read these first. They decide how much weight the "
                        "results below can carry &mdash; a volcano plot looks equally convincing "
                        "whether or not the run was any good.</div>")
        for _, fn, path in sorted(figs[sec]):
            cap = caps.get(fn) or os.path.splitext(fn)[0].replace("_", " ")
            cls = " class='qc'" if sec == "Quality control" else ""
            try:
                uri = data_uri(path)
            except Exception as e:
                body.append(f"<p class='sub'>[could not embed {html.escape(fn)}: {html.escape(str(e))}]</p>")
                continue
            body.append(f"<figure{cls}><img src='{uri}' alt='{html.escape(cap)}'>"
                        f"<figcaption>{md_inline(cap)}</figcaption></figure>")

    for path, title, anchor in ((a.quality, "Sample quality notes", "quality"),
                                (a.audit, "Audit &amp; caveats", "audit"),
                                (a.report, "Analysis", "analysis")):
        if path and os.path.exists(path):
            toc.append(f'<li><a href="#{anchor}">{title}</a></li>')
            body.append(f'<h2 id="{anchor}">{title}</h2>')
            body.append(md_to_html(open(path, encoding="utf-8", errors="replace").read()))

    if not body:
        sys.exit("[make_analysis_html] nothing to render — check --session/--report/--figures")

    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(a.title)}</title><style>{CSS}</style></head>
<body><button class="toggle" id="tt">Dark</button><div class="wrap">
<h1>{html.escape(a.title)}</h1>
<p class="sub">Self-contained &mdash; every figure is embedded, so this one file is the whole
report. No network needed; copy it anywhere and double-click to open.</p>
<div class="toc"><strong>Contents</strong><ul>{''.join(toc)}</ul></div>
{''.join(body)}
</div><script>{JS}</script></body></html>"""

    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write(doc)
    nfig = sum(len(v) for v in figs.values())
    print(json.dumps({"wrote": a.out, "bytes": os.path.getsize(a.out),
                      "figures_embedded": nfig, "contrasts": len(de),
                      "self_contained": True}, indent=2))


if __name__ == "__main__":
    main()
