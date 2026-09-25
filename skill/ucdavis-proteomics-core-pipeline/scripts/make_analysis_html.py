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

FIGURES SIT WHERE THE TEXT DISCUSSES THEM
-----------------------------------------
Each `![alt](figures/x.png)` in AI_Analysis_Report.md becomes an embedded, numbered
figure AT THAT POSITION -- the alt text is its title, figures.json's caption (when
there is one) its caption. The report used to print those references as literal
Markdown (28 of them on Silva08172026, 2026-09-24) while dumping every image in
figures/ into galleries at the top, a stale pca_original_labels.png included.
  * an image the report does not reference is NOT embedded; one warning line names it;
  * a referenced image that is missing becomes a visible "figure missing: <file>" note
    (and a warning) -- never raw Markdown;
  * only paths inside the session are embedded -- no URLs, no ../ escapes;
  * the galleries appear only when there is no report at all (the quick pre-analysis
    page), and then list figures.json's figures only.
The tool's own "Sample quality notes" / "Audit & caveats" sections are left out when
the report already has that section, so no heading appears twice.

Usage
  python3 make_analysis_html.py --session <session dir> --out report.html
  python3 make_analysis_html.py --report AI_Analysis_Report.md --figures ./figures \\
      --tables ./tables --out report.html [--title "..."] [--quality SAMPLE_QUALITY.md]
"""
import argparse, base64, csv, html, json, mimetypes, os, re, sys, urllib.parse

# Galleries (no report only): QC first, then overview, then per-contrast results.
FIGURE_ORDER = [
    ("qc_detected_vs_inferred", "Quality control"),
    ("qc_protein_counts", "Quality control"),
    ("pca", "Overview"),
    ("heatmap_top", "Overview"),
    ("volcano", "Differential expression"),
    ("pvalue", "Differential expression"),
]
SECTION_ORDER = ["Quality control", "Overview", "Differential expression", "Other figures"]

IMAGE_EXT = (".png", ".jpg", ".jpeg", ".svg", ".webp")
# Markdown images ![alt](path "title") -> groups 1 (alt), 2 (path); inline <img src="path">
# -> group 3 (path).
IMG_REF = re.compile(r"""!\[([^\]]*)\]\(\s*<?([^)\s>]+)>?(?:\s+["'][^)]*["'])?\s*\)"""
                     r"""|<img\b[^>]*?\bsrc\s*=\s*["']([^"']+)["'][^>]*>""", re.I)
_ALT_ATTR = re.compile(r"""\balt\s*=\s*["']([^"']*)["']""", re.I)
# The tool's own sections, and the report headings that already cover them.
SUPERSEDED_BY = {"quality": {"sample quality notes", "sample quality", "data quality notes"},
                 "audit": {"audit & caveats", "audit and caveats", "audit"}}


def _ref(m):
    """(alt, path) of an IMG_REF match."""
    if m.group(2) is not None:
        return m.group(1), m.group(2)
    a = _ALT_ATTR.search(m.group(0))
    return (a.group(1) if a else ""), m.group(3)


def _clean(ref):
    return urllib.parse.unquote(ref).split("?")[0].split("#")[0]


def _external(ref):
    return "://" in ref or ref.startswith(("data:", "//", "mailto:"))


def report_figures(md_text):
    """Basenames of the images a Markdown report references, in order of first mention."""
    out = []
    for m in IMG_REF.finditer(md_text or ""):
        ref = _ref(m)[1]
        if _external(ref):
            continue                      # remote or inline: not a file in figures/
        name = os.path.basename(_clean(ref))
        if name.lower().endswith(IMAGE_EXT) and name not in out:
            out.append(name)
    return out


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


def _inside(path, root):
    root = os.path.realpath(root)
    return path == root or path.startswith(root.rstrip(os.sep) + os.sep)


class Figures:
    """Turns image references into numbered, embedded figures, and keeps the ledger the
    warnings are written from. One instance per page, so numbering runs through it."""

    def __init__(self, base_dir, root, figures_dir=None, captions=None):
        self.base, self.root, self.figdir = base_dir, root, figures_dir
        self.caps = captions or {}
        self.n, self.by_path = 0, {}
        self.embedded, self.missing, self.rejected = [], [], []

    def _locate(self, ref):
        rel = _clean(ref)
        cand = [os.path.realpath(rel if os.path.isabs(rel) else os.path.join(self.base, rel))]
        if self.figdir:           # --figures given apart from --report: look there by name
            cand.append(os.path.realpath(os.path.join(self.figdir, os.path.basename(rel))))
        inside = [c for c in cand if _inside(c, self.root) or
                  (self.figdir and _inside(c, self.figdir))]
        return next((c for c in inside if os.path.isfile(c)), None), bool(inside)

    def render(self, alt, ref, section=None):
        name = os.path.basename(_clean(ref)) or ref
        if _external(ref):
            self.rejected.append(ref)
            return self.note(f"figure not embedded (not a file in this session): {ref}")
        path, allowed = self._locate(ref)
        if not allowed:
            self.rejected.append(ref)
            return self.note(f"figure not embedded (outside the session folder): {ref}")
        if path is None:
            self.missing.append(name)
            return self.note(f"figure missing: {name}")
        if path in self.by_path:
            n = self.by_path[path]
            return f'<p class="figref">(See <a href="#fig-{n}">Figure {n}</a>.)</p>'
        try:
            uri = data_uri(path)
        except OSError as e:
            self.missing.append(name)
            return self.note(f"figure missing: {name} (unreadable: {e})")
        self.n += 1
        self.by_path[path] = self.n
        self.embedded.append(os.path.basename(path))
        return self.figure(self.n, uri, alt, self.caps.get(os.path.basename(path)), section)

    # presentation -- overridden by report_style when the page is styled
    @staticmethod
    def note(text):
        return f'<div class="figmissing" role="note">{html.escape(text)}</div>'

    @staticmethod
    def figure(n, uri, alt, caption, section=None):
        cls = " class='qc'" if section == "Quality control" else ""
        title = md_inline(alt) if alt else ""
        cap = (" &mdash; " + md_inline(caption)) if caption and alt else md_inline(caption or "")
        return (f'<figure id="fig-{n}"{cls}><img src="{uri}" alt="{html.escape(alt or caption or "")}">'
                f"<figcaption><strong>Figure {n}.</strong> {title}{cap}</figcaption></figure>")


def md_inline(t, figs=None):
    """Inline Markdown. Image references become figures (via `figs`) and never reach the
    page as text; without `figs` they are dropped to their alt text."""
    out, last = [], 0
    for m in IMG_REF.finditer(t or ""):
        out.append(_inline_text(t[last:m.start()]))
        alt, ref = _ref(m)
        out.append(figs.render(alt, ref) if figs else html.escape(alt))
        last = m.end()
    out.append(_inline_text((t or "")[last:]))
    return "".join(out)


def _inline_text(t):
    t = html.escape(t)
    t = re.sub(r"`([^`]+)`", r"<code>\1</code>", t)
    t = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", t)
    t = re.sub(r"(?<![*\w])\*([^*]+)\*(?!\w)", r"<em>\1</em>", t)
    return t


def _para(text, figs):
    """A paragraph, split around its images: text runs stay <p>, each image becomes a
    block figure at its position."""
    out, buf, last = [], [], 0
    for m in IMG_REF.finditer(text):
        buf.append(text[last:m.start()])
        chunk = "".join(buf).strip()
        if chunk:
            out.append(f"<p>{md_inline(chunk)}</p>")
        buf = []
        out.append(figs.render(*_ref(m)) if figs else html.escape(_ref(m)[0]))
        last = m.end()
    rest = ("".join(buf) + text[last:]).strip()
    if rest:
        out.append(f"<p>{md_inline(rest)}</p>")
    return "\n".join(out)


def anchor(title, used):
    base = re.sub(r"[^a-z0-9]+", "-", html.unescape(re.sub(r"<[^>]+>", "", title)).lower()).strip("-") or "section"
    a, i = base, 2
    while a in used:
        a, i = f"{base}-{i}", i + 1
    used.add(a)
    return a


def md_to_html(md, figs=None, used=None, drop_h1=False):
    """Enough Markdown for the report the model writes: headings, lists, tables,
    fenced code, blockquotes, images. Not a general converter -- deliberately small so it
    has no dependencies to install on a cluster.
    -> (html, [(level, id, text)], h1 text or None). With drop_h1 the first level-1
    heading is returned as the title instead of rendered."""
    used = used if used is not None else set()
    out, heads, h1 = [], [], None
    lines = md.splitlines()
    i, n = 0, len(lines)
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
            lvl, text = len(m.group(1)), m.group(2).strip()
            i += 1
            if lvl == 1 and drop_h1 and h1 is None:
                h1 = text
                continue
            aid = anchor(text, used)
            heads.append((lvl, aid, text))
            out.append(f'<h{lvl} id="{aid}">{md_inline(text)}</h{lvl}>')
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
            out.append(table_html(head, body, figs))
            continue
        if re.match(r"^\s*[-*+]\s+", ln):
            items = []
            while i < n and re.match(r"^\s*[-*+]\s+", lines[i]):
                items.append(re.sub(r"^\s*[-*+]\s+", "", lines[i])); i += 1
            out.append("<ul>" + "".join(f"<li>{md_inline(x, figs)}</li>" for x in items) + "</ul>")
            continue
        if re.match(r"^\s*\d+[.)]\s+", ln):
            items = []
            while i < n and re.match(r"^\s*\d+[.)]\s+", lines[i]):
                items.append(re.sub(r"^\s*\d+[.)]\s+", "", lines[i])); i += 1
            out.append("<ol>" + "".join(f"<li>{md_inline(x, figs)}</li>" for x in items) + "</ol>")
            continue
        if ln.strip().startswith(">"):
            q = []
            while i < n and lines[i].strip().startswith(">"):
                q.append(re.sub(r"^\s*>\s?", "", lines[i])); i += 1
            out.append(f"<blockquote>{md_inline(' '.join(q), figs)}</blockquote>")
            continue
        if not ln.strip():
            i += 1
            continue
        para = []
        while i < n and lines[i].strip() and not re.match(r"^(#{1,6}\s|```|\s*[-*+]\s|\s*\d+[.)]\s|\s*>)", lines[i]) \
                and not lines[i].strip().startswith("|"):
            para.append(lines[i]); i += 1
        out.append(_para(" ".join(para), figs))
    return "\n".join(out), heads, h1


def table_html(head, body, figs=None):
    t = ['<div class="tablewrap"><table><thead><tr>']
    t += [f"<th>{md_inline(c)}</th>" for c in head]
    t.append("</tr></thead><tbody>")
    for r in body:
        t.append("<tr>" + "".join(f"<td>{md_inline(c, figs)}</td>" for c in r) + "</tr>")
    t.append("</tbody></table></div>")
    return "".join(t)


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
.figmissing{margin:1.2rem 0;padding:.7rem 1rem;border:1px dashed var(--warn);border-radius:8px;color:var(--warn);font-size:.92rem}
.figref{color:var(--mut);font-size:.92rem}
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


def norm_title(t):
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", "", t or "")).strip().lower())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--session", help="session dir; infers report/figures/tables under output/")
    ap.add_argument("--report", help="analysis report markdown")
    ap.add_argument("--figures", help="figures dir")
    ap.add_argument("--tables", help="DE tables dir")
    ap.add_argument("--quality", help="SAMPLE_QUALITY.md")
    ap.add_argument("--audit", help="AUDIT.md")
    ap.add_argument("--title", help="page title (default: the report's own # heading)")
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
    has_report = bool(a.report and os.path.exists(a.report))

    caps, listed = {}, []
    if a.figures and os.path.exists(os.path.join(a.figures, "figures.json")):
        try:
            fj = json.load(open(os.path.join(a.figures, "figures.json")))
            entries = fj if isinstance(fj, list) else fj.get("figures", [])
            for e in entries:
                if isinstance(e, dict):
                    k = e.get("file") or e.get("filename") or e.get("name")
                    if k:
                        listed.append(os.path.basename(k))
                        caps[os.path.basename(k)] = e.get("caption") or e.get("title") or ""
        except Exception as e:
            print(f"[make_analysis_html] WARNING: figures.json unreadable ({e}); captions and "
                  f"its figure list are not used", file=sys.stderr)
    available = (sorted(fn for fn in os.listdir(a.figures) if fn.lower().endswith(IMAGE_EXT))
                 if a.figures and os.path.isdir(a.figures) else [])

    # Images resolve relative to the report's folder and must stay inside the session.
    base = os.path.dirname(os.path.abspath(a.report)) if has_report else (a.figures or ".")
    root = os.path.abspath(a.session) if a.session else base
    figs = Figures(base, root, a.figures, caps)

    used, body, toc = set(), [], []

    def sect(title, content):
        aid = anchor(title, used)
        toc.append(f'<li><a href="#{aid}">{html.escape(title)}</a></li>')
        body.append(f'<h2 id="{aid}">{html.escape(title)}</h2>')
        body.append(content)

    report_html, heads, h1 = "", [], None
    if has_report:
        with open(a.report, encoding="utf-8", errors="replace") as fh:
            md = fh.read()
        # Rendered first: its headings decide which tool sections are redundant, and its
        # ids are reserved in `used` so the tool's sections can never collide with them.
        report_html, heads, h1 = md_to_html(md, figs, used=used, drop_h1=True)
    report_h2 = {norm_title(t) for lvl, _, t in heads if lvl == 2}

    de = de_summary(a.tables, a.adjp, a.logfc)
    if de:
        t = ["<p class='sub'>Counted directly from the DE tables at adjusted "
             f"p &lt; {a.adjp} and |log2 fold change| &ge; {a.logfc}.</p>",
             "<div class='tablewrap'><table><thead><tr><th>Contrast</th><th>Proteins tested</th>"
             "<th>Higher</th><th>Lower</th><th>Total changed</th></tr></thead><tbody>"]
        for r in de:
            t.append(f"<tr><td>{html.escape(r['contrast'])}</td><td>{r['tested']:,}</td>"
                     f"<td>{r['up']:,}</td><td>{r['down']:,}</td><td>{r['up']+r['down']:,}</td></tr>")
        t.append("</tbody></table></div>")
        sect("Results at a glance", "".join(t))

    left_out = []
    if not has_report:
        # The quick pre-analysis page: galleries of figures.json's figures only.
        gal = {}
        for fn in listed:
            sec, rank = classify(fn)
            gal.setdefault(sec, []).append((rank, fn))
        for sec in SECTION_ORDER:
            if sec not in gal:
                continue
            parts = []
            if sec == "Quality control":
                parts.append("<div class='banner'>Read these first. They decide how much weight "
                             "the results below can carry &mdash; a volcano plot looks equally "
                             "convincing whether or not the run was any good.</div>")
            for _, fn in sorted(gal[sec]):
                parts.append(figs.render("", fn, section=sec))
            sect(sec, "".join(parts))
        left_out = [f for f in available if f not in set(listed)]
        why = "not listed in figures.json" if listed else "no report and no figures.json"
    else:
        shown = set(figs.embedded)
        referenced = set(report_figures(md))
        left_out = [f for f in available if f not in shown and f not in referenced]
        why = f"not referenced in {os.path.basename(a.report)}"

    for path, title, key in ((a.quality, "Sample quality notes", "quality"),
                             (a.audit, "Audit & caveats", "audit")):
        if not (path and os.path.exists(path)):
            continue
        if report_h2 & (SUPERSEDED_BY[key] | {norm_title(title)}):
            continue                    # the report has its own -- never show it twice
        with open(path, encoding="utf-8", errors="replace") as fh:
            sect(title, md_to_html(fh.read(), None, used=used, drop_h1=True)[0])

    if has_report:
        for lvl, aid, text in heads:
            if lvl == 2:
                toc.append(f'<li><a href="#{aid}">{md_inline(text)}</a></li>')
        body.append(report_html)

    if not body:
        sys.exit("[make_analysis_html] nothing to render — check --session/--report/--figures")

    if left_out:
        print(f"[make_analysis_html] WARNING: {len(left_out)} image(s) in {a.figures} are "
              f"{why} and were NOT embedded: {', '.join(left_out)}", file=sys.stderr)
    if figs.missing:
        print(f"[make_analysis_html] WARNING: {len(figs.missing)} referenced image(s) are "
              f"missing and are shown as a 'figure missing' note: {', '.join(figs.missing)}",
              file=sys.stderr)
    if figs.rejected:
        print(f"[make_analysis_html] WARNING: {len(figs.rejected)} image reference(s) point "
              f"outside the session and were not embedded: {', '.join(figs.rejected)}",
              file=sys.stderr)

    title = a.title or (html.unescape(re.sub(r"[*`]", "", h1)) if h1 else "Proteomics Analysis Report")
    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(title)}</title><style>{CSS}</style></head>
<body><button class="toggle" id="tt">Dark</button><div class="wrap">
<h1>{html.escape(title)}</h1>
<p class="sub">Self-contained &mdash; every figure is embedded, so this one file is the whole
report. No network needed; copy it anywhere and double-click to open.</p>
<div class="toc"><strong>Contents</strong><ul>{''.join(toc)}</ul></div>
{''.join(body)}
</div><script>{JS}</script></body></html>"""

    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write(doc)
    print(json.dumps({"wrote": a.out, "bytes": os.path.getsize(a.out),
                      "figures_embedded": figs.n,
                      "figure_list": os.path.basename(a.report) if has_report else "figures.json",
                      "figures_not_embedded": left_out, "figures_missing": figs.missing,
                      "figures_rejected": figs.rejected,
                      "contrasts": len(de), "self_contained": True}, indent=2))


if __name__ == "__main__":
    main()
