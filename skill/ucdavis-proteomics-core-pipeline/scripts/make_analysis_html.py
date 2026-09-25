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
import argparse, base64, csv, datetime, html, json, mimetypes, os, re, sys, urllib.parse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import report_style as rs  # noqa: E402  -- the ONE look shared by the skill's HTML pages

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

    @staticmethod
    def note(text):
        return rs.note(text)

    @staticmethod
    def figure(n, uri, alt, caption, section=None):
        return rs.figure_card(n, uri, alt or caption or "", md_inline(alt) if alt else "",
                              md_inline(caption) if caption else "")


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
            # keep the Markdown's own number: a list broken up by nested bullets must go on
            # 2, 3, ... instead of restarting at 1 after every break
            start = int(re.match(r"^\s*(\d+)", ln).group(1))
            items = []
            while i < n and re.match(r"^\s*\d+[.)]\s+", lines[i]):
                items.append(re.sub(r"^\s*\d+[.)]\s+", "", lines[i])); i += 1
            out.append((f'<ol start="{start}">' if start != 1 else "<ol>")
                       + "".join(f"<li>{md_inline(x, figs)}</li>" for x in items) + "</ol>")
            continue
        if ln.strip().startswith(">"):
            q = []
            while i < n and lines[i].strip().startswith(">"):
                q.append(re.sub(r"^\s*>\s?", "", lines[i])); i += 1
            out.append(rs.callout("info", f"<p>{md_inline(' '.join(q), figs)}</p>"))
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
    return rs.table([md_inline(c) for c in head], [[md_inline(c, figs) for c in r] for r in body])


def significance_rule(tables_dir, default_adjp=0.05):
    """-> (adjp, source). The cutoff the DE actually applied, from de_provenance.json --
    significance is the adjusted p-value alone there (run_de.R: "adj.P.Val < adjp (BH); no
    fold-change filter"), so it is here too."""
    try:
        with open(os.path.join(tables_dir or "", "de_provenance.json")) as fh:
            adjp = json.load(fh).get("adjp")
        if isinstance(adjp, (int, float)):
            return float(adjp), "de_provenance.json"
    except (OSError, ValueError):
        pass
    return default_adjp, "--adjp"


def de_summary(tables_dir, adjp=0.05):
    """Count significant proteins per contrast: adj.P.Val < adjp ONLY, split by the sign of
    the fold change. Read from the DE CSVs rather than re-stating whatever the prose claimed
    -- if the two disagree, the reader can see it. It used to also require |log2FC| >= 1,
    the volcano reference line, re-imposing a fold-change filter the DE never applied."""
    rows = []
    if not tables_dir or not os.path.isdir(tables_dir):
        return rows
    for fn in sorted(os.listdir(tables_dir)):
        # The method is one lowercase word (dpc, maxlfq): \w+ also ate the contrast's first
        # word, so "DE_dpc_Old_JPH3.Old_IgG.csv" read as method "dpc_Old", contrast "JPH3...".
        m = re.match(r"^DE_([a-z0-9]+)_(.+)\.csv$", fn)
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
                    if p < adjp:
                        up += lf > 0
                        dn += lf < 0
        except (OSError, csv.Error) as e:
            print(f"[make_analysis_html] WARNING: {fn} unreadable ({e}); left out of the "
                  f"results summary", file=sys.stderr)
            continue
        rows.append({"contrast": m.group(2).replace(".", " vs "), "method": m.group(1),
                     "tested": tot, "up": up, "down": dn, "file": fn})
    return rows


ENGINE_LABEL = {"diann": "DIA-NN", "sage": "Sage", "fragpipe": "FragPipe", "radiant": "Radiant",
                "alphadia": "AlphaDIA"}


def norm_title(t):
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", "", t or "")).strip().lower())


def _load(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError, TypeError):
        return {}


def session_facts(a, prov, submission=None):
    """The header band's study facts, each read from a record of the run -- a fact no record
    holds is left out, never filled in (architectural rule 2). `submission` is the attached
    CoreOmics record's label (submission_report.label)."""
    s = os.path.abspath(a.session) if a.session else None
    man = _load(os.path.join(s, "input", "wf", "workflow.manifest.json")) if s else {}
    fmeta = _load(os.path.join(s, "input", "search.fasta.meta.json")) if s else {}
    sprov = _load(os.path.join(s, "output", "search", "search_provenance.json")) if s else {}
    org = fmeta.get("organism") or man.get("organism")
    tax = fmeta.get("taxid") or man.get("organism_taxid")
    eng = (sprov.get("engine") or (man.get("engine") or {}).get("name") or "").lower()
    ver = sprov.get("version") or (man.get("engine") or {}).get("version")
    return [("Submission", submission),
            ("Organism", f"{org} (taxid {tax})" if org and tax else org),
            ("Instrument", ", ".join(man.get("instruments") or []) or None),
            ("Search", " ".join(x for x in (ENGINE_LABEL.get(eng, eng), str(ver) if ver else "") if x) or None),
            ("DE", prov.get("display_label")),
            ("Report generated", datetime.date.today().isoformat())]


def _make_names(c):
    """R's make.names -- how run_de.R turned a contrast into its DE_<method>_<name>.csv."""
    n = re.sub(r"[^A-Za-z0-9._]", ".", c)
    return n if re.match(r"^([A-Za-z]|\.(?!\d))", n) else "X" + n


def glance(de, prov, adjp, adjp_src, tables):
    """"Results at a glance": design tiles, a tile per contrast, and the inferred-values caveat."""
    label = {_make_names(c): c for c in prov.get("contrasts") or []}
    tiles = []
    if isinstance(prov.get("n_samples"), int):
        tiles.append((prov["n_samples"], "samples", None, "key"))
    if isinstance(prov.get("groups"), dict):
        tiles.append((len(prov["groups"]), "groups", None, "key"))
    if de:
        tiles.append((len(de), "contrasts", None, "key"))
        tiles.append((max(r["tested"] for r in de), "proteins tested", None, "key"))
    out = [rs.stat_tiles(tiles)] if tiles else []
    if de:
        order = list(label)                       # the run's own contrast order
        rows = sorted(de, key=lambda r: (order.index(r["file"][len(f"DE_{r['method']}_"):-4])
                                         if r["file"][len(f"DE_{r['method']}_"):-4] in order
                                         else len(order), r["file"]))
        per = []
        for r in rows:
            raw = r["file"][len(f"DE_{r['method']}_"):-4]
            name = label.get(raw, r["contrast"]).replace("-", " vs ").replace("_", " ")
            per.append((r["up"] + r["down"], name,
                        f"&#9650;&thinsp;{r['up']:,} up &nbsp;&#9660;&thinsp;{r['down']:,} down",
                        None))
        out.append(rs.stat_tiles(per, heading=f"Significant proteins per contrast "
                                              f"(adj. p < {adjp:g})"))
        out.append(f"<p class='lead'>Counted directly from the DE tables: significant = adjusted "
                   f"p &lt; {adjp:g} (Benjamini&ndash;Hochberg; {adjp_src}), the only rule the DE "
                   f"applied &mdash; no fold-change filter. Up / down = sign of the fold "
                   f"change.</p>")
    qc = os.path.join(tables or "", "QC_detected_vs_inferred.csv")
    if os.path.exists(qc):
        try:
            with open(qc, newline="") as fh:
                pct = [float(r["PctInferred"]) for r in csv.DictReader(fh)]
        except (OSError, KeyError, ValueError) as e:
            pct = []
            print(f"[make_analysis_html] WARNING: {qc} unreadable ({e}); the inferred-values "
                  f"note is left out", file=sys.stderr)
        if pct:
            pct.sort()
            med = pct[len(pct) // 2]
            dm = os.path.exists(os.path.join(tables, "Detection_Matrix.csv"))
            out.append(rs.callout(
                "warning" if pct[-1] >= 50 else "info",
                f"<p>limpa's detection-probability model gives every protein a value in every "
                f"sample. Where no precursor of a protein was observed, that value is a model "
                f"estimate, not a measurement: {pct[0]:.0f}&ndash;{pct[-1]:.0f}% of each sample's "
                f"protein values (median {med:.0f}%) are inferred here "
                f"(<code>QC_detected_vs_inferred.csv</code>). Weigh a large fold change carried by "
                f"inferred values accordingly"
                + (" &mdash; <code>Detection_Matrix.csv</code> marks every value." if dm else ".")
                + "</p>",
                title="Some values are inferred, not measured"))
    return "".join(out)


def severity(title, content):
    """Callout severity for the sections a PI must not miss; None for ordinary sections."""
    t = norm_title(title)
    txt = html.unescape(re.sub(r"<[^>]+>", " ", content))
    crit = "⛔" in txt or re.search(r"\bFAIL\b|\bcritical\b|CONFOUNDED", txt, re.I)
    if (t in SUPERSEDED_BY["audit"] or "audit" in t or t in SUPERSEDED_BY["quality"]
            or "quality notes" in t or "expert review" in t or "caveat" in t):
        return "critical" if crit else "warning"
    return None


def split_sections(report_html):
    """-> (content before the first h2, [(anchor, title_html, content)])."""
    parts = re.split(r'<h2 id="([^"]+)">(.*?)</h2>', report_html)
    pre, secs = parts[0], []
    for i in range(1, len(parts), 3):
        secs.append((parts[i], parts[i + 1], parts[i + 2]))
    return pre, secs


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
    ap.add_argument("--submission", help="CoreOmics submission to show (a record file or fetch's "
                                         "folder); default: the one attached to --session")
    ap.add_argument("--adjp", type=float, default=0.05,
                    help="only when the tables carry no de_provenance.json (its adjp wins)")
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
    prov = _load(os.path.join(a.tables, "de_provenance.json")) if a.tables else {}

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

    used, sections = set(), []          # sections: (anchor, title_html, content, kind)
    md, report_html, heads, h1 = "", "", [], None
    if has_report:
        with open(a.report, encoding="utf-8", errors="replace") as fh:
            md = fh.read()
        # Rendered first: its headings decide which tool sections are redundant, and its
        # ids are reserved in `used` so the tool's sections can never collide with them.
        report_html, heads, h1 = md_to_html(md, figs, used=used, drop_h1=True)
    report_h2 = {norm_title(t) for lvl, _, t in heads if lvl == 2}

    adjp, adjp_src = significance_rule(a.tables, a.adjp)
    de = de_summary(a.tables, adjp)
    g = glance(de, prov, adjp, adjp_src, a.tables)
    if g:
        sections.append((anchor("Results at a glance", used), "Results at a glance", g, None))

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
                parts.append(rs.callout("info", "<p>Read these first. They decide how much weight "
                                        "the results below can carry &mdash; a volcano plot looks "
                                        "equally convincing whether or not the run was any good.</p>"))
            for _, fn in sorted(gal[sec]):
                parts.append(figs.render("", fn, section=sec))
            sections.append((anchor(sec, used), html.escape(sec), "".join(parts), None))
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
            content = md_to_html(fh.read(), None, used=used, drop_h1=True)[0]
        sections.append((anchor(title, used), html.escape(title), content,
                         severity(title, content)))

    pre = ""
    if has_report:
        pre, rsecs = split_sections(report_html)
        for aid, title_html, content in rsecs:
            sections.append((aid, title_html, content, severity(title_html, content)))

    if not sections and not pre.strip():
        sys.exit("[make_analysis_html] nothing to render — check --session/--report/--figures")

    # The CoreOmics submission this run answers goes first: ONE record, attached to the session
    # (submission_report.py attach) or named by --submission. The record and its rendering
    # (allowlisted: never a contact or billing field) live in submission_report.py, and its
    # label is the header's Submission fact -- never a number read from other text.
    import submission_report
    try:
        rec, _ = submission_report.resolve(a.submission, a.session)
    except submission_report.RecordError:
        rec = None                      # html_section says why, in the Submission section
    submission = submission_report.html_section(a.submission, a.session)
    if submission:
        sections.insert(0, (anchor("Submission", used), "Submission", submission, None))

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
    # A single short paragraph before the first section is the report's own standfirst.
    subtitle = None
    m = re.fullmatch(r"\s*<p>(.*?)</p>\s*", pre, re.S)
    if m and len(m.group(1)) < 600:
        subtitle, pre = m.group(1), ""
    body = (f'<div class="lead">{pre}</div>' if pre.strip() else "") + "".join(
        rs.section(aid, t, c, k) for aid, t, c, k in sections)
    doc = rs.page(title, body, toc=[(aid, t) for aid, t, _, _ in sections],
                  facts=session_facts(a, prov, submission_report.label(rec) if rec else None),
                  subtitle=subtitle,
                  footer=(f"Self-contained: all {figs.n} figure(s) are embedded, so this one file is "
                          f"the whole report &mdash; no network needed; copy it anywhere and "
                          f"double-click to open. Generated by the UC Davis Proteomics Core "
                          f"pipeline skill (make_analysis_html.py). Click a figure to enlarge it."))

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
