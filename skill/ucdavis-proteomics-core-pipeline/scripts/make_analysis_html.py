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
import argparse, base64, csv, datetime, html, io, json, mimetypes, os, re, sys, urllib.parse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import report_style as rs  # noqa: E402  -- the ONE look shared by the skill's HTML pages
import html_to_pdf  # noqa: E402  -- the PDF: the same page printed by a headless browser
import make_podcast       # noqa: E402  -- an optional audio discussion keeps its Listen card
from session_docs import IP_CONTROL_NAME  # noqa: E402  -- which groups are pull-down controls

# Galleries (no report only): QC first, then overview, then per-contrast results, and the
# p-value calibration check last, as an appendix.
APPENDIX = "Appendix: p-value calibration"
FIGURE_ORDER = [
    ("qc_detected_vs_inferred", "Quality control"),
    ("qc_protein_counts", "Quality control"),
    ("pca", "Overview"),
    ("heatmap_top", "Overview"),
    ("volcano", "Differential expression"),
    ("violin_top", "Differential expression"),     # make_figures.R's top-protein plots
    ("qc_pvalue_panel", APPENDIX),                 # every contrast's p-values, one panel (2.8.0)
    ("pvalue", APPENDIX),                          # one per contrast (sessions before 2.8.0)
]
SECTION_ORDER = ["Quality control", "Overview", "Differential expression", "Other figures",
                 APPENDIX]
# Figures the page adds as the appendix itself when the report does not embed them: a
# calibration check, not a finding, so the narrative never has to interleave it.
APPENDIX_FIGURES = ("qc_pvalue_panel",)

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


def read_figures_json(figures_dir):
    """make_figures.R's figures.json -> {"figures": [{file, type, caption}], "failed": [{file,
    type, reason}], "found", "error"}. Two shapes: the object {"figures": [...], "failed":
    [...]} (2.8.0+) and the bare list of figures older sessions wrote. The ONE reader, for this
    page and for analysis_prompt.py's brief: a schema change cannot leave one of them silently
    figure-less. `error` says why a figures.json present could not be used."""
    out = {"figures": [], "failed": [], "found": False, "error": None}
    path = os.path.join(figures_dir or "", "figures.json")
    if not (figures_dir and os.path.exists(path)):
        return out
    out["found"] = True
    try:
        fj = json.loads(read_text(path))
    except (OSError, ValueError) as e:
        out["error"] = f"unreadable ({e})"
        return out
    if isinstance(fj, list):
        figs, failed = fj, []
    elif isinstance(fj, dict) and isinstance(fj.get("figures"), list):
        figs, failed = fj.get("figures", []), fj.get("failed") or []
    else:
        out["error"] = ("neither a list of figures nor an object with a \"figures\" list "
                        f"(got {type(fj).__name__})")
        return out
    for e in figs:
        if not isinstance(e, dict):
            continue
        k = e.get("file") or e.get("filename") or e.get("name")
        if k:
            out["figures"].append({"file": os.path.basename(k), "type": e.get("type") or "",
                                   "caption": e.get("caption") or e.get("title") or ""})
    for e in failed if isinstance(failed, list) else []:
        if isinstance(e, str):
            e = {"file": e}
        if not isinstance(e, dict):
            continue
        k = e.get("file") or e.get("name")
        if k:
            out["failed"].append({"file": os.path.basename(k), "type": e.get("type") or "",
                                  "reason": str(e.get("reason") or e.get("error") or "")})
    return out


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
    """Turns image references into numbered figures and keeps the ledger the warnings are
    written from. ONE instance per page, shared by both renderers: the HTML and the Markdown
    twin resolve every reference through it, so they give the same figure the same number,
    caption and data summary, and report the same missing / rejected files."""

    def __init__(self, base_dir, root, figures_dir=None, captions=None, summarize=None,
                 suppress=(), current=None, failed=None):
        self.base, self.root, self.figdir = base_dir, root, figures_dir
        self.caps = captions or {}
        self.summarize = summarize or (lambda name: None)
        self.suppress = tuple(suppress)          # figure-name prefixes that say nothing here
        # This run's figures.json: only a figure it lists is embedded. None = no figures.json
        # (nothing to check against). `failed`: name -> why make_figures.R could not draw it.
        self.current = None if current is None else set(current)
        self.failed = dict(failed or {})
        self.n, self.entries = 0, {}
        self.embedded, self.missing, self.rejected, self.suppressed = [], [], [], []
        self.stale = []

    def _locate(self, ref):
        rel = _clean(ref)
        cand = [os.path.realpath(rel if os.path.isabs(rel) else os.path.join(self.base, rel))]
        if self.figdir:           # --figures given apart from --report: look there by name
            cand.append(os.path.realpath(os.path.join(self.figdir, os.path.basename(rel))))
        inside = [c for c in cand if _inside(c, self.root) or
                  (self.figdir and _inside(c, self.figdir))]
        return next((c for c in inside if os.path.isfile(c) and os.access(c, os.R_OK)), None), \
            bool(inside)

    def resolve(self, alt, ref):
        """-> entry dict. The first sight of a file numbers it; later sights return it."""
        name = os.path.basename(_clean(ref)) or ref
        if self.suppress and name.startswith(self.suppress):
            if name not in self.suppressed:
                self.suppressed.append(name)
            return {"status": "suppressed", "name": name}
        if _external(ref):
            key, e = ("x", ref), {"status": "rejected", "name": name,
                                  "text": f"figure not embedded (not a file in this session): {ref}"}
        else:
            path, allowed = self._locate(ref)
            if not allowed:
                key, e = ("x", ref), {"status": "rejected", "name": name,
                                      "text": f"figure not embedded (outside the session folder): {ref}"}
            elif name in self.failed:
                key, e = ("m", name), {"status": "missing", "name": name,
                                       "text": f"figure could not be drawn in this run: {name}"
                                               + (f" ({self.failed[name]})" if self.failed[name]
                                                  else "")}
            elif path is None:
                # A report written for the per-contrast p-value histograms, rendered against a
                # run that draws them as one panel: say where they went.
                moved = (name.startswith("pvalue_") and self.current is not None and
                         any(c.startswith(APPENDIX_FIGURES) for c in self.current))
                key, e = ("m", name), {"status": "missing", "name": name,
                                       "text": f"figure missing: {name}"
                                               + (f" (this run draws every contrast's p-values "
                                                  f"in one panel: see {APPENDIX})" if moved
                                                  else "")}
            elif self.current is not None and os.path.basename(path) not in self.current:
                # On disk but not in this run's figures.json: left over from an earlier run,
                # so it may show data this report no longer describes.
                key, e = ("s", name), {"status": "stale", "name": name,
                                       "text": f"figure not shown: {name} is not in this run's "
                                               f"figures.json (a file left from an earlier run)"}
            else:
                key = ("f", path)
                if key not in self.entries:
                    self.n += 1
                    nm = os.path.basename(path)
                    e = {"status": "ok", "n": self.n, "path": path, "name": nm, "alt": alt,
                         "caption": self.caps.get(nm) or "", "summary": self.summarize(nm)}
                    self.embedded.append(nm)
        if key not in self.entries:
            self.entries[key] = e
            if e["status"] == "missing":
                self.missing.append(name)
            elif e["status"] == "rejected":
                self.rejected.append(ref)
            elif e["status"] == "stale":
                self.stale.append(name)
        return self.entries[key]

    def html(self, alt, ref, emitted):
        e = self.resolve(alt, ref)
        if e["status"] == "suppressed":
            return ""
        if e["status"] != "ok":
            return rs.note(e["text"])
        n = e["n"]
        if n in emitted:
            return f'<p class="figref">(See <a href="#fig-{n}">Figure {n}</a>.)</p>'
        emitted.add(n)
        return rs.figure_card(n, data_uri(e["path"]), alt or e["caption"],
                              md_inline(alt) if alt else "",
                              md_inline(e["caption"]) if e["caption"] else "",
                              summary_html=md_inline(e["summary"]) if e["summary"] else "")

    def md(self, alt, ref, emitted, md_dir):
        """The Markdown twin of a figure: NotebookLM and the like read text, not images, so
        the caption and a data summary taken from the tables travel with the reference."""
        e = self.resolve(alt, ref)
        if e["status"] == "suppressed":
            return ""
        if e["status"] != "ok":
            return f"> **Note:** {e['text']}"
        n = e["n"]
        if n in emitted:
            return f"(See Figure {n}.)"
        emitted.add(n)
        title = (alt or e["caption"] or e["name"]).strip().rstrip(".")
        # realpath on both sides: macOS /var and /tmp are symlinks into /private
        rel = os.path.relpath(e["path"], os.path.realpath(md_dir)).replace(os.sep, "/")
        out = [f"![Figure {n}. {title}]({rel})", "", f"**Figure {n}. {title}.**"
               + (f" {e['caption']}" if e["caption"] and e["caption"] != alt else "")]
        if e["summary"]:
            out += ["", f"*Data in this figure:* {e['summary']}"]
        return "\n".join(out)


class _HtmlFigs:
    """Adapter for md_to_html: one render pass, so repeats become "See Figure N"."""

    def __init__(self, figs, emitted):
        self.figs, self.emitted = figs, emitted

    def render(self, alt, ref, section=None):
        return self.figs.html(alt, ref, self.emitted)


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
    fold-change filter"), so it is here too. Source "--adjp" means NOTHING recorded it: the
    page then says so and tags the value (architectural rule 2), never "the rule the DE
    applied"."""
    adjp = load_record(os.path.join(tables_dir or "", "de_provenance.json")).get("adjp")
    if isinstance(adjp, (int, float)) and not isinstance(adjp, bool):
        return float(adjp), "de_provenance.json"
    return default_adjp, "--adjp"


DEFAULT_TAG = "(DEFAULT — not user-confirmed)"


def rule_text(adjp, src):
    """The significance sentence, as recorded -- or tagged when nothing recorded it."""
    if src == "de_provenance.json":
        return (f"Significant = adjusted p < {adjp:g} (Benjamini–Hochberg within each contrast; "
                f"de_provenance.json), the only rule the DE applied — no fold-change filter. "
                f"Up / down = sign of the fold change. {LOGFC_DIRECTION}")
    return (f"Significant here = adjusted p < {adjp:g} {DEFAULT_TAG}: no de_provenance.json "
            f"records the rule the DE applied, so this page cannot say which cutoff — or whether "
            f"a fold-change filter — it used. Up / down = sign of the fold change. "
            f"{LOGFC_DIRECTION}")


ENGINE_LABEL = {"diann": "DIA-NN", "sage": "Sage", "fragpipe": "FragPipe", "radiant": "Radiant",
                "alphadia": "AlphaDIA"}


def norm_title(t):
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", "", t or "")).strip().lower())


# Records that exist but could not be read, {path: why}: said on stderr when it happens and on
# the page (Results at a glance) -- never a silent {}. A de_provenance.json the platform's
# locale could not decode once dropped the contaminant-database warning without a word.
_UNREADABLE = {}


def _unreadable(path, why):
    if path not in _UNREADABLE:
        _UNREADABLE[path] = why
        print(f"[{os.path.basename(sys.argv[0] or 'make_analysis_html.py')}] WARNING: {path} "
              f"{why}", file=sys.stderr)


def read_text(path):
    """A record's text, decoded as UTF-8 (a BOM dropped) whatever the platform's locale: the
    skill writes UTF-8, which Windows' cp1252 default mis-decodes or rejects. A file that is not
    valid UTF-8 is still read, bad bytes replaced, and recorded as unreadable."""
    with open(path, "rb") as fh:
        raw = fh.read()
    try:
        return raw.decode("utf-8-sig")
    except UnicodeDecodeError as e:
        _unreadable(path, f"is not UTF-8 (byte {e.start}): read with the bad bytes replaced")
        return raw.decode("utf-8-sig", errors="replace")


def csv_text(path):
    """read_text() as a file for csv.reader / csv.DictReader."""
    return io.StringIO(read_text(path), newline="")


def load_record(path):
    """A JSON record -> dict. {} when there is no such file (nothing was recorded); a file that
    exists but cannot be read, parsed, or is not a JSON object gives {} AND is recorded as
    unreadable (stderr + the page), so what it holds is said to be missing, not absent."""
    if not path or not os.path.exists(path):
        return {}
    try:
        v = json.loads(read_text(path))
    except (OSError, ValueError) as e:
        _unreadable(path, f"could not be read ({type(e).__name__}: {e})")
        return {}
    if not isinstance(v, dict):
        _unreadable(path, f"is not a JSON object ({type(v).__name__})")
        return {}
    return v


def unreadable_note():
    """The fixed callout for records that could not be read, or None."""
    if not _UNREADABLE:
        return None
    return {"kind": "warning",
            "title": f"{len(_UNREADABLE)} record{'' if len(_UNREADABLE) == 1 else 's'} of this "
                     f"run could not be read",
            "text": "; ".join(f"`{os.path.basename(p)}` {why}" for p, why in _UNREADABLE.items())
                    + ". What this page would show from them (a recorded cutoff, a contaminant "
                      "warning, a study fact) may be missing here — missing from the page, not "
                      "from the run. Re-save the file as UTF-8 and render again."}


def session_facts(a, prov, submission=None):
    """The header band's study facts, each read from a record of the run -- a fact no record
    holds is left out, never filled in (architectural rule 2). `submission` is the attached
    CoreOmics record's label (submission_report.label)."""
    s = os.path.abspath(a.session) if a.session else None
    man = load_record(os.path.join(s, "input", "wf", "workflow.manifest.json")) if s else {}
    fmeta = load_record(os.path.join(s, "input", "search.fasta.meta.json")) if s else {}
    sprov = load_record(os.path.join(s, "output", "search", "search_provenance.json")) if s else {}
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


def contrast_label(c):
    """THE display form of a contrast, used everywhere a contrast is named (tiles, tables,
    figure data, the analysis brief): "Old_JPH3-Old_IgG" -> "Old JPH3 vs Old IgG"."""
    return c.replace("-", " vs ").replace("_", " ")


# Which way a fold change points, said once and shown wherever log2FC values are.
LOGFC_DIRECTION = "Positive log2FC = higher in the first-named group of the contrast."


def _make_names(c):
    """R's make.names -- how run_de.R turned a contrast into its DE_<method>_<name>.csv."""
    n = re.sub(r"[^A-Za-z0-9._]", ".", c)
    return n if re.match(r"^([A-Za-z]|\.(?!\d))", n) else "X" + n


class Tables:
    """The DE tables, read once and shared by every place that quotes a number from them --
    the glance, the figure summaries, the top-protein tables -- so all agree."""

    def __init__(self, tables_dir, prov, adjp, src):
        self.dir, self.prov, self.adjp, self.src = tables_dir, prov, adjp, src
        self._rows = {}
        self.label = {_make_names(c): c for c in prov.get("contrasts") or []}
        self.files = {}                        # make.names contrast -> DE csv
        if tables_dir and os.path.isdir(tables_dir):
            for fn in sorted(os.listdir(tables_dir)):
                m = re.match(r"^DE_([a-z0-9]+)_(.+)\.csv$", fn)
                if m:
                    self.files[m.group(2)] = fn
        order = list(self.label)
        self.contrasts = sorted(self.files, key=lambda c: (order.index(c) if c in order
                                                           else len(order), c))

    def display(self, c):
        return contrast_label(self.label.get(c, c.replace(".", "-")))

    def rows(self, c):
        if c not in self._rows:
            out = []
            for r in csv.DictReader(csv_text(os.path.join(self.dir, self.files[c]))):
                try:
                    r["_p"] = float(r.get("adj.P.Val") or r.get("padj") or r.get("FDR"))
                    r["_lfc"] = float(r.get("logFC") or r.get("log2FoldChange"))
                except (TypeError, ValueError):
                    r["_p"], r["_lfc"] = None, None
                try:
                    r["_raw_p"] = float(r.get("P.Value") or r.get("pvalue"))
                except (TypeError, ValueError):
                    r["_raw_p"] = None
                out.append(r)
            self._rows[c] = out
        return self._rows[c]

    def counts(self, c):
        """(tested, up, down) at adj.P < adjp ONLY -- the only significance rule the DE applied.
        The page used to also require |log2FC| >= 1, the volcano reference line, re-imposing a
        fold-change filter the DE never applied."""
        rows = self.rows(c)
        sig = [r for r in rows if r["_p"] is not None and r["_p"] < self.adjp]
        return len(rows), sum(r["_lfc"] > 0 for r in sig), sum(r["_lfc"] < 0 for r in sig)

    def top(self, c, k):
        """The k best by adj.P, ties broken by the raw p-value (topTable's order, as the brief
        and make_figures.R rank them) -- never by |log2FC|: BH ties are common, and the largest
        fold changes there are presence calls, not the strongest evidence. A row without a raw
        p-value keeps its place in the table among its ties (sorted() is stable)."""
        rows = [r for r in self.rows(c) if r["_p"] is not None]
        return sorted(rows, key=lambda r: (r["_p"], r["_raw_p"] if r["_raw_p"] is not None
                                           else float("inf")))[:k]


def gene_label(r):
    g = (r.get("Genes") or "").split(";")[0].strip()
    return g or (r.get("Protein.Group") or "?").split(";")[0]


def fmt_p(p):
    return f"{p:.3g}" if p >= 1e-3 else f"{p:.2e}"


# Figures that say nothing when the expression matrix is complete by construction: every
# bar of "proteins quantified per sample" is the same number (Silva08172026: 30 bars of
# 6,112). make_figures.R no longer draws it then, but an older figures/ folder -- and the
# narrative written against it -- still can; Brett flagged it twice (2026-09-25).
SUPPRESS_WHEN_COMPLETE = ("qc_protein_counts",)
_COUNTS_PARAGRAPH = re.compile(r"^\s*\*\*\s*Proteins?\s+(?:quantified\s+)?per\s+sample\b", re.I)


def matrix_complete(tables_dir, prov):
    """-> (complete, why). make_figures.R's test -- every sample has the same number of
    non-missing proteins in Expression_Matrix.csv -- and, when that file is not at hand,
    the pipeline itself: a DPC-Quant (limpa) matrix is complete by construction."""
    em = os.path.join(tables_dir or "", "Expression_Matrix.csv")
    if os.path.exists(em):
        try:
            with csv_text(em) as fh:
                rd = csv.reader(fh)
                head = next(rd)
                cols = [i for i, h in enumerate(head)
                        if h not in ("Protein.Group", "Genes", "Protein.Names")]
                n = [0] * len(cols)
                for rec in rd:
                    for k, i in enumerate(cols):
                        if i < len(rec) and rec[i] not in ("", "NA", "NaN"):
                            n[k] += 1
            if n:
                return (len(set(n)) == 1,
                        f"every sample has {n[0]:,} proteins in Expression_Matrix.csv"
                        if len(set(n)) == 1 else None)
        except (OSError, StopIteration, csv.Error):
            pass
    if (prov.get("pipeline_id") or prov.get("method")) == "dpc":
        return True, "a DPC-Quant (limpa) matrix gives every protein a value in every sample"
    return False, None


def drop_suppressed(md, prefixes):
    """The report text without references to suppressed figures, and without a paragraph
    right after one that only describes that plot (it opens "**Proteins per sample.**" or
    "**Proteins quantified per sample**"). Anything else is left exactly as written.
    -> (text, [suppressed figure names], paragraphs dropped)."""
    if not prefixes:
        return md, [], 0
    out, after, dropping, names, n_par = [], False, False, [], 0
    for ln in md.splitlines():
        if dropping:
            if ln.strip():
                continue
            dropping = False
        hits = [m for m in IMG_REF.finditer(ln)
                if os.path.basename(_clean(_ref(m)[1])).startswith(prefixes)]
        if hits:
            for m in reversed(hits):
                nm = os.path.basename(_clean(_ref(m)[1]))
                if nm not in names:
                    names.append(nm)
                ln = ln[:m.start()] + ln[m.end():]
            after = True
            if not ln.strip():
                continue
        elif after and ln.strip():
            after = False
            if _COUNTS_PARAGRAPH.match(ln):
                dropping, n_par = True, n_par + 1
                continue
        out.append(ln)
    return "\n".join(out) + ("\n" if md.endswith("\n") else ""), names, n_par


def raw_pvalues(tables, c):
    out = []
    for r in tables.rows(c):
        try:
            out.append(float(r.get("P.Value")))
        except (TypeError, ValueError):
            pass
    return out


def figure_summary(name, tables, qc_path, em_path):
    """A short factual summary of the data a figure draws, taken from the tables -- never a
    number that is not in them. None for figures without one (PCA, heatmap)."""
    stem = os.path.splitext(name)[0]
    if stem.startswith(APPENDIX_FIGURES) and tables and tables.contrasts:
        parts = []
        for c in tables.contrasts:
            raw = raw_pvalues(tables, c)
            if raw:
                parts.append(f"{tables.display(c)} {100 * sum(x < 0.05 for x in raw) / len(raw):.0f}%")
        if parts:
            return ("Share of raw p-values below 0.05 per contrast (about 5% when nothing "
                    "changed): " + ", ".join(parts) + ".")
    for prefix in ("volcano_", "pvalue_"):
        c = stem[len(prefix):] if stem.startswith(prefix) else None
        if c and tables and c in tables.files:
            tested, up, dn = tables.counts(c)
            if prefix == "volcano_":
                top = "; ".join(f"{gene_label(r)} (log2FC {r['_lfc']:+.2f}, adj.P {fmt_p(r['_p'])})"
                                for r in tables.top(c, 5))
                return (f"{tables.display(c)}: {up + dn:,} of {tested:,} proteins significant at "
                        f"adj.P < {tables.adjp:g} ({up:,} up, {dn:,} down; no fold-change "
                        f"filter). Top 5 by adj.P: {top}.")
            raw = raw_pvalues(tables, c)
            if raw:
                below = sum(x < 0.05 for x in raw)
                return (f"{tables.display(c)}: {len(raw):,} p-values, {below:,} "
                        f"({100 * below / len(raw):.0f}%) below 0.05; {up + dn:,} significant "
                        f"at adj.P < {tables.adjp:g}.")
    if stem.startswith("qc_detected_vs_inferred") and qc_path and os.path.exists(qc_path):
        q = list(csv.DictReader(csv_text(qc_path)))
        if q:
            det = sorted(int(float(r["Detected"])) for r in q)
            pct = sorted(float(r["PctInferred"]) for r in q)
            return (f"{len(q)} samples; proteins detected per sample {det[0]:,}–{det[-1]:,} "
                    f"(median {det[len(det) // 2]:,}) of {int(float(q[0]['Total'])):,}; "
                    f"inferred {pct[0]:.0f}–{pct[-1]:.0f}% per sample.")
    if stem.startswith("qc_protein_counts") and em_path and os.path.exists(em_path):
        with csv_text(em_path) as fh:
            rd = csv.reader(fh)
            head = next(rd)
            cols = [i for i, h in enumerate(head) if h not in ("Protein.Group", "Genes", "Protein.Names")]
            n = [0] * len(cols)
            for rec in rd:
                for k, i in enumerate(cols):
                    if i < len(rec) and rec[i] not in ("", "NA", "NaN"):
                        n[k] += 1
        if n:
            return (f"proteins with a value per sample: {min(n):,}–{max(n):,} across {len(n)} "
                    f"samples" + (" (identical: the matrix is complete by construction)"
                                  if min(n) == max(n) else "") + ".")
    return None


# Rows a reader must not take for the sample's biology: the antibody's own chains and
# common-contaminant entries. Human IMGT symbols (IGHG1, IGKC, IGLV1-40, JCHAIN) and their
# mouse forms (Ighg2c, Igkc, Iglv1) alike; IgLON (Iglon5) and IGF are NOT Ig chains.
_IG_CHAIN = re.compile(r"^(IGH[GAMDEVJ]|IGK[CVJ]|IGL[CVJ]|JCHAIN\b)", re.I)


def background_flag(gene, protein):
    """"Ig chain" / "contaminant" / None for one DE row."""
    if any(t.strip().startswith("Cont_") for t in (protein or "").split(";")):
        return "contaminant"
    if _IG_CHAIN.match((gene or "").split(";")[0].strip()):
        return "Ig chain"
    return None


def detected_columns(row, contrast):
    """ "k/n A, k/n B" from the DE table's own Detected_<group> columns (run_de.R: k of n
    samples in the group with the protein measured), or None when the table has none."""
    parts = [g.strip() for g in (contrast or "").split("-")]
    if len(parts) != 2:
        return None
    vals = [row.get(f"Detected_{g}") for g in parts]
    if not all(v not in (None, "", "NA") for v in vals):
        return None
    return ", ".join(f"{v} {contrast_label(g)}" for v, g in zip(vals, parts))


def detection_note_fn(tables_dir, prov, conditions):
    """-> f(protein, contrast) giving "detected 3/3 Old_Kv21, 0/3 Old_IgG" from
    Detection_Matrix.csv, or None when there is no matrix. The word follows the record:
    "detected" when 0 means inferred (dpc), "quantified" when 0 means missing (maxlfq)."""
    path = os.path.join(tables_dir or "", "Detection_Matrix.csv")
    if not os.path.exists(path):
        return None
    word = "quantified" if (prov.get("detection_matrix") or {}).get("zero_means") == "missing" \
        else "detected"
    rd = csv.reader(csv_text(path))
    head = next(rd)
    det = {rec[0]: {head[i]: rec[i] for i in range(1, len(rec))} for rec in rd}
    groups = {}
    if conditions and os.path.exists(conditions):
        for r in csv.DictReader(csv_text(conditions)):
            groups.setdefault((r.get("Group") or "").strip(), []).append(
                (r.get("File.Name") or "").strip())

    def on(v):
        try:
            return float(v) > 0
        except (TypeError, ValueError):
            return False

    def note(protein, contrast, row=None):
        own = detected_columns(row or {}, contrast)
        if own:                           # the DE table's own per-group counts win
            ev = ((row or {}).get("Evidence") or "").strip()   # run_de.R's own category
            return f"{word} {own}" + (f" — {ev}" if ev and ev != "NA" else "")
        d = det.get(protein)
        if d is None:
            return "not recorded"
        parts = [g.strip() for g in contrast.split("-")] if contrast else []
        if len(parts) == 2 and all(g in groups for g in parts):
            return f"{word} " + ", ".join(
                f"{sum(on(d.get(s)) for s in groups[g])}/{len(groups[g])} {contrast_label(g)}"
                for g in parts)
        return f"{word} in {sum(on(v) for v in d.values())}/{len(d)} samples"
    return note


def build_page(a, prov, tables, figs, md_text, used):
    """Assemble the page ONCE, as data. render_html() and render_md() draw this same list,
    so the HTML report and its Markdown twin cannot say different things (rule 3).
    -> (title, subtitle_md, sections); section = {anchor, title, kind, blocks}, block =
    ("md", text) | ("glance", data) | ("gallery", section, [names]) | ("top", data) |
    ("submission", html, markdown) -- the CoreOmics record, main() puts it first."""
    sections, h1, pre = [], None, ""
    report_secs = []
    if md_text is not None:
        h1, pre, report_secs = split_md_sections(md_text)
    report_h2 = {norm_title(t) for t, _ in report_secs}
    referenced = set(report_figures(md_text)) if md_text is not None else set()
    glance_anchor = anchor("Results at a glance", used)
    # Every record the page reads is read before the glance is assembled, so its "could not be
    # read" callout lists them all (a record read later would be on stderr only).
    extras = []
    for path, title, key in ((a.quality, "Sample quality notes", "quality"),
                             (a.audit, "Audit & caveats", "audit")):
        if not (path and os.path.exists(path)):
            continue
        if report_h2 & (SUPERSEDED_BY[key] | {norm_title(title)}):
            continue                    # the report has its own -- never show it twice
        extras.append((title, strip_h1(read_text(path))))
    top = top_data(tables, a, prov)
    g = glance_data(prov, tables, a.tables, a.session)
    for n in (not_drawn_note(figs.failed, referenced), unreadable_note()):
        if n:
            g["notes"].append(n)
    if g["tiles"] or g["contrasts"] or g["notes"]:
        sections.append({"anchor": glance_anchor,
                         "title": "Results at a glance", "kind": None, "blocks": [("glance", g)]})
    if md_text is None:
        gal = {}
        for fn in a._listed:
            if figs.suppress and fn.startswith(figs.suppress):
                figs.resolve("", fn)               # recorded as suppressed
                continue
            sec, rank = classify(fn)
            gal.setdefault(sec, []).append((rank, fn))
        for sec in SECTION_ORDER:
            if sec in gal:
                sections.append({"anchor": anchor(sec, used), "title": sec, "kind": None,
                                 "blocks": [("gallery", sec, [fn for _, fn in sorted(gal[sec])])]})
    for title, body in extras:
        sections.append({"anchor": anchor(title, used), "title": title,
                         "kind": severity(title, body), "blocks": [("md", body)]})
    for title, body in report_secs:
        sections.append({"anchor": anchor(title, used), "title": title,
                         "kind": severity(title, body), "blocks": [("md", body)]})
    if top:
        sections.append({"anchor": anchor("Top proteins per contrast", used),
                         "title": "Top proteins per contrast", "kind": None,
                         "blocks": [("top", top)]})
    if md_text is not None:
        app = [fn for fn in a._listed if fn.startswith(APPENDIX_FIGURES) and fn not in referenced]
        if app:
            sections.append({"anchor": anchor(APPENDIX, used), "title": APPENDIX, "kind": None,
                             "blocks": [("gallery", APPENDIX, app)]})
    subtitle = None
    if pre.strip() and len(pre.strip()) < 600 and "\n\n" not in pre.strip():
        subtitle, pre = pre.strip(), ""
    if pre.strip():                      # anything else before the first section stays first
        sections.insert(0, {"anchor": anchor("Introduction", used), "title": "Introduction",
                            "kind": None, "blocks": [("md", pre)]})
    return h1, subtitle, sections


def strip_h1(text):
    return re.sub(r"\A\s*#\s+[^\n]*\n", "", text, count=1)


def split_md_sections(md):
    """-> (h1 text, text before the first ## heading, [(h2 title, body)]), fences respected."""
    h1, pre, secs, cur, fence = None, [], [], None, False
    for ln in md.splitlines():
        if ln.startswith("```"):
            fence = not fence
        m = None if fence else re.match(r"^(#{1,2})\s+(.*)", ln)
        if m and len(m.group(1)) == 1 and h1 is None and cur is None:
            h1 = m.group(2).strip()
            continue
        if m and len(m.group(1)) == 2:
            cur = [m.group(2).strip(), []]
            secs.append(cur)
            continue
        (cur[1] if cur else pre).append(ln)
    return h1, "\n".join(pre), [(t, "\n".join(b)) for t, b in secs]


def glance_data(prov, tables, tables_dir, session=None):
    tiles = []
    if isinstance(prov.get("n_samples"), int):
        tiles.append((prov["n_samples"], "samples"))
    if isinstance(prov.get("groups"), dict):
        tiles.append((len(prov["groups"]), "groups"))
    rows = []
    for c in tables.contrasts:
        tested, up, dn = tables.counts(c)
        rows.append({"contrast": tables.display(c), "tested": tested, "up": up, "down": dn})
    if rows:
        tiles.append((len(rows), "contrasts"))
        tiles.append((max(r["tested"] for r in rows), "proteins tested"))
    notes = [n for n in (inferred_note(prov, tables_dir), database_note(prov, tables_dir, session))
             if n]
    return {"tiles": tiles, "contrasts": rows, "adjp": tables.adjp, "adjp_src": tables.src,
            "recorded": tables.src == "de_provenance.json", "notes": notes,
            "rule": rule_text(tables.adjp, tables.src)}


def not_drawn_note(failed, referenced=()):
    """The fixed "not drawn in this run" callout: every figure make_figures.R's figures.json
    lists under `failed`, with its reason, so a missing PCA or top-protein plot is never just
    absent. A figure the report references already carries its own note where it would be."""
    items = [(n, why) for n, why in (failed or {}).items() if n not in set(referenced)]
    if not items:
        return None
    return {"kind": "warning",
            "title": f"{len(items)} figure{'' if len(items) == 1 else 's'} could not be drawn "
                     f"in this run",
            "text": "; ".join(f"`{n}`" + (f": {why.rstrip('.')}" if why else "")
                              for n, why in items)
                    + " (make_figures.R's `failed` list in `figures/figures.json`)."}


def inferred_note(prov, tables_dir):
    """The fixed "inferred, not measured" callout. Only for a DPC-Quant run (pipeline_id dpc):
    it is the detection-probability model that supplies inferred values, and its own record
    (missing_policy) says how -- the page does not restate the policy in its own words."""
    if (prov.get("pipeline_id") or prov.get("method")) != "dpc":
        return None
    qc = os.path.join(tables_dir or "", "QC_detected_vs_inferred.csv")
    if not os.path.exists(qc):
        return None
    try:
        pct = sorted(float(r["PctInferred"]) for r in csv.DictReader(csv_text(qc)))
    except (OSError, KeyError, ValueError) as e:
        print(f"[make_analysis_html] WARNING: {qc} unreadable ({e}); the inferred-values "
              f"note is left out", file=sys.stderr)
        return None
    if not pct:
        return None
    dm = os.path.exists(os.path.join(tables_dir, "Detection_Matrix.csv"))
    policy = (prov.get("missing_policy") or "").strip().rstrip(".")
    return {"kind": "warning" if pct[-1] >= 50 else "info",
            "title": "Some values are inferred, not measured",
            "text": ((f"{policy}. " if policy else
                      f"Missing-value policy: not recorded in de_provenance.json. ")
                     + f"Where no precursor of a protein was observed in a sample, its value there "
                       f"is a model estimate, not a measurement: {pct[0]:.0f}–{pct[-1]:.0f}% of "
                       f"each sample's protein values (median {pct[len(pct) // 2]:.0f}%) are "
                       f"inferred in this run (`QC_detected_vs_inferred.csv`). A fold change "
                       f"where one group was never measured is a detection event, not a measured "
                       f"magnitude"
                     + (" — `Detection_Matrix.csv` marks every value." if dm else "."))}


def database_note(prov, tables_dir, session=None):
    """The fixed contaminant-database caveat: run_de.R's record says real proteins may have been
    removed as contaminants (database_risk). Shown by the page itself -- it must never depend
    on the report writer remembering it -- and naming the proteins: the ones fetch_fasta.py can
    show are identical to target proteins when the searched FASTA is still readable, else every
    contaminant group the filter removed (some of which are then real)."""
    cont = prov.get("contaminants") if isinstance(prov.get("contaminants"), dict) else {}
    if cont.get("database_risk") is not True:
        return None
    removed = []
    table = cont.get("removed_table")
    path = os.path.join(tables_dir or "", table) if table else None
    if path and os.path.exists(path):
        removed = [r for r in csv.DictReader(csv_text(path))
                   if (r.get("Contaminant.Group") or "").upper() == "TRUE"]
    named, how, extra = [], "", ""
    meta_path = next((m for m in (cont.get("fasta_meta"),
                                  os.path.join(session, "input", "search.fasta.meta.json")
                                  if session else None) if m and os.path.exists(m)), None)
    if meta_path and removed:
        try:
            import fetch_fasta
            meta = load_record(meta_path)     # an unreadable one is said, on the page too
            if not meta:
                raise ValueError(f"{os.path.basename(meta_path)} could not be read")
            tc = fetch_fasta.target_contaminants(meta)
            named = fetch_fasta.seen_only_as_cont(
                tc["kept_as_contaminant"],
                [({t.strip().upper() for t in (r.get("Protein.Group") or "").split(";")},
                  r.get("Genes") or "?") for r in removed])
            if named:
                how = "identical to target proteins, re-checked in the searched FASTA"
            elif tc.get("legacy_note"):
                extra = "the searched FASTA could not be re-checked"
        except Exception as e:                           # said, not swallowed
            extra = f"the FASTA check could not run: {type(e).__name__}: {e}"
    if not named and removed:
        named = sorted({(r.get("Genes") or r.get("Protein.Group") or "?").split(";")[0]
                        for r in removed}, key=str.lower)
        how = ("every contaminant group the filter removed; which of them are real cannot be "
               "told from here" + (f" ({extra})" if extra else ""))
    text = str(cont.get("database_note") or "").strip()
    text = text[:1].upper() + text[1:]
    if named:
        text += f" Proteins affected ({how}; {len(named)}): {', '.join(named)}."
    return {"kind": "warning", "title": "Real proteins may have been removed as contaminants",
            "text": text}


TOP_K = 20


def top_data(tables, a, prov):
    """Top TOP_K proteins by adj.P per contrast, with a measured-vs-inferred note."""
    if not tables.contrasts:
        return []
    cond = os.path.join(a.session, "input", "conditions.csv") if a.session else None
    note = detection_note_fn(a.tables, prov, cond)
    word = "quantified" if (prov.get("detection_matrix") or {}).get("zero_means") == "missing" \
        else "detected"
    out = []
    for c in tables.contrasts:
        raw = tables.label.get(c)
        rows = []
        for r in tables.top(c, TOP_K):
            det = note(r.get("Protein.Group"), raw, r) if note else None
            if det is None and detected_columns(r, raw):
                det = f"{word} {detected_columns(r, raw)}"
            rows.append({"protein": (r.get("Protein.Group") or "?"), "gene": gene_label(r),
                         "lfc": r["_lfc"], "p": r["_p"], "sig": r["_p"] < tables.adjp,
                         "det": det, "flag": background_flag(gene_label(r), r.get("Protein.Group"))})
        out.append({"contrast": tables.display(c), "rows": rows, "counts": tables.counts(c),
                    "adjp": tables.adjp})
    return out


def severity(title, content):
    """Callout severity for the sections a PI must not miss; None for ordinary sections."""
    t = norm_title(title)
    txt = html.unescape(re.sub(r"<[^>]+>", " ", content))
    crit = "⛔" in txt or re.search(r"\bFAIL\b|\bcritical\b|CONFOUNDED", txt, re.I)
    if (t in SUPERSEDED_BY["audit"] or "audit" in t or t in SUPERSEDED_BY["quality"]
            or "quality notes" in t or "expert review" in t or "caveat" in t):
        return "critical" if crit else "warning"
    return None


# ------------------------------------------------------------------ renderer 1: HTML
def render_html(title, subtitle, sections, figs, facts, used):
    emitted = set()
    hf = _HtmlFigs(figs, emitted)
    body = []
    for sec in sections:
        parts = []
        for b in sec["blocks"]:
            if b[0] == "md":
                parts.append(md_to_html(b[1], hf, used=used)[0])
            elif b[0] == "glance":
                parts.append(glance_html(b[1]))
            elif b[0] == "gallery":
                if b[1] == "Quality control":
                    parts.append(rs.callout("info", "<p>Read these first. They decide how much "
                                            "weight the results below can carry &mdash; a volcano "
                                            "plot looks equally convincing whether or not the run "
                                            "was any good.</p>"))
                parts += [figs.html("", fn, emitted) for fn in b[2]]
            elif b[0] == "top":
                parts.append(top_html(b[1]))
            elif b[0] == "submission":
                parts.append(b[1])
        body.append(rs.section(sec["anchor"], md_inline(sec["title"]), "".join(parts), sec["kind"]))
    return rs.page(title, "".join(body),
                   toc=[(s["anchor"], md_inline(s["title"])) for s in sections],
                   facts=facts, subtitle=md_inline(subtitle) if subtitle else None,
                   footer=(f"Self-contained: all {figs.n} figure(s) are embedded, so this one file "
                           f"is the whole report &mdash; no network needed; copy it anywhere and "
                           f"double-click to open. A plain-text twin, Analysis_Report.md, and a "
                           f"PDF, Analysis_Report.pdf, carry the same content. Generated by the UC Davis "
                           f"Proteomics Core pipeline skill (make_analysis_html.py). Click a figure "
                           f"to enlarge it."))


def glance_html(g):
    out = [rs.stat_tiles([(v, lab, None, "key") for v, lab in g["tiles"]])] if g["tiles"] else []
    if g["contrasts"]:
        out.append(rs.stat_tiles(
            [(r["up"] + r["down"], r["contrast"],
              f"&#9650;&thinsp;{r['up']:,} up &nbsp;&#9660;&thinsp;{r['down']:,} down", None)
             for r in g["contrasts"]],
            heading=f"Significant proteins per contrast (adj. p < {g['adjp']:g}"
                    + ("" if g["recorded"] else ", default — not recorded") + ")"))
        out.append(f"<p class='lead'>Counted directly from the DE tables. {md_inline(g['rule'])}</p>")
    for i in g["notes"]:
        out.append(rs.callout(i["kind"], f"<p>{md_inline(i['text'])}</p>", title=i["title"]))
    return "".join(out)


def _top_lead():
    return (f"The {TOP_K} proteins with the smallest adjusted p per contrast, from the DE tables "
            f"(the full lists are the DE_*.csv files). {LOGFC_DIRECTION} A row flagged Ig chain "
            f"(the antibody's own chains) or contaminant is not the sample's biology.")


def _ns_divider(t):
    return f"not significant below this line (adj.P ≥ {t['adjp']:g})"


def top_html(top):
    out = [f"<p class='lead'>{md_inline(_top_lead())}</p>"]
    for t in top:
        tested, up, dn = t["counts"]
        det = any(r["det"] for r in t["rows"])
        flag = any(r["flag"] for r in t["rows"])
        head = (["Protein", "Gene", "log2FC", "adj.P"] + (["Measured"] if det else [])
                + (["Note"] if flag else []))
        body, prev_sig = [], True
        for r in t["rows"]:
            if prev_sig and not r["sig"]:
                body.append({"divider": _ns_divider(t)})
            prev_sig = r["sig"]
            body.append([rs.esc(r["protein"]), rs.esc(r["gene"]), f"{r['lfc']:+.2f}", fmt_p(r["p"])]
                        + ([rs.esc(r["det"] or "")] if det else [])
                        + ([f"<b>{rs.esc(r['flag'])}</b>" if r["flag"] else ""] if flag else []))
        # open: a closed <details> does not print, and the PDF must carry these tables
        out.append(f"<details open><summary><strong>{rs.esc(t['contrast'])}</strong> &mdash; "
                   f"{up + dn:,} significant ({up:,} up, {dn:,} down)</summary>"
                   f"{rs.table([rs.esc(h) for h in head], body)}</details>")
    return "".join(out)


# ------------------------------------------------------------------ renderer 2: Markdown
_TAGS = re.compile(r"</?(?:br|img|div|span|p|b|i|em|strong|sup|sub|code|a|table|thead|tbody|"
                   r"tr|td|th|details|summary|font|u|hr|section|figure|figcaption)\b[^>]*>", re.I)


def md_expand(text, figs, emitted, md_dir):
    """The report's Markdown with every image reference replaced by its text-bearing figure
    block (image link + caption + data summary), and any HTML tag stripped."""
    out, fence = [], False
    for ln in text.splitlines():
        if ln.startswith("```"):
            fence = not fence
            out.append(ln)
            continue
        if fence or not IMG_REF.search(ln):
            out.append(_TAGS.sub("", ln) if not fence else ln)
            continue
        if ln.lstrip().startswith("|"):          # inside a table: a reference, not a block
            out.append(_TAGS.sub("", IMG_REF.sub(
                lambda m: _fig_word(figs.resolve(*_ref(m))), ln)))
            continue
        # A figure is a block: blank lines around it, so neither the text before nor the
        # line after (often the next figure, or a list) runs into it.
        last = 0
        for m in IMG_REF.finditer(ln):
            text = _TAGS.sub("", ln[last:m.start()]).strip()
            if text:
                out += ["", text]
            out += ["", figs.md(*_ref(m), emitted, md_dir), ""]
            last = m.end()
        text = _TAGS.sub("", ln[last:]).strip()
        if text:
            out += [text, ""]
    return "\n".join(out)


def _fig_word(e):
    if e["status"] == "suppressed":
        return ""
    return f"Figure {e['n']}" if e["status"] == "ok" else f"[{e['text']}]"


def quote(text):
    return "\n".join(("> " + ln) if ln.strip() else ">" for ln in text.strip("\n").splitlines())


def render_md(title, subtitle, sections, figs, facts, md_dir):
    emitted = set()
    L = [f"# {title}", ""]
    fl = " · ".join(f"**{k}:** {v}" for k, v in facts if v not in (None, ""))
    if fl:
        L += [fl, ""]
    if subtitle:
        L += [_TAGS.sub("", subtitle), ""]
    L += ["*This is the plain-text twin of Analysis_Report.html, made from the same sections, for "
          "NotebookLM or other AI notebooks: each figure's caption and the numbers it shows are "
          "written out as text.*", ""]
    for sec in sections:
        L += [f"## {sec['title']}", ""]
        parts = []
        for b in sec["blocks"]:
            if b[0] == "md":
                parts.append(md_expand(b[1], figs, emitted, md_dir).strip("\n"))
            elif b[0] == "glance":
                parts.append(glance_md(b[1]))
            elif b[0] == "gallery":
                parts += [figs.md("", fn, emitted, md_dir) for fn in b[2]]
            elif b[0] == "top":
                parts.append(top_md(b[1]))
            elif b[0] == "submission":
                parts.append(b[2])
        body = "\n\n".join(p for p in parts if p.strip())
        if sec["kind"]:
            body = quote(f"**{rs.CALLOUT_KINDS[sec['kind']][1]}.**\n\n{body}")
        L += [body, ""]
    doc = "\n".join(L)
    return re.sub(r"\n{3,}", "\n\n", doc).rstrip() + "\n"


def glance_md(g):
    out = []
    if g["tiles"]:
        out.append(" · ".join(f"**{v:,}** {lab}" for v, lab in g["tiles"]))
    if g["contrasts"]:
        rows = ["| Contrast | Significant | Up | Down | Tested |", "|---|---:|---:|---:|---:|"]
        rows += [f"| {r['contrast']} | {r['up'] + r['down']:,} | {r['up']:,} | {r['down']:,} | "
                 f"{r['tested']:,} |" for r in g["contrasts"]]
        out.append("\n".join(rows))
        out.append(g["rule"])
    for i in g["notes"]:
        out.append(quote(f"**{rs.CALLOUT_KINDS[i['kind']][1]}:** {i['title']}. {i['text']}"))
    return "\n\n".join(out)


def top_md(top):
    out = [_top_lead()]
    for t in top:
        tested, up, dn = t["counts"]
        det = any(r["det"] for r in t["rows"])
        flag = any(r["flag"] for r in t["rows"])
        ncol = 4 + det + flag
        rows = [f"### {t['contrast']}", "",
                f"{up + dn:,} significant ({up:,} up, {dn:,} down) of {tested:,} tested.", "",
                "| Protein | Gene | log2FC | adj.P |" + (" Measured |" if det else "")
                + (" Note |" if flag else ""),
                "|---|---|---:|---:|" + ("---|" if det else "") + ("---|" if flag else "")]
        prev_sig = True
        for r in t["rows"]:
            if prev_sig and not r["sig"]:
                rows.append(f"| *{_ns_divider(t)}* |" + " |" * (ncol - 1))
            prev_sig = r["sig"]
            rows.append(f"| {r['protein']} | {r['gene']} | {r['lfc']:+.2f} | {fmt_p(r['p'])} |"
                        + (f" {r['det'] or ''} |" if det else "")
                        + (f" **{r['flag']}** |" if r["flag"] else " |" if flag else ""))
        out.append("\n".join(rows))
    return "\n\n".join(out)


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
    ap.add_argument("--md-out", help="the Markdown twin (default: --out with .md, e.g. "
                                     "output/Analysis_Report.md)")
    ap.add_argument("--no-pdf", action="store_true",
                    help="skip the PDF (default: print the HTML to --out with .pdf with a "
                         "headless Chrome/Chromium/Edge when one is installed)")
    a = ap.parse_args()
    _UNREADABLE.clear()                  # one page per run

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
    prov = load_record(os.path.join(a.tables, "de_provenance.json")) if a.tables else {}
    md_out = a.md_out or (os.path.splitext(a.out)[0] + ".md")
    if os.path.abspath(md_out) == os.path.abspath(a.report or ""):
        sys.exit("[make_analysis_html] --md-out would overwrite the report it is made from")

    fjd = read_figures_json(a.figures)
    if fjd["error"]:
        print(f"[make_analysis_html] WARNING: figures.json {fjd['error']}; captions and its "
              f"figure list are not used", file=sys.stderr)
    have_fj = fjd["found"] and not fjd["error"]
    a._listed = [f["file"] for f in fjd["figures"]]
    caps = {f["file"]: f["caption"] for f in fjd["figures"]}
    # make_figures.R's `failed` list: figures it tried and could not draw this run
    failed = {f["file"]: f["reason"] for f in fjd["failed"]}
    available = (sorted(fn for fn in os.listdir(a.figures) if fn.lower().endswith(IMAGE_EXT))
                 if a.figures and os.path.isdir(a.figures) else [])

    adjp, adjp_src = significance_rule(a.tables, a.adjp)
    tables = Tables(a.tables, prov, adjp, adjp_src)
    qc = os.path.join(a.tables or "", "QC_detected_vs_inferred.csv")
    em = os.path.join(a.tables or "", "Expression_Matrix.csv")
    # Images resolve relative to the report's folder and must stay inside the session.
    base = os.path.dirname(os.path.abspath(a.report)) if has_report else (a.figures or ".")
    root = os.path.abspath(a.session) if a.session else base
    complete, complete_why = matrix_complete(a.tables, prov)
    figs = Figures(base, root, a.figures, caps,
                   summarize=lambda nm: figure_summary(nm, tables, qc, em),
                   suppress=SUPPRESS_WHEN_COMPLETE if complete else (),
                   current=a._listed if have_fj else None, failed=failed)

    md_text = md_orig = None
    n_par = 0
    if has_report:
        md_text = md_orig = make_podcast.strip_block(read_text(a.report))  # its Listen line: the card is below
        # The one source both renderers draw from, so the HTML, .md and PDF all lose it.
        md_text, gone, n_par = drop_suppressed(md_text, figs.suppress)
        figs.suppressed += [g for g in gone if g not in figs.suppressed]
    # The CoreOmics submission this run answers goes first: ONE record, attached to the session
    # (submission_report.py attach) or named by --submission. The record and its rendering
    # (allowlisted: never a contact or billing field) live in submission_report.py, and its
    # label is the header's Submission fact -- never a number read from other text.
    import submission_report
    sub = submission_report.report_section(a.submission, a.session)
    # The session's records (manifest, FASTA meta, search provenance) and the QC table a figure
    # summary reads are read now, before the page is assembled: see build_page().
    facts = session_facts(a, prov, sub["label"] if sub else None)
    if os.path.exists(qc):
        read_text(qc)
    used = set()
    sub_anchor = anchor("Submission", used) if sub else None
    h1, subtitle, sections = build_page(a, prov, tables, figs, md_text, used)
    if not sections and not subtitle and md_text is None:
        sys.exit("[make_analysis_html] nothing to render — check --session/--report/--figures")
    title = a.title or (re.sub(r"[*`]", "", h1) if h1 else "Proteomics Analysis Report")
    if sub:
        sections.insert(0, {"anchor": sub_anchor, "title": "Submission",
                            "kind": None, "blocks": [("submission", sub["html"], sub["md"])]})

    doc = render_html(title, subtitle, sections, figs, facts, used)
    md_doc = render_md(title, subtitle, sections, figs, facts,
                       os.path.dirname(os.path.abspath(md_out)))

    if has_report:
        referenced = set(report_figures(md_orig))
        left_out = [f for f in available if f not in set(figs.embedded) and f not in referenced]
        why = f"not referenced in {os.path.basename(a.report)}"
    else:
        left_out = [f for f in available if f not in set(a._listed)
                    and not (figs.suppress and f.startswith(figs.suppress))]
        why = "not listed in figures.json" if a._listed else "no report and no figures.json"
    if left_out:
        print(f"[make_analysis_html] WARNING: {len(left_out)} image(s) in {a.figures} are "
              f"{why} and were NOT embedded: {', '.join(left_out)}", file=sys.stderr)
    if figs.suppressed:
        print(f"[make_analysis_html] left out {', '.join(figs.suppressed)}"
              + (f" and the report paragraph describing it" if has_report and n_par else "")
              + f": {complete_why}, so every bar is the same number; "
              f"qc_detected_vs_inferred.png is the per-sample depth view", file=sys.stderr)
    if figs.missing:
        print(f"[make_analysis_html] WARNING: {len(figs.missing)} referenced image(s) are "
              f"missing and are shown as a 'figure missing' note: {', '.join(figs.missing)}",
              file=sys.stderr)
    if figs.rejected:
        print(f"[make_analysis_html] WARNING: {len(figs.rejected)} image reference(s) point "
              f"outside the session and were not embedded: {', '.join(figs.rejected)}",
              file=sys.stderr)
    if figs.stale:
        print(f"[make_analysis_html] WARNING: {len(figs.stale)} referenced image(s) are not in "
              f"this run's figures.json (left from an earlier run) and are shown as a note, not "
              f"embedded: {', '.join(figs.stale)}", file=sys.stderr)

    # The optional podcast's Listen card and line (make_podcast.py): added before the files are
    # written, so the HTML, its .md twin and the PDF printed from the HTML all carry it.
    outdir = os.path.dirname(os.path.abspath(a.out))
    doc = make_podcast.add_listen_card(doc, outdir)
    md_doc = make_podcast.add_listen_md(md_doc, outdir)
    for path, text in ((a.out, doc), (md_out, md_doc)):
        os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
    # The PDF: the same HTML through its print stylesheet, so it carries the figures -- for
    # NotebookLM (which reads a PDF's images too) and for printing. Never fatal: without a
    # browser (HIVE) it says so, and finalize retries on the laptop.
    pdf_out = os.path.splitext(a.out)[0] + ".pdf"
    if a.no_pdf:
        pdf_status, pdf_note = "INFO", "--no-pdf was given"
    else:
        # print_report(): a failed print never leaves an older PDF looking current
        pdf_status, pdf_note = html_to_pdf.print_report(a.out, pdf_out)
    pdf_ok = pdf_status == "OK"
    print(f"[make_analysis_html] {'PDF: ' + pdf_out + ' (' + pdf_note + ')' if pdf_ok else pdf_status + ': no PDF -- ' + pdf_note}",
          file=sys.stderr)
    print(json.dumps({"wrote": a.out, "bytes": os.path.getsize(a.out),
                      "markdown_twin": md_out, "markdown_bytes": os.path.getsize(md_out),
                      "pdf": pdf_out if pdf_ok else None, "pdf_note": pdf_note,
                      "figures_embedded": figs.n,
                      "figure_list": os.path.basename(a.report) if has_report else "figures.json",
                      "figures_not_embedded": left_out, "figures_missing": figs.missing,
                      "figures_rejected": figs.rejected, "figures_suppressed": figs.suppressed,
                      "figures_stale": figs.stale, "pdf_status": pdf_status,
                      "records_unreadable": dict(_UNREADABLE),
                      "contrasts": len(tables.contrasts), "self_contained": True}, indent=2))


if __name__ == "__main__":
    main()
