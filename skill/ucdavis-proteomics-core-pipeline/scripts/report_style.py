"""
report_style.py -- the ONE look for every HTML page the skill writes.

Analysis_Report.html (make_analysis_html.py) is a page(); README.html and HOW_TO_SUBMIT.html
(make_deposit.md_to_html) are a document(), so a restyle happens here once, not in three copies.

Everything is inline -- no fonts, CSS or JS fetched from anywhere -- so a page opens by
double-click from a zip or a network share, on Windows too, with no network. Stdlib only.

  page(title, body, toc=..., facts=..., subtitle=...)   the whole document
  document(title, body)                                  a plain Markdown page, no script
  header_band(title, facts, subtitle)                    study facts under the title
  stat_tiles([(value, label, sub, kind)])                "results at a glance"
  callout(kind, body, title=None)                        info / warning / critical
  figure_card(n, uri, title_html, caption_html)          numbered, click to enlarge
  note(text)                                             a visible "figure missing" note
  table(head_cells, body_rows)                           zebra, sticky header, numbers right
  section(anchor, title_html, content, kind=None)        an h2 section, optionally a callout

Design: system fonts, a ~72ch reading column, :root colour tokens with a real dark mode
(figures stay on white cards so plots remain legible), a sticky contents rail on wide screens
that folds above the text on narrow ones, no horizontal page scroll at phone width (wide
tables scroll inside their own box), and a print stylesheet (no rail, figures unsplit,
callout colours kept).
"""
import html
import re

CALLOUT_KINDS = {"info": ("ℹ", "Note"), "warning": ("⚠", "Warning"), "critical": ("⛔", "Critical")}

CSS = r"""
:root{
  --bg:#f6f7f9;--surface:#ffffff;--surface-2:#f0f2f5;--fg:#17191c;--muted:#5d6570;--line:#dfe3e8;
  --accent:#2463c9;--accent-weak:#e7eefb;--header-bg:#15253f;--header-fg:#f4f7fb;--header-muted:#b8c4d6;
  --info:#2463c9;--info-bg:#eaf1fc;--warning:#a8590a;--warning-bg:#fdf3e4;--critical:#b42318;--critical-bg:#fdecea;
  --zebra:#f7f8fa;--shadow:0 1px 2px rgba(16,24,40,.06),0 1px 3px rgba(16,24,40,.08);
  --radius:12px;--measure:72ch;
  --font:system-ui,-apple-system,"Segoe UI",Roboto,"Helvetica Neue",Arial,"Noto Sans",sans-serif;
  --mono:ui-monospace,SFMono-Regular,Menlo,Consolas,"Liberation Mono",monospace;
  color-scheme:light;
}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
  --bg:#0f1318;--surface:#171c23;--surface-2:#1e252e;--fg:#e6e9ee;--muted:#9aa5b3;--line:#2a323d;
  --accent:#7aa7f5;--accent-weak:#1c2a42;--header-bg:#0b1628;--header-fg:#eef2f8;--header-muted:#9fb0c8;
  --info:#7aa7f5;--info-bg:#15223a;--warning:#f0a04b;--warning-bg:#2b2014;--critical:#f2837a;--critical-bg:#321716;
  --zebra:#1a2029;--shadow:0 1px 2px rgba(0,0,0,.4);color-scheme:dark;}}
:root[data-theme="dark"]{
  --bg:#0f1318;--surface:#171c23;--surface-2:#1e252e;--fg:#e6e9ee;--muted:#9aa5b3;--line:#2a323d;
  --accent:#7aa7f5;--accent-weak:#1c2a42;--header-bg:#0b1628;--header-fg:#eef2f8;--header-muted:#9fb0c8;
  --info:#7aa7f5;--info-bg:#15223a;--warning:#f0a04b;--warning-bg:#2b2014;--critical:#f2837a;--critical-bg:#321716;
  --zebra:#1a2029;--shadow:0 1px 2px rgba(0,0,0,.4);color-scheme:dark;}
*,*::before,*::after{box-sizing:border-box}
html{-webkit-text-size-adjust:100%;scroll-padding-top:1rem}
body{margin:0;background:var(--bg);color:var(--fg);font:16px/1.65 var(--font);overflow-wrap:break-word}
a{color:var(--accent)}
.skip{position:absolute;left:-999px}.skip:focus{left:1rem;top:1rem;z-index:20;background:var(--surface);padding:.5rem}
/* header band */
.band{background:var(--header-bg);color:var(--header-fg);padding:2rem 1rem 1.6rem}
.band-in{max-width:78rem;margin:0 auto;position:relative}
.band h1{font-size:clamp(1.45rem,1.1rem + 1.6vw,2.15rem);line-height:1.2;margin:0 0 .35rem;font-weight:700;letter-spacing:-.01em;padding-right:6.5rem}
.band .subtitle{color:var(--header-muted);margin:0 0 1rem;max-width:var(--measure)}
.facts{display:flex;flex-wrap:wrap;gap:.5rem .6rem;margin:0;padding:0;list-style:none}
.facts li{background:rgba(255,255,255,.08);border:1px solid rgba(255,255,255,.14);border-radius:999px;padding:.22rem .75rem;font-size:.86rem}
.facts b{color:var(--header-muted);font-weight:600;margin-right:.35rem}
.toggle{position:absolute;top:0;right:0;background:rgba(255,255,255,.1);color:var(--header-fg);border:1px solid rgba(255,255,255,.25);border-radius:8px;padding:.35rem .7rem;font:inherit;font-size:.82rem;cursor:pointer}
/* layout: contents rail + reading column */
.layout{max-width:78rem;margin:0 auto;padding:1.5rem 1rem 4rem;display:block}
.toc{background:var(--surface);border:1px solid var(--line);border-radius:var(--radius);padding:.2rem .9rem;margin:0 0 1.5rem;box-shadow:var(--shadow)}
.toc summary{cursor:pointer;font-weight:650;padding:.55rem 0;list-style:none}
.toc summary::-webkit-details-marker{display:none}
.toc summary::after{content:" \25BE";color:var(--muted)}
.toc ol{margin:0 0 .7rem;padding:0;list-style:none;counter-reset:s}
.toc li{margin:0}
.toc a{display:block;padding:.28rem .5rem;border-radius:6px;color:var(--fg);text-decoration:none;font-size:.9rem;line-height:1.35}
.toc a:hover{background:var(--accent-weak);color:var(--accent)}
main{min-width:0}
@media (min-width:1100px){
  .layout{display:grid;grid-template-columns:15.5rem minmax(0,1fr);gap:2.5rem;align-items:start}
  .toc{position:sticky;top:1rem;max-height:calc(100vh - 2rem);overflow:auto;margin:0}
  .toc summary{pointer-events:none}.toc summary::after{content:""}
}
/* typography */
main h2{font-size:1.5rem;line-height:1.25;margin:2.8rem 0 1rem;padding-bottom:.4rem;border-bottom:1px solid var(--line);letter-spacing:-.005em}
main h3{font-size:1.15rem;margin:1.9rem 0 .5rem}
main h4{font-size:1rem;margin:1.4rem 0 .4rem}
main p,main li,main blockquote{max-width:var(--measure)}
main section:first-child h2{margin-top:.5rem}
code{font-family:var(--mono);font-size:.86em;background:var(--surface-2);border:1px solid var(--line);border-radius:5px;padding:.08em .35em;overflow-wrap:anywhere}
pre{font-family:var(--mono);background:var(--surface-2);border:1px solid var(--line);border-radius:10px;padding:.9rem 1rem;overflow-x:auto;font-size:.86rem}
pre code{background:none;border:0;padding:0}
.lead{color:var(--muted);max-width:var(--measure)}
.muted{color:var(--muted)}
/* stat tiles */
.tiles{display:grid;grid-template-columns:repeat(auto-fill,minmax(10.5rem,1fr));gap:.75rem;margin:1rem 0 1.2rem}
.tile{background:var(--surface);border:1px solid var(--line);border-radius:var(--radius);padding:.85rem 1rem;box-shadow:var(--shadow);min-width:0}
.tile .v{font-size:1.75rem;font-weight:700;line-height:1.1;font-variant-numeric:tabular-nums;letter-spacing:-.01em}
.tile .l{color:var(--muted);font-size:.83rem;margin-top:.25rem;line-height:1.3}
.tile .s{font-size:.8rem;margin-top:.35rem;font-variant-numeric:tabular-nums;color:var(--muted)}
.tile.key{border-top:3px solid var(--accent)}
.tiles-h{font-size:.78rem;text-transform:uppercase;letter-spacing:.06em;color:var(--muted);margin:1.2rem 0 -.3rem;font-weight:650}
/* callouts */
.callout{--c:var(--info);--cb:var(--info-bg);background:var(--cb);border:1px solid var(--line);border-left:5px solid var(--c);border-radius:10px;padding:.9rem 1.1rem;margin:1.1rem 0;max-width:calc(var(--measure) + 4rem)}
.callout.warning{--c:var(--warning);--cb:var(--warning-bg)}
.callout.critical{--c:var(--critical);--cb:var(--critical-bg)}
.callout > .ch{display:flex;align-items:center;gap:.45rem;font-weight:700;color:var(--c);margin:0 0 .35rem;font-size:.95rem}
.callout > .ch .ic{font-size:1.05rem;line-height:1}
.callout > .ch .kind{font-size:.72rem;text-transform:uppercase;letter-spacing:.07em;border:1px solid var(--c);border-radius:999px;padding:.02rem .45rem}
.callout p:last-child,.callout ul:last-child{margin-bottom:0}
.callout p:first-of-type{margin-top:.2rem}
section.callout-sec > .callout{max-width:none}
/* figures */
figure.fig{margin:1.6rem 0;background:var(--surface);border:1px solid var(--line);border-radius:var(--radius);padding:.8rem;box-shadow:var(--shadow);max-width:52rem;break-inside:avoid}
figure.fig .imgbox{display:block;width:100%;background:#fff;border-radius:8px;padding:.35rem;border:0;cursor:zoom-in}
figure.fig img{display:block;width:100%;height:auto;border-radius:4px}
figure.fig figcaption{font-size:.9rem;color:var(--muted);margin:.7rem .2rem .1rem;line-height:1.5}
figure.fig .fignum{font-weight:700;color:var(--fg);margin-right:.35rem}
figure.fig .figt{color:var(--fg);font-weight:600}
figure.fig .figdata{display:block;margin-top:.45rem;padding-top:.45rem;border-top:1px dashed var(--line);font-variant-numeric:tabular-nums}
details{margin:.6rem 0}details > summary{cursor:pointer;padding:.35rem 0}
.figmissing{margin:1.2rem 0;padding:.75rem 1rem;border:1.5px dashed var(--warning);border-radius:10px;color:var(--warning);background:var(--warning-bg);font-size:.92rem;max-width:52rem}
.figref{color:var(--muted);font-size:.92rem}
/* lightbox */
.lb{position:fixed;inset:0;background:rgba(8,12,20,.88);display:none;align-items:center;justify-content:center;padding:1.2rem;z-index:50;cursor:zoom-out}
.lb.open{display:flex;flex-direction:column}
.lb img{max-width:100%;max-height:85vh;background:#fff;border-radius:8px;padding:.4rem}
.lb p{color:#e8edf5;max-width:60rem;text-align:center;font-size:.92rem;margin:.8rem 0 0}
/* tables */
.tablewrap{overflow:auto;max-height:70vh;margin:1rem 0 1.3rem;border:1px solid var(--line);border-radius:10px;background:var(--surface);box-shadow:var(--shadow)}
table{border-collapse:separate;border-spacing:0;width:100%;font-size:.9rem;font-variant-numeric:tabular-nums}
th,td{text-align:left;padding:.5rem .75rem;border-bottom:1px solid var(--line);vertical-align:top}
thead th{position:sticky;top:0;z-index:1;background:var(--surface-2);font-size:.78rem;text-transform:uppercase;letter-spacing:.04em;color:var(--muted);font-weight:650;white-space:nowrap}
tbody tr:nth-child(even) td{background:var(--zebra)}
tbody tr:last-child td{border-bottom:0}
td.num,th.num{text-align:right;white-space:nowrap}
tr.divider td{background:var(--surface-2) !important;color:var(--muted);font-size:.78rem;text-transform:uppercase;letter-spacing:.05em;text-align:center;border-top:2px solid var(--line)}
.scrollhint{font-size:.8rem;color:var(--muted);margin:.8rem 0 -.6rem;text-align:right}
blockquote{margin:1rem 0;padding:.6rem 1rem;border-left:4px solid var(--info);background:var(--info-bg);border-radius:8px}
footer.foot{max-width:78rem;margin:0 auto;padding:0 1rem 2.5rem;color:var(--muted);font-size:.82rem}
/* the CoreOmics submission (submission_report.render_html): the form as a definition list */
.subm .sub,.subm dt,.subm .blank{color:var(--muted)}
.subm dl{display:grid;grid-template-columns:max-content 1fr;gap:.3rem 1.2rem;margin:.8rem 0;background:var(--surface);border:1px solid var(--line);border-left:3px solid var(--accent);border-radius:10px;padding:.9rem 1.2rem;max-width:calc(var(--measure) + 4rem)}
.subm dt{font-size:.85rem}
.subm dd{margin:0;white-space:pre-line}
.subm blockquote{white-space:pre-line}
.subm ol.sheet{columns:17rem;column-gap:1.6rem;padding-left:0;list-style:none;font-size:.9rem}
.subm ol.sheet li{break-inside:avoid}
.subm ol.sheet code{margin-right:.5rem}
@media (max-width:40rem){.subm dl{grid-template-columns:1fr}}
@media (max-width:600px){
  .band{padding:1.4rem 1rem 1.2rem}
  .tile .v{font-size:1.45rem}
  .tiles{grid-template-columns:repeat(auto-fill,minmax(8.5rem,1fr))}
  figure.fig{padding:.5rem}
}
@media print{
  /* Light tokens for paper whatever the screen's theme: (0,2,0) and later in the sheet, so
     they beat both dark rules -- a dark-mode reader's printout is not grey on white. */
  :root:not([data-print]),:root[data-theme]{
    --bg:#fff;--surface:#fff;--surface-2:#f2f4f7;--fg:#111;--muted:#444;--line:#cfd5dc;
    --accent:#1d4fa3;--accent-weak:#eef3fb;--header-bg:#fff;--header-fg:#111;--header-muted:#444;
    --info:#1d4fa3;--info-bg:#eef3fb;--warning:#8a4a06;--warning-bg:#fdf3e4;--critical:#a01d12;
    --critical-bg:#fdecea;--zebra:#f6f7f9;--shadow:none;color-scheme:light}
  @page{margin:14mm 12mm}
  body{background:#fff;color:#000;font-size:11pt}
  .toc,.toggle,.lb,.skip,.scrollhint{display:none !important}
  .layout{display:block;padding:0;max-width:none}
  .band{background:#fff;color:#000;border-bottom:2px solid #000;padding:0 0 .8rem;margin-bottom:.8rem}
  .band h1{padding-right:0}
  .band .subtitle,.facts b{color:#333}
  .facts li{border-color:#999;background:#fff}
  figure.fig,.callout,.tile,tr,.figmissing{break-inside:avoid}
  figure.fig{max-width:none;box-shadow:none}
  figure.fig img{max-height:22cm;width:auto;max-width:100%;margin:0 auto}
  main p,main li{orphans:3;widows:3}
  h2,h3,details > summary{break-after:avoid}
  .tablewrap{max-height:none;overflow:visible;box-shadow:none}
  thead th{position:static}
  .callout,.figmissing,.tile.key,tbody tr:nth-child(even) td{-webkit-print-color-adjust:exact;print-color-adjust:exact}
  a{color:#000;text-decoration:none}
}
"""


# A plain document page (document()): Markdown rendered as header / nav / section /
# aside.callout / table by make_deposit.md_to_html -- no header band, contents rail or script.
# Everything else (tokens, dark mode, type, code, callouts, tables, print) is CSS above.
DOC_CSS = r"""
main.doc{max-width:60rem;margin:0 auto;padding:1.5rem 1rem 3rem}
main.doc > header{border-bottom:3px solid var(--accent);padding-bottom:.4rem;margin-bottom:1.2rem}
main.doc h1{font-size:clamp(1.45rem,1.1rem + 1.6vw,2.15rem);line-height:1.2;margin:.4rem 0 .5rem;letter-spacing:-.01em}
main.doc > nav{background:var(--surface);border:1px solid var(--line);border-radius:var(--radius);padding:0 1.1rem .4rem;margin:1.2rem 0;box-shadow:var(--shadow)}
main.doc > nav h2{border-bottom:0;margin-top:.9rem}
main.doc table{display:block;overflow-x:auto;width:max-content;max-width:100%;margin:1rem 0 1.3rem;border:1px solid var(--line);border-radius:10px;background:var(--surface)}
main.doc td code{word-break:break-all}
@media print{main.doc{max-width:none;padding:0}main.doc > nav{box-shadow:none}main.doc table{display:table;overflow:visible}}
"""

JS = r"""
(function(){
 var root=document.documentElement,b=document.getElementById('theme');
 function cur(){var a=root.getAttribute('data-theme');if(a)return a;
  return window.matchMedia&&matchMedia('(prefers-color-scheme: dark)').matches?'dark':'light';}
 function set(t){root.setAttribute('data-theme',t);if(b)b.textContent=t==='dark'?'☀ Light':'☾ Dark';}
 if(b){set(cur());b.addEventListener('click',function(){set(cur()==='dark'?'light':'dark');});}
 // a table wider than the screen scrolls in its own box: say so, or it just looks cut off
 function hints(){var w=document.querySelectorAll('.tablewrap');for(var i=0;i<w.length;i++){
  var t=w[i],h=t.previousElementSibling,has=h&&h.className==='scrollhint';
  if(t.scrollWidth>t.clientWidth+4){if(!has){h=document.createElement('div');h.className='scrollhint';
   h.textContent='\u21c6 wider than the screen: swipe the table sideways';t.parentNode.insertBefore(h,t);}}
  else if(has){h.parentNode.removeChild(h);}}}
 hints();window.addEventListener('resize',hints);
 var toc=document.querySelector('details.toc');
 if(toc&&window.matchMedia&&!matchMedia('(min-width: 1100px)').matches)toc.removeAttribute('open');
 var lb=document.getElementById('lb'),li,lc;
 if(lb){li=document.createElement('img');li.alt='';lc=document.createElement('p');lb.appendChild(li);lb.appendChild(lc);}
 function close(){if(lb){lb.classList.remove('open');li.removeAttribute('src');}}
 document.addEventListener('click',function(e){
  var bx=e.target.closest&&e.target.closest('.imgbox');
  if(bx&&lb){var im=bx.querySelector('img'),fc=bx.parentNode.querySelector('figcaption');
   li.src=im.src;li.alt=im.alt;lc.textContent=fc?fc.textContent:'';lb.classList.add('open');e.preventDefault();return;}
  if(lb&&lb.classList.contains('open')&&lb.contains(e.target))close();});
 document.addEventListener('keydown',function(e){if(e.key==='Escape')close();});
})();
"""

_NUM = re.compile(r"^[<>≤≥~±+\-−]?\s*[\d.,]+(?:e[+-]?\d+)?\s*%?$", re.I)


def esc(s):
    return html.escape(str(s if s is not None else ""))


def is_number(text):
    """True for a cell that reads as a number ("6,112", "0.05", "−1.3", "54%", "<0.001")."""
    t = re.sub(r"<[^>]+>", "", text or "").strip()
    return bool(t) and bool(_NUM.match(t))


def header_band(title, facts=(), subtitle=None):
    items = "".join(f"<li><b>{esc(k)}</b>{esc(v)}</li>" for k, v in facts if v not in (None, ""))
    return (f'<header class="band"><div class="band-in">'
            f'<button class="toggle" id="theme" type="button">Dark</button>'
            f"<h1>{esc(title)}</h1>"
            + (f'<p class="subtitle">{subtitle}</p>' if subtitle else "")
            + (f'<ul class="facts">{items}</ul>' if items else "")
            + "</div></header>")


def stat_tiles(tiles, heading=None):
    """tiles: [(value, label, sub_html or None, "key" or None)]"""
    out = [f'<div class="tiles-h">{esc(heading)}</div>'] if heading else []
    out.append('<div class="tiles">')
    for value, label, sub, kind in tiles:
        v = f"{value:,}" if isinstance(value, int) else esc(value)
        out.append(f'<div class="tile{" " + kind if kind else ""}"><div class="v">{v}</div>'
                   f'<div class="l">{esc(label)}</div>'
                   + (f'<div class="s">{sub}</div>' if sub else "") + "</div>")
    out.append("</div>")
    return "".join(out)


def callout(kind, body_html, title=None):
    kind = kind if kind in CALLOUT_KINDS else "info"
    icon, word = CALLOUT_KINDS[kind]
    head = (f'<div class="ch"><span class="ic" aria-hidden="true">{icon}</span>'
            f'<span class="kind">{word}</span>{" " + esc(title) if title else ""}</div>')
    return f'<div class="callout {kind}" role="note">{head}{body_html}</div>'


def figure_card(n, uri, alt, title_html="", caption_html="", summary_html=""):
    sep = " &mdash; " if title_html and caption_html else ""
    return (f'<figure class="fig" id="fig-{n}">'
            f'<button class="imgbox" type="button" aria-label="Enlarge figure {n}">'
            f'<img src="{uri}" alt="{esc(alt)}"></button>'
            f'<figcaption><span class="fignum">Figure {n}.</span>'
            f'<span class="figt">{title_html}</span>{sep}{caption_html}'
            + (f'<span class="figdata"><b>Data:</b> {summary_html}</span>' if summary_html else "")
            + "</figcaption></figure>")


def note(text):
    return f'<div class="figmissing" role="note">&#9888; {esc(text)}</div>'


def table(head_cells, body_rows):
    """head_cells / body_rows hold already-rendered inline HTML. A column whose body cells
    all read as numbers is right-aligned, header included. A body row given as
    {"divider": text} is a full-width divider line (e.g. "not significant below")."""
    ncol = len(head_cells)
    cells = [r for r in body_rows if not isinstance(r, dict)]
    numeric = [bool(cells) and all(is_number(r[j]) for r in cells if j < len(r) and
                                   re.sub(r"<[^>]+>", "", r[j]).strip() not in ("", "—", "-", "NA"))
               and any(j < len(r) and is_number(r[j]) for r in cells) for j in range(ncol)]
    cls = lambda j: ' class="num"' if j < ncol and numeric[j] else ""   # noqa: E731
    t = ['<div class="tablewrap"><table><thead><tr>']
    t += [f"<th{cls(j)}>{c}</th>" for j, c in enumerate(head_cells)]
    t.append("</tr></thead><tbody>")
    for r in body_rows:
        if isinstance(r, dict):
            t.append(f'<tr class="divider"><td colspan="{ncol}">{esc(r["divider"])}</td></tr>')
            continue
        t.append("<tr>" + "".join(f"<td{cls(j)}>{c}</td>" for j, c in enumerate(r)) + "</tr>")
    t.append("</tbody></table></div>")
    return "".join(t)


def section(anchor, title_html, content, kind=None):
    """An h2 section. With `kind`, its content sits in a callout of that severity -- for the
    sections a PI must not miss (data-quality notes, audit, expert review)."""
    body = callout(kind, content) if kind else content
    # not inlined into the f-string: a quote of the f-string's own kind inside {} needs
    # Python 3.12, and macOS's /usr/bin/python3 is 3.9
    cls = ' class="callout-sec"' if kind else ""
    return (f'<section id="sec-{anchor}"{cls}>'
            f'<h2 id="{anchor}">{title_html}</h2>{body}</section>')


def page(title, body, toc=(), facts=(), subtitle=None, footer=None):
    """The whole self-contained document. toc: [(anchor, label_html)]."""
    nav = ""
    if toc:
        nav = ('<details class="toc" open><summary>Contents</summary><ol>'
               + "".join(f'<li><a href="#{a}">{lab}</a></li>' for a, lab in toc)
               + "</ol></details>")
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="light dark">
<title>{esc(title)}</title><style>{CSS}</style></head>
<body><a class="skip" href="#main">Skip to content</a>
{header_band(title, facts, subtitle)}
<div class="layout">{nav}<main id="main">{body}</main></div>
{f'<footer class="foot">{footer}</footer>' if footer else ''}
<div class="lb" id="lb" role="dialog" aria-label="Enlarged figure"></div>
<script>{JS}</script></body></html>"""


def document(title, body):
    """A plain document page: README.html, HOW_TO_SUBMIT.html. `body` is Markdown already
    rendered to header / nav / section / aside.callout / table (make_deposit.md_to_html). The
    same CSS as page() plus DOC_CSS; no script, so dark mode follows the system setting."""
    return (f"<!DOCTYPE html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
            f"<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">"
            f"<meta name=\"color-scheme\" content=\"light dark\">"
            f"<title>{esc(title)}</title><style>{CSS}{DOC_CSS}</style></head><body>\n"
            f"<main class=\"doc\" id=\"main\">\n{body}\n</main>\n</body></html>\n")
