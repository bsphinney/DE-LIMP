#!/usr/bin/env python3
"""
make_comparison_report.py -- build the cross-tool comparison report, interactive
HTML + markdown, from the outputs the comparison scripts already produce.

WHY
---
compare_searches.py and compare_analyses.R emit correct numbers in CSVs. Nobody reads
CSVs. This assembles them into the report structure DE-LIMP's Run Comparator prompt
defines (docs/AI_PROMPTS.md section 1) so a cross-tool bake-off lands as a document a
collaborator can act on:

  1 Factual Observations   2 Sources of Disagreement   3 Case for Run A
  4 Case for Run B         5 Settings Audit            6 Concordant Proteins
  7 Synthesis              8 Recommended Follow-ups

It writes the FACTUAL half — every number, table and chart — and leaves the
interpretive sections as explicit TODO blocks for the model to fill from the data.
That split is deliberate: the numbers must not be invented, and the judgement must not
be templated.

The HTML is standalone (no CDN, no fonts, inline SVG), theme-aware, and interactive:
metric toggles, hover tooltips, sortable tables, a 3x3 concordance matrix.

Usage
  python3 make_comparison_report.py --out <dir> \\
      --search-comparison <dir from compare_searches.py> \\
      [--de-comparison <dir from compare_analyses.R>] \\
      [--qc "LabelA:<sessionA>/output/tables/QC_detected_vs_inferred.csv"] \\
      [--qc "LabelB:<sessionB>/output/tables/QC_detected_vs_inferred.csv"] \\
      [--instrument "timsTOF HT"] [--title "..."]

Then: python3 to_docx.py --in <out>/COMPARISON_REPORT.md --out <out>/COMPARISON_REPORT.docx
"""
import argparse, csv, html, json, os, sys

PALETTE = {  # validated with dataviz/scripts/validate_palette.js, light + dark
    "a_light": "#2a78d6", "b_light": "#eb6834",
    "a_dark": "#3987e5", "b_dark": "#d95926",
}


def read_csv(p):
    if not p or not os.path.exists(p):
        return []
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh))


def num(x, d=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--search-comparison", required=True,
                    help="directory written by compare_searches.py")
    ap.add_argument("--de-comparison", help="directory written by compare_analyses.R")
    ap.add_argument("--qc", action="append", default=[],
                    help='"Label:path/to/QC_detected_vs_inferred.csv" (repeatable)')
    ap.add_argument("--instrument", default="")
    ap.add_argument("--title", default="Cross-tool search comparison")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    searches = read_csv(os.path.join(a.search_comparison, "search_summary.csv"))
    if not searches:
        sys.exit(f"no search_summary.csv in {a.search_comparison} — run compare_searches.py first")
    conc = read_csv(os.path.join(a.de_comparison, "concordance_summary.csv")) if a.de_comparison else []
    universe = read_csv(os.path.join(a.de_comparison, "protein_universe.csv")) if a.de_comparison else []

    qc = {}
    for spec in a.qc:
        lab, _, path = spec.partition(":")
        rows = read_csv(path)
        if rows:
            qc[lab] = rows

    # 3x3 matrices, one per shared contrast
    mats = {}
    if a.de_comparison:
        for fn in sorted(os.listdir(a.de_comparison)):
            if fn.startswith("concordance_") and fn.endswith(".csv") and "summary" not in fn:
                rows = read_csv(os.path.join(a.de_comparison, fn))
                if rows:
                    ct = fn.replace("concordance_", "").replace(".csv", "").split("__")[-1]
                    mats[ct] = [[int(num(r.get(k, 0))) for k in ("Up", "Down", "NS")] for r in rows]

    labels = [s["search"] for s in searches]
    A = labels[0]
    B = labels[1] if len(labels) > 1 else None

    data = {
        "title": a.title, "instrument": a.instrument,
        "searches": searches, "concordance": conc, "universe": universe,
        "qc": qc, "matrices": mats, "labels": labels,
    }

    # ---------------------------------------------------------------- markdown
    md = [f"# {a.title}", ""]
    if a.instrument:
        md += [f"**Instrument:** {a.instrument}", ""]
    md += ["> Both pipelines are established and peer-reviewed; neither is inherently superior.",
           "> Evaluate on the data.", ""]

    md += ["## 1. Factual observations", "",
           "| Search | Runs | Protein groups | Per-run min/med/max | In every run | Median CV |",
           "|---|---:|---:|---|---:|---:|"]
    for s in searches:
        md.append(f"| {s['search']} | {s['runs']} | {s['protein_groups']} | "
                  f"{s['per_run_min']} / {s['per_run_median']} / {s['per_run_max']} | "
                  f"{s['complete_in_all_runs']} ({s['completeness_pct']}%) | {s['median_cv_pct']}% |")
    md.append("")
    md.append("Depth alone does not decide which search is better — read protein groups, "
              "how many survive in *every* run, and the CV together.")
    md.append("")

    if conc:
        md += ["### Differential expression — same DE pipeline, search as the only variable", "",
               "| Contrast | " + " | ".join(f"{c} sig" for c in (A, B)) +
               " | Co-significant | Direction concordance | logFC r |",
               "|---|---:|---:|---:|---:|---:|"]
        for c in conc:
            md.append(f"| {c['contrast']} | {c['sig_A']} | {c['sig_B']} | {c['sig_both']} | "
                      f"{100*num(c['direction_concordance']):.1f}% | {num(c['logFC_correlation']):.3f} |")
        md += ["", "<!-- TODO(model): does one engine give more DE proteins? State whether the "
                   "answer is consistent across contrasts or flips, and quantify. Compare the "
                   "logFC correlation against the search-level intensity correlation — ratios "
                   "amplify disagreement, so a high intensity r with a lower logFC r means "
                   "direction is trustworthy and magnitude is engine-dependent. -->", ""]
    for ct, m in mats.items():
        md += [f"#### 3x3 concordance — {ct}", "",
               f"| {A} \\ {B} | Up | Down | NS |", "|---|---:|---:|---:|"]
        for lab, row in zip(("Up", "Down", "NS"), m):
            md.append(f"| **{lab}** | " + " | ".join(f"{v:,}" for v in row) + " |")
        md.append("")

    if qc:
        md += ["### QC panels", "",
               "| Run | " + " | ".join(f"{k} % detected" for k in qc) + " |",
               "|---|" + "---:|" * len(qc)]
        keys = list(qc)
        runs = {r["Sample"]: r for r in qc[keys[0]]}
        for samp in runs:
            cells = []
            for k in keys:
                row = next((r for r in qc[k] if r["Sample"] == samp), None)
                cells.append(f"{num(row['PctDetected']):.1f}" if row else "—")
            md.append(f"| {samp} | " + " | ".join(cells) + " |")
        md += ["", "<!-- TODO(model): note where each engine's detected fraction is higher and "
                   "whether the advantage tracks sample quality. -->", ""]

    for n, head, todo in [
        (2, "Sources of disagreement",
         "For each source: what the difference is, how many proteins it plausibly affects "
         "(systematic vs protein-specific), and the expected direction. Ground claims in the "
         "literature. Name any engine-VERSION difference explicitly — it confounds everything."),
        (3, f"Case for {A}",
         f"Argue the strongest case for trusting {A} for THIS experiment, citing numbers above."),
        (4, f"Case for {B or 'Run B'}",
         f"Argue the strongest case for trusting {B or 'Run B'} for THIS experiment, citing numbers above."),
        (5, "Settings audit",
         "One row per differing setting: could it explain an observed discrepancy? Flag anything "
         "that looks misconfigured for this experiment type. Be specific about which discrepancy."),
        (6, "Concordant proteins",
         "Characterise the shared significant set — pathways, compartments, known markers. Then "
         "warn that protein-GROUP identifiers differ between tools, so some 'unique' proteins are "
         "grouping artefacts, not detection differences."),
        (7, "Synthesis",
         "Weigh sections 3-4. If one pipeline is clearly more appropriate for this design, say why. "
         "If the evidence does not clearly favour one, SAY SO — do not force a recommendation."),
        (8, "Recommended follow-ups",
         "Two concrete follow-ups, one probing each pipeline (not both aimed at the same tool). "
         "Name specific proteins, settings or validations."),
    ]:
        md += [f"## {n}. {head}", "", f"<!-- TODO(model): {todo} -->", ""]

    md_path = os.path.join(a.out, "COMPARISON_REPORT.md")
    open(md_path, "w").write("\n".join(md) + "\n")

    # ---------------------------------------------------------------- html
    html_path = os.path.join(a.out, "COMPARISON_REPORT.html")
    open(html_path, "w").write(build_html(data))

    print(json.dumps({
        "markdown": md_path, "html": html_path,
        "searches": labels, "contrasts": [c["contrast"] for c in conc],
        "qc_panels": list(qc),
        "next": ["Fill every TODO(model) block in the markdown from the data — do not invent numbers.",
                 f"python3 to_docx.py --in {md_path} --out {os.path.join(a.out,'COMPARISON_REPORT.docx')}"],
    }, indent=2))


def build_html(d):
    esc = html.escape
    j = json.dumps
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{esc(d['title'])}</title>
<style>
:root{{color-scheme:light;--s0:#f4f3ef;--s1:#fcfcfb;--s2:#eeece6;--line:#dcd9d0;--lineS:#c4c0b4;
 --t1:#0b0b0b;--t2:#52514e;--t3:#7d7b74;--a:{PALETTE['a_light']};--b:{PALETTE['b_light']};--warn:#eda100;--good:#1baf7a}}
@media (prefers-color-scheme:dark){{:root:where(:not([data-theme="light"])){{color-scheme:dark;
 --s0:#121211;--s1:#1a1a19;--s2:#232322;--line:#33322f;--lineS:#4a4945;
 --t1:#fff;--t2:#c3c2b7;--t3:#928f85;--a:{PALETTE['a_dark']};--b:{PALETTE['b_dark']};--warn:#c98500;--good:#199e70}}}}
:root[data-theme="dark"]{{color-scheme:dark;--s0:#121211;--s1:#1a1a19;--s2:#232322;--line:#33322f;--lineS:#4a4945;
 --t1:#fff;--t2:#c3c2b7;--t3:#928f85;--a:{PALETTE['a_dark']};--b:{PALETTE['b_dark']};--warn:#c98500;--good:#199e70}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--s0);color:var(--t1);line-height:1.6;
 font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}}
.wrap{{max-width:1080px;margin:0 auto;padding:40px 22px 80px}}
header{{border-bottom:2px solid var(--t1);padding-bottom:18px;margin-bottom:24px}}
h1{{font-size:1.8rem;margin:0 0 6px;letter-spacing:-.015em}}h2{{font-size:1.2rem;margin:40px 0 10px}}
.sub{{color:var(--t2);margin:0}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:14px;margin:20px 0}}
.card{{background:var(--s1);border:1px solid var(--line);border-radius:12px;padding:16px 18px}}
.lab{{font-size:.73rem;text-transform:uppercase;letter-spacing:.07em;color:var(--t3);font-weight:700}}
.n{{font-size:1.85rem;font-weight:700;margin:4px 0 2px;font-variant-numeric:tabular-nums}}
.note{{font-size:.83rem;color:var(--t2)}}
.sw{{display:inline-block;width:10px;height:10px;border-radius:3px;margin-right:6px}}
.pa{{background:var(--a)}}.pb{{background:var(--b)}}
.panel{{background:var(--s1);border:1px solid var(--line);border-radius:12px;padding:18px 20px;margin:16px 0}}
.ctl{{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin-bottom:12px}}
button.seg{{background:var(--s2);border:1px solid var(--lineS);color:var(--t2);padding:6px 13px;
 border-radius:999px;font-size:.85rem;cursor:pointer;font-family:inherit}}
button.seg[aria-pressed="true"]{{background:var(--t1);color:var(--s1);border-color:var(--t1)}}
.leg{{display:flex;gap:15px;font-size:.84rem;color:var(--t2);margin-left:auto}}
svg{{display:block;width:100%;height:auto;overflow:visible}}.grid line{{stroke:var(--line)}}
.bar:hover{{opacity:.78;cursor:pointer}}.dot{{stroke:var(--s1);stroke-width:2}}
.tip{{position:fixed;pointer-events:none;opacity:0;transition:opacity .1s;background:var(--t1);
 color:var(--s1);padding:8px 11px;border-radius:8px;font-size:.8rem;z-index:50;max-width:260px}}
table{{border-collapse:collapse;width:100%;font-size:.9rem;margin:.6rem 0}}
th,td{{border-bottom:1px solid var(--line);padding:8px 10px;text-align:left}}
th{{color:var(--t2);font-size:.77rem;text-transform:uppercase;letter-spacing:.05em;cursor:pointer}}
td.num,th.num{{text-align:right;font-variant-numeric:tabular-nums}}
.tw{{overflow-x:auto}}
.todo{{background:var(--s2);border-left:3px solid var(--warn);border-radius:8px;padding:12px 15px;
 margin:14px 0;font-size:.9rem;color:var(--t2)}}
footer{{margin-top:48px;padding-top:16px;border-top:1px solid var(--line);color:var(--t3);font-size:.83rem}}
.tb{{position:fixed;top:14px;right:14px;z-index:60;background:var(--s1);border:1px solid var(--lineS);
 color:var(--t2);border-radius:999px;padding:6px 13px;font-size:.8rem;cursor:pointer}}
</style></head><body>
<button class="tb" id="tb">◐ theme</button>
<div class="wrap">
<header><h1>{esc(d['title'])}</h1>
<p class="sub">{esc(d['instrument'])}</p></header>

<div class="todo"><strong>Both pipelines are established and peer-reviewed; neither is inherently
superior.</strong> Sections 2–8 of the markdown carry <code>TODO(model)</code> blocks to be filled
from the data below. The numbers here are generated, not written.</div>

<div class="cards" id="cards"></div>

<h2>1. Factual observations</h2>
<div class="panel">
  <div class="ctl" id="mctl"></div>
  <svg id="bars" viewBox="0 0 900 400" role="img" aria-label="Per-search identification summary"></svg>
</div>

<div class="panel" id="concpanel" hidden>
  <div class="ctl" id="cctl"><strong style="font-size:.92rem">3&times;3 DE concordance</strong></div>
  <svg id="conc" viewBox="0 0 900 290" role="img" aria-label="Concordance matrix"></svg>
  <p class="note" id="cnote"></p>
</div>

<div class="panel" id="qcpanel" hidden>
  <div class="ctl"><strong style="font-size:.92rem">Detected vs DPC-inferred, per sample</strong>
    <div class="leg" id="qcleg"></div></div>
  <svg id="qc" viewBox="0 0 900 340" role="img" aria-label="Percent detected per sample"></svg>
</div>

<div class="tw"><table id="tbl"><thead></thead><tbody></tbody></table></div>

<footer>Generated by <code>make_comparison_report.py</code> ·
structure follows DE-LIMP's Run Comparator (<code>docs/AI_PROMPTS.md</code> §1) ·
palette validated with the dataviz six-checks validator</footer>
</div>
<div class="tip" id="tip"></div>
<script>
const D = {j(d)};
const NS="http://www.w3.org/2000/svg", tip=document.getElementById("tip");
const fmt=n=>Number(n).toLocaleString();
const el=(t,a={{}})=>{{const e=document.createElementNS(NS,t);for(const k in a)e.setAttribute(k,a[k]);return e}};
function T(h,ev){{tip.innerHTML=h;tip.style.opacity=1;const r=tip.getBoundingClientRect();
 let x=ev.clientX+14,y=ev.clientY-10;
 if(x+r.width>innerWidth-10)x=ev.clientX-r.width-14;
 if(y+r.height>innerHeight-10)y=innerHeight-r.height-10;
 tip.style.left=x+"px";tip.style.top=y+"px"}}
const H=()=>tip.style.opacity=0;
const COL=i=>i===0?"var(--a)":"var(--b)";

/* cards */
(function(){{
  const c=document.getElementById("cards");
  D.searches.forEach((s,i)=>{{
    c.insertAdjacentHTML("beforeend",
     `<div class="card"><div class="lab"><span class="sw ${{i===0?'pa':'pb'}}"></span>${{s.search}}</div>
      <div class="n">${{fmt(s.protein_groups)}}</div>
      <div class="note">protein groups · ${{s.median_cv_pct}}% median CV</div></div>`);
  }});
  D.concordance.forEach(k=>{{
    c.insertAdjacentHTML("beforeend",
     `<div class="card"><div class="lab">${{k.contrast}} · significant</div>
      <div class="n">${{fmt(k.sig_A)}} <span style="font-size:1rem;color:var(--t3)">vs ${{fmt(k.sig_B)}}</span></div>
      <div class="note">${{fmt(k.sig_both)}} co-significant · r ${{Number(k.logFC_correlation).toFixed(3)}}</div></div>`);
  }});
}})();

/* metric bars */
const METRICS=[["protein_groups","Protein groups"],["complete_in_all_runs","In every run"],
               ["median_cv_pct","Median CV %"]];
let mi=0;
(function(){{const w=document.getElementById("mctl");
  METRICS.forEach(([k,l],i)=>{{const b=document.createElement("button");
    b.className="seg";b.textContent=l;b.setAttribute("aria-pressed",i===0);
    b.onclick=()=>{{mi=i;[...w.querySelectorAll("button")].forEach((o,jj)=>o.setAttribute("aria-pressed",jj===i));bars()}};
    w.appendChild(b)}});
  const lg=document.createElement("div");lg.className="leg";
  lg.innerHTML=D.searches.map((s,i)=>`<span><span class="sw ${{i===0?'pa':'pb'}}"></span>${{s.search}}</span>`).join("");
  w.appendChild(lg)}})();
function bars(){{
  const svg=document.getElementById("bars");svg.innerHTML="";
  const [key,label]=METRICS[mi];
  const W=900,L=210,R=70,Tp=20,B=44,ih=D.searches.length*72;
  const vals=D.searches.map(s=>Number(s[key])||0), mx=Math.max(...vals)*1.12;
  svg.setAttribute("viewBox",`0 0 ${{W}} ${{Tp+ih+B}}`);
  const x=v=>L+(v/mx)*(W-L-R);
  for(let i=0;i<=4;i++){{const gx=x(mx*i/4);
    svg.appendChild(el("line",{{x1:gx,x2:gx,y1:Tp,y2:Tp+ih,class:"grid","stroke-width":1}}));
    const t=el("text",{{x:gx,y:Tp+ih+18,"text-anchor":"middle",fill:"var(--t3)","font-size":"11"}});
    t.textContent=(mx*i/4).toFixed(mx<100?1:0);svg.appendChild(t)}}
  D.searches.forEach((s,i)=>{{
    const cy=Tp+72*i+36, v=Number(s[key])||0;
    const lb=el("text",{{x:L-14,y:cy+4,"text-anchor":"end","font-size":"12",fill:"var(--t2)"}});
    lb.textContent=s.search;svg.appendChild(lb);
    const r=el("rect",{{x:L,y:cy-13,width:Math.max(1,x(v)-L),height:26,rx:4,fill:COL(i),class:"bar"}});
    r.addEventListener("mousemove",e=>T(`<b>${{s.search}}</b><br>${{label}}: <b>${{fmt(v)}}</b><br>`+
      `runs ${{s.runs}} · per-run ${{s.per_run_min}}–${{s.per_run_max}}`,e));
    r.addEventListener("mouseleave",H);svg.appendChild(r);
    const vt=el("text",{{x:x(v)+10,y:cy+4,"font-size":"12",fill:"var(--t2)"}});
    vt.textContent=fmt(v);svg.appendChild(vt)}});
  const at=el("text",{{x:(L+W-R)/2,y:Tp+ih+38,"text-anchor":"middle",fill:"var(--t2)","font-size":"11"}});
  at.textContent=label;svg.appendChild(at)}}

/* 3x3 concordance */
const CT=Object.keys(D.matrices||{{}});let ck=CT[0];
if(CT.length){{
  document.getElementById("concpanel").hidden=false;
  const w=document.getElementById("cctl");
  CT.forEach((k,i)=>{{const b=document.createElement("button");b.className="seg";b.textContent=k;
    b.setAttribute("aria-pressed",i===0);
    b.onclick=()=>{{ck=k;[...w.querySelectorAll("button")].forEach(o=>o.setAttribute("aria-pressed",o.textContent===k));conc()}};
    w.appendChild(b)}});
}}
function conc(){{
  if(!CT.length)return;
  const svg=document.getElementById("conc");svg.innerHTML="";
  const m=D.matrices[ck],labs=["Up","Down","NS"],L=200,Tp=44,cw=118,ch=62;
  const mx=Math.max(...m.flat());
  labs.forEach((l,jj)=>{{const t=el("text",{{x:L+cw*jj+cw/2,y:Tp-14,"text-anchor":"middle",
    "font-size":"12",fill:"var(--t2)","font-weight":"600"}});
    t.textContent=(D.labels[1]||"B")+" "+l;svg.appendChild(t)}});
  m.forEach((row,i)=>{{
    const t=el("text",{{x:L-14,y:Tp+ch*i+ch/2+4,"text-anchor":"end","font-size":"12",
      fill:"var(--t2)","font-weight":"600"}});
    t.textContent=(D.labels[0]||"A")+" "+labs[i];svg.appendChild(t);
    row.forEach((v,jj)=>{{
      const agree=i===jj, flip=(i===0&&jj===1)||(i===1&&jj===0);
      const f=Math.pow(v/mx,0.45);
      const fill = flip ? `color-mix(in oklab, var(--warn) ${{Math.max(18,f*100)}}%, var(--s1))`
        : agree ? `color-mix(in oklab, var(--a) ${{Math.max(10,f*88)}}%, var(--s1))`
        : `color-mix(in oklab, var(--t3) ${{Math.max(7,f*42)}}%, var(--s1))`;
      const r=el("rect",{{x:L+cw*jj+1,y:Tp+ch*i+1,width:cw-3,height:ch-3,rx:6,fill:fill,
        stroke:"var(--line)"}});
      r.style.cursor="pointer";
      const mean=agree?"both agree":(flip?"OPPOSITE DIRECTION":(i===2?"B only":"A only"));
      r.addEventListener("mousemove",e=>T(`<b>${{fmt(v)}}</b> protein groups<br>${{mean}}`,e));
      r.addEventListener("mouseleave",H);svg.appendChild(r);
      const tx=el("text",{{x:L+cw*jj+cw/2,y:Tp+ch*i+ch/2+5,"text-anchor":"middle","font-size":"14",
        "font-weight":(agree||flip)?"700":"500",fill:f>0.62?"var(--s1)":"var(--t1)"}});
      tx.textContent=fmt(v);svg.appendChild(tx)}})}});
  const k=(D.concordance||[]).find(c=>ck.includes(c.contrast)||c.contrast.includes(ck));
  document.getElementById("cnote").textContent = k
    ? `${{k.contrast}} — direction concordance ${{(100*Number(k.direction_concordance)).toFixed(1)}}%, logFC r ${{Number(k.logFC_correlation).toFixed(3)}}. Off-diagonal NS cells are detection differences; the Up/Down corners are genuine direction disagreements.`
    : "";
}}

/* QC dumbbell */
const QK=Object.keys(D.qc||{{}});
function qcplot(){{
  if(QK.length<1)return;
  document.getElementById("qcpanel").hidden=false;
  document.getElementById("qcleg").innerHTML=
    QK.map((k,i)=>`<span><span class="sw ${{i===0?'pa':'pb'}}"></span>${{k}}</span>`).join("");
  const rows=D.qc[QK[0]].map(r=>r.Sample);
  const svg=document.getElementById("qc");svg.innerHTML="";
  const W=900,L=210,R=60,Tp=16,B=44,ih=Math.max(120,rows.length*32);
  svg.setAttribute("viewBox",`0 0 ${{W}} ${{Tp+ih+B}}`);
  const lo=Math.min(50,...QK.flatMap(k=>D.qc[k].map(r=>+r.PctDetected)))-4, hi=100;
  const x=v=>L+((v-lo)/(hi-lo))*(W-L-R), band=ih/rows.length;
  for(let v=Math.ceil(lo/10)*10;v<=100;v+=10){{const gx=x(v);
    svg.appendChild(el("line",{{x1:gx,x2:gx,y1:Tp,y2:Tp+ih,class:"grid","stroke-width":1}}));
    const t=el("text",{{x:gx,y:Tp+ih+18,"text-anchor":"middle",fill:"var(--t3)","font-size":"11"}});
    t.textContent=v+"%";svg.appendChild(t)}}
  rows.forEach((samp,i)=>{{
    const cy=Tp+band*i+band/2;
    const short=samp.replace(/^.*?DIA_/,"");
    const lb=el("text",{{x:L-12,y:cy+4,"text-anchor":"end","font-size":"11",fill:"var(--t2)"}});
    lb.textContent=short.length>26?short.slice(0,25)+"…":short;svg.appendChild(lb);
    const vs=QK.map(k=>{{const r=D.qc[k].find(r=>r.Sample===samp);return r?+r.PctDetected:null}});
    const ok=vs.filter(v=>v!==null);
    if(ok.length>1) svg.appendChild(el("line",{{x1:x(Math.min(...ok)),x2:x(Math.max(...ok)),
      y1:cy,y2:cy,stroke:"var(--lineS)","stroke-width":2}}));
    vs.forEach((v,k)=>{{if(v===null)return;
      const c=el("circle",{{cx:x(v),cy:cy,r:6,fill:COL(k),class:"dot"}});
      c.addEventListener("mousemove",e=>T(`<b>${{short}}</b><br>${{QK[k]}}: <b>${{v.toFixed(1)}}%</b> detected`+
        `<br>${{(100-v).toFixed(1)}}% inferred by DPC`,e));
      c.addEventListener("mouseleave",H);svg.appendChild(c)}})}});
  const at=el("text",{{x:(L+W-R)/2,y:Tp+ih+38,"text-anchor":"middle",fill:"var(--t2)","font-size":"11"}});
  at.textContent="% of protein matrix actually measured (not DPC-inferred)";svg.appendChild(at)}}

/* sortable table */
let sk="search",sa=true;
function tbl(){{
  const cols=Object.keys(D.searches[0]);
  const th=document.querySelector("#tbl thead");
  th.innerHTML="<tr>"+cols.map(c=>`<th class="${{isNaN(D.searches[0][c])?'':'num'}}" data-k="${{c}}">${{c.replace(/_/g," ")}}</th>`).join("")+"</tr>";
  th.querySelectorAll("th").forEach(h=>h.onclick=()=>{{const k=h.dataset.k;
    if(sk===k)sa=!sa;else{{sk=k;sa=true}}tbl()}});
  const rows=[...D.searches].sort((p,q)=>{{const A=p[sk],Bv=q[sk];
    const v=(!isNaN(A)&&!isNaN(Bv))?A-Bv:String(A).localeCompare(String(Bv));return sa?v:-v}});
  document.querySelector("#tbl tbody").innerHTML=rows.map(r=>"<tr>"+
    cols.map(c=>`<td class="${{isNaN(r[c])?'':'num'}}">${{r[c]}}</td>`).join("")+"</tr>").join("");
}}

document.getElementById("tb").onclick=()=>{{
  const cur=document.documentElement.getAttribute("data-theme");
  const nx=cur==="dark"?"light":(cur==="light"?"dark":
    (matchMedia("(prefers-color-scheme: dark)").matches?"light":"dark"));
  document.documentElement.setAttribute("data-theme",nx);
  bars();conc();qcplot()}};

bars();conc();qcplot();tbl();
</script></body></html>
"""


if __name__ == "__main__":
    main()
