"""
CrackPipe — Proteomics Data Quality Explorer
A humorous educational tool inspired by FragPipe.
All data is synthetically generated. Not affiliated with FragPipe or Nesvilab.
"""

import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import base64

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="CrackPipe — Proteomics Data Quality Explorer",
    page_icon="🧪",
    layout="wide",
)

# ---------------------------------------------------------------------------
# Inline SVG logo — cracked Erlenmeyer flask
# ---------------------------------------------------------------------------
LOGO_SVG = """
<svg viewBox="0 0 120 120" xmlns="http://www.w3.org/2000/svg" width="{size}" height="{size}">
  <defs>
    <linearGradient id="flask_grad_{id}" x1="0" y1="0" x2="1" y2="1">
      <stop offset="0%" stop-color="#38bdf8"/>
      <stop offset="100%" stop-color="#0284c7"/>
    </linearGradient>
    <linearGradient id="liquid_grad_{id}" x1="0" y1="0" x2="0" y2="1">
      <stop offset="0%" stop-color="#f59e0b"/>
      <stop offset="100%" stop-color="#ef4444"/>
    </linearGradient>
    <linearGradient id="crack_grad_{id}" x1="0" y1="0" x2="1" y2="1">
      <stop offset="0%" stop-color="#fbbf24"/>
      <stop offset="100%" stop-color="#f59e0b"/>
    </linearGradient>
  </defs>
  <!-- Background circle -->
  <circle cx="60" cy="60" r="58" fill="#1a2332"/>
  <!-- Spectral lines -->
  <line x1="15" y1="50" x2="25" y2="50" stroke="#38bdf8" stroke-width="1" opacity="0.3"/>
  <line x1="15" y1="55" x2="22" y2="55" stroke="#38bdf8" stroke-width="1" opacity="0.2"/>
  <line x1="15" y1="60" x2="28" y2="60" stroke="#38bdf8" stroke-width="1" opacity="0.25"/>
  <line x1="95" y1="50" x2="105" y2="50" stroke="#38bdf8" stroke-width="1" opacity="0.3"/>
  <line x1="98" y1="55" x2="105" y2="55" stroke="#38bdf8" stroke-width="1" opacity="0.2"/>
  <line x1="92" y1="60" x2="105" y2="60" stroke="#38bdf8" stroke-width="1" opacity="0.25"/>
  <!-- Flask body -->
  <path d="M50 28 L50 55 L32 90 Q30 95 35 97 L85 97 Q90 95 88 90 L70 55 L70 28"
        fill="none" stroke="url(#flask_grad_{id})" stroke-width="2.5" stroke-linejoin="round"/>
  <!-- Flask neck -->
  <rect x="48" y="22" width="24" height="8" rx="2" fill="none"
        stroke="url(#flask_grad_{id})" stroke-width="2"/>
  <!-- Liquid -->
  <path d="M38 82 Q40 78 45 80 Q50 82 55 79 Q60 76 65 80 Q70 78 75 80 L82 82 L86 90
        Q88 95 83 97 L37 97 Q32 95 34 90 Z"
        fill="url(#liquid_grad_{id})" opacity="0.85"/>
  <!-- Bubbles -->
  <circle cx="50" cy="88" r="2" fill="white" opacity="0.4"/>
  <circle cx="62" cy="85" r="1.5" fill="white" opacity="0.3"/>
  <circle cx="55" cy="92" r="1" fill="white" opacity="0.35"/>
  <circle cx="68" cy="90" r="1.8" fill="white" opacity="0.25"/>
  <!-- Crack lines -->
  <path d="M52 45 L48 55 L53 62 L47 72" stroke="url(#crack_grad_{id})"
        stroke-width="2" fill="none" stroke-linecap="round"/>
  <path d="M68 50 L72 58 L66 65 L73 75" stroke="url(#crack_grad_{id})"
        stroke-width="1.8" fill="none" stroke-linecap="round"/>
  <path d="M55 52 L58 48" stroke="#fbbf24" stroke-width="1.2" opacity="0.7"/>
  <path d="M65 55 L62 52" stroke="#fbbf24" stroke-width="1.2" opacity="0.7"/>
  <!-- Lightning bolt -->
  <polygon points="58,15 54,26 59,24 55,35 62,22 57,24"
           fill="#fbbf24" opacity="0.9"/>
  <!-- Sparkles -->
  <circle cx="44" cy="48" r="1.2" fill="#fbbf24" opacity="0.6"/>
  <circle cx="76" cy="52" r="1" fill="#fbbf24" opacity="0.5"/>
  <circle cx="42" cy="65" r="0.8" fill="#fbbf24" opacity="0.4"/>
  <circle cx="78" cy="68" r="1.1" fill="#fbbf24" opacity="0.5"/>
</svg>
"""


def logo_html(size: int = 80, uid: str = "hdr") -> str:
    svg = LOGO_SVG.replace("{size}", str(size)).replace("{id}", uid)
    b64 = base64.b64encode(svg.encode("utf-8")).decode("utf-8")
    return f'<img src="data:image/svg+xml;base64,{b64}" width="{size}" height="{size}" alt="CrackPipe logo"/>'


# ---------------------------------------------------------------------------
# CSS
# ---------------------------------------------------------------------------
CSS = """
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@500;700&family=IBM+Plex+Sans:wght@400;500;600&display=swap');

/* Global */
html, body, [class*="css"] { font-family: 'IBM Plex Sans', sans-serif; color: #1e293b; }
h1,h2,h3,h4 { font-family: 'IBM Plex Mono', monospace; }

/* Sidebar */
section[data-testid="stSidebar"] {
    background: #1e293b !important;
}
section[data-testid="stSidebar"] * {
    color: #e2e8f0 !important;
}
section[data-testid="stSidebar"] label[data-testid="stWidgetLabel"] {
    color: #94a3b8 !important;
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.75rem !important;
    text-transform: uppercase;
    letter-spacing: 0.05em;
}
section[data-testid="stSidebar"] .stRadio > div {
    gap: 0.15rem;
}
section[data-testid="stSidebar"] .stRadio > div label {
    padding: 0.35rem 0.5rem;
    border-radius: 6px;
    transition: background 0.15s;
    font-size: 0.92rem !important;
}
section[data-testid="stSidebar"] .stRadio > div label:hover {
    background: rgba(56, 189, 248, 0.12);
}

/* Header bar */
.fp-header {
    background: linear-gradient(135deg, #0f172a, #1e293b);
    border-radius: 10px;
    padding: 1rem 1.5rem;
    display: flex;
    align-items: center;
    gap: 1rem;
    margin-bottom: 1.2rem;
}
.fp-header .title {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 1.6rem;
    font-weight: 700;
    color: #38bdf8;
    margin: 0;
}
.fp-header .subtitle {
    color: #94a3b8;
    font-size: 0.85rem;
    margin: 0;
}
.fp-header .version {
    background: #0284c7;
    color: white;
    font-size: 0.7rem;
    padding: 2px 8px;
    border-radius: 999px;
    font-family: 'IBM Plex Mono', monospace;
    margin-left: 0.5rem;
}

/* Cards */
.fp-card { background:#fff; border-radius:8px; padding:1rem 1.2rem; margin:0.7rem 0;
           border:1px solid #e2e8f0; border-left:4px solid #e2e8f0; }
.fp-card-danger { border-left-color: #ef4444; }
.fp-card-ok { border-left-color: #22c55e; }
.fp-card-warn { border-left-color: #f59e0b; }
.fp-card-info { border-left-color: #0284c7; }

/* Severity badges */
.sev-badge {
    display:inline-block; font-size:0.72rem; font-weight:600; padding:2px 10px;
    border-radius:999px; font-family:'IBM Plex Mono',monospace; letter-spacing:0.04em;
}
.sev-high { background:#fef2f2; color:#dc2626; border:1px solid #fca5a5; }
.sev-med  { background:#fffbeb; color:#d97706; border:1px solid #fcd34d; }
.sev-low  { background:#eff6ff; color:#2563eb; border:1px solid #93c5fd; }

/* Status bar */
.fp-status {
    font-family:'IBM Plex Mono',monospace; font-size:0.75rem; color:#64748b;
    background:#f8fafc; border:1px solid #e2e8f0; border-radius:6px;
    padding:0.4rem 1rem; margin-top:1.5rem; text-align:center;
}

/* Module category cards on home */
.cat-card {
    background: #f8fafc; border:1px solid #e2e8f0; border-radius:8px;
    padding: 0.9rem 1rem; margin-bottom:0.5rem;
}
.cat-card h4 { margin:0 0 0.4rem 0; font-size:0.95rem; color:#0284c7; }
.cat-card p  { margin:0; font-size:0.85rem; color:#64748b; }

/* Plot background helper */
.plot-col-title {
    text-align:center; font-family:'IBM Plex Mono',monospace; font-size:0.9rem;
    font-weight:600; margin-bottom:0.3rem;
}
"""

st.markdown(f"<style>{CSS}</style>", unsafe_allow_html=True)


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def header_bar():
    st.markdown(
        f"""<div class="fp-header">
            {logo_html(64, "hdr")}
            <div>
                <p class="title">CrackPipe <span class="version">v1.0.0</span></p>
                <p class="subtitle">Proteomics Data Quality Explorer</p>
            </div>
        </div>""",
        unsafe_allow_html=True,
    )


def severity_badge(level: str) -> str:
    cls = {"HIGH": "sev-high", "MEDIUM": "sev-med", "LOW": "sev-low"}[level]
    return f'<span class="sev-badge {cls}">SEVERITY: {level}</span>'


def card(html: str, kind: str = "info"):
    st.markdown(f'<div class="fp-card fp-card-{kind}">{html}</div>', unsafe_allow_html=True)


PLOTLY_LAYOUT = dict(
    template="plotly_white",
    paper_bgcolor="#fafbfc",
    plot_bgcolor="#fafbfc",
    font=dict(family="IBM Plex Sans, sans-serif", color="#1e293b"),
    margin=dict(l=50, r=30, t=40, b=40),
)


def styled_fig(fig, **kw):
    fig.update_layout(**{**PLOTLY_LAYOUT, **kw})
    return fig


# ---------------------------------------------------------------------------
# NAV ITEMS
# ---------------------------------------------------------------------------
NAV_ITEMS = [
    "🏠  Home",
    "👻  Chimeric Spectra",
    "⚡  Spray Instability",
    "📏  Mass Cal Drift",
    "🎲  Batch Effects",
    "🕳️  Missing Values",
    "📉  TMT Compression",
    "🧤  Contamination",
    "✂️  Incomplete Digestion",
    "🔄  LC Carryover",
    "📡  DIA Pitfalls",
    "🔍  Missing Proteins",
    "🧲  AP-MS Normalization",
    "🩻  Western vs MS",
    "🎰  P-Value Shopping",
    "🍒  Cherry-Picking Data",
    "📐  Multiple Testing",
    "🧠  Quiz",
]


# ===================================================================
# MODULE: Home
# ===================================================================
def page_home():
    header_bar()

    card(
        "<h3 style='margin-top:0'>Welcome to CrackPipe!</h3>"
        "<p>An educational tool that teaches proteomics researchers how to spot "
        "common data quality problems in mass spectrometry-based proteomics.</p>"
        "<p style='color:#64748b;font-size:0.85rem;'>Inspired by "
        "<a href='https://github.com/Nesvilab/FragPipe' target='_blank'>FragPipe</a> — "
        "but instead of processing data, we show what happens when things go "
        "<b>terribly wrong</b>.</p>",
        "info",
    )

    cols = st.columns(2)
    categories = [
        ("🔬 Spectra & MS", "Chimeric Spectra · Spray Instability · Mass Cal Drift"),
        ("📊 Quantification", "Batch Effects · Missing Values · TMT Compression · AP-MS Norm"),
        ("🧪 Sample Prep & LC", "Contamination · Incomplete Digestion · LC Carryover"),
        ("🔎 Search & Validation", "DIA Pitfalls · Missing Proteins · Western vs MS"),
        ("📈 Statistics & Reproducibility", "P-Value Shopping · Cherry-Picking · Multiple Testing"),
        ("🧠 Quiz", "15 questions to test your diagnostic skills"),
    ]
    for i, (title, desc) in enumerate(categories):
        with cols[i % 2]:
            st.markdown(
                f'<div class="cat-card"><h4>{title}</h4><p>{desc}</p></div>',
                unsafe_allow_html=True,
            )

    st.markdown(
        '<div class="fp-status">STATUS: Ready &nbsp;|&nbsp; Modules: 16 &nbsp;|&nbsp; '
        "Quiz questions: 15</div>",
        unsafe_allow_html=True,
    )


# ===================================================================
# MODULE: Chimeric Spectra
# ===================================================================
def page_chimeric():
    header_bar()
    st.markdown("### 👻 Chimeric Spectra")
    st.markdown(severity_badge("HIGH"), unsafe_allow_html=True)

    card(
        "<b>The problem:</b> When two peptides co-elute and get co-isolated in the same "
        "MS2 scan, you get a <i>chimeric</i> (mixed) spectrum. The search engine can't "
        "explain all the peaks with a single peptide, leading to poor scores and false IDs.",
        "danger",
    )

    rng = np.random.RandomState(42)

    # --- Good spectrum (single peptide) ---
    residues = rng.uniform(57, 186, 12)
    residues = residues / residues.sum() * 1400
    b_ions = np.cumsum(residues) + 1.0078
    y_ions = 1400 - b_ions + 18.015 + 1.0078
    good_mz = np.concatenate([b_ions, y_ions])
    good_int = rng.uniform(30, 100, len(good_mz))
    good_labels = [f"b{i+1}" for i in range(len(b_ions))] + [f"y{i+1}" for i in range(len(y_ions))]
    good_colors = ["#0284c7"] * len(good_mz)
    # noise
    noise_mz = rng.uniform(100, 1400, 15)
    noise_int = rng.uniform(2, 12, 15)
    good_mz = np.concatenate([good_mz, noise_mz])
    good_int = np.concatenate([good_int, noise_int])
    good_labels += [""] * len(noise_mz)
    good_colors += ["#94a3b8"] * len(noise_mz)

    # --- Bad spectrum (chimeric — two peptides merged) ---
    residues2 = rng.uniform(57, 186, 10)
    residues2 = residues2 / residues2.sum() * 1650
    b2 = np.cumsum(residues2) + 1.0078
    y2 = 1650 - b2 + 18.015 + 1.0078
    bad_mz = np.concatenate([good_mz, b2, y2])
    bad_int = np.concatenate(
        [good_int * rng.uniform(0.3, 0.7, len(good_int)),
         rng.uniform(20, 80, len(b2)),
         rng.uniform(20, 80, len(y2))]
    )

    # Helper: build stem plot traces (vertical lines from baseline to peak)
    def add_stems(fig, mzs, intensities, color, name=None):
        stem_x, stem_y = [], []
        for m, i in zip(mzs, intensities):
            stem_x.extend([m, m, None])
            stem_y.extend([0, i, None])
        fig.add_trace(go.Scatter(
            x=stem_x, y=stem_y, mode="lines",
            line=dict(color=color, width=1.5),
            showlegend=False, name=name or "",
            hoverinfo="skip",
        ))

    col1, col2 = st.columns(2)
    with col1:
        st.markdown('<p class="plot-col-title">Good Data ✅</p>', unsafe_allow_html=True)
        fig = go.Figure()
        # Matched ions (b/y) as blue stems
        matched = [(m, i) for m, i, lab in zip(good_mz, good_int, good_labels) if lab]
        if matched:
            add_stems(fig, [p[0] for p in matched], [p[1] for p in matched], "#0284c7")
        # Noise as gray stems
        noise = [(m, i) for m, i, lab in zip(good_mz, good_int, good_labels) if not lab]
        if noise:
            add_stems(fig, [p[0] for p in noise], [p[1] for p in noise], "#94a3b8")
        # Ion labels
        for m, i, lab in zip(good_mz, good_int, good_labels):
            if lab:
                fig.add_annotation(x=m, y=i, text=lab, showarrow=False,
                                   yshift=8, font=dict(size=8, color="#0284c7"))
        styled_fig(fig, title="Clean MS2 — single peptide",
                   xaxis_title="m/z", yaxis_title="Intensity", height=420)
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.markdown('<p class="plot-col-title">Bad Data ❌</p>', unsafe_allow_html=True)
        fig2 = go.Figure()
        # Peptide 1 ions (attenuated original) in orange
        n_pep1 = len(good_mz)
        add_stems(fig2, bad_mz[:n_pep1], bad_int[:n_pep1], "#ef4444")
        # Peptide 2 ions (interfering) in darker red
        add_stems(fig2, bad_mz[n_pep1:], bad_int[n_pep1:], "#b91c1c")
        styled_fig(fig2, title="Chimeric MS2 — two peptides mixed",
                   xaxis_title="m/z", yaxis_title="Intensity", height=420)
        st.plotly_chart(fig2, use_container_width=True)

    with st.expander("📖 How to spot & fix chimeric spectra"):
        st.markdown("""
**Detection:**
- Poor fragment ion coverage (< 50% of expected ions matched)
- Many intense unexplained peaks remaining after the best PSM match
- Two apparent b-ion or y-ion series at different mass offsets
- Check MS1 scan for multiple precursor ions within the isolation window

**Fixes:**
- Use **narrower isolation windows** (1.2 Da instead of 2.0 Da)
- Switch to **DIA with deconvolution** software (DIA-NN, Spectronaut)
- Filter PSMs by **explained intensity** (% of total intensity assigned to the match)
- Use **chimera-aware search engines** (MSFragger, Comet with multi-match)
- Check **precursor purity** metrics in your search results
""")


# ===================================================================
# MODULE: Spray Instability
# ===================================================================
def page_spray():
    header_bar()
    st.markdown("### ⚡ Spray Instability")
    st.markdown(severity_badge("HIGH"), unsafe_allow_html=True)

    card(
        "<b>The problem:</b> Electrospray instability causes Total Ion Current (TIC) "
        "dropouts, meaning you lose peptides in large chunks of your gradient. Entire "
        "regions of your chromatogram may have zero signal.",
        "danger",
    )

    rng = np.random.RandomState(7)
    t = np.linspace(0, 120, 2400)
    base = 1e8 * (0.3 + 0.7 * np.exp(-((t - 50) ** 2) / 800))

    # Good: stable
    stable = base + rng.normal(0, base * 0.03)
    stable = np.clip(stable, 0, None)

    # Bad: dropouts + oscillation + crash
    unstable = base.copy()
    dropout_mask = rng.random(len(t))
    unstable[dropout_mask < 0.03] *= rng.uniform(0, 0.2, (dropout_mask < 0.03).sum())
    unstable += base * 0.08 * np.sin(t * 2)
    crash_start = rng.randint(800, 1200)
    unstable[crash_start:crash_start + 150] *= rng.uniform(0, 0.15, 150)
    unstable += rng.normal(0, base * 0.05)
    unstable = np.clip(unstable, 0, None)

    col1, col2 = st.columns(2)
    with col1:
        st.markdown('<p class="plot-col-title">Good Data ✅</p>', unsafe_allow_html=True)
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=t, y=stable, fill="tozeroy",
                                 line=dict(color="#0284c7", width=1),
                                 fillcolor="rgba(2,132,199,0.25)", showlegend=False))
        styled_fig(fig, title="Stable TIC", xaxis_title="Retention Time (min)",
                   yaxis_title="TIC", height=400)
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.markdown('<p class="plot-col-title">Bad Data ❌</p>', unsafe_allow_html=True)
        fig2 = go.Figure()
        fig2.add_trace(go.Scatter(x=t, y=unstable, fill="tozeroy",
                                  line=dict(color="#ef4444", width=1),
                                  fillcolor="rgba(239,68,68,0.25)", showlegend=False))
        styled_fig(fig2, title="Unstable TIC — dropouts & crash",
                   xaxis_title="Retention Time (min)", yaxis_title="TIC", height=400)
        st.plotly_chart(fig2, use_container_width=True)

    with st.expander("📖 How to spot & fix spray instability"):
        st.markdown("""
**Causes:**
- Clogged or damaged emitter tip
- Air bubbles in the LC lines
- Wrong spray voltage (too high or too low)
- Dirty MS inlet or capillary
- High salt / buffer concentration in sample

**Fixes:**
- Inspect and replace the emitter tip
- Adjust spray voltage ± 200 V and observe TIC stability
- Run a blank gradient to isolate the cause (column vs spray)
- Desalt samples properly (C18 cleanup, SPE)
- Purge LC lines and degas solvents
- Check grounding and electrical connections
""")


# ===================================================================
# MODULE: Mass Calibration Drift
# ===================================================================
def page_masscal():
    header_bar()
    st.markdown("### 📏 Mass Calibration Drift")
    st.markdown(severity_badge("MEDIUM"), unsafe_allow_html=True)

    card(
        "<b>The problem:</b> Systematic mass error that changes over the LC-MS run. "
        "Early peptides may be well-calibrated, but by the end of the gradient the "
        "mass error has drifted several ppm, causing missed identifications.",
        "warn",
    )

    rng = np.random.RandomState(123)
    n = 500
    rt = rng.uniform(2, 110, n)

    # Good
    good_err = rng.normal(0, 1.5, n)

    # Bad: drift + sinusoid
    bad_err = 3 * np.sin(rt / 30) + np.linspace(-1, 4, n)[np.argsort(rt).argsort()] + rng.normal(0, 1.0, n)

    col1, col2 = st.columns(2)
    with col1:
        st.markdown('<p class="plot-col-title">Good Data ✅</p>', unsafe_allow_html=True)
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=rt, y=good_err, mode="markers",
                                 marker=dict(size=4, color="#0284c7", opacity=0.6),
                                 showlegend=False))
        fig.add_hline(y=5, line_dash="dash", line_color="#94a3b8", annotation_text="+5 ppm")
        fig.add_hline(y=-5, line_dash="dash", line_color="#94a3b8", annotation_text="-5 ppm")
        fig.add_hline(y=0, line_color="#64748b", line_width=1)
        styled_fig(fig, title="Centered mass error", xaxis_title="RT (min)",
                   yaxis_title="Mass error (ppm)", height=420)
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.markdown('<p class="plot-col-title">Bad Data ❌</p>', unsafe_allow_html=True)
        fig2 = go.Figure()
        fig2.add_trace(go.Scatter(x=rt, y=bad_err, mode="markers",
                                  marker=dict(size=4, color="#ef4444", opacity=0.6),
                                  showlegend=False))
        # Trend line
        rt_sorted = np.sort(rt)
        trend = 3 * np.sin(rt_sorted / 30) + np.linspace(-1, 4, n)
        fig2.add_trace(go.Scatter(x=rt_sorted, y=trend, mode="lines",
                                  line=dict(color="#f59e0b", width=3),
                                  name="Drift trend"))
        fig2.add_hline(y=5, line_dash="dash", line_color="#94a3b8", annotation_text="+5 ppm")
        fig2.add_hline(y=-5, line_dash="dash", line_color="#94a3b8", annotation_text="-5 ppm")
        fig2.add_hline(y=0, line_color="#64748b", line_width=1)
        styled_fig(fig2, title="Drifting mass error", xaxis_title="RT (min)",
                   yaxis_title="Mass error (ppm)", height=420)
        st.plotly_chart(fig2, use_container_width=True)

    with st.expander("📖 How to spot & fix mass calibration drift"):
        st.markdown("""
**What to look for:**
- A centered, symmetric cloud around 0 ppm = **good**
- A visible trend or drift = **bad**
- Errors > ±5 ppm on an Orbitrap are concerning
- Sinusoidal patterns suggest temperature fluctuations

**Fixes:**
- **Recalibrate** the instrument (external calibration)
- Use **lock-mass** (internal calibrant ions) during acquisition
- Apply **software recalibration** post-acquisition (MSFragger mass-offset search, MaxQuant first-search recalibration)
- Check instrument temperature stability and vacuum
""")


# ===================================================================
# MODULE: Batch Effects
# ===================================================================
def page_batch():
    header_bar()
    st.markdown("### 🎲 Batch Effects")
    st.markdown(severity_badge("HIGH"), unsafe_allow_html=True)

    card(
        "<b>The problem:</b> When PCA separates your samples by processing batch "
        "instead of biological condition, your biology is buried under technical "
        "variation. This is one of the most common and destructive problems in proteomics.",
        "danger",
    )

    rng = np.random.RandomState(99)
    # 3 conditions × 2 batches × 4 reps = 24 points
    conditions = ["Control", "Treatment A", "Treatment B"] * 8
    batches = (["Batch 1"] * 4 + ["Batch 2"] * 4) * 3
    cond_arr = np.array(conditions)
    batch_arr = np.array(batches)

    # Good PCA: PC1 = condition, PC2 = noise
    cond_means = {"Control": -3, "Treatment A": 1, "Treatment B": 3}
    good_pc1 = np.array([cond_means[c] for c in conditions]) + rng.normal(0, 0.6, 24)
    good_pc2 = rng.normal(0, 1.5, 24)

    # Bad PCA: PC1 = batch, PC2 = condition (buried)
    batch_means = {"Batch 1": -3, "Batch 2": 3}
    bad_pc1 = np.array([batch_means[b] for b in batches]) + rng.normal(0, 0.7, 24)
    bad_pc2 = np.array([cond_means[c] for c in conditions]) * 0.5 + rng.normal(0, 0.8, 24)

    symbols = {"Batch 1": "circle", "Batch 2": "diamond"}
    colors = {"Control": "#0284c7", "Treatment A": "#22c55e", "Treatment B": "#f59e0b"}

    col1, col2 = st.columns(2)
    with col1:
        st.markdown('<p class="plot-col-title">Good Data ✅</p>', unsafe_allow_html=True)
        fig = go.Figure()
        for cond in ["Control", "Treatment A", "Treatment B"]:
            for batch in ["Batch 1", "Batch 2"]:
                mask = (cond_arr == cond) & (batch_arr == batch)
                fig.add_trace(go.Scatter(
                    x=good_pc1[mask], y=good_pc2[mask], mode="markers",
                    marker=dict(size=10, color=colors[cond],
                                symbol=symbols[batch], line=dict(width=1, color="white")),
                    name=f"{cond} / {batch}",
                ))
        styled_fig(fig, title="PCA — biology drives PC1",
                   xaxis_title="PC1 (45% var)", yaxis_title="PC2 (12% var)",
                   height=450, legend=dict(font=dict(size=10)))
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.markdown('<p class="plot-col-title">Bad Data ❌</p>', unsafe_allow_html=True)
        fig2 = go.Figure()
        for cond in ["Control", "Treatment A", "Treatment B"]:
            for batch in ["Batch 1", "Batch 2"]:
                mask = (cond_arr == cond) & (batch_arr == batch)
                fig2.add_trace(go.Scatter(
                    x=bad_pc1[mask], y=bad_pc2[mask], mode="markers",
                    marker=dict(size=10, color=colors[cond],
                                symbol=symbols[batch], line=dict(width=1, color="white")),
                    name=f"{cond} / {batch}",
                ))
        styled_fig(fig2, title="PCA — batch drives PC1!",
                   xaxis_title="PC1 (52% var — BATCH)", yaxis_title="PC2 (11% var)",
                   height=450, legend=dict(font=dict(size=10)))
        st.plotly_chart(fig2, use_container_width=True)

    with st.expander("📖 How to spot & fix batch effects"):
        st.markdown("""
**Prevention (before the experiment):**
- **Randomize** sample processing across batches
- Include **all conditions in every batch**
- Use the **same reagent lots** across batches
- Include **QC samples** (pooled reference) in each batch

**Correction (after the fact):**
- **ComBat** (sva package) for batch correction
- **limma::removeBatchEffect** for visualization (but include batch as a covariate in the model for DE testing, don't remove it from the data)
- Include **batch as a covariate** in your linear model: `~ condition + batch`
- Never correct batch effects by simply subtracting means — use proper statistical methods
""")


# ===================================================================
# MODULE: Missing Values
# ===================================================================
def page_missing():
    header_bar()
    st.markdown("### 🕳️ Missing Values")
    st.markdown(severity_badge("MEDIUM"), unsafe_allow_html=True)

    card(
        "<b>The problem:</b> Missing values in proteomics are <i>not random</i>. "
        "Low-abundance proteins are more likely to be missing (MNAR = Missing Not At Random). "
        "Naive imputation (e.g., replacing with the mean) can bias your results.",
        "warn",
    )

    rng = np.random.RandomState(55)
    n_prot, n_samp = 80, 10
    base_intensity = np.sort(rng.normal(22, 3, n_prot))  # sorted: low at top
    data = np.tile(base_intensity, (n_samp, 1)).T + rng.normal(0, 0.5, (n_prot, n_samp))

    # MCAR: 15% random
    mcar_mask = rng.random((n_prot, n_samp)) < 0.15
    mcar_display = np.where(mcar_mask, np.nan, 1)

    # MNAR: low-abundance proteins much more likely missing
    mnar_mask = np.zeros((n_prot, n_samp), dtype=bool)
    for i in range(n_prot):
        if base_intensity[i] < np.percentile(base_intensity, 30):
            mnar_mask[i] = rng.random(n_samp) < 0.70
        else:
            mnar_mask[i] = rng.random(n_samp) < 0.05
    mnar_display = np.where(mnar_mask, np.nan, 1)

    col1, col2 = st.columns(2)
    with col1:
        st.markdown('<p class="plot-col-title">MCAR (Random) ✅</p>', unsafe_allow_html=True)
        fig = go.Figure(data=go.Heatmap(
            z=mcar_display, colorscale=[[0, "#e2e8f0"], [1, "#0284c7"]],
            showscale=False, xgap=1, ygap=1,
        ))
        styled_fig(fig, title="Missing values scattered randomly",
                   xaxis_title="Samples", yaxis_title="Proteins (sorted by abundance ↓)",
                   height=500)
        fig.update_yaxes(autorange="reversed")
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.markdown('<p class="plot-col-title">MNAR (Abundance-dependent) ❌</p>', unsafe_allow_html=True)
        fig2 = go.Figure(data=go.Heatmap(
            z=mnar_display, colorscale=[[0, "#fef2f2"], [1, "#ef4444"]],
            showscale=False, xgap=1, ygap=1,
        ))
        styled_fig(fig2, title="Low-abundance proteins systematically missing",
                   xaxis_title="Samples", yaxis_title="Proteins (sorted by abundance ↓)",
                   height=500)
        fig2.update_yaxes(autorange="reversed")
        st.plotly_chart(fig2, use_container_width=True)

    with st.expander("📖 How to diagnose & handle missing values"):
        st.markdown("""
**Diagnosis:**
- Plot **missingness vs. intensity** — if low-abundance proteins are preferentially missing, it's MNAR
- Check missingness pattern across conditions — is it condition-specific?

**Imputation strategies:**
- **MCAR** (truly random): kNN imputation, mean/median imputation are acceptable
- **MNAR** (abundance-dependent): Use MinProb, MinDet, QRILC, or Perseus-style left-shifted Gaussian
- **Or don't impute** — use methods that handle missing values natively (MSstats, limma with proper NA handling)

**Validation:**
- Compare results with and without imputation
- Sensitivity analysis: try multiple imputation methods and check consistency
- Be suspicious of any protein that becomes significant only with one imputation method
""")


# ===================================================================
# MODULE: TMT Ratio Compression
# ===================================================================
def page_tmt():
    header_bar()
    st.markdown("### 📉 TMT Ratio Compression")
    st.markdown(severity_badge("MEDIUM"), unsafe_allow_html=True)

    card(
        "<b>The problem:</b> Co-isolation interference in TMT experiments compresses "
        "fold changes toward 1:1. A protein with a true 8-fold change may appear as "
        "only 3-fold in your TMT data. This causes massive underestimation of effect sizes "
        "and loss of statistical power.",
        "warn",
    )

    rng = np.random.RandomState(77)
    n = 200
    true_fc = rng.uniform(-3, 3, n)
    lf_fc = true_fc + rng.normal(0, 0.3, n)
    tmt_fc = true_fc * 0.4 + rng.normal(0, 0.2, n)

    fig = make_subplots(rows=1, cols=2, subplot_titles=("Label-Free (accurate)", "TMT (compressed)"))

    fig.add_trace(go.Scatter(
        x=true_fc, y=lf_fc, mode="markers",
        marker=dict(size=5, color="#0284c7", opacity=0.6), showlegend=False,
    ), row=1, col=1)
    fig.add_trace(go.Scatter(
        x=[-3, 3], y=[-3, 3], mode="lines",
        line=dict(dash="dash", color="#64748b"), showlegend=False,
    ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=true_fc, y=tmt_fc, mode="markers",
        marker=dict(size=5, color="#ef4444", opacity=0.6), showlegend=False,
    ), row=1, col=2)
    fig.add_trace(go.Scatter(
        x=[-3, 3], y=[-3, 3], mode="lines",
        line=dict(dash="dash", color="#64748b"), showlegend=False,
    ), row=1, col=2)

    styled_fig(fig, height=450)
    fig.update_xaxes(title_text="True log2 FC", row=1, col=1)
    fig.update_xaxes(title_text="True log2 FC", row=1, col=2)
    fig.update_yaxes(title_text="Measured log2 FC", row=1, col=1)
    fig.update_yaxes(title_text="Measured log2 FC", row=1, col=2)
    st.plotly_chart(fig, use_container_width=True)

    with st.expander("📖 How to spot & fix TMT ratio compression"):
        st.markdown("""
**Why it happens:**
- The MS2 isolation window captures co-eluting peptides alongside your target
- These co-isolated peptides contribute **~1:1 reporter ion ratios** (they're unchanged between conditions)
- The measured ratio is a weighted average of the true ratio + interference ratios → **compressed toward 1:1**

**Fixes:**
- **SPS-MS3** (Synchronous Precursor Selection) — isolate MS2 fragments and re-fragment for cleaner reporter ions
- **Narrower isolation windows** — reduces co-isolation but lowers sensitivity
- **FAIMS** — gas-phase separation reduces co-eluting species
- **Computational correction** — IRS (Internal Reference Scaling), Complement Reporter Ion approach
- **Always validate** key findings with an orthogonal method (Western blot, PRM, label-free)
""")


# ===================================================================
# MODULE: Contamination
# ===================================================================
def page_contamination():
    header_bar()
    st.markdown("### 🧤 Contamination")
    st.markdown(severity_badge("HIGH"), unsafe_allow_html=True)

    card(
        "<b>The problem:</b> Keratin from skin/hair and polymer (PEG) contamination from "
        "plastics are the most common contaminants in proteomics. They consume your MS "
        "acquisition time and eat into your FDR budget.",
        "danger",
    )

    tab1, tab2, tab3 = st.tabs(["🧤 Keratin", "🧪 Polymer (PEG)", "✅ Clean"])

    rng = np.random.RandomState(33)

    with tab1:
        keratins = ["KRT1", "KRT10", "KRT2", "KRT9", "KRT14", "KRT5", "KRT6A", "KRT16"]
        real_prots = [f"PROT_{i}" for i in range(1, 23)]
        all_names = keratins + real_prots
        abundances = np.sort(rng.lognormal(15, 2, 30))[::-1]
        abundances[:8] *= rng.uniform(2, 5, 8)
        abundances = np.sort(abundances)[::-1]
        colors = ["#ef4444"] * 8 + ["#0284c7"] * 22

        fig = go.Figure(go.Bar(
            x=all_names, y=abundances, marker_color=colors,
        ))
        styled_fig(fig, title="Top 30 proteins — keratin-contaminated sample",
                   xaxis_title="Protein", yaxis_title="Abundance (log scale)",
                   height=450, yaxis_type="log")
        fig.update_xaxes(tickangle=45)
        st.plotly_chart(fig, use_container_width=True)

    with tab2:
        peg_mz = [400 + i * 44.026 for i in range(20)]
        envelope = np.exp(-((np.arange(20) - 10) ** 2) / 30) * 1e6
        peptide_mz = rng.uniform(450, 850, 12)
        peptide_int = rng.uniform(5e4, 3e5, 12)

        fig2 = go.Figure()
        fig2.add_trace(go.Bar(x=peg_mz, y=envelope, width=2, marker_color="#ef4444",
                              name="PEG series (Δ44 Da)"))
        fig2.add_trace(go.Bar(x=peptide_mz, y=peptide_int, width=2,
                              marker_color="#94a3b8", name="Peptide peaks"))
        # Annotate spacing
        fig2.add_annotation(x=peg_mz[5], y=envelope[5] * 1.2, text="Δ44 Da spacing",
                            showarrow=True, arrowhead=2, font=dict(color="#ef4444", size=11))
        styled_fig(fig2, title="Polymer (PEG) contamination — regular 44 Da spacing",
                   xaxis_title="m/z", yaxis_title="Intensity", height=420)
        st.plotly_chart(fig2, use_container_width=True)

    with tab3:
        clean_names = [f"TARGET_{i}" for i in range(1, 31)]
        clean_abund = np.sort(rng.lognormal(15, 2, 30))[::-1]
        fig3 = go.Figure(go.Bar(x=clean_names, y=clean_abund, marker_color="#0284c7"))
        styled_fig(fig3, title="Top 30 proteins — clean sample",
                   xaxis_title="Protein", yaxis_title="Abundance (log scale)",
                   height=450, yaxis_type="log")
        fig3.update_xaxes(tickangle=45)
        st.plotly_chart(fig3, use_container_width=True)

    with st.expander("📖 How to prevent contamination"):
        st.markdown("""
**Keratin prevention:**
- Work under a **laminar flow hood** or in a clean space
- Always wear **nitrile gloves** (not latex)
- Use **keratin-free reagents** and clean tubes with methanol
- Include the **cRAP contaminant database** in your search to prevent keratins from consuming your FDR budget
- Low-bind tubes, pre-washed with LC-MS grade methanol

**Polymer (PEG) prevention:**
- **Avoid plastic** wherever possible — use glass vials for solvents
- Use **LC-MS grade solvents** (not HPLC grade)
- **Low-bind polypropylene** tubes if plastic is necessary
- Watch for **detergent contamination** (Tween, Triton, PEG) from upstream sample prep
- The telltale sign: perfectly regular 44.026 Da spacing = PEG
""")


# ===================================================================
# MODULE: Incomplete Digestion
# ===================================================================
def page_digestion():
    header_bar()
    st.markdown("### ✂️ Incomplete Digestion")
    st.markdown(severity_badge("MEDIUM"), unsafe_allow_html=True)

    card(
        "<b>The problem:</b> When trypsin digestion is incomplete, too many peptides "
        "retain missed cleavage sites (internal Lys/Arg not cut). This reduces "
        "quantification reproducibility and shifts the peptide population.",
        "warn",
    )

    rng = np.random.RandomState(88)
    good_mc = [82, 14, 4]
    bad_mc = [38, 35, 27]
    labels = ["0 MC", "1 MC", "2+ MC"]
    mc_colors = ["#0284c7", "#f59e0b", "#ef4444"]

    col1, col2 = st.columns(2)
    with col1:
        st.markdown('<p class="plot-col-title">Good Digestion ✅</p>', unsafe_allow_html=True)
        fig = go.Figure(go.Pie(
            labels=labels, values=good_mc, hole=0.5,
            marker=dict(colors=mc_colors),
            textinfo="label+percent", textfont=dict(size=13),
        ))
        styled_fig(fig, title="82% zero missed cleavages", height=400)
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.markdown('<p class="plot-col-title">Bad Digestion ❌</p>', unsafe_allow_html=True)
        fig2 = go.Figure(go.Pie(
            labels=labels, values=bad_mc, hole=0.5,
            marker=dict(colors=mc_colors),
            textinfo="label+percent", textfont=dict(size=13),
        ))
        styled_fig(fig2, title="Only 38% zero missed cleavages!", height=400)
        st.plotly_chart(fig2, use_container_width=True)

    with st.expander("📖 How to improve digestion quality"):
        st.markdown("""
**Key factors for good trypsin digestion:**
- **Trypsin:protein ratio** of 1:50 to 1:20 (w/w)
- **Full denaturation** (8 M urea, or heat)
- **Reduction** (DTT or TCEP) and **alkylation** (IAA or CAA)
- **16–18 hours at 37°C**, pH 7.5–8.5
- Fresh trypsin (avoid freeze-thaw cycles)

**When to worry:** < 70% zero-MC rate indicates a problem

**Why it matters for quantification:**
- Missed cleavage peptides are longer → different ionization efficiency
- Variability in MC rate between samples introduces systematic bias
- Some missed cleavage peptides overlap with other proteins → ambiguous quantification
""")


# ===================================================================
# MODULE: LC Carryover
# ===================================================================
def page_carryover():
    header_bar()
    st.markdown("### 🔄 LC Carryover")
    st.markdown(severity_badge("MEDIUM"), unsafe_allow_html=True)

    card(
        "<b>The problem:</b> Hydrophobic peptides stick to the LC column and elute as "
        "ghost peaks in subsequent runs. If your blank run shows peaks, those peptides "
        "are contaminating the next sample.",
        "warn",
    )

    rng = np.random.RandomState(44)
    t = np.linspace(0, 90, 1800)

    # Run 1: several peaks
    peaks_rt = [25, 35, 45, 55, 60, 70]
    peaks_h = [8e6, 6e6, 7e6, 9e6, 5e6, 1e7]
    run1 = np.zeros_like(t)
    for rt, h in zip(peaks_rt, peaks_h):
        run1 += h * np.exp(-((t - rt) ** 2) / 2)
    run1 += rng.normal(0, 1e5, len(t))
    run1 = np.clip(run1, 0, None)

    # Run 2 (blank): only late peaks at 1-5% intensity
    run2 = np.zeros_like(t)
    for rt, h in zip([55, 60, 70], [9e6, 5e6, 1e7]):
        run2 += h * rng.uniform(0.01, 0.05) * np.exp(-((t - rt) ** 2) / 2)
    run2 += rng.normal(0, 5e4, len(t))
    run2 = np.clip(run2, 0, None)

    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.08,
                        subplot_titles=("Run 1: Sample", "Run 2: Blank (should be empty!)"))

    fig.add_trace(go.Scatter(x=t, y=run1, fill="tozeroy",
                             line=dict(color="#0284c7", width=1),
                             fillcolor="rgba(2,132,199,0.25)", showlegend=False), row=1, col=1)
    fig.add_trace(go.Scatter(x=t, y=run2, fill="tozeroy",
                             line=dict(color="#ef4444", width=1),
                             fillcolor="rgba(239,68,68,0.25)", showlegend=False), row=2, col=1)
    # Carryover annotations
    for rt in [55, 60, 70]:
        fig.add_annotation(x=rt, y=run2[np.argmin(np.abs(t - rt))] * 1.5,
                           text="⬆ Carryover!", font=dict(color="#ef4444", size=10),
                           showarrow=False, row=2, col=1)

    styled_fig(fig, height=550)
    fig.update_xaxes(title_text="Retention Time (min)", row=2, col=1)
    fig.update_yaxes(title_text="Intensity", row=1, col=1)
    fig.update_yaxes(title_text="Intensity", row=2, col=1)
    st.plotly_chart(fig, use_container_width=True)

    with st.expander("📖 How to prevent & detect LC carryover"):
        st.markdown("""
**Prevention:**
- Run **1–2 blank injections** between samples (especially after high-abundance samples)
- Use **aggressive column washes** (high organic solvent, extended wash gradients)
- Consider a **trap column** with back-flush capability
- **Randomize** sample injection order (don't run all treatments consecutively)

**Detection:**
- **Always include blank runs** in your acquisition sequence
- Check blanks for identifiable peptides — even 1 ID is a warning
- **> 0.1% carryover** may affect quantification in sensitive experiments
- Look for peaks at late retention times (hydrophobic peptides stick most)
""")


# ===================================================================
# STUB modules for 10-17 (to be completed in phase 2)
# ===================================================================
def _stub(title, emoji):
    header_bar()
    st.markdown(f"### {emoji} {title}")
    card(
        f"<b>{title}</b> — This module is coming soon! "
        "Check back for the full interactive demonstration.",
        "info",
    )


def page_dia():
    _stub("DIA Pitfalls", "📡")

def page_missing_proteins():
    _stub("Missing Proteins", "🔍")

def page_apms():
    _stub("AP-MS Normalization", "🧲")

def page_western():
    _stub("Western vs MS", "🩻")

def page_pvalue_shopping():
    _stub("P-Value Shopping", "🎰")

def page_cherrypick():
    _stub("Cherry-Picking Data", "🍒")

def page_multiple_testing():
    _stub("Multiple Testing", "📐")

def page_quiz():
    _stub("Quiz", "🧠")


# ===================================================================
# MAIN
# ===================================================================
def main():
    with st.sidebar:
        st.markdown(
            f'<div style="text-align:center;padding:0.5rem 0;">{logo_html(48, "sb")}</div>',
            unsafe_allow_html=True,
        )
        st.markdown(
            '<p style="text-align:center;font-family:IBM Plex Mono,monospace;'
            'font-size:1.1rem;font-weight:700;color:#38bdf8;margin:0;">CrackPipe</p>',
            unsafe_allow_html=True,
        )
        st.markdown("---")
        page = st.radio("NAVIGATION", NAV_ITEMS, label_visibility="collapsed")
        st.markdown("---")
        st.markdown(
            '<p style="font-size:0.75rem;color:#64748b;text-align:center;">'
            'Inspired by <a href="https://github.com/Nesvilab/FragPipe" '
            'style="color:#38bdf8;">FragPipe</a><br>Not affiliated with Nesvilab</p>',
            unsafe_allow_html=True,
        )

    routes = {
        NAV_ITEMS[0]: page_home,
        NAV_ITEMS[1]: page_chimeric,
        NAV_ITEMS[2]: page_spray,
        NAV_ITEMS[3]: page_masscal,
        NAV_ITEMS[4]: page_batch,
        NAV_ITEMS[5]: page_missing,
        NAV_ITEMS[6]: page_tmt,
        NAV_ITEMS[7]: page_contamination,
        NAV_ITEMS[8]: page_digestion,
        NAV_ITEMS[9]: page_carryover,
        NAV_ITEMS[10]: page_dia,
        NAV_ITEMS[11]: page_missing_proteins,
        NAV_ITEMS[12]: page_apms,
        NAV_ITEMS[13]: page_western,
        NAV_ITEMS[14]: page_pvalue_shopping,
        NAV_ITEMS[15]: page_cherrypick,
        NAV_ITEMS[16]: page_multiple_testing,
        NAV_ITEMS[17]: page_quiz,
    }

    routes.get(page, page_home)()


if __name__ == "__main__":
    main()
