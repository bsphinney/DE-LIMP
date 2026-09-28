# Hsp104 and Ssa1 pulldowns in heat-shocked yeast, Early vs Late (synthetic example)

*Synthetic test fixture for make_podcast.py: every value and finding here is invented. There is
no such study. Submission 0421.*

## Results at a glance

**24** samples · **6** groups · **5** contrasts · **3,412** proteins tested. 2 baits (Hsp104,
Ssa1) plus IgG controls, Early vs Late heat shock at 42 °C, n = 4 per group.

## How the data were measured

timsTOF HT, dia-PASEF, 44 minute gradient. DIA-NN 2.6.1, library-free search with
match-between-runs, 1% precursor FDR: 41,236 precursors mapped to 3,412 protein groups. DE with
limpa (DPC-Quant) and limma; significant = adjusted p < 0.05 (Benjamini–Hochberg), no
fold-change filter. 18–41% of the values per sample are inferred (median 27%), not measured;
PropObs gives the observed fraction per protein. Empirical Bayes moderates the variances: limma
borrows variance information across thousands of proteins.

## Key findings

Hsp104 is the top hit in its own pulldown (log2FC 9.84, adj.P 3.1e-12), hundreds of times more
than in IgG. Ssa1 tops its pulldown
(log2FC 8.27, adj.P 7.4e-11). Early Hsp104 vs IgG: 412 significant (367 up). Late Hsp104 vs IgG:
158 significant. Sis1 and Ydj1 co-purify with Ssa1. Hsp26 and Hsp42 rise in Late Hsp104
pulldowns (log2FC 2.35 and 1.92); Hsp42 was detected in fewer than half of the runs.

## Caveats

The Late IgG runs are thinner: 2,104 and 2,388 proteins detected against a mean of 2,961 for the
Early IgG runs. Hsp104 itself falls by −1.1 in Late and its partners by a median of −0.8
(n = 96): mostly bait recovery, not rewiring.

## Data Quality Notes

Keratins (KRT1, KRT10) appear in 3 of 24 runs and were removed as contaminants before
quantification.

## Files

README.html, Analysis_Report.html, the DE_dpc_*.csv tables (with the Detected_<group> columns, detected in k of n runs per group, and the Evidence column), methods.md.
