# De novo → homology species ID (SKILL.md step 6e)

For samples whose organism is unknown or has no proteome. The method and every measurement are
written up in `docs/denovo_decoy_method.html` in the skill's repository, which is not part of the
installed plugin, so the settings that decide how the route is run, and the numbers behind each,
are all here.

> **Exploratory, not a species identification.** The result is the closest sequenced relative
> in nr that the de novo peptides point to. It has no error control: no minimum evidence, no
> decoy tally for the taxon. On the HeLa entrapment (true answer: human) the peptides accepted at
> the model's nominal 1% were only 59.8% correct, and best-hit attribution put 29.5% of them in
> primates. A species that is not in nr comes out as its nearest sequenced relative. Every output
> says so (`status: exploratory` in `<out>_species_call.json` and in the FDR script's record).
> Tell the user this in plain words. Do not use it for forensic or client conclusions until the
> route is validated on samples of known species.

## The route and the five measured settings

For samples where the organism is unknown or has no proteome — forensic hair,
feather, fossil, wildlife — the route is Casanovo de novo → DIAMOND vs nr → LCA.
Four scripts, in order:

```
python3 scripts/denovo_decoy_gen.py <mzML> <decoy.mgf> 0.8 1337     # 1. decoy spectra
python3 scripts/denovo_search_cmd.py --query peptides.fasta \
        --db $NR --out hits.tsv                                     # 2. emits the DIAMOND cmd
python3 scripts/denovo_homology_fdr.py --target-hits ... --decoy-hits ... \
        --target-conf real.mztab --decoy-conf decoy.mztab \
        --n-target N --n-decoy N --out accepted.tsv                 # 3. confidence tier
python3 scripts/denovo_lca.py hits.tsv out_prefix                   # 4. closest relative
```

**Five things here are measured, and each is a way to get it silently wrong.**
- **FRAC ≥ 0.8 for the decoy.** At 0.5 the decoy leaks — Casanovo reconstructs
  abundant real peptides from half-kept peaks, and the null stops being a null
  (356 decoy peptides hitting nr vs 36 at 0.8).
- **`--max-target-seqs 25`, never 1.** The species call is LCA over *all* hits;
  with one hit per peptide LCA silently becomes best-hit, which mis-assigns ~36%.
  `denovo_lca.py` warns if the input averages < 2 hits/peptide.
- **The EXTRA short-peptide DIAMOND settings**, and the Riffle build — stock
  DIAMOND finds almost nothing on 7–30 aa queries. EXTRA nearly doubles confirmed
  homology (9,801 → 19,671) and the gain is concentrated at ≤15 aa (168 → 1,661).
  Costs ~5×.
- **Linear model, calibrated on E-value.** Gradient boosting overfits the decoy
  and scores *worse than raw bitscore* (9,554 vs 10,404). Raw bitscore is
  length-confounded and starves short peptides.
- **Five features, not thirteen.** On ground truth, more features accept more
  peptides at *lower* precision and recover *fewer* correct ones (12-feature:
  2,532 accepted / 55.5% / 1,405 correct, vs 5-feature: 2,400 / 59.8% / 1,436).

**The model's cut is not an FDR on the sequences.** `--fdr` is a target-decoy level for chance
homology (decoy spectra whose de novo peptides still hit nr). It says nothing about whether a
de novo sequence is right: at a nominal 0.01 the accepted peptides were 59.8% correct. Never
describe it as "1% FDR" for the peptides or the species.

**Do not gate the species call on the FDR model.** Filtering lifts peptide
precision 42.4% → 59.8% but moves primate attribution only 29.5% → 32.9%, leaves
mammal flat, and throws away 715 correct peptides. Species error comes from
conserved peptides matching many taxa — that is LCA's job. Run LCA over the
unfiltered hits and report the model's score as a confidence tier.

De novo sequence error is the binding error mode: 45% raw, 21% at conf ≥ 0.95,
12% at conf ≥ 0.99. Say which threshold you used.

→ method + all measurements: `docs/denovo_decoy_method.html` in the skill's repository.
