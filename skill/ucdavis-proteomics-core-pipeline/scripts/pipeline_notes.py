#!/usr/bin/env python3
"""
pipeline_notes.py  --  What the pipeline is doing right now, in plain language, plus a
rotating proteomics note to read while a multi-hour search grinds away.

watch_run.sh calls this each poll and folds the result into its JSON, so the orchestrator
can tell the user something true and specific instead of "still running" for four hours.

Every FACT below carries a `source`. Nothing here is invented: if a claim could not be
attributed it did not go in. Prefer qualitative statements over numbers -- a number the
user cannot check is worse than no number (DE-LIMP architectural rule #2).

Usage:
  python3 pipeline_notes.py --stage step2 --index 7
  python3 pipeline_notes.py --stage de --index 0 --no-fact
"""
import argparse, json

# --- what is happening, and why it takes as long as it does -------------------
# Keys match the stages watch_run.sh infers from the output directory.
STAGES = {
    "step1": ("Predicting a spectral library from your FASTA",
              "A deep-learning model predicts the retention time, ion mobility and fragment "
              "pattern of every peptide your proteins could produce. One job, no raw data "
              "read yet -- this is why nothing looks like it is touching your files."),
    "step2": ("First pass: searching each file against the predicted library",
              "Each raw file is searched independently as its own cluster task, which is "
              "where the parallel chain earns its time back. Each finished file leaves a "
              ".quant behind, so progress here is real files completed, not an estimate."),
    "step3": ("Building an empirical library from what was actually seen",
              "The first-pass results are pooled into a library measured from your own data "
              "rather than predicted. This replaces match-between-runs, and it is why the "
              "second pass identifies more than the first."),
    "step4": ("Final pass: re-searching every file against the empirical library",
              "The same per-file fan-out as the first pass, now against the better library. "
              "Files that failed the first pass are skipped rather than blocking the run."),
    "step5": ("Assembling the cross-run report",
              "Per-file results are combined into one protein-group matrix, with quantities "
              "made comparable across runs. This produces report.parquet, the file the "
              "statistics step reads."),
    "single": ("Searching your files",
               "One job is working through the raw data: extracting fragment traces, matching "
               "them to peptides, and scoring each match against decoys to control the false "
               "discovery rate."),
    "de": ("Differential expression",
           "Protein abundances are modelled per group and moderated across proteins, so a "
           "protein measured in few replicates borrows variance information from the rest."),
}

# --- notes. Each: (text, source, stages it suits best or None for any) --------
FACTS = [
    ("Trypsin cuts after lysine and arginine, but not when proline follows. Nearly every "
     "search you will ever run assumes those rules.",
     "Keil, Specificity of Proteolysis (1992)", None),
    ("DIA-NN is short for 'Data-Independent Acquisition by Neural Networks' — the neural "
     "networks are what turn a FASTA into a usable spectral library without you measuring one.",
     "Demichev et al. 2020, Nat Methods 17:41-44", ("step1", "step2")),
    ("A predicted library is a hypothesis about your peptides; an empirical one is a "
     "measurement of them. Searching twice — predicted, then empirical — is why the second "
     "pass finds more.",
     "Demichev et al. 2020, Nat Methods 17:41-44", ("step3", "step4")),
    ("An Orbitrap has no magnet. Ions orbit a spindle-shaped central electrode, and their "
     "oscillation frequency scales with the inverse square root of m/z — mass is measured as "
     "a frequency, which is why the resolution is so high.",
     "Makarov 2000, Anal Chem 72:1156-1162", None),
    ("On a timsTOF, ions are held stationary against a flowing gas by an electric field, then "
     "released in order of size and shape. That extra separation is the 'TIMS' in the name.",
     "Meier et al. 2015, J Proteome Res 14:5378-5387", None),
    ("Two peptides from the same protein can differ enormously in how well they ionize. That "
     "is why label-free quantification compares a peptide to itself across runs, and never to "
     "another peptide.",
     "Cox et al. 2014, MaxLFQ, Mol Cell Proteomics 13:2513-2526", ("step5", "de")),
    ("Shared peptides mean the protein you detected is sometimes ambiguous. Software reports "
     "'protein groups' rather than proteins precisely because the evidence cannot always tell "
     "family members apart — this is the protein inference problem.",
     "Nesvizhskii & Aebersold 2005, Mol Cell Proteomics 4:1419-1440", ("step5", "de")),
    ("Controlling the false discovery rate at 1% for peptides does not give you 1% at the "
     "protein level. The error rates are separate, and protein-level FDR is the stricter one "
     "to satisfy.",
     "Nesvizhskii 2010, J Proteomics 73:2092-2123", ("step5", "de")),
    ("Keratin is the classic contaminant — skin, hair and dust find their way into almost "
     "every sample. It is only a problem when it is unevenly distributed across your groups.",
     "Frankenfield et al. 2022, J Proteome Res 21:2104-2113 (the contaminant "
     "library appended to the search)", None),
    ("The moderated t-test used here borrows variance information across all proteins, which "
     "is what makes small-replicate experiments analysable at all.",
     "Smyth 2004, Stat Appl Genet Mol Biol 3:Article3", ("de",)),
    ("The Benjamini-Hochberg procedure does not tell you which discoveries are wrong. It "
     "promises that, on average, only a set fraction of the ones you report will be.",
     "Benjamini & Hochberg 1995, J R Stat Soc B 57:289-300", ("de",)),
    ("Decoy sequences exist to be wrong on purpose. How often a search prefers a decoy over a "
     "real peptide is what calibrates the false discovery rate.",
     "Elias & Gygi 2007, Nat Methods 4:207-214", ("step2", "step4")),
    ("Monoisotopic mass uses only the lightest isotope of each element; average mass uses the "
     "natural mixture. Mixing them up shifts a peptide by about one dalton per 1000 — enough "
     "to lose the identification entirely.",
     "IUPAC isotopic composition conventions", None),
    ("In data-independent acquisition the instrument stops choosing which peptides to "
     "fragment and simply fragments everything in a window. Nothing is missed because it was "
     "not selected — the difficulty moves from acquisition into the software.",
     "Venable et al. 2004, Nat Methods 1:39-45", ("step2", "step4", "single")),
    ("Thousands of human proteins still have little or no direct protein-level evidence, "
     "despite their genes being known for decades. They are tracked as 'missing proteins'.",
     "neXtProt protein-existence (PE) levels, nextprot.org", None),
    ("Retention time is a strong second opinion. A peptide eluting far from where it should "
     "is treated as suspicious even when its fragments match well.",
     "Escher et al. 2012, Proteomics 12:1111-1121 (iRT)", ("step1", "step2")),
    ("Cysteine is usually reduced and capped with iodoacetamide before digestion, which is "
     "why searches carry a fixed +57 Da carbamidomethyl modification. Forgetting it costs "
     "most of your cysteine-containing peptides.",
     "Standard bottom-up sample preparation", None),
    ("A protein's abundance in a cell can span roughly ten orders of magnitude. No single "
     "mass-spec run sees all of it at once — depth is always a compromise with time.",
     "Beck et al. 2011, Mol Syst Biol 7:549", None),
]


def pick(stage, index, prefer_stage=True):
    """Deterministic rotation: same (stage, index) always yields the same note, and a
    caller stepping index forward each poll walks the whole list without repeating.

    Stage-relevant notes come first, then everything else -- rather than cycling only the
    matching ones, which for a stage with four facts would repeat every four minutes at a
    60-second poll. The full list is always reachable, so a long search keeps varying."""
    matched = [f for f in FACTS if f[2] and stage in f[2]] if prefer_stage else []
    pool = matched + [f for f in FACTS if f not in matched]
    return pool[index % len(pool)]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage", default="single", help="step1..step5 | single | de")
    ap.add_argument("--index", type=int, default=0, help="poll counter; rotates the note")
    ap.add_argument("--no-fact", action="store_true")
    a = ap.parse_args()

    headline, why = STAGES.get(a.stage, STAGES["single"])
    out = {"stage": a.stage, "doing": headline, "why": why}
    if not a.no_fact:
        text, source, _ = pick(a.stage, a.index)
        out["note"] = text
        out["note_source"] = source
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
