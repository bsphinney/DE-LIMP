#!/usr/bin/env python3
"""The sample-sheet check reads every sheet exactly as skill 2.13.2 did; only who is told changes.

submission_report's sheet check takes a trailing replicate index off the sheet's condition names
to find its groups, and compares them with the groups the analysis used ("conditions differ").
Read that way, ATG5 and ATG7 are replicates of "ATG" -- and so are the sample IDs Parental12 and
Parental51: on names alone the two cannot be told apart. Three reviews found every rule that
changed the groups reading some sheets worse than 2.13.2 did, so the groups stay 2.13.2's
exactly. What 2.14.1 changes is the note: a part of it that rests on a group made by taking a
number glued to the names off them (no separator, ')', '#', rep / replicate or ID letter) is
withheld from the report and said to staff instead (submission_report.sheet_cautions).

This test holds today's code to 2.13.2's on a corpus of sheet shapes: the reviews' examples,
then stems x separators x numbering (per condition, on across conditions, sample IDs) x levels
(L/R, R/S, A/B) x mixed formats x slips x case x duplicated names, generated with a fixed seed.
The groups the analysis used are each sample's true condition, known by construction -- never
computed with the code under test. fixtures/replicate_parity_2_13_2.json holds each sheet with
2.13.2's sheet groups, its "every sample has its own condition" note and its "conditions
differ" note, made by 2.13.2's own code (dd932cc), so the test needs no git. Groups and the
"every sample" note must be identical, always. The "conditions differ" note must be 2.13.2's
less exactly the parts about a glued-number group, and those parts must be in the staff
caution -- which is there only then.

Regenerate the fixture (new shapes; 2.13.2's code is what must answer):
  git -C <repo> archive dd932cc skill/ucdavis-proteomics-core-pipeline/scripts | tar -x -C <dir>
  python3 tests/test_replicate_parity.py --write-fixture <dir>/<that scripts path>

Stdlib only.
"""
import csv
import importlib
import json
import os
import random
import re
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
FIXTURE = os.path.join(HERE, "fixtures", "replicate_parity_2_13_2.json")
RECORD = os.path.join(HERE, "fixtures", "coreomics_prot0756", "submission.json")
SEED, N_GENERATED = 20261020, 520

def twice(names):
    """Every name listed twice, in place (an IP-MS sheet's IP and input rows)."""
    return [x for x in names for _ in (0, 1)]


# The reviews' sheets: (condition names, the analysed group of each sample -- by hand).
REVIEWS = [
    ([f"Cis_{x}_{n}" for x in "RS" for n in (1, 2, 3)], ["Cis_R"] * 3 + ["Cis_S"] * 3),
    ([f"Cis_{x}{n}" for x in "RS" for n in (1, 2, 3)], ["Cis_R"] * 3 + ["Cis_S"] * 3),
    ([f"Brain_{x}-{n}" for x in "RL" for n in (1, 2, 3)], ["Brain_R"] * 3 + ["Brain_L"] * 3),
    ([f"Brain_{x} {n}" for x in "RL" for n in (1, 2, 3)], ["Brain_R"] * 3 + ["Brain_L"] * 3),
    ([f"Brain_{x}{n}" for x in "LR" for n in (1, 2, 3)], ["Brain_L"] * 3 + ["Brain_R"] * 3),
    (["Cis"] * 3 + ["Cis_R_1", "Cis_R_2", "Cis_R_3"], ["Cis"] * 3 + ["Cis_R"] * 3),
    (["A549"] * 3 + ["A549_R1", "A549_R2", "A549_R3"], ["A549"] * 3 + ["A549_R"] * 3),
    (["Cis_S"] * 3 + ["Cis_R1", "Cis_R2", "Cis_R3"] + ["Cis"] * 3,
     ["Cis_S"] * 3 + ["Cis_R"] * 3 + ["Cis"] * 3),
    (["Brain_L"] * 3 + ["Brain_R1", "Brain_R2", "Brain_R3"], ["Brain_L"] * 3 + ["Brain_R"] * 3),
    (["WT1", "WT2", "WT3", "KO", "KO", "KO"], ["WT"] * 3 + ["KO"] * 3),
    (["Ctrl"] * 3 + ["Treated1", "Treated2", "Treated3"], ["Ctrl"] * 3 + ["Treated"] * 3),
    (["HSP90"] * 3 + ["WT1", "WT2", "WT3"], ["HSP90"] * 3 + ["WT"] * 3),
    (["WT1", "WT2", "WT2", "WT3", "KO1", "KO2", "KO3"], ["WT"] * 4 + ["KO"] * 3),
    (["WT_1", "WT_2", "WT_2", "WT_3", "KO1", "KO2", "KO3"], ["WT"] * 4 + ["KO"] * 3),
    (["WT_1", "WT_2", "WT_2", "WT_3", "KO_1", "KO_2", "KO_3"], ["WT"] * 4 + ["KO"] * 3),
    (["WT_1", "WT_1", "WT_2", "WT_2", "WT_3", "KO_1", "KO_2", "KO_3"], ["WT"] * 5 + ["KO"] * 3),
    (["Mutant-1"] * 3 + ["Mutant-2"] * 3 + ["Mutant-3", "Mutant-4"],
     ["Mutant_1"] * 3 + ["Mutant_2"] * 3 + ["Mutant_3", "Mutant_4"]),
    (["Mutant-1"] * 3 + ["Mutant-2"] * 3 + ["Mutant-3"], ["Mutant_1"] * 3 + ["Mutant_2"] * 3
     + ["Mutant_3"]),
    (["Mutant-1"] * 3 + ["Mutant-2"] * 3 + ["WT"] * 3,
     ["Mutant_1"] * 3 + ["Mutant_2"] * 3 + ["WT"] * 3),
    ([f"{g}{n}" for g in ("WT", "KO") for n in (2, 3, 4)], ["WT"] * 3 + ["KO"] * 3),
    (["WT2", "WT3", "WT4", "KO", "KO", "KO"], ["WT"] * 3 + ["KO"] * 3),
    ([f"{g}{n}" for g in ("WT", "KO") for n in (1, 2, 3)], ["WT"] * 3 + ["KO"] * 3),
    ([f"WT{n}" for n in (1, 2, 3)] + [f"KO_{n}" for n in (1, 2, 3)], ["WT"] * 3 + ["KO"] * 3),
    (["wt1", "wt2", "WT3", "KO1", "KO2", "KO3"], ["WT"] * 3 + ["KO"] * 3),
    ([f"{g}{n}" for g in ("WT", "KO") for n in (1, 2, 4)], ["WT"] * 3 + ["KO"] * 3),
    (["ATG5"] * 3 + ["ATG7"] * 3 + ["WT"] * 3, ["ATG5"] * 3 + ["ATG7"] * 3 + ["WT"] * 3),
    (["ATG5"] * 3 + ["ATG7"] * 3 + ["HSP90"] * 3, ["ATG5"] * 3 + ["ATG7"] * 3 + ["HSP90"] * 3),
    (["ATG5"] * 3 + ["ATG7"] * 3, ["ATG5"] * 3 + ["ATG7"] * 3),
    (["HSP70"] * 3 + ["HSP90"] * 3 + ["WT"] * 3, ["HSP70"] * 3 + ["HSP90"] * 3 + ["WT"] * 3),
    (["ATG5", "ATG7", "CD45", "HSP90"], ["ATG5", "ATG7", "CD45", "HSP90"]),
    (["ATG5", "ATG7", "WT_1", "WT_2"], ["ATG5", "ATG7", "WT", "WT"]),
    (["HSP70", "HSP90", "WT1", "WT2"], ["HSP70", "HSP90", "WT", "WT"]),
    (["Kv2.1", "Kv2.2", "WT_1", "WT_2"], ["Kv2.1", "Kv2.2", "WT", "WT"]),
    (["WT1"] * 2 + ["WT2"] * 2 + ["KO1"] * 2 + ["KO2"] * 2, ["WT"] * 4 + ["KO"] * 4),
    (["WT_1"] * 2 + ["WT_2"] * 2 + ["KO_1"] * 2 + ["KO_2"] * 2, ["WT"] * 4 + ["KO"] * 4),
    ([f"WT_{n}" for n in (1, 2, 3)] * 2 + [f"KO_{n}" for n in (1, 2, 3)] * 2,
     ["WT"] * 6 + ["KO"] * 6),
    (["T0"] * 3 + ["T24"] * 3 + ["Ctrl"] * 3, ["T0"] * 3 + ["T24"] * 3 + ["Ctrl"] * 3),
    ([f"{g}_{n}" for g in ("T0", "T24") for n in (1, 2)], ["T0"] * 2 + ["T24"] * 2),
    (["Day1", "Day2", "Day3", "Day4"], ["Day1", "Day2", "Day3", "Day4"]),
    (["Day1", "Day2", "Day3", "Ctl1", "Ctl2", "Ctl3"], ["Day"] * 3 + ["Ctl"] * 3),
    (["sample 1", "sample 2", "sample 3"], ["sample 1", "sample 2", "sample 3"]),
    (["Old - A - Rep 1", "Old - A - Rep 2", "Old - B - Rep 1", "Old - B - Rep 2"],
     ["Old_A"] * 2 + ["Old_B"] * 2),
    (["WT"] * 3 + ["WT4"] + ["KO"] * 3, ["WT"] * 4 + ["KO"] * 3),
    (["Strep_1", "Strep_2", "SamplePrep_1", "SamplePrep_2"], ["Strep"] * 2 + ["SamplePrep"] * 2),
    (["ATG5", "ATG5", "ATG7", "ATG8"], ["ATG5", "ATG5", "ATG7", "ATG8"]),
    ([f"{g}_{n}" for g in ("ATG5", "ATG7") for n in (1, 2, 3)], ["ATG5"] * 3 + ["ATG7"] * 3),
    ([f"{g}-{n}" for g in ("HSP70", "HSP90") for n in (1, 2, 3)], ["HSP70"] * 3 + ["HSP90"] * 3),
    ([f"{g}{n}" for g in ("WT", "KO") for n in (1, 2, 3)] + ["Pool"],
     ["WT"] * 3 + ["KO"] * 3 + ["Pool"]),
    (["WT1", "WT2", "WT3", "KO1", "KO2", "KO4"], ["WT"] * 3 + ["KO"] * 3),
    ([f"{g} \u2013 {n}" for g in ("WT", "KO") for n in (1, 2)], ["WT"] * 2 + ["KO"] * 2),
    ([f"{g}.{n}" for g in ("WT", "KO") for n in (1, 2)], ["WT"] * 2 + ["KO"] * 2),
    (["Mutant-1"] * 2 + ["WT_1", "WT_2"], ["Mutant_1"] * 2 + ["WT"] * 2),
    (["WT_1", "WT_2", "WT_2", "WT_3"], ["WT"] * 4),
    (["WT1"] * 2 + ["WT2"] * 2 + ["KO"] * 2, ["WT"] * 4 + ["KO"] * 2),
    (["ATG5"] * 3 + ["ATG7"] * 3 + ["ATG8"] + ["WT"] * 2,
     ["ATG5"] * 3 + ["ATG7"] * 3 + ["ATG8"] + ["WT"] * 2),
    (["ATG5"] * 3 + ["ATG7"] * 3 + ["WT_1", "WT_2", "WT_3"],
     ["ATG5"] * 3 + ["ATG7"] * 3 + ["WT"] * 3),
    ([f"Dox(+)_{n}" for n in (1, 2, 3)] + [f"Dox(-)_{n}" for n in (1, 2, 3)],
     ["Dox(+)"] * 3 + ["Dox(-)"] * 3),
    (["WT_Rep1", "WT_Rep2", "KO replicate 1", "KO replicate 2"], ["WT"] * 2 + ["KO"] * 2),
    (["WTrep1", "WTrep2", "KORep1", "KORep2"], ["WT"] * 2 + ["KO"] * 2),
    # every replicate listed twice (IP and input rows), with a dropped replicate or from 2
    ([x for x in ("WT_rep1", "WT_rep2", "WT_rep4") for _ in (0, 1)]
     + [x for x in ("KO_rep1", "KO_rep2", "KO_rep3") for _ in (0, 1)], ["WT"] * 6 + ["KO"] * 6),
    ([x for x in ("IP_rep1", "IP_rep2", "IP_rep4") for _ in (0, 1)] + ["IgG"] * 6,
     ["IP"] * 6 + ["IgG"] * 6),
    (["WT1"] * 2 + ["WT3"] * 2 + ["KO1"] * 2 + ["KO2"] * 2, ["WT"] * 4 + ["KO"] * 4),
    ([x for x in ("WT_2", "WT_3", "WT_4") for _ in (0, 1)]
     + [x for x in ("KO_2", "KO_3", "KO_4") for _ in (0, 1)], ["WT"] * 6 + ["KO"] * 6),
    (["CD4"] * 2 + ["CD8"] * 2 + ["WT"] * 2, ["CD4"] * 2 + ["CD8"] * 2 + ["WT"] * 2),
    # only the replicate word says so: replicates 5 and 7, each listed twice
    ([x for x in ("WT_rep5", "WT_rep7", "KO_rep5", "KO_rep7") for _ in (0, 1)],
     ["WT"] * 4 + ["KO"] * 4),
    # samples numbered on across conditions, or by sample / animal ID, each listed twice
    (twice(["Ctrl_1", "Ctrl_2", "Ctrl_3", "LPS_4", "LPS_5", "LPS_6"]), ["Ctrl"] * 6 + ["LPS"] * 6),
    (twice(["Bait_IP 1", "Bait_IP 2", "Bait_IP 3", "IgG 4", "IgG 5", "IgG 6"]),
     ["Bait_IP"] * 6 + ["IgG"] * 6),
    (twice(["WT1", "WT2", "WT3", "KO4", "KO5", "KO6"]), ["WT"] * 6 + ["KO"] * 6),
    (twice([f"{g} {n}" for g, ns in (("WT", (1, 2, 3)), ("KO", (4, 5, 6)), ("DKO", (7, 8, 9)))
            for n in ns]), ["WT"] * 6 + ["KO"] * 6 + ["DKO"] * 6),
    (twice(["KO_R1", "KO_R2", "KO_R3", "WT_R4", "WT_R5", "WT_R6"]), ["KO"] * 6 + ["WT"] * 6),
    (twice(["WT_M12", "WT_M15", "KO_M13", "KO_M14"]), ["WT"] * 4 + ["KO"] * 4),
    (twice(["WT_1", "WT_2", "KO_3", "KO_4"]), ["WT"] * 4 + ["KO"] * 4),
    (twice(["WT1", "WT2", "KO3", "KO4"]), ["WT"] * 4 + ["KO"] * 4),
    (twice(["Parental_12", "Parental_51"]) + ["Dox(+)"] * 4, ["Parental"] * 4 + ["Dox(+)"] * 4),
    (twice(["FLAG_IP-4", "FLAG_IP-5"]) + ["Ctrl"] * 6, ["FLAG_IP"] * 4 + ["Ctrl"] * 6),
    (twice(["Sham_R5", "Sham_R6"]) + ["T48"] * 8, ["Sham"] * 4 + ["T48"] * 8),
    (twice(["Mock_S5", "Mock_S34", "Mock_S59"]) + ["Cis_S"] * 6, ["Mock"] * 6 + ["Cis_S"] * 6),
    (["siNT"] * 6 + twice(["sg24", "sg25"]), ["siNT"] * 6 + ["sg"] * 4),
    (twice(["Bait15", "Bait17", "Normal.15", "Normal.18"]), ["Bait"] * 4 + ["Normal"] * 4),
    # the boundaries of replicate numbering and of the markers, each listed twice
    (twice(["WT3", "WT4", "WT5"]) + ["KO"] * 6, ["WT"] * 6 + ["KO"] * 6),
    (twice(["WT1", "WT5"]) + ["KO"] * 4, ["WT"] * 4 + ["KO"] * 4),
    (twice(["WT#5", "WT#7"]) + ["KO"] * 4, ["WT"] * 4 + ["KO"] * 4),
    (twice(["WT_REP5", "WT_REP7"]) + ["KO"] * 4, ["WT"] * 4 + ["KO"] * 4),
    (twice(["WT replicate 5", "WT replicate 7"]) + ["KO"] * 4, ["WT"] * 4 + ["KO"] * 4),
    (twice(["Old - A - R5", "Old - A - R7"]) + ["Old - B"] * 4, ["Old_A"] * 4 + ["Old_B"] * 4),
    (twice(["WT bio rep 5", "WT bio rep 7"]) + ["KO"] * 4, ["WT"] * 4 + ["KO"] * 4),
    # one condition numbered by glued sample IDs, others not: the twin of ATG5 / ATG7 (r9)
    (twice(["Parental12", "Parental51"]) + ["Dox(+)"] * 4, ["Parental"] * 4 + ["Dox(+)"] * 4),
    (twice(["Parental12", "Parental51"]) + ["Dox(+)"] * 4,
     twice(["Parental12", "Parental51"]) + ["Dox(+)"] * 4),
    (twice(["WT3", "WT4", "WT6"]) + ["KO"] * 6, ["WT"] * 6 + ["KO"] * 6),
    (twice(["WT3", "WT4", "WT6"]) + ["KO"] * 6, twice(["WT3", "WT4", "WT6"]) + ["KO"] * 6),
    (twice(["Bait15", "Bait17"]) + ["Ctrl"] * 6, ["Bait"] * 4 + ["Ctrl"] * 6),
    (twice(["Tumor12", "Tumor35", "Tumor40"]) + ["Normal"] * 6, ["Tumor"] * 6 + ["Normal"] * 6),
    (twice(["WT12", "WT40"]) + ["KO"] * 4 + ["Het"] * 4, ["WT"] * 4 + ["KO"] * 4 + ["Het"] * 4),
    (["WT"] * 2 + twice(["WT5", "WT9"]) + ["KO"] * 6, ["WT"] * 6 + ["KO"] * 6),
    (twice(["WT05", "WT07"]) + ["KO"] * 4, ["WT"] * 4 + ["KO"] * 4),
    (twice(["Mock1", "Mock8"]) + ["Infected"] * 4, ["Mock"] * 4 + ["Infected"] * 4),
    (twice(["DMSO3", "DMSO5"]) + ["Drug"] * 4, ["DMSO"] * 4 + ["Drug"] * 4),
    (twice(["Sample7", "Sample9"]) + ["Blank"] * 4, ["Sample"] * 4 + ["Blank"] * 4),
    (twice(["CtrlS12", "CtrlS40"]) + ["LPS"] * 4, ["CtrlS"] * 4 + ["LPS"] * 4),
    (twice(["IP12", "IP15"]) + ["IgG"] * 4, ["IP"] * 4 + ["IgG"] * 4),
    (twice(["ATG5", "ATG7"]) + ["WT"] * 4, twice(["ATG5", "ATG7"]) + ["WT"] * 4),
    (["HSP90"] + twice(["HSP701", "HSP702", "HSP703", "HSP704"]), ["HSP90"] + ["HSP70"] * 8),
    # r9's boundaries
    (twice(["WT2", "WT5"]) + ["KO"] * 4, ["WT"] * 4 + ["KO"] * 4),
    (twice(["WT2", "WT7"]) + ["KO"] * 4, ["WT"] * 4 + ["KO"] * 4),
    (twice(["WT11", "WT13"]) + ["KO"] * 4, ["WT"] * 4 + ["KO"] * 4),
    (twice(["WT_replicate5", "WT_replicate7"]) + ["KO"] * 4, ["WT"] * 4 + ["KO"] * 4),
    (twice(["WT.S5", "WT.S7"]) + ["KO"] * 4, ["WT"] * 4 + ["KO"] * 4),
]
STEMS = ("WT", "KO", "Ctrl", "Control", "Treated", "ATG5", "ATG7", "HSP70", "HSP90", "CD45",
         "IgG1", "Cis", "Brain", "Mutant", "Mutant-1", "Mutant-2", "T0", "T24", "Day", "A549",
         "Kv2.1", "Strep", "Prep", "Old - IgG", "Bait_IP", "sample", "WT1", "WT2", "Dox(+)")
LEVELS = (("L", "R"), ("R", "S"), ("A", "B"))
STYLES = ("{b}", "{b}_{n}", "{b}-{n}", "{b} {n}", "{b}{n}", "{b}_rep{n}", "{b} rep {n}",
          "{b}_R{n}", "{b}#{n}", "{b}.{n}", "{b} - Rep {n}", "{b}-rep{n}", "{b}_Rep{n}",
          "{b} replicate {n}", "{b}_REP{n}", "{b} \u2013 {n}", "{b}rep{n}", "{b}Rep{n}",
          "{b}_S{n}", "{b}_M{n}")
NUMBERING = ((1, 2, 3), (0, 1, 2), (2, 3, 4), (1, 2, 4), (1, 2), (3, 4, 5), (1, 5), (4, 5, 6))


def corpus():
    """[(names, analysed)]: REVIEWS, then N_GENERATED generated sheets (fixed seed, no repeats)."""
    rnd = random.Random(SEED)
    sheets, seen = [], set()
    for names, analysed in REVIEWS:
        sheets.append((list(names), list(analysed)))
        seen.add(tuple(names))
    while len(sheets) < len(REVIEWS) + N_GENERATED:
        if rnd.random() < 0.35:
            stem, sep = rnd.choice(STEMS), rnd.choice(("_", "-", " ", ""))
            a, b = rnd.choice(LEVELS)
            conds = [f"{stem}{sep}{a}", f"{stem}{sep}{b}"] + (
                [rnd.choice(STEMS)] if rnd.random() < 0.4 else [])
        else:
            conds = rnd.sample(STEMS, rnd.choice((2, 2, 3)))
        conds = list(dict.fromkeys(conds))
        mixed, style = rnd.random() < 0.3, rnd.choice(STYLES)
        # how the samples are numbered: per condition, on across the conditions, or by an ID
        scheme, last = rnd.random(), 0
        names, analysed = [], []
        for c in conds:
            st, nums = (rnd.choice(STYLES) if mixed else style), rnd.choice(NUMBERING)
            if scheme < 0.2:
                nums, last = tuple(range(last + 1, last + 1 + len(nums))), last + len(nums)
            elif scheme < 0.35:
                nums = tuple(sorted(rnd.sample(range(1, 60), len(nums))))
            reps = [c] * len(nums) if st == "{b}" else [st.format(b=c, n=n) for n in nums]
            names += reps
            analysed += [c] * len(reps)
        r = rnd.random()
        if r < 0.15:                                   # a slip: one name twice
            i = rnd.randrange(len(names))
            names.insert(i, names[i])
            analysed.insert(i, analysed[i])
        elif r < 0.3:                                  # one name in another case
            i = rnd.randrange(len(names))
            names[i] = (names[i].lower() if rnd.random() < 0.5 else names[i].upper())
        elif r < 0.5:                                  # every name twice (technical duplicates)
            names = [x for x in names for _ in (0, 1)]
            analysed = [x for x in analysed for _ in (0, 1)]
        if tuple(names) not in seen:
            seen.add(tuple(names))
            sheets.append((names, analysed))
    return sheets


def outcome(sr, names, analysed, tmp):
    """(the sheet's groups in sample order, "every sample has its own condition" note or not,
    the "conditions differ" note or None, the staff caution or None) from submission_report
    module `sr`, end to end: a
    session whose conditions.csv holds `analysed`, each run assigned to its sample."""
    with open(RECORD, encoding="utf-8") as fh:
        r = json.load(fh)
    r["submission_data"]["samples"] = [
        {"unique_id": f"SMP{100 + i}", "sample_name": f"sample {i}", "condition_name": c}
        for i, c in enumerate(names, 1)]
    rec = sr.sanitize(r)
    groups = sr.sheet_groups(rec)
    s = os.path.join(tmp, "2026-10-20_parity")
    os.makedirs(os.path.join(s, "input"))
    os.makedirs(os.path.join(s, "output"))
    runs = {x["unique_id"]: f"10202026_DIA_{x['unique_id']}_1_{i}"
            for i, x in enumerate(rec["samples"], 1)}
    with open(os.path.join(s, "input", "raw_files.txt"), "w") as fh:
        fh.write("# raw\n" + "".join(f"/data/OCT26/{x}.d\n" for x in runs.values()))
    with open(os.path.join(s, "input", "conditions.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["File.Name", "Group"])
        for x, g in zip(rec["samples"], analysed):
            w.writerow([runs[x["unique_id"]], g])
    with open(os.path.join(s, "input", "locate_decisions.json"), "w") as fh:
        json.dump({"assignments": [{"sample": u, "file": f"{x}.d"} for u, x in runs.items()]}, fh)
    sr.attach(s, rec)
    notes = {n["id"]: n["text"] for n in sr.quality_notes(sr.load(s), s)}
    cautions = getattr(sr, "sheet_cautions", lambda r, sess: [])(sr.load(s), s)  # not in 2.13.2
    return ([groups[x["unique_id"]] for x in rec["samples"]], "sheet_conditions_unique" in notes,
            notes.get("conditions_differ"), " ".join(c["text"] for c in cautions) or None)


def _key(s):
    return " ".join(str(s).split()).casefold()


def glued_groups(names, groups):
    """The sheet groups (as 2.13.2 labels them) made by taking a number glued straight to the
    name's letters off it: no separator, ')', '#', rep / replicate or single ID letter after a
    separator before the number (the rule, as specified)."""
    out = set()
    for c, g in zip(names, groups):
        c = c.strip()
        if c != g and c.startswith(g) and c[len(g):].isdigit() and g[-1:].isalpha():
            word = re.search(r"[A-Za-z]*$", g).group(0)
            before = g[:len(g) - len(word)]
            if not re.search(r"rep(?:licate)?$", word, re.I) and not (
                    len(word) == 1 and before[-1:] in (" ", "_", ".", "-")):
                out.add(g)
    return out


def bits(note):
    """The parts of a "Conditions analysed vs the sample sheet: a; b; c." note."""
    if not note:
        return []
    return note.split(": ", 1)[1].rstrip(".").split("; ")


def names_a_glued_group(bit, glued):
    """Does this part of the note name a group made from glued numbers?"""
    if bit.startswith("the sheet's “"):
        return bit.split("“", 1)[1].split("”", 1)[0] in glued
    if " pools the sheet's " in bit:
        listed = re.sub(r" and \d+ more$", "", bit.split(" pools the sheet's ", 1)[1])
        return any(x in glued for x in listed.split(", "))
    return False


def write_fixture(scripts):
    """2.13.2's answers, from its own scripts directory `scripts` -> FIXTURE."""
    sys.path.insert(0, os.path.abspath(scripts))
    sr = importlib.import_module("submission_report")
    assert os.path.dirname(os.path.abspath(sr.__file__)) == os.path.abspath(scripts), sr.__file__
    rows = []
    for names, analysed in corpus():
        with tempfile.TemporaryDirectory() as tmp:
            g, unique, differ, _caution = outcome(sr, names, analysed, tmp)
        rows.append(json.dumps({"names": names, "analysed": analysed, "groups": g,
                                "unique": unique, "differ": differ}, ensure_ascii=False))
    with open(FIXTURE, "w", encoding="utf-8") as fh:
        fh.write('{"made_by": "skill 2.13.2 (dd932cc) submission_report.py", "seed": %d,\n'
                 ' "sheets": [\n%s\n]}\n' % (SEED, ",\n".join(rows)))
    return len(rows)


class SameAsTwoThirteenTwo(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, SCRIPTS)
        cls.sr = importlib.import_module("submission_report")
        with open(FIXTURE, encoding="utf-8") as fh:
            cls.base = json.load(fh)["sheets"]

    def test_the_fixture_is_this_corpus(self):
        self.assertGreater(len(self.base), 400)
        self.assertEqual([(b["names"], b["analysed"]) for b in self.base], corpus())

    def test_every_sheet_reads_as_2_13_2_and_glued_parts_go_to_staff(self):
        withheld = 0
        for b in self.base:
            names = b["names"]
            with tempfile.TemporaryDirectory() as tmp:
                groups, unique, differ, caution = outcome(self.sr, names, b["analysed"], tmp)
            self.assertEqual((groups, unique), (b["groups"], b["unique"]), names)
            glued = glued_groups(names, b["groups"])
            to_staff = [x for x in bits(b["differ"]) if names_a_glued_group(x, glued)]
            self.assertEqual(bits(differ), [x for x in bits(b["differ"]) if x not in to_staff],
                             names)
            if to_staff:
                withheld += 1
                for x in to_staff:
                    self.assertIn(x, caution or "", names)
            else:
                self.assertIsNone(caution, names)
        self.assertGreater(withheld, 0)


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--write-fixture":
        print(f"{write_fixture(sys.argv[2])} sheets -> {FIXTURE}")
    else:
        unittest.main()
