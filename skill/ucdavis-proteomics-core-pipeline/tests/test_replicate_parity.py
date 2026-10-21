#!/usr/bin/env python3
"""The sample-sheet check reads every sheet as skill 2.13.2 did, except the one case 2.14 changed.

submission_report's sheet check takes a trailing replicate index off the sheet's condition names
to find its groups, and compares them with the groups the analysis used ("conditions differ").
2.14 added one guard (sheet_groups, _own_conditions_apart): two names that several samples carry
each are never merged into one base (ATG5 x3 and ATG7 x3 are not "ATG"). Three reviews of a
broader rule kept finding sheets that 2.13.2 read right and the new rule read wrong, so this test
holds today's code to 2.13.2's on a corpus of sheet shapes: the reviews' examples, then stems x
separators x numbering x levels (L/R, R/S, A/B) x mixed formats x slips x case x duplicated
names, generated with a fixed seed. The groups the analysis used are each sample's true
condition, known by construction -- never computed with the code under test.

fixtures/replicate_parity_2_13_2.json holds each sheet with 2.13.2's sheet groups, its "every
sample has its own condition" note and its "conditions differ" note, made by 2.13.2's own code
(dd932cc), so the test needs no git. The only difference allowed: a sheet where 2.13.2's groups
merged two names that several samples carry each -- there today's groups may only split 2.13.2's,
and those names must be apart.

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
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
FIXTURE = os.path.join(HERE, "fixtures", "replicate_parity_2_13_2.json")
RECORD = os.path.join(HERE, "fixtures", "coreomics_prot0756", "submission.json")
SEED, N_GENERATED = 20261020, 420

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
]
STEMS = ("WT", "KO", "Ctrl", "Control", "Treated", "ATG5", "ATG7", "HSP70", "HSP90", "CD45",
         "IgG1", "Cis", "Brain", "Mutant", "Mutant-1", "Mutant-2", "T0", "T24", "Day", "A549",
         "Kv2.1", "Strep", "Prep", "Old - IgG", "Bait_IP", "sample", "WT1", "WT2")
LEVELS = (("L", "R"), ("R", "S"), ("A", "B"))
STYLES = ("{b}", "{b}_{n}", "{b}-{n}", "{b} {n}", "{b}{n}", "{b}_rep{n}", "{b} rep {n}",
          "{b}_R{n}", "{b}#{n}", "{b}.{n}", "{b} - Rep {n}", "{b}-rep{n}")
NUMBERING = ((1, 2, 3), (0, 1, 2), (2, 3, 4), (1, 2, 4), (1, 2))


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
        names, analysed = [], []
        for c in conds:
            st, nums = (rnd.choice(STYLES) if mixed else style), rnd.choice(NUMBERING)
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
            names[i] = names[i].lower()
        elif r < 0.4:                                  # every name twice (technical duplicates)
            names = [x for x in names for _ in (0, 1)]
            analysed = [x for x in analysed for _ in (0, 1)]
        if tuple(names) not in seen:
            seen.add(tuple(names))
            sheets.append((names, analysed))
    return sheets


def outcome(sr, names, analysed, tmp):
    """(the sheet's groups in sample order, "every sample has its own condition" note or not,
    the "conditions differ" note or None) from submission_report module `sr`, end to end: a
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
    return ([groups[x["unique_id"]] for x in rec["samples"]], "sheet_conditions_unique" in notes,
            notes.get("conditions_differ"))


def _key(s):
    return " ".join(str(s).split()).casefold()


def merged_shared_names(names, groups):
    """The names several samples carry that `groups` puts in one group with another such name."""
    counts = {}
    for c in names:
        counts[_key(c)] = counts.get(_key(c), 0) + 1
    by_group = {}
    for c, g in zip(names, groups):
        if counts[_key(c)] > 1:
            by_group.setdefault(_key(g), set()).add(_key(c))
    return {c for cs in by_group.values() if len(cs) > 1 for c in cs}


def write_fixture(scripts):
    """2.13.2's answers, from its own scripts directory `scripts` -> FIXTURE."""
    sys.path.insert(0, os.path.abspath(scripts))
    sr = importlib.import_module("submission_report")
    assert os.path.dirname(os.path.abspath(sr.__file__)) == os.path.abspath(scripts), sr.__file__
    rows = []
    for names, analysed in corpus():
        with tempfile.TemporaryDirectory() as tmp:
            g, unique, differ = outcome(sr, names, analysed, tmp)
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

    def test_every_sheet_reads_as_2_13_2_but_the_guards_case(self):
        guarded = 0
        for b in self.base:
            names = b["names"]
            with tempfile.TemporaryDirectory() as tmp:
                groups, unique, differ = outcome(self.sr, names, b["analysed"], tmp)
            merged = merged_shared_names(names, b["groups"])
            if not merged:
                self.assertEqual((groups, unique, differ), (b["groups"], b["unique"], b["differ"]),
                                 names)
                continue
            guarded += 1
            # the guard's case: today's groups only split 2.13.2's, and keep those names apart
            within = {}
            for g, old in zip(groups, b["groups"]):
                within.setdefault(_key(g), set()).add(_key(old))
            self.assertTrue(all(len(v) == 1 for v in within.values()), names)
            self.assertEqual(merged_shared_names(names, groups), set(), names)
            self.assertEqual(unique, b["unique"], names)
        self.assertGreater(guarded, 0)


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--write-fixture":
        print(f"{write_fixture(sys.argv[2])} sheets -> {FIXTURE}")
    else:
        unittest.main()
