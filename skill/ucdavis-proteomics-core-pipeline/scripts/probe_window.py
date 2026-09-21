#!/usr/bin/env python3
"""probe_window.py -- measure the scan-window radius DIA-NN infers for an acquisition.

WHY THIS EXISTS
---------------
The 5-step parallel chain reuses .quant files across steps. DIA-NN warns:

    WARNING: combining reuse of .quant files with automatic optimisation of mass
    accuracies or scan window will lead to results that are different from those of
    the original analysis that produced the .quant files and is strongly not
    recommended

Mass accuracy is pinned by estimate_params.py -- except an Orbitrap level with no documented
DIA-NN value, which is measured here too (MASS ACCURACY, below). The scan window was NOT: the
flag was omitted, so DIA-NN optimised it PER FILE. On a real 18-file poplar run that produced a
radius of 7 for seventeen files and 8 for one, which the chain then stitched together.

The radius is a property of the ACQUISITION SCHEME (cycle time vs chromatographic peak
width), not of the individual sample -- but it must be MEASURED, not guessed, because it
depends on the gradient and method.

WHICH RUNS
----------
Not simply the first one. DIA-NN's README ("Changing default settings") says its automatic
optimisation "is inherently noisy: even replicate injections may not produce identical
results, and therefore the analysis results will depend on which run is first in the list",
and recommends running "on several representative runs". Measured on HIVE (DIA-NN 2.7.0,
2026-09-16): one-file probes of one timsTOF cohort gave 10, 11 or 14 depending on the run,
while each run's own radius was reproducible. So, given the whole cohort (--raw ... /
--raw-list):

  * never probe a run whose data DIA-NN cannot fully read (BRUKER INDEX, below);
  * never probe a run under HALF the size of a typical run OF ITS OWN KIND, where "typical" is
    the median of the LARGER half of those runs (see MIN_FRACTION) -- blanks, washes, failed
    acquisitions. Per kind, because size_bytes is not one quantity: a .d's is the INDEXED bytes
    of its analysis.tdf_bin, a .raw's is the file. One floor across both called every 0.9 GB
    Orbitrap .raw a blank beside 8 GB timsTOF .d and pinned a cohort of two instruments from one
    of them. A cohort holding more than one kind is WARNED about instead: one radius is not valid
    for two acquisition schemes (MIXED COHORTS, below);
  * rank what is left by how much was acquired: when every run is a .d with a readable
    Frames.Time, by acquisition time, with the same half-of-typical floor on it (a short wash the
    size of a real run); otherwise by size;
  * probe the median run first, then the lower and upper quartile runs (when >= 3 remain);
  * a run that logs no radius is REPLACED by the remaining run nearest the median IN RANK
    POSITION, and the radius pinned is the MEDIAN of the radii measured. Every probe is recorded;
  * warn when the median run is under half the largest run left, and when the measured radii
    differ by more than 2. The floor is measured from the cohort itself, so a cohort that is
    MOSTLY blanks or washes drags it down with them and excludes nothing -- nothing here can tell
    a wash from a sample, so it says so rather than pretending.

BRUKER INDEX
------------
A .d is only as long as its analysis.tdf frame index says. 342 .d on HIVE have a truncated
index beside a complete analysis.tdf_bin: the copier left a stale mid-acquisition -wal beside
the finished tdf, and a READ-WRITE sqlite open checkpointed those pages into it. DIA-NN then
reads the run silently short -- in the pilot, 121 cycles and 80 precursors -- and a bytes rule
picked exactly such a run (tdf_bin full size, index 0.7%). So a .d is never probed when:

  * its analysis.tdf header is in WAL mode (bytes 18-19 = 2,2; every truncated file on HIVE is,
    intact ones are 1,1),
  * a non-empty analysis.tdf-wal or analysis.tdf-journal sits beside it (the next read-write
    open by anything would rewrite it), or
  * the last frame block (its TimsId offset + the uint32 block size stored there) ends before
    99.9% of analysis.tdf_bin, or past its end.

The header and side files are checked with plain byte reads, and the tdf is only opened with
`file:<tdf>?mode=ro&immutable=1` (tdf_uri). `mode=ro` alone still reads a stale WAL and leaves
-wal/-shm files beside it. Such runs are named in the job log: they are still SEARCHED by the
chain, and someone should look at them.

There is one last resort. A WAL-mode header is a warning sign, not damage in itself: a run can
carry it with a complete index and nothing stale beside it, and whole cohorts do. When NOTHING
else in the cohort is probeable, those runs -- and only those: complete index, no -wal/-journal,
nothing else wrong -- are probed, loudly. The alternative is that step 1b measures nothing while
steps 2-5 go on to SEARCH the very same files, which leaves them in DependencyNeverSatisfied over
data the chain was about to read anyway. Naming runs with --raw overrides none of these checks;
it only narrows the cohort they are applied to.

MIXED COHORTS
-------------
A radius is a property of the acquisition scheme, so a cohort spanning two instruments has no one
right answer. The probe cannot pick for you: it warns, keeps every kind's size floor measured
against runs of its own kind (so no kind is thrown out as a "blank" for being on another scale),
and pins the median of whatever it measured. Run the chain once per instrument.

THERMO .RAW AND OTHER FILES
---------------------------
Ranked by file size only. Reading a .raw's acquisition length or checking its integrity needs
Thermo's RawFileReader (.NET), the same library DIA-NN loads, so it is not done here. A .raw
that DIA-NN cannot open logs no radius and is replaced by the next representative run, with its
log in window.json. The one failure that is NOT replaced is DIA-NN's missing-.NET error: no
other .raw can succeed in that environment, so probing stops at once and says how to fix it.

PROBES RUN ONE AT A TIME
------------------------
Measured on HIVE (16-CPU allocation, full mouse predicted library, the same three runs each
way): timsTOF .d took 1083 s one after another at 16 threads against 1314 s with three
concurrent probes at 5 threads each; Orbitrap .raw 289 s against 250 s. Neither order wins, but
concurrent probes need three probes' memory at once (~10 GB each on timsTOF) and triple the
simultaneous reads of the raw data, and one at a time is what lets a failed run be replaced.
Each probe is announced on stderr -- the job log -- and DIA-NN's output is copied there as it
arrives (and to the probe's probe.log): a silent 18-minute job log is past watch_run.sh's
15-minute stall rule, whose playbook is to cancel the job.

TIME
----
--timeout is per probe (default 3600 s): a run that hangs is cut there, so the next run still
gets its turn. --budget is for ALL probes together: no probe starts after it, and a running one
is cut at it, so the probes cannot outlive the job's wall clock and the evidence is always
written. The chain's step 1b passes its own wall clock less 10 minutes.

STOPPING DIA-NN
---------------
--diann is not always the engine itself: acquire_tools.sh records `apptainer exec ... <sif>
/diann-*/diann-linux` where there is no native build, and sites wrap binaries in scripts. So
DIA-NN runs in its own session (process group), writes to probe.log rather than a pipe a
grandchild could hold open, and every stop -- radius found (SIGTERM, then SIGKILL after 30 s),
timeout or budget (SIGKILL), SIGTERM/SIGINT/SIGHUP to this script -- goes to the whole group.
Measured on HIVE under apptainer 1.5.3: a 5 s timeout returned after 5.3 s with 0 container
processes left (a probe that signalled only its child returned after 41.2 s and left 2).

MASS ACCURACY (--measure window mass-acc [--ms1-ppm X | --ms2-ppm Y])
--------------------------------------------------------------------
For an Orbitrap with a level outside DIA-NN's README table -- resolution unknown, or outside
30k-240k, e.g. the 15k MS2 both pilot Orbitraps acquire -- estimate_params.py omits
--mass-acc/--mass-acc-ms1 and plans `measure_with_diann`. The measurement is adapted from the
README's item 6 of "Changing default settings" ("run DIA-NN on several representative runs (best
to use any suitable empirical library, as this is the quickest) with Unrelated runs option
checked and review the 'Averaged recommended settings for this experiment' values reported at the
end of the log") -- adapted, not followed: it uses the predicted library, one DIA-NN per run
instead of one "Unrelated runs" search, and pins the median of what each run printed instead of
DIA-NN's averaged line. The same DIA-NN run per probe logs the values, after the radius and
before the search proper (DIA-NN 2.7.0 on HIVE, srun job 23522741, Exploris 480 120k/15k, full
mouse predicted library, 32 threads):

    DIA-NN will automatically optimise the mass accuracy for the first run of the experiment, ...
    [1:34] Scan window radius set to 7
    [1:35] Recommended MS1 mass accuracy setting: 4.1 ppm
    [2:27] Optimised mass accuracy: 14 ppm
    [2:45] Searching decoys                     (the whole run: 5:13)

So the probe stops at the later of the lines it needs. A probe counts only when its run logged
EVERYTHING asked; one that did not is replaced like a run with no radius, and what it did log is
recorded but never pinned (the radius and the mass accuracy then describe the same runs). The
MEDIAN (high) of each measured level is pinned, as DIA-NN printed it, never rounded (see
pin_mass_acc). A level with a documented value is passed as --ms1-ppm / --ms2-ppm and pinned as
given; its per-run values are still recorded. (Measuring the documented 120k MS1 instead pinned
4.2 ppm, and every DIA-NN pass then warned that it "deviates significantly from the value
recommended (7 ppm) for the Orbitrap resolution of this run (120000)".)

A MEASUREMENT IS NOT AUTOMATICALLY A VALUE. Measuring the wrong thing is still measuring, and
DIA-NN will happily settle on a tolerance for a search pointed at the wrong FASTA or the wrong
species, or for a batch whose calibration is out. The median over three runs does not catch that
-- all three move together -- so what the runs said must also be plausible before it is pinned
for the whole cohort:

  * each measured median must lie inside MASS_ACC_BAND (MS2 3-30 ppm, MS1 1.5-25 ppm: the
    derivation, all of it from DIA-NN's own numbers, is at the constant);
  * the per-run values of a measured level must span no more than MASS_ACC_MAX_SPREAD (50%) of
    that median -- one acquisition method has one tolerance, or none;
  * at least MASS_ACC_MIN_RUNS (2) runs must have logged everything asked; one run is DIA-NN's
    own first-run auto mode, which this probe exists to replace.

Failing any of them is a PROBE FAILURE, handled exactly like a run that logged nothing: nothing
is written to the cfg or massacc.txt, the evidence JSON keeps every per-run value and the reason
under `mass_acc_refused`, and the exit status is non-zero. The band is a plausibility check, not
a recommendation -- it does not decide between a measured 14 ppm and DIA-NN's own 25 ppm for a
15k MS2, and both are inside it.

THE SOP IS A FLOOR (estimate_params.SOP_MASS_ACC, MS2 20 ppm / MS1 7 ppm). A measured level that
passes the checks above is pinned at max(measured, SOP): the measurement is used only where it
is WIDER than the SOP. That is where the probe earns its keep -- an instrument that genuinely
needs a wider window than the SOP is exactly what nothing else would have caught -- while a
measured value TIGHTER than the SOP buys nothing and costs identifications: on the one cohort
benchmarked (references/diann_parallel.md) the measured 14/7 gave 18,476 precursors against
19,592 at the SOP's 20/7, same runs, same library, same FDR.

The floor NEVER rescues a refused measurement: the band and the floor answer different
questions, and a measured 0.4 ppm is a probe failure, not a 20. A level given as
--ms1-ppm/--ms2-ppm from DIA-NN's resolution table is pinned as given and never floored -- DIA-NN
reads the run's resolution itself and warns when a pass deviates from its own tier.

The measurement is never erased by the floor. `mass_acc` records `measured_ms2_ppm` /
`measured_ms1_ppm` (what the runs said, for a documented level too), `pinned_ms2_ppm` /
`pinned_ms1_ppm` (what the search will use), `floored` per level, and `sop_floor`; the job log
says both numbers when a floor applies.

The flags must not set EITHER mass-accuracy flag: DIA-NN 2.7.0 fixes both levels when one is
given ("automatic optimisation will not be performed as at least one of MS1/MS2 mass accuracies
is user-provided"; `--mass-acc-ms1 7` alone -> "Mass accuracy will be fixed to 2e-05 (MS2) and
7e-06 (MS1)", HIVE srun job 23528991). And a pinned run still prints "Recommended MS1 mass
accuracy setting". So values are accepted only from a run that ANNOUNCED automatic optimisation;
a run whose settings block ends without that line is stopped at once, and no other run is tried:
they all get the same flags. The two ways that happens are reported SEPARATELY, because their
fixes are opposite. DIA-NN saying it is fixing the tolerance (FIXED_ACC_RE) means the flags carry
--mass-acc/--mass-acc-ms1: remove them. The settings block simply ending without the announcement
means either that, worded in a way this probe does not know, or -- when the flags carry neither
flag, which the message checks and says -- that the ANNOUNCEMENT's wording changed between DIA-NN
releases: AUTO_ACC_RE and FIXED_ACC_RE were written against 2.7.0's exact words, and a reworded or
reordered settings block would land here with nothing wrong with the flags at all.

USAGE
-----
    probe_window.py --diann <diann cmd> --raw <run> [<run> ...] | --raw-list <file> \\
                    --fasta <f.fasta> --lib <predicted.speclib> [--threads 16] \\
                    [--timeout 3600] [--budget S] [--max-probes 3] [--max-failures 3] \\
                    [--measure window [mass-acc]] [--ms1-ppm X | --ms2-ppm Y] \\
                    [--workdir DIR] [--write-cfg CFG] [-- <DIA-NN flags of the real search>]

The library must already exist -- run this after step 1 of the chain, not before. Everything
after a bare `--` goes to DIA-NN verbatim (the chain passes the cfg flags this way, so the probe
runs under exactly the flags steps 2-5 run). --diann is exec'd without a shell, so the .NET 8
exports Thermo .raw needs must be in THIS script's environment (ensure_dotnet8.sh, next to
this script, prints DOTNET_ROOT). Prints JSON: the pinned radius and/or mass accuracy and every
probe's evidence -- also on failure, with window_radius and mass_acc null and a non-zero exit
status.
"""
import argparse, json, os, re, shlex, signal, sqlite3, statistics, struct, subprocess, sys
import tempfile, time
from urllib.parse import quote

HERE = os.path.dirname(os.path.abspath(__file__))

sys.path.insert(0, HERE)
# The facility's SOP tolerance, defined once, next to DIA-NN's own tables. It is the FLOOR under
# a measured level (MASS ACCURACY, above): never re-typed here, so the floor tracks the SOP.
from estimate_params import SOP_MASS_ACC, SOP_MASS_ACC_SOURCE      # noqa: E402

# DIA-NN prints e.g. "Scan window radius set to 7". Match loosely (case-insensitive,
# tolerant of the leading '[m:ss]' timestamp) but require the integer.
WINDOW_RE = re.compile(r"window\s+radius\s+set\s+to\s+(\d+)", re.I)

# Mass accuracy, when --mass-acc/--mass-acc-ms1 are omitted (DIA-NN's auto mode). Verbatim from
# DIA-NN 2.7.0 on HIVE (srun job 23522741; Exploris 480, 120k/15k, full mouse predicted library):
#
#     [0:36] Calibrating with mass accuracies 25 (MS1), 25 (MS2)
#     [1:34] Scan window radius set to 7
#     [1:35] Recommended MS1 mass accuracy setting: 4.1 ppm
#     [2:27] Optimised mass accuracy: 14 ppm
#     [2:45] Searching decoys
#
# "Optimised mass accuracy" is the MS2 tolerance DIA-NN settles on; the MS1 line is its
# recommendation. Both are printed during calibration, after the radius and before the search
# proper, so one DIA-NN per run yields all three and is stopped there.
MS2_ACC_RE = re.compile(r"Optimised\s+mass\s+accuracy:\s*([0-9]*\.?[0-9]+)\s*ppm", re.I)
MS1_ACC_RE = re.compile(r"Recommended\s+MS1\s+mass\s+accuracy\s+setting:\s*([0-9]*\.?[0-9]+)\s*ppm",
                        re.I)
# Printed at startup when the flags FIX mass accuracy; DIA-NN then optimises nothing, so a probe
# asked to measure it would wait for lines that never come. Seen with --mass-acc 20
# --mass-acc-ms1 7: "Mass accuracy will be fixed to 2e-05 (MS2) and 7e-06 (MS1)"; with ONE of the
# two flags, DIA-NN 2.7.0 first warns "note the mass accuracy settings used by DIA-NN, automatic
# optimisation will not be performed as at least one of MS1/MS2 mass accuracies is user-provided".
FIXED_ACC_RE = re.compile(r"Mass\s+accuracy\s+will\s+be\s+fixed\s+to|"
                          r"automatic\s+optimisation\s+will\s+not\s+be\s+performed", re.I)
# ...but a list of what means "fixed" is not enough: a pinned run STILL prints "Recommended MS1
# mass accuracy setting" (HIVE, --mass-acc 20 --mass-acc-ms1 7: "[1:50] Recommended MS1 mass
# accuracy setting: 4.3 ppm"), so a reworded notice would let a pinned run's recommendation pass
# as a measurement. Mass accuracy is therefore read only from a run that announced automatic
# optimisation among its settings -- "DIA-NN will automatically optimise the mass accuracy for
# the first run of the experiment" -- and the settings block ends at "N files will be processed".
AUTO_ACC_RE = re.compile(r"automatically\s+optimise\s+the\s+mass\s+accuracy", re.I)
SETTINGS_END_RE = re.compile(r"^\s*\d+\s+files?\s+will\s+be\s+processed", re.I)

# PLAUSIBLE MASS ACCURACY, ppm, per level. A measured median outside its band is NOT pinned:
# pin_mass_acc() rejects it and the probe fails like a run that logged nothing -- loudly, with
# the evidence -- rather than writing a number into massacc.txt that every pass would search at.
#
# Every number below is DIA-NN's own, none is ours:
#   * its README's Orbitrap resolution table spans 4 ppm (240k) to 15 ppm (30k), MS1 and MS2;
#   * it begins every analysis at 25 ppm on both levels -- "Calibrating with mass accuracies
#     25 (MS1), 25 (MS2)" (DIA-NN 2.7.0, HIVE srun job 23522741);
#   * its runtime value for a 15,000-resolution MS2 -- the very case this probe exists for -- is
#     also 25 ppm ("deviates significantly from the value recommended (25 ppm) for the Orbitrap
#     resolution of this run (15000)").
# So the MS2 ceiling is DIA-NN's own widest number plus 20% (25 -> 30), and the floor sits below
# its tightest documented tier (240k -> 4). MS1 is never the wider level, so its ceiling is that
# calibration window itself and its floor is set below half the tightest tier (measured MS1 on
# HIVE: 4.1-4.3 ppm at 120k, where the table says 7).
#
# This is a PLAUSIBILITY band, not a recommendation: 14 ppm and 25 ppm both sit inside it, and
# the band deliberately does not settle which of them a 15k MS2 should use. It catches a number
# that cannot be an Orbitrap mass accuracy at all -- the wrong FASTA, the wrong species, a
# systematically miscalibrated batch -- which the median over three runs cannot, because all
# three move together.
MASS_ACC_BAND = {"ms2_ppm": (3.0, 30.0), "ms1_ppm": (1.5, 25.0)}

# ...and the probed runs must agree to within this fraction of the median, measured as
# (max - min) / median. The radius and the mass accuracy are properties of the ACQUISITION
# METHOD, so one number describes the cohort or no number does.
#
# The threshold is set from the widest REAL disagreement this branch has measured, not from a
# round number. On HIVE (DIA-NN 2.7.0, three Exploris 480 runs of one cohort, one library) the
# same three runs gave MS2 14, 17, 14 ppm with the window left to DIA-NN -- 21% of the median --
# and 12, 17, 14 with `--window 7` pinned, which is 36%. That second set is not noise: it
# reproduced exactly in two separate sruns (23524590, 23526415), pinning the window simply
# changes the conditions DIA-NN optimises under, and probing under a pinned window is a
# supported mode. A 35% rule would refuse the cohort this branch was validated on, so the rule
# sits at 50%: wide enough for every per-run disagreement seen from correctly-configured runs,
# tight enough that one run at 30 ppm beside two at 14 (114%) cannot reach a cfg. The BAND above
# is the primary check; this one catches runs that do not describe one method at all.
MASS_ACC_MAX_SPREAD = 0.50

# Mass accuracy is pinned from the median of at least this many runs. One run is exactly DIA-NN's
# own first-run auto mode, "use this mode for preliminary analyses only", which measuring
# representative runs exists to replace -- and a lone value agrees with itself vacuously.
MASS_ACC_MIN_RUNS = 2

# The level each measured key belongs to, and the DIA-NN log line it comes from -- for messages.
LEVEL_OF = {"ms2_ppm": "MS2", "ms1_ppm": "MS1"}
LOG_LINE_OF = {"ms2_ppm": "Optimised mass accuracy",
               "ms1_ppm": "Recommended MS1 mass accuracy setting"}

# What a probe can measure, and the values each yields per run.
MEASURES = {"window": ("radius",), "mass-acc": ("ms2_ppm", "ms1_ppm")}

# Printed by DIA-NN when it has no usable .NET 8 runtime -- and it still exits 0. Measured on
# HIVE with DIA-NN 2.7.0 (srun job 23509324). Not a property of the run: nothing else will read.
DOTNET_RE = re.compile(r"cannot read \.raw files", re.I)

# A run under this fraction of a typical run OF ITS OWN KIND (see _size_scale) is treated as a
# blank, wash or failed injection and never probed. It only decides which runs are PROBED; every
# run is still searched.
#
# "Typical" is the median of the LARGER half of the runs, not of all of them: with a blank
# injected after every sample the blanks are the majority, the overall median is a blank's
# size, and a floor measured from it excludes nothing (review reproduction: 6 samples and 7
# blanks -> probes blank3, blank6 and s2). Measured sizes on HIVE, raw_data/Lumos1/noi25:
# washes 120-240 MB and 60-min washes 527-567 MB, beside 80 90-min DIA runs of 0.93-2.85 GB
# (median 2.12 GB; median of the larger half 2.42 GB). The 60-min washes are 25-27% of that DIA
# median, so a 25% floor let them in. At 50% of the larger half's median they are out; 2 of the
# 80 DIA runs (38-40%) are not probed. They are still searched. The rule keeps blanks out while
# they are at most 60-75% of the list; on a folder that is 86% washes it fails, so hand the
# search the cohort, not the folder.
MIN_FRACTION = 0.5

# A .d whose last indexed frame block ends before this fraction of analysis.tdf_bin is treated
# as a truncated index (the pilot's was 0.7%).
INDEX_COVERAGE_MIN = 0.999

SELECTION_RULE = (
    "never probe a Bruker .d whose analysis.tdf is in WAL mode, has a non-empty -wal/-journal "
    "beside it, or whose frame index ends before 99.9% of analysis.tdf_bin -- unless the WAL "
    "header is the only objection and no other run in the cohort can be probed at all, when "
    "those runs are probed with a WARNING; drop runs under "
    f"{MIN_FRACTION:.0%} of the median size of the larger half of the runs OF THEIR OWN KIND "
    "(directories, .d among them, by their bytes -- for a TDF .d the indexed bytes of "
    "analysis.tdf_bin; files by their own bytes), and warn when a cohort holds more than one "
    "kind, because one radius is not valid for two acquisition schemes; when every run left is "
    "a .d with a readable Frames.Time, rank by acquisition time and drop runs under "
    f"{MIN_FRACTION:.0%} of its larger-half median too, "
    "otherwise rank by size; probe the median run, then the lower and upper quartile runs; "
    "replace a run that logs no radius (or, when mass accuracy is measured too, not all of what "
    "it measures) with the remaining run nearest the median in rank position; pin the median of "
    "what the runs that logged everything measured")

POLL_S = 0.2          # how often DIA-NN's log is read; nothing against minutes of calibration
STOP_GRACE_S = 30     # SIGTERM -> SIGKILL, once a radius is in hand
_HARD = getattr(signal, "SIGKILL", signal.SIGTERM)

_LIVE = set()    # running DIA-NN processes, so a signal to the probe can stop them
_STOP = []       # signals received; once set, no probe starts and running ones are stopped


# ------------------------------------------------------------------------------------------
# Bruker .d integrity
# ------------------------------------------------------------------------------------------
def tdf_uri(tdf):
    """The ONLY way this script opens an analysis.tdf. immutable=1: SQLite neither reads a -wal
    nor creates -wal/-shm, and cannot write. `mode=ro` alone still reads (and depends on) a
    stale WAL -- the state that truncated 342 tdfs on HIVE when something opened them
    read-write. Percent-encoded, so a path with a space, `#` or `?` is still that path."""
    return "file:" + quote(os.path.abspath(tdf)) + "?mode=ro&immutable=1"


def bruker_tdf_status(dotd, ignore_wal_header=False):
    """What a Bruker TDF .d's frame index says, and whether it can be trusted.

    Returns None when `dotd` has neither analysis.tdf nor analysis.tdf_bin (not a TDF run).
    Otherwise {time_s, n_frames, indexed_bytes, tdf_bin_bytes, index_coverage, wal_header,
    time_problem, problems}: time_s is max(Frames.Time) in seconds, indexed_bytes where the last
    indexed frame block ends in analysis.tdf_bin, and `problems` is empty only for a run safe to
    probe. The tdf is not opened at all when its header or side files already disqualify it.

    `time_problem` is NOT a problem with the data: a Frames table this reader cannot get a time
    out of (no Time column, or NULL for every frame) is a schema it does not know how to RANK by,
    not a damaged run. Its spectra are all there, it just falls back to ranking by size.

    `ignore_wal_header=True` reads the index even though the header says WAL mode. Only
    select_representative()'s last resort passes it, and only when nothing else is probeable:
    the open is immutable, so it reads the file as it sits on disk and cannot write to it."""
    d = dotd.rstrip("/") or dotd
    tdf, tdf_bin = os.path.join(d, "analysis.tdf"), os.path.join(d, "analysis.tdf_bin")
    has_tdf, has_bin = os.path.isfile(tdf), os.path.isfile(tdf_bin)
    if not (has_tdf or has_bin):
        return None
    st = {"time_s": None, "n_frames": None, "indexed_bytes": None,
          "tdf_bin_bytes": os.path.getsize(tdf_bin) if has_bin else None,
          "index_coverage": None, "wal_header": False, "time_problem": None, "problems": []}
    problems = st["problems"]
    if not has_tdf:
        problems.append("analysis.tdf_bin without analysis.tdf (no frame index to read it by)")
        return st
    if not has_bin:
        problems.append("analysis.tdf without analysis.tdf_bin (no spectra)")

    # Plain bytes, never SQLite: nothing here may touch the file's journal state.
    with open(tdf, "rb") as fh:
        head = fh.read(100)
    if len(head) < 100 or not head.startswith(b"SQLite format 3\x00"):
        problems.append("analysis.tdf has no SQLite header")
        return st
    if head[18] == 2 or head[19] == 2:
        st["wal_header"] = True
        if not ignore_wal_header:
            problems.append(f"analysis.tdf is in WAL mode (header bytes 18-19 = {head[18]},"
                            f"{head[19]}; every truncated tdf found on HIVE is, intact ones "
                            "are 1,1)")
    for side in ("-wal", "-journal"):
        path = tdf + side
        if os.path.isfile(path) and os.path.getsize(path) > 0:
            problems.append(f"non-empty analysis.tdf{side} ({os.path.getsize(path)} bytes) beside "
                            "it -- a read-write open would rewrite (and may truncate) the index")
    if problems:
        return st

    # Two queries, because they answer different questions. MAX(TimsId)/COUNT(*) is the frame
    # INDEX: a Frames table that cannot answer it is a damaged run. MAX(Time) is only how this
    # cohort is RANKED; a Frames table without a Time column (or with NULL in it) is a schema
    # this reader does not know, not damage -- asking for both at once reported one as the other.
    t_max = None
    try:
        con = sqlite3.connect(tdf_uri(tdf), uri=True)
        try:
            tims_max, n = con.execute("SELECT MAX(TimsId), COUNT(*) FROM Frames").fetchone()
            try:
                t_max = con.execute("SELECT MAX(Time) FROM Frames").fetchone()[0]
            except sqlite3.Error as e:
                st["time_problem"] = (f"analysis.tdf Frames has no readable Time column ({e}) -- "
                                      "the spectra are there, but this run cannot be ranked by "
                                      "acquisition time")
        finally:
            con.close()
    except sqlite3.Error as e:
        problems.append(f"analysis.tdf frame index cannot be read ({e})")
        return st
    if not n or tims_max is None:
        problems.append("analysis.tdf indexes no frames")
        return st
    st["n_frames"] = int(n)
    if t_max is None:
        st["time_problem"] = st["time_problem"] or (
            "analysis.tdf Frames.Time is NULL for every frame -- the spectra are there, but this "
            "run cannot be ranked by acquisition time")
    else:
        st["time_s"] = float(t_max)

    size = st["tdf_bin_bytes"]
    with open(tdf_bin, "rb") as fh:
        fh.seek(int(tims_max))
        raw = fh.read(4)
    if len(raw) < 4:
        problems.append(f"the last indexed frame starts at byte {tims_max}, past the end of "
                        f"analysis.tdf_bin ({size} bytes) -- tdf_bin is truncated")
        return st
    end = int(tims_max) + struct.unpack("<I", raw)[0]
    st["indexed_bytes"] = end
    st["index_coverage"] = end / size if size else 0.0
    if end > size:
        problems.append(f"the last indexed frame block ends at byte {end}, past the end of "
                        f"analysis.tdf_bin ({size} bytes) -- tdf_bin is truncated")
    elif end < INDEX_COVERAGE_MIN * size:
        span = "" if t_max is None else f", {t_max / 60:.1f} min"
        problems.append(f"the frame index covers only {st['index_coverage']:.1%} of "
                        f"analysis.tdf_bin ({n} frames{span}) -- a truncated "
                        "index; DIA-NN would read only that part of the run")
    return st


# ------------------------------------------------------------------------------------------
# Which runs
# ------------------------------------------------------------------------------------------
def _dir_bytes(p):
    """Every file under `p`, at any depth. An Agilent .d keeps its data in AcqData/ and a Waters
    .raw directory in _FUNC*.DAT: counting only the top level made every such run tie at ~0 bytes,
    so the ranking fell back to the path tie-break -- alphabetical, i.e. first-file order again."""
    total = 0
    for root, _dirs, files in os.walk(p):
        for f in files:
            fp = os.path.join(root, f)
            if os.path.isfile(fp):          # not a broken symlink
                total += os.path.getsize(fp)
    return total


def measure_run(path, ignore_wal_header=False):
    """{file, format, readable, size_bytes, time_s, time_problem, problems} for one input run.

    A Bruker TDF .d: time_s from its frame index, size_bytes the indexed part of
    analysis.tdf_bin, problems from bruker_tdf_status(). Any other directory: every file under
    it. A file (.raw, .mzML, ...): its size, time_s None. `format` is what those bytes are
    (see _size_scale); `ignore_wal_header` goes to bruker_tdf_status()."""
    p = path.rstrip("/") or path
    m = {"file": path, "format": None, "readable": False, "size_bytes": None, "time_s": None,
         "time_problem": None, "problems": []}
    try:
        if os.path.isdir(p):
            st = bruker_tdf_status(p, ignore_wal_header=ignore_wal_header)
            if st is None:
                m["format"] = "directory"
                m["size_bytes"] = _dir_bytes(p)
            else:
                m["format"] = "bruker_tdf"
                m["time_s"], m["problems"] = st["time_s"], list(st["problems"])
                m["time_problem"] = st["time_problem"]
                m["size_bytes"] = (st["indexed_bytes"] if st["indexed_bytes"] is not None
                                   else st["tdf_bin_bytes"] or 0)
                m["index_coverage"] = st["index_coverage"]
                m["wal_header"] = st["wal_header"]
            m["readable"] = True
        elif os.path.isfile(p):
            m["format"] = os.path.splitext(p)[1].lower() or "file"
            m["size_bytes"], m["readable"] = os.path.getsize(p), True
    except OSError as e:
        m["readable"], m["error"] = False, str(e)
    return m


def _size_scale(m):
    """Which runs a run's size_bytes may be compared with -- its "kind", for the size floor.

    A .raw's size is a FILE's bytes; a .d's is a directory's, and for a TDF .d the INDEXED bytes
    of analysis.tdf_bin. Measuring one floor over both compares quantities that are not the same
    (verified: 3 x 0.9 GB Orbitrap .raw beside 3 x 8 GB timsTOF .d excluded every .raw as "likely
    a blank", flipped the ranking to time because only .d were left, and pinned a two-instrument
    cohort from one instrument).

    All directories share one scale, TDF or not: a FAILED Bruker acquisition is a .d with neither
    analysis.tdf nor analysis.tdf_bin in it, and it is exactly the size floor against the real .d
    beside it that has to drop it (HIVE, garg cohort) -- put it in a partition of its own and it
    survives, has no acquisition time, and drags the whole cohort back to size ranking."""
    return "directory" if m["format"] in ("bruker_tdf", "directory") else m["format"]


_SCALE_NAMES = {".raw": "Thermo .raw", ".mzml": ".mzML", ".wiff": "SCIEX .wiff"}


def _scale_name(scale, runs):
    """What to call a size scale in a message, given the runs on it."""
    if scale == "directory":
        return ("Bruker .d" if any(m["format"] == "bruker_tdf" for m in runs)
                else "directory of files")
    return _SCALE_NAMES.get(scale, scale)


def _mixed_cohort_warning(names):
    """One radius for two acquisition schemes is not a measurement, it is an average of two."""
    joined = names[0] if len(names) == 1 else ", ".join(names[:-1]) + " and " + names[-1]
    return ("cohort mixes " + joined + " -- one scan-window radius is not valid for "
            + ("both" if len(names) == 2 else "all of them")
            + ". The radius is a property of the acquisition scheme, and these runs were not all "
            "acquired on one. Each kind's size floor is measured against runs of its own kind "
            "here, so none of them is dropped as a 'blank' merely for being on another byte "
            "scale -- but the radius pinned comes from whichever runs rank in the middle of the "
            "mixture. Run the chain once per instrument.")


class NoProbeableRun(ValueError):
    """No input run can be probed. `selection` holds what was found, for the evidence JSON."""

    def __init__(self, msg, selection):
        super().__init__(msg)
        self.selection = selection


def _typical(values):
    """Median of the larger half -- see MIN_FRACTION."""
    v = sorted(values)
    return statistics.median(v[len(v) // 2:])


def _wal_header_only_runs(measured, sel):
    """The last resort: runs excluded ONLY because their analysis.tdf header says WAL mode.

    A WAL-mode header is how every one of the 342 truncated tdfs on HIVE presents, so it is
    refused while anything else can be probed. It is not itself damage, though: a run can carry
    it with a complete frame index and nothing stale beside it, and a whole cohort can (6 of 6 in
    the review's reproduction). Refusing such a cohort outright fails step 1b and leaves steps
    2-5 in DependencyNeverSatisfied -- over runs those steps would have SEARCHED regardless. So
    when nothing else is probeable, these are probed and named loudly.

    Each candidate is measured again with that one objection suppressed; the immutable open reads
    the file as it sits on disk and cannot write to it. Everything else -- a short index, a stale
    -wal, an unreadable Frames table -- is still reported, so an empty problems list the second
    time means the header really was the only objection. Mutates `sel`; returns the runs."""
    recovered = []
    for m in measured:
        if not (m["readable"] and m["problems"] and m.get("wal_header")):
            continue
        again = measure_run(m["file"], ignore_wal_header=True)
        if again["readable"] and not again["problems"]:
            recovered.append(again)
            sel["probed_despite_wal_header"].append(
                {"file": m["file"], "problems": m["problems"],
                 "index_coverage": again.get("index_coverage"),
                 "why": "no run in this cohort is probeable and the WAL-mode header is this "
                        "run's only objection: its frame index covers all of analysis.tdf_bin "
                        "and nothing stale sits beside it, so the data is intact. The chain "
                        "searches it either way -- look at it."})
    if recovered:
        named = {x["file"] for x in sel["probed_despite_wal_header"]}
        sel["excluded_damaged"] = [x for x in sel["excluded_damaged"] if x["file"] not in named]
    return recovered


def select_representative(paths, max_probes=3, min_fraction=MIN_FRACTION):
    """Choose the runs to probe. Deterministic and independent of input order (ties are broken
    by path), because "which run is first" is exactly the dependence being removed.

    Returns {chosen, reserves, rank_by, reference_size_bytes, min_size_bytes,
    size_reference_by_kind, mixed_cohort, median_much_smaller, reference_time_s, min_time_s,
    no_acquisition_time, excluded_small, excluded_damaged, probed_despite_wal_header, unreadable,
    n_inputs, n_eligible, rule}. `chosen` is in probing order: median, lower quartile, upper
    quartile (positions n//2, floor(m/4), ceil(3m/4) of the m+1 ranked eligible runs -- distinct
    for every n >= 3; on an even count the upper-middle run, the less likely to be a weak
    injection). `reserves` are the other eligible runs, nearest the median by RANK POSITION
    first -- |i - n//2| in the ranking, not |value - the median value|, which differ wherever the
    ranked values are unevenly spaced -- the larger run on a tie.
    Raises NoProbeableRun (a ValueError) when nothing is eligible."""
    measured = [measure_run(p) for p in paths]
    usable = [m for m in measured if m["readable"] and not m["problems"]]
    sel = {"chosen": [], "reserves": [], "rank_by": None,
           "reference_size_bytes": None, "min_size_bytes": None,
           "size_reference_by_kind": {}, "mixed_cohort": None, "median_much_smaller": None,
           "reference_time_s": None, "min_time_s": None, "no_acquisition_time": [],
           "excluded_small": [],
           "excluded_damaged": [{"file": m["file"], "problems": m["problems"]}
                                for m in measured if m["readable"] and m["problems"]],
           "probed_despite_wal_header": [],
           "unreadable": [m["file"] for m in measured if not m["readable"]],
           "n_inputs": len(paths), "n_eligible": 0, "rule": SELECTION_RULE}
    if not usable:
        usable = _wal_header_only_runs(measured, sel)
    if not usable:
        raise NoProbeableRun(
            f"none of the {len(paths)} input runs can be probed ({len(sel['unreadable'])} "
            f"unreadable, {len(sel['excluded_damaged'])} with a damaged Bruker index), so there "
            "is nothing to measure the scan window on: " + ", ".join(paths[:5])
            + ". --raw does not override these checks -- it narrows the cohort they are applied "
              "to, and a run named there still has to clear them. Repair the runs above (a "
              "truncated index is fixed by copying the .d from the instrument again) or give "
              "step 1b a cohort with at least one intact run in it.", sel)

    def small(m, why):
        sel["excluded_small"].append({"file": m["file"], "size_bytes": m["size_bytes"],
                                      "time_s": m["time_s"], "why": why})

    # The size floor first, over every usable run: it is what removes a failed acquisition with
    # no data at all -- a .d with neither analysis.tdf nor tdf_bin has no acquisition time, and
    # deciding the ranking before it was gone made a whole timsTOF cohort rank by bytes (HIVE,
    # garg cohort). Then, when every run left has a frame index, the time floor and time ranking.
    #
    # Per KIND (see _size_scale): a .d's size_bytes and a .raw's are not the same quantity, and
    # one floor over both threw out every .raw of a two-instrument cohort as a blank.
    kinds = {}
    for m in usable:
        kinds.setdefault(_size_scale(m), []).append(m)
    refs = {k: _typical([m["size_bytes"] for m in g]) for k, g in kinds.items()}
    sel["size_reference_by_kind"] = {
        _scale_name(k, g): {"n_runs": len(g), "reference_size_bytes": refs[k],
                            "min_size_bytes": min_fraction * refs[k]}
        for k, g in kinds.items()}
    if len(kinds) == 1:
        only = list(kinds)[0]
        sel["reference_size_bytes"] = refs[only]
        sel["min_size_bytes"] = min_fraction * refs[only]
    else:
        sel["mixed_cohort"] = _mixed_cohort_warning(
            sorted(_scale_name(k, g) for k, g in kinds.items()))
    eligible = []
    for m in usable:
        ref_size = refs[_size_scale(m)]
        if m["size_bytes"] < min_fraction * ref_size:
            small(m, f"size {m['size_bytes'] / 1e9:.2f} GB is under {min_fraction:.0%} of a "
                     f"typical {_scale_name(_size_scale(m), kinds[_size_scale(m)])} run's "
                     f"{ref_size / 1e9:.2f} GB")
        else:
            eligible.append(m)
    # A .d whose Frames table this reader cannot get a time out of is not damaged; it just
    # cannot be ranked by acquisition time, which drops the whole cohort back to size ranking.
    sel["no_acquisition_time"] = [{"file": m["file"], "why": m["time_problem"]}
                                  for m in eligible if m.get("time_problem")]
    by_time = bool(eligible) and all(m["time_s"] is not None for m in eligible)
    sel["rank_by"] = "time" if by_time else "size"
    if by_time:
        ref_time = _typical([m["time_s"] for m in eligible])
        sel["reference_time_s"], sel["min_time_s"] = ref_time, min_fraction * ref_time
        timed, eligible = eligible, []
        for m in timed:
            if m["time_s"] < sel["min_time_s"]:
                small(m, f"acquisition time {m['time_s'] / 60:.1f} min is under "
                         f"{min_fraction:.0%} of a typical run's {ref_time / 60:.1f} min")
            else:
                eligible.append(m)
    # Not a way out: the largest run of every kind clears both floors on its own (its size is at
    # least the median of the larger half of its kind, its time at least the median of the larger
    # half of the times left, and min_fraction <= 1), so `eligible` cannot be empty here. Kept as
    # an assert because a min_fraction above 1 -- which no caller passes -- would break that.
    assert eligible, "no run clears the size and acquisition-time floors"

    if by_time:
        eligible.sort(key=lambda m: (m["time_s"], m["size_bytes"], m["file"]))
    else:
        eligible.sort(key=lambda m: (m["size_bytes"], m["file"]))
    n, last = len(eligible), len(eligible) - 1
    mid = n // 2
    # The floor is measured from the cohort itself, so a cohort that is MOSTLY blanks or washes
    # drags it down with them and excludes none of them: 12 washes and 2 samples picks three
    # washes, and 9 one-byte files beside a 10 GB run picks three one-byte files. Nothing here
    # can tell a wash from a sample -- so say what the shape of the cohort looks like.
    biggest = max(eligible, key=lambda m: m["size_bytes"])
    if eligible[mid]["size_bytes"] < 0.5 * biggest["size_bytes"]:
        sel["median_much_smaller"] = (
            f"the median run {_base(eligible[mid]['file'])} "
            f"({eligible[mid]['size_bytes'] / 1e9:.2f} GB) is under half the largest run left "
            f"({_base(biggest['file'])}, {biggest['size_bytes'] / 1e9:.2f} GB). If most of this "
            "cohort is blanks, washes or failed injections, the floor was measured from THEM and "
            "kept them in, and the radius is about to be measured on one. Hand step 1b the "
            "sample runs, not the folder.")
    if max_probes >= 3 and n >= 3:
        picks = [("median", mid), ("lower_quartile", last // 4),
                 ("upper_quartile", -(-3 * last // 4))]
    else:
        picks = [("median", mid)]
    taken = {i for _, i in picks}
    # "Nearest the median" is by RANK POSITION, |i - mid|, not by how near the value is to the
    # median value; on a ranked list the two agree except where the values are spaced unevenly.
    rest = sorted((i for i in range(n) if i not in taken), key=lambda i: (abs(i - mid), -i))

    def entry(m, role):
        return {"file": m["file"], "size_bytes": m["size_bytes"], "time_s": m["time_s"],
                "role": role}

    sel["chosen"] = [entry(eligible[i], role) for role, i in picks]
    sel["reserves"] = [entry(eligible[i], "reserve") for i in rest]
    sel["n_eligible"] = n
    return sel


# ------------------------------------------------------------------------------------------
# Running DIA-NN
# ------------------------------------------------------------------------------------------
def _signal_group(p, sig):
    """Send `sig` to DIA-NN and everything it started (its process group). Where there are no
    process groups (Windows) only the direct child can be reached."""
    try:
        if hasattr(os, "killpg"):
            os.killpg(p.pid, sig)
        elif p.poll() is None:
            p.kill()
    except OSError:          # the whole group is already gone
        pass


def _group_alive(pgid):
    try:
        os.killpg(pgid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


def _end_group(p, grace=STOP_GRACE_S):
    """SIGTERM DIA-NN's whole group, then SIGKILL whatever is left after `grace` s. A wrapper
    that forks (a bash script without `exec`; `apptainer exec`) leaves the engine as a
    grandchild, which signalling only the direct child never reaches. A process that moves
    itself into a NEW session escapes this."""
    _signal_group(p, signal.SIGTERM)
    if not hasattr(os, "killpg"):
        try:
            p.wait(timeout=grace)
        except subprocess.TimeoutExpired:
            p.kill()
        return
    deadline = time.time() + grace
    while time.time() < deadline:
        p.poll()                           # reap our child so it stops counting as a member
        if not _group_alive(p.pid):
            return
        time.sleep(0.1)
    _signal_group(p, _HARD)
    p.wait()


def _kill_now(p):
    """Timeout, budget or a signal: the whole group, at once."""
    _signal_group(p, _HARD)
    try:
        p.wait(timeout=STOP_GRACE_S)
    except subprocess.TimeoutExpired:
        pass


def _on_signal(signum, _frame):
    """SIGTERM (scancel, a SLURM time limit), SIGINT, SIGHUP. DIA-NN is in its own process group,
    so a signal to the job's group does not reach it; stop it here. The probes then return and
    main() exits with 128 + signum."""
    _STOP.append(signum)
    for p in list(_LIVE):
        _signal_group(p, _HARD)


def _say(text):
    """One line into the job log (stderr), now -- not when the job ends."""
    sys.stderr.write(text + "\n")
    sys.stderr.flush()


def read_log_line(line):
    """{radius | ms2_ppm | ms1_ppm: value} for whatever one DIA-NN log line reports.
    Non-positive values are ignored: DIA-NN rejects a 0 radius, and 0 ppm matches nothing."""
    out = {}
    m = WINDOW_RE.search(line)
    if m and int(m.group(1)) > 0:
        out["radius"] = int(m.group(1))
    for key, rx in (("ms2_ppm", MS2_ACC_RE), ("ms1_ppm", MS1_ACC_RE)):
        m = rx.search(line)
        if m and float(m.group(1)) > 0:
            out[key] = float(m.group(1))
    return out


class LogReader:
    """What one DIA-NN run's log says about the values in `measure`, fed a line at a time.

    `found` maps radius / ms2_ppm / ms1_ppm to the first value DIA-NN logged. When mass accuracy
    was asked for and the run is plainly not optimising it, `found` holds ONE of two flags --
    they are different faults with different fixes, so they are never merged:

      * `fixed_mass_acc`    -- DIA-NN SAID it is fixing the tolerance (FIXED_ACC_RE). The
                               DIA-NN flags carry --mass-acc/--mass-acc-ms1; drop them.
      * `no_auto_announcement` -- the settings block ended (SETTINGS_END_RE) without
                               "DIA-NN will automatically optimise the mass accuracy". Values are
                               read only from a run that announced it, because a run with a FIXED
                               tolerance still prints "Recommended MS1 mass accuracy setting".
                               The flags may fix it under wording this probe does not know, OR
                               the announcement itself may be worded differently in this DIA-NN
                               release -- AUTO_ACC_RE and FIXED_ACC_RE were written against
                               2.7.0. main() reports which by looking at the flags it passed.

    `skip` names levels given as documented: still recorded when logged, never waited for."""

    def __init__(self, measure, skip=()):
        self.measure = list(measure)
        self.wanted = [k for m in self.measure for k in MEASURES[m] if k not in skip]
        self.found = {}
        self.auto = False

    def feed(self, line):
        for k, v in read_log_line(line).items():
            if k == "radius":
                # Only when asked for. With --window N pinned, DIA-NN 2.7.0 echoes "Scan window
                # radius set to N" among its startup settings: the cfg's value, not a measurement.
                if "window" in self.measure:
                    self.found.setdefault(k, v)
            elif "mass-acc" in self.measure and self.auto:
                self.found.setdefault(k, v)       # the first value DIA-NN settles on
        if "mass-acc" in self.measure and not self.auto and not self.not_optimising():
            if AUTO_ACC_RE.search(line):
                self.auto = True
            elif FIXED_ACC_RE.search(line):
                self.found["fixed_mass_acc"] = True
            elif SETTINGS_END_RE.search(line):
                self.found["no_auto_announcement"] = True

    def not_optimising(self):
        """This run will never log a measured mass accuracy -- for either reason."""
        return bool(self.found.get("fixed_mass_acc")
                    or self.found.get("no_auto_announcement"))

    def missing(self):
        return [k for k in self.wanted if k not in self.found]

    def done(self):
        return not self.missing() or self.not_optimising()


def ppm_text(v):
    """14.0 -> '14', 4.1 -> '4.1': DIA-NN's own spelling, so the pinned flag reads as printed."""
    return f"{v:g}"


def band_text(key):
    """'3-30 ppm' -- the plausible band for one level, for a message."""
    lo, hi = MASS_ACC_BAND[key]
    return f"{lo:g}-{hi:g} ppm"


def describe_missing(keys):
    """Human names for the values a run did not log."""
    names = {"radius": "a scan-window radius",
             "ms2_ppm": "an 'Optimised mass accuracy' (MS2)",
             "ms1_ppm": "a 'Recommended MS1 mass accuracy setting'"}
    return " or ".join(names[k] for k in keys) or "anything"


def pin_mass_acc(probes, documented=None):
    """The mass accuracy to pin from the per-run values: the MEDIAN of each measured level, taken
    high, and the documented value of a level in `documented` ({"ms1_ppm": 7}).

    The median rather than the mean for the same reason the radius uses it: one atypical run of
    three cannot move it. No rounding, so a measured number is always one DIA-NN itself printed
    for one of the runs -- DIA-NN already rounds what it prints (MS2 "14 ppm", MS1 "4.1 ppm"),
    and rounding again, e.g. up to 0.5 ppm, would pin a tolerance no run produced. median_HIGH,
    not low: on an even number of runs (a replacement makes four) the two are different numbers,
    and a tolerance that is too tight loses identifications while one that is too wide only costs
    specificity -- 14 and 20 must pin 20, not 14. A documented level is pinned as documented and
    its per-run values are kept only as evidence.

    `probes` are the runs that logged everything (main() never passes a failed one). Returns None
    unless every probe logged every measured level.

    The result carries `rejected`: the reasons this measurement must NOT be pinned, empty when
    there are none. A measurement is refused when a measured median falls outside MASS_ACC_BAND,
    or when the per-run spread of a measured level exceeds MASS_ACC_MAX_SPREAD of that median.
    `pin_as` is then None, so no caller can splice a rejected number onto a DIA-NN command line.
    A `documented` level is not band-checked here: main() checks --ms1-ppm/--ms2-ppm as it parses
    them, before any DIA-NN starts.

    THE MEASURED VALUE AND THE PINNED VALUE CAN DIFFER, so the result records both. A measured
    level that survives the checks is floored at the facility's SOP tolerance (SOP_MASS_ACC): the
    pinned value is max(measured, SOP). `measured_ms2_ppm` / `measured_ms1_ppm` are always what
    the runs said -- for a documented level too -- `ms2_ppm` / `ms1_ppm` are what the search will
    use, and `floored` says per level which of the two `pin_as` carries."""
    documented = documented or {}
    if not probes:
        return None
    out = {"sources": {}, "documented": dict(documented), "rejected": [], "spread": {},
           "sop_floor": dict(SOP_MASS_ACC), "sop_floor_source": SOP_MASS_ACC_SOURCE,
           "floored": {}}
    for key in ("ms2_ppm", "ms1_ppm"):
        level, line = LEVEL_OF[key], LOG_LINE_OF[key]
        per_run = [p.get(key) for p in probes]
        out[key.replace("_ppm", "_per_run")] = per_run
        # What the runs measured, kept for EVERY level -- including one pinned from DIA-NN's
        # table or raised to the SOP floor. Neither may erase the measurement: a reader has to
        # be able to see that the instrument measured 14 and the search ran at 20.
        out["measured_" + key] = (statistics.median_high(per_run)
                                  if per_run and None not in per_run else None)
        out["floored"][key] = False
        if key in documented:
            out[key] = float(documented[key])
            measured_note = ("" if out["measured_" + key] is None else
                             f"; the runs measured {ppm_text(out['measured_' + key])} ppm "
                             "(recorded, not used)")
            out["sources"][key] = (f"given, not measured: {ppm_text(out[key])} ppm, passed to the "
                                   f"probe as --{key.replace('_ppm', '-ppm')} (a value from "
                                   "DIA-NN's resolution table: a documented tier, or interpolated "
                                   "between tiers) and pinned as given" + measured_note)
            continue
        if None in per_run:
            return None
        out[key] = out["measured_" + key]
        out["sources"][key] = (f"measured: median (high) of DIA-NN's per-run '{line}' ({level}), "
                               "as printed -- not rounded")
    measured = [k for k in ("ms2_ppm", "ms1_ppm") if k not in documented]
    # A measured value is only ever as good as what DIA-NN was pointed at. The wrong FASTA, the
    # wrong species or a batch calibrated against the wrong lock mass moves EVERY run the same
    # way, so the median absorbs it and `agree` says the runs agreed. Only a magnitude check
    # catches that, and only a spread check catches runs that do not describe one method.
    for key in measured:
        per_run = out[key.replace("_ppm", "_per_run")]
        value, lo, hi = out["measured_" + key], *MASS_ACC_BAND[key]
        printed = ", ".join(ppm_text(v) for v in per_run)
        bad = []
        if not lo <= value <= hi:
            bad.append(f"{LEVEL_OF[key]} {ppm_text(value)} ppm is outside {band_text(key)}, the "
                       f"plausible band for an Orbitrap (per-run: {printed})")
        spread = (max(per_run) - min(per_run)) / value if value else None
        out["spread"][key] = spread
        if spread is not None and spread > MASS_ACC_MAX_SPREAD:
            bad.append(f"the probed runs disagree on {LEVEL_OF[key]}: {printed} ppm spans "
                       f"{spread:.0%} of the median ({ppm_text(value)} ppm), over the "
                       f"{MASS_ACC_MAX_SPREAD:.0%} one acquisition method is allowed")
        out["rejected"] += bad
        if bad:
            # A refused level is NOT floored. The band and the floor answer different questions:
            # the band says this cannot be a mass accuracy at all, and the SOP is not a fallback
            # for a measurement that made no sense -- 0.4 ppm is a probe failure, not a 20.
            continue

        # THE FLOOR. The probe earns its keep on an instrument that genuinely needs a WIDER
        # window than the SOP; a measured value TIGHTER than the SOP buys nothing and costs
        # identifications. On the one cohort benchmarked (references/diann_parallel.md) the
        # measured 14/7 gave 18,476 precursors against 19,592 at the SOP's 20/7, on the same
        # runs, library and FDR. Flooring keeps the win and drops the loss. The measurement is
        # still made, still recorded, and still has to be plausible to get this far.
        floor = float(SOP_MASS_ACC[key])
        if value < floor:
            out[key] = floor
            out["floored"][key] = True
            out["sources"][key] = (
                f"measured {ppm_text(value)} ppm (median (high) of DIA-NN's per-run "
                f"'{LOG_LINE_OF[key]}'), "
                f"pinned at {ppm_text(floor)} ppm: the measurement is TIGHTER than the SOP "
                f"floor, and a tolerance tighter than the SOP costs identifications without "
                f"buying anything. Floor from {SOP_MASS_ACC_SOURCE}")
        else:
            out["sources"][key] += (f"; at or above the {ppm_text(floor)} ppm SOP floor, so the "
                                    "measurement is what is pinned")
        # The floor can only widen a value, and a mis-set SOP is the one way that could put an
        # implausible number on a command line. Check what will actually be pinned.
        lo, hi = MASS_ACC_BAND[key]
        if not lo <= out[key] <= hi:
            out["rejected"].append(
                f"the {ppm_text(floor)} ppm SOP floor for {LEVEL_OF[key]} is outside "
                f"{band_text(key)}: {SOP_MASS_ACC_SOURCE} needs fixing, not this cohort")
    # Spelled out, beside the historical `ms2_ppm`/`ms1_ppm` keys that hold the same numbers:
    # after the floor, "measured" and "pinned" are different questions and a reader of this file
    # must not have to know which key means which.
    out["pinned_ms2_ppm"], out["pinned_ms1_ppm"] = out["ms2_ppm"], out["ms1_ppm"]
    out.update(
        pin_as=(None if out["rejected"] else
                f"--mass-acc {ppm_text(out['ms2_ppm'])} "
                f"--mass-acc-ms1 {ppm_text(out['ms1_ppm'])}"),
        agree=all(len(set(out[k.replace("_ppm", "_per_run")])) == 1 for k in measured),
        rule=("median (high) of each measured level over the probed runs, as DIA-NN printed it, "
              f"then floored at the SOP ({', '.join(LEVEL_OF[k] + ' ' + ppm_text(v) for k, v in SOP_MASS_ACC.items())} "
              "ppm) so a measured tolerance is used only where it is WIDER; a level given as "
              "--ms1-ppm/--ms2-ppm is pinned as given and never floored; refused, and never "
              "floored, when a measured median is outside "
              f"{'; '.join(LEVEL_OF[k] + ' ' + band_text(k) for k in LEVEL_OF)} or the per-run "
              f"spread of a measured level exceeds {MASS_ACC_MAX_SPREAD:.0%} of it"))
    return out


def run_probe(diann, raw, fasta, lib, threads, timeout, extra="", extra_args=(), workdir=None,
              tag="", measure=("window",), skip=()):
    """Run DIA-NN on one run until it has logged everything in `measure` (bar the `skip` levels).

    Returns {radius, found, missing, lines, timed_out, log, environmental}. `found` is
    LogReader.found -- what was logged, plus `fixed_mass_acc` or `no_auto_announcement` when the
    run is not optimising the mass accuracy the probe was asked to measure (DIA-NN is stopped at
    once: it would never log it). `missing` lists what was asked and not logged; the probe
    succeeded only when it is empty. `environmental` means the failure is not about this run
    (DIA-NN cannot read .raw here, could not be started, or no run will optimise mass accuracy
    under these flags and this DIA-NN build -- every run gets the same ones).

    `workdir` gives DIA-NN its own --temp and --out and holds probe.log. Without them a probe
    that never logs a radius runs on, writes its report into the working directory (in the
    chain, <out> -- the path step 5's check reads) and, if the search completes, its .quant
    NEXT TO THE RAW FILE. They go before --threads, so the flags of the real search stay the
    tail of DIA-NN's argv, as in steps 2-5."""
    workdir = workdir or tempfile.mkdtemp(prefix="probe_window_")
    os.makedirs(workdir, exist_ok=True)
    flags = shlex.split(extra) + list(extra_args)
    cmd = shlex.split(diann) + ["--f", raw, "--fasta", fasta, "--lib", lib]
    if "--temp" not in flags:
        cmd += ["--temp", workdir]
    if "--out" not in flags:
        cmd += ["--out", os.path.join(workdir, "report.parquet")]
    cmd += ["--threads", str(threads)] + flags
    log_path = os.path.join(workdir, "probe.log")
    reader = LogReader(measure, skip)
    res = {"radius": None, "found": reader.found, "missing": reader.missing(), "lines": [],
           "timed_out": False, "log": log_path, "environmental": False}
    deadline = time.time() + max(0.0, timeout)

    def note(msg):
        res["lines"].append(msg)
        _say(tag + msg)
        with open(log_path, "a") as fh:
            fh.write(msg + "\n")

    if _STOP or time.time() >= deadline:
        open(log_path, "w").close()
        res["timed_out"] = not _STOP
        note("[probe_window] stopped by a signal before this probe started" if _STOP else
             "[probe_window] timeout: no time left to start this probe")
        return res

    try:
        with open(log_path, "wb") as sink:
            p = subprocess.Popen(cmd, stdout=sink, stderr=subprocess.STDOUT,
                                 stdin=subprocess.DEVNULL, start_new_session=True)
    except OSError as e:
        res["environmental"] = True
        note(f"[probe_window] could not start DIA-NN ({cmd[0]}): {e}")
        return res
    _LIVE.add(p)
    if _STOP:                       # a signal arrived while it was starting
        _signal_group(p, _HARD)
    try:
        with open(log_path, "rb") as src:
            pending = b""
            while not reader.done():
                exited = p.poll() is not None   # before reading: the read then has it all
                pending += src.read()
                if exited and pending and not pending.endswith(b"\n"):
                    pending += b"\n"              # DIA-NN is gone: its last line is complete
                *complete, pending = pending.split(b"\n")
                for raw_line in complete:
                    line = raw_line.decode("utf-8", "replace").rstrip()
                    res["lines"].append(line)
                    _say(tag + line)
                    reader.feed(line)
                    if reader.done():
                        break              # got it all -- no need to finish the search
                if reader.done() or exited or _STOP:
                    break
                if time.time() >= deadline:
                    res["timed_out"] = True
                    break
                time.sleep(POLL_S)
            if not reader.done() and pending.strip():
                res["lines"].append(pending.decode("utf-8", "replace").rstrip())
                _say(tag + res["lines"][-1])
    finally:
        if p.poll() is None and not (res["timed_out"] or _STOP):
            _end_group(p)                   # values in hand: let it stop gracefully
        else:
            _kill_now(p)                    # timeout, signal, or the leader already gone
        _LIVE.discard(p)
    res["radius"], res["missing"] = reader.found.get("radius"), reader.missing()
    if res["timed_out"]:
        note(f"[probe_window] timeout: DIA-NN stopped after {timeout:.0f} s without logging "
             + describe_missing(res["missing"]))
    res["environmental"] = (res["environmental"] or reader.not_optimising()
                            or any(DOTNET_RE.search(ln) for ln in res["lines"]))
    return res


def probe(diann, raw, fasta, lib, threads, timeout, extra="", extra_args=(), workdir=None,
          tag=""):
    """(radius | None, log lines) for one run -- see run_probe()."""
    r = run_probe(diann, raw, fasta, lib, threads, timeout, extra, extra_args, workdir, tag)
    return r["radius"], r["lines"]


# ------------------------------------------------------------------------------------------
def _gb(b):
    return "?" if b is None else f"{b / 1e9:.2f} GB"


def _min(t):
    return "" if t is None else f", {t / 60:.1f} min"


def _base(path):
    return os.path.basename(path.rstrip("/"))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diann", required=True,
                    help="DIA-NN command (binary path or 'apptainer exec ... diann-linux')")
    ap.add_argument("--raw", nargs="+", default=[],
                    help="the cohort's runs (all of them -- representative runs are chosen here)")
    ap.add_argument("--raw-list", help="file with one run path per line (handles spaces in paths)")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--lib", required=True, help="predicted spectral library (step 1 output)")
    ap.add_argument("--threads", type=int, default=16, help="DIA-NN threads per probe")
    ap.add_argument("--timeout", type=int, default=3600,
                    help="seconds for ONE probe (default 3600); a run that hangs is cut here "
                         "and replaced")
    ap.add_argument("--budget", type=int, default=None,
                    help="seconds for ALL probes together (default: no limit). The chain passes "
                         "its step's wall clock less a margin, so the evidence is always written")
    ap.add_argument("--max-probes", type=int, choices=(1, 3), default=3,
                    help="radii to measure: 3 (default) = median + quartile runs; 1 = the median "
                         f"run. 1 cannot measure mass accuracy: that is pinned from the median of "
                         f"at least {MASS_ACC_MIN_RUNS} runs")
    ap.add_argument("--max-failures", type=int, default=3,
                    help="stop after this many runs logged no radius (default 3)")
    ap.add_argument("--workdir", help="per-probe DIA-NN --temp/--out and logs "
                    "(default: a fresh temporary directory)")
    ap.add_argument("--extra", default="", help="extra DIA-NN flags to match the real search, "
                    "as ONE shlex-quoted string (for hand use; the chain passes them after --)")
    ap.add_argument("--measure", nargs="+", choices=tuple(MEASURES), default=["window"],
                    help="what to measure (default: window). `mass-acc` needs the DIA-NN flags "
                         "WITHOUT --mass-acc/--mass-acc-ms1 -- DIA-NN optimises only what is "
                         "omitted -- and is read from the same DIA-NN run as the radius")
    ap.add_argument("--ms1-ppm", type=float, help="with --measure mass-acc: MS1 has a "
                    "documented DIA-NN value -- pin it as given instead of the measured median "
                    "(the per-run values are still recorded)")
    ap.add_argument("--ms2-ppm", type=float, help="with --measure mass-acc: the same for MS2")
    ap.add_argument("--write-cfg", help="append the measured '--window N' and/or "
                    "'--mass-acc X --mass-acc-ms1 Y' to this cfg file on success")
    # Everything after a bare `--` goes to DIA-NN verbatim, as separate arguments, so bash quotes
    # and expands them exactly as it does for steps 2-5.
    argv = sys.argv[1:]
    extra_args = []
    if "--" in argv:
        k = argv.index("--")
        argv, extra_args = argv[:k], argv[k + 1:]
    a = ap.parse_args(argv)

    measure = [m for m in MEASURES if m in a.measure]       # canonical order, no repeats
    documented = {k: v for k, v in (("ms1_ppm", a.ms1_ppm), ("ms2_ppm", a.ms2_ppm))
                  if v is not None}
    if documented and "mass-acc" not in measure:
        sys.exit("--ms1-ppm/--ms2-ppm pin a documented level of a mass-accuracy measurement; "
                 "add --measure mass-acc, or drop them")
    if any(not v > 0 for v in documented.values()):
        sys.exit("--ms1-ppm/--ms2-ppm must be positive ppm values")
    # A documented level is pinned as given and never measured, so pin_mass_acc() cannot check
    # it. Check it here, before any DIA-NN starts: every number that can reach `pin_as` -- and
    # so massacc.txt, params.resolved.cfg and a DIA-NN command line -- is inside the band.
    for key, value in sorted(documented.items()):
        lo, hi = MASS_ACC_BAND[key]
        if not lo <= value <= hi:
            sys.exit(f"--{key.replace('_ppm', '-ppm')} {value:g} is outside {band_text(key)}, "
                     f"the plausible band for an Orbitrap {LEVEL_OF[key]} mass accuracy. "
                     "DIA-NN's own Orbitrap table spans 4-15 ppm and it calibrates from 25 ppm; "
                     "a value outside that is not a documented tier. Check the resolution the "
                     "value came from.")
    if len(documented) == 2:
        sys.exit("--ms1-ppm and --ms2-ppm both given: both levels are documented, so there is "
                 "nothing to measure -- pin them in the cfg instead of running DIA-NN")
    if "mass-acc" in measure and a.max_probes < MASS_ACC_MIN_RUNS:
        sys.exit(f"--max-probes {a.max_probes} cannot measure mass accuracy: it is pinned from "
                 f"the median of at least {MASS_ACC_MIN_RUNS} representative runs. One run is "
                 "exactly DIA-NN's own first-run auto mode -- 'use this mode for preliminary "
                 "analyses only' -- which measuring representative runs exists to replace, and "
                 "a lone value agrees with itself whatever it is. Use --max-probes 3, or "
                 "measure only the scan window (--measure window).")
    try:
        flag_words = shlex.split(a.extra) + extra_args
    except ValueError as e:
        sys.exit(f"--extra is not a valid shell word list: {e}")
    if "window" in measure and "--window" in flag_words:
        # DIA-NN 2.7.0 echoes a given --window N as "Scan window radius set to N" among its
        # startup settings (HIVE, 2026-09-16), so this would report the flags' value as measured.
        sys.exit("the DIA-NN flags (--extra, or after --) pin --window, so there is no radius to "
                 "measure: DIA-NN only echoes the given value. Drop --window from them, or "
                 "measure only mass-acc (--measure mass-acc).")

    started = time.time()
    deadline = started + a.budget if a.budget is not None else None
    raws = list(a.raw)
    if a.raw_list:
        with open(a.raw_list) as fh:
            raws += [ln.strip() for ln in fh if ln.strip()]
    if not raws:
        sys.exit("no runs given (--raw or --raw-list)")
    for f, what in ((a.fasta, "fasta"), (a.lib, "library")):
        if not os.path.exists(f):
            sys.exit(f"{what} not found: {f}")

    def report(result):
        print(json.dumps(result, indent=2))
        sys.stdout.flush()

    # Before anything slow. select_representative() stats every run in the cohort -- thousands
    # of files, over the cluster's storage -- and a SLURM kill (or scancel) during it used to end
    # the probe on the default disposition, leaving no window.json at all. With the handlers in
    # place the signal is recorded, selection finishes, and the probe loop stops at once and
    # still writes the evidence.
    for name in ("SIGTERM", "SIGINT", "SIGHUP"):
        if hasattr(signal, name):
            signal.signal(getattr(signal, name), _on_signal)

    try:
        sel = select_representative(raws, max_probes=a.max_probes)
    except NoProbeableRun as e:
        report({"measured": measure, "window_radius": None, "pin_as": None, "radii": [],
                "radii_agree": None, "mass_acc": None, "mass_acc_refused": None,
                "incomplete": None, "failed": [],
                "stopped_because": "no_probeable_run", "probes": [], "selection": e.selection})
        sys.exit(str(e))

    for f in sel["unreadable"]:
        _say(f"[probe_window] WARNING: cannot read {f}; not considered")
    for x in sel["excluded_damaged"]:
        _say(f"[probe_window] WARNING: not probing {_base(x['file'])}: "
             + "; ".join(x["problems"]) + " -- it is still SEARCHED by the chain; check it")
    for x in sel["probed_despite_wal_header"]:
        _say(f"[probe_window] WARNING: probing {_base(x['file'])} anyway: "
             + "; ".join(x["problems"]) + " -- " + x["why"])
    for x in sel["no_acquisition_time"]:
        _say(f"[probe_window] {_base(x['file'])}: {x['why']}; this cohort is ranked by size")
    if sel["mixed_cohort"]:
        _say("[probe_window] WARNING: " + sel["mixed_cohort"])
    for x in sel["excluded_small"]:
        _say(f"[probe_window] not probing {_base(x['file'])}: {x['why']} -- likely a blank, "
             "wash or failed injection")
    _say(f"[probe_window] {sel['n_eligible']} of {sel['n_inputs']} runs eligible, ranked by "
         + ("acquisition time" if sel["rank_by"] == "time" else "size")
         + f"; measuring {len(sel['chosen'])}, {len(sel['reserves'])} in reserve")
    if sel["median_much_smaller"]:
        _say("[probe_window] WARNING: " + sel["median_much_smaller"])

    workdir = a.workdir or tempfile.mkdtemp(prefix="probe_window_")
    targets, reserves = sel["chosen"], list(sel["reserves"])
    queue = [dict(c) for c in targets]
    probes, tails, failures, stopped, k = [], {}, 0, None, 0
    while queue:
        if _STOP:
            stopped = "signal"
            break
        if deadline is not None and time.time() >= deadline:
            stopped = "budget"
            _say(f"[probe_window] budget: {a.budget} s spent; no further probe started")
            break
        c = queue.pop(0)
        k += 1
        base = _base(c["file"])
        if c["role"] == "reserve":
            label = f"probe {k}: {base} (reserve for {_base(c['replaces'])}"
        else:
            label = f"probe {k}/{len(targets)}: {base} ({c['role']}"
        _say(f"[probe_window] {label}, {_gb(c['size_bytes'])}{_min(c['time_s'])}), "
             f"{a.threads} threads" + (", measuring " + " + ".join(measure)
                                       if measure != ["window"] else ""))
        limit = a.timeout
        budget_cut = deadline is not None and deadline - time.time() < a.timeout
        if budget_cut:
            limit = max(0.0, deadline - time.time())
        t0 = time.time()
        wd = os.path.join(workdir, f"probe{k}_{base}")
        r = run_probe(a.diann, c["file"], a.fasta, a.lib, a.threads, limit, a.extra,
                      extra_args, workdir=wd, tag=f"[probe {k}] ", measure=measure,
                      skip=tuple(documented))
        secs = round(time.time() - t0, 1)
        rec = dict(c, radius=r["radius"], seconds=secs, threads=a.threads,
                   timed_out=r["timed_out"], missing=r["missing"], log=r["log"])
        if "mass-acc" in measure:
            rec.update(ms2_ppm=r["found"].get("ms2_ppm"), ms1_ppm=r["found"].get("ms1_ppm"),
                       mass_acc_fixed_by_flags=bool(r["found"].get("fixed_mass_acc")),
                       mass_acc_no_auto_announcement=bool(
                           r["found"].get("no_auto_announcement")))
        probes.append(rec)
        if not r["missing"]:
            _say(f"[probe_window] probe {k}: {base} -> {_values(rec, measure)} in {secs} s")
            continue
        tails[k] = r["lines"]
        failures += 1
        _say(f"[probe_window] probe {k}: {base} -> NO {_names(r['missing'])}"
             + (f" (logged {_values(rec, measure)})" if _values(rec, measure) else "")
             + (" (timed out)" if r["timed_out"] else "") + f" in {secs} s")
        if _STOP:
            stopped = "signal"
            break
        if r["environmental"]:
            stopped = "environment"
            break
        if r["timed_out"] and budget_cut:
            stopped = "budget"
            _say(f"[probe_window] budget: {a.budget} s spent during this probe; no further "
                 "probe started")
            break
        if failures >= a.max_failures:
            stopped = "max_failures"
            break
        if reserves:
            nxt = dict(reserves.pop(0), replaces=c["file"])
            queue.append(nxt)
            _say(f"[probe_window] {base} gave no {_names(r['missing'])}; {_base(nxt['file'])}, "
                 "the next run nearest the median, takes its place")

    # Only runs that logged EVERYTHING asked count, for every value: a run with a radius but no
    # mass accuracy would otherwise put the radius and the mass accuracy on different run sets.
    good = [p for p in probes if not p["missing"]]
    if stopped is None:
        stopped = "measured" if len(good) >= len(targets) else "no_more_runs"
    pin = bool(good) and stopped not in ("environment", "signal")
    # Refusals that are NOT "no run answered": runs DID answer, and what they said must not be
    # pinned. They fail the probe exactly as a run that logged nothing does -- loudly, with the
    # evidence written -- because a wrong number pinned cohort-wide is worse than no number.
    refused = []
    if pin and "mass-acc" in measure and len(good) < MASS_ACC_MIN_RUNS:
        # `stopped_because: budget` after one run leaves one good probe, whose value "agrees"
        # with itself. That is DIA-NN's first-run auto mode wearing the evidence of a measurement.
        refused.append(
            f"only {len(good)} run logged everything asked, and mass accuracy is pinned from "
            f"the median of at least {MASS_ACC_MIN_RUNS} representative runs")
        pin = False
    radii = [p["radius"] for p in good] if "window" in measure else []
    radius = statistics.median_low(radii) if pin and radii else None
    mass_acc = pin_mass_acc(good, documented) if pin and "mass-acc" in measure else None
    if mass_acc:
        refused += mass_acc["rejected"]
    pinned = pin and not refused and (radius is not None or "window" not in measure) and \
        (mass_acc is not None or "mass-acc" not in measure)
    what = _what(measure)
    result = {
        "measured": measure,
        "window_radius": radius if pinned else None,
        "pin_as": f"--window {radius}" if pinned and radius is not None else None,
        "radii": radii,
        # One radius corroborates nothing: with n <= 2 runs there is only ever one probe, and
        # `true` there claimed an agreement that was never measured.
        "radii_agree": (len(set(radii)) == 1) if len(radii) > 1 else None,
        # null unless mass accuracy was asked for AND a run logged all of it
        "mass_acc": mass_acc if pinned else None,
        # why a measurement that WAS made is not being pinned (band, spread, too few runs);
        # null when there is nothing to refuse. `mass_acc` above is null whenever this is not.
        "mass_acc_refused": refused or None,
        "incomplete": (len(good) < len(targets)) if pinned else None,
        "failed": [_base(p["file"]) for p in probes if p["missing"]],
        "stopped_because": stopped,
        "seconds": round(time.time() - started, 1),
        "probes": probes,
        "selection": dict({x: v for x, v in sel.items() if x not in ("chosen", "reserves")},
                          planned=[c["file"] for c in targets],
                          reserves=[c["file"] for c in sel["reserves"]]),
        "note": ("Pin this in the cfg for EVERY step of the parallel chain. The radius and the "
                 "mass accuracies are properties of the acquisition method and instrument, so "
                 "they are valid for all files acquired the same way -- but re-probe for a "
                 "different gradient, cycle time, resolution or instrument."),
    }

    _say("[probe_window] per-run radii:" if measure == ["window"] else
         f"[probe_window] per-run {what}:")
    for p in probes:
        got = _values(p, measure)
        _say(f"[probe_window]   {p['role']:<15} {_base(p['file'])}  {_gb(p['size_bytes'])}"
             f"{_min(p['time_s'])}  "
             + (got if not p["missing"] else
                "NO " + _names(p["missing"]) + (f" (logged {got})" if got else "")
                + (" (timed out)" if p["timed_out"] else "")))

    if len(radii) > 1 and max(radii) - min(radii) > 2:
        _say("[probe_window] WARNING: the measured radii disagree by more than 2 ("
             + ", ".join(str(r) for r in radii) + ")"
             + (f"; {radius} is pinned, their median" if radius is not None else "")
             + " -- DIA-NN's optimisation is noisy, but a spread this wide usually means these "
               "runs are not one acquisition scheme. Look at the per-run evidence before one "
               "radius is used for the whole cohort.")

    # The floor in the job log, not only in the JSON: whoever reads the log has to see that the
    # instrument measured one number and the search will run at another.
    if pinned and mass_acc:
        for key in ("ms2_ppm", "ms1_ppm"):
            if mass_acc["floored"].get(key):
                _say(f"[probe_window] {LEVEL_OF[key]}: measured "
                     f"{ppm_text(mass_acc['measured_' + key])} ppm, PINNED "
                     f"{ppm_text(mass_acc[key])} ppm -- the SOP floor. A measured tolerance is "
                     "used only where it is WIDER than the SOP; tighter than the SOP costs "
                     "identifications and buys nothing. Both numbers are in the JSON.")

    if stopped == "signal":
        report(result)
        sys.exit(128 + _STOP[0])

    if not pinned:
        report(result)                     # the evidence first: window.json is what is read
        for i, p in enumerate(probes, 1):
            if p["missing"]:
                _say(f"\n--- {p['file']}: no {_names(p['missing'])}; DIA-NN log tail "
                     f"(full log {p['log']}) ---\n" + "\n".join(tails.get(i, [])[-25:]))
        if any(DOTNET_RE.search(ln) for t in tails.values() for ln in t):
            _say("\nDIA-NN could not open Thermo .raw: there is no .NET 8 runtime in this "
                 "environment. Export DOTNET_ROOT before running this probe:\n"
                 f'    export DOTNET_ROOT="$(bash {shlex.quote(os.path.join(HERE, "ensure_dotnet8.sh"))} '
                 '| tail -1)"; export PATH="$DOTNET_ROOT:$PATH"\n'
                 "(diann_parallel.py puts that export into every generated step.)")
        # Two different faults, never merged into one message: the second one's fix is not to
        # touch the flags, and saying "remove both" when neither is there sends the reader after
        # something that is already correct.
        if any(p.get("mass_acc_fixed_by_flags") for p in probes):
            _say("\nDIA-NN said it is FIXING the mass accuracy this probe was asked to measure. "
                 "DIA-NN 2.7.0 fixes BOTH levels when either --mass-acc or --mass-acc-ms1 is "
                 "given, so remove both from the DIA-NN flags (pass a documented level as "
                 "--ms1-ppm/--ms2-ppm instead), or do not ask for --measure mass-acc. Every run "
                 "gets the same flags, so no other run was tried.")
        if any(p.get("mass_acc_no_auto_announcement") for p in probes):
            ma_given = [f for f in ("--mass-acc", "--mass-acc-ms1") if f in flag_words]
            _say("\nDIA-NN's settings block ended ('N files will be processed') without "
                 "announcing 'DIA-NN will automatically optimise the mass accuracy', so this run "
                 "would never have logged a measured value and was stopped there. The "
                 "announcement is required because a run with a FIXED tolerance still prints "
                 "'Recommended MS1 mass accuracy setting', which would otherwise pass as a "
                 "measurement.\n"
                 + ("The DIA-NN flags DO carry " + " and ".join(ma_given) + ", worded in a way "
                    "this probe does not recognise as a fixed-accuracy notice: remove them "
                    "(pass a documented level as --ms1-ppm/--ms2-ppm instead).\n" if ma_given else
                    "The DIA-NN flags carry NEITHER --mass-acc NOR --mass-acc-ms1, so nothing "
                    "here fixes the tolerance -- the flags are right, and changing them will not "
                    "help. The likeliest cause is that this DIA-NN release words the "
                    "announcement differently: probe_window.py's AUTO_ACC_RE and FIXED_ACC_RE "
                    "were written against DIA-NN 2.7.0's exact wording. Compare the log tail "
                    "above with them and update the pattern if the wording has changed.\n")
                 + "Every run gets the same flags, so no other run was tried.")
        if refused:
            sys.exit(f"Refusing to pin the {what}: " + "; ".join(refused)
                     + ". Every probe's values are in the JSON above. Do NOT widen the band or "
                     "pin one of these numbers by hand: a mass accuracy this far from what "
                     "DIA-NN itself works with normally means the search was pointed at the "
                     "wrong FASTA or the wrong species, or the batch is miscalibrated -- and the "
                     "median over three runs cannot see that, because all three move together. "
                     "Check the library, the FASTA and the instrument, then re-run.")
        why = {"environment": "a failure no other run can fix -- DIA-NN could not start, "
                              "cannot read .raw in this environment, or no run will optimise the "
                              "mass accuracy under these flags and this DIA-NN build (above) -- "
                              "so no other run was tried",
               "max_failures": f"{failures} runs did not log it",
               "budget": f"the {a.budget} s budget ran out",
               "no_more_runs": "every eligible run was tried"}.get(stopped, stopped)
        sys.exit(f"Could not read the {what}: " + why + " ("
                 + ", ".join(result["failed"]) + "). Do NOT guess a value, and do not pin one "
                 "from what a failed run did log -- a value inconsistent across files is exactly "
                 "the defect this probe exists to prevent. Fix the cause above and re-run. --raw "
                 "is not a way past it: it narrows the cohort the checks are applied to, and "
                 "every run named there still has to clear them, so a run left out as damaged or "
                 "as a blank is left out however it is named.")

    if a.write_cfg:
        with open(a.write_cfg, "a") as fh:
            fh.write("\n")                        # the cfg may not end in a newline
            if radius is not None:
                fh.write(f"--window {radius}\n")
            if mass_acc:
                fh.write(f"--mass-acc {ppm_text(mass_acc['ms2_ppm'])}\n"
                         f"--mass-acc-ms1 {ppm_text(mass_acc['ms1_ppm'])}\n")
    report(result)
    if result["incomplete"]:
        pins = "; ".join(x for x in (result["pin_as"], mass_acc and mass_acc["pin_as"]) if x)
        _say(f"[probe_window] WARNING: pinned {pins} from {len(good)} of the "
             f"{len(targets)} runs planned ({stopped}; nothing usable from "
             + ", ".join(result["failed"]) + "). It is the median of the runs that answered; "
             "the per-run evidence is in the JSON.")


def _what(measure):
    """'scan-window radius', 'mass accuracy', or both -- for messages."""
    return " and ".join({"window": "scan-window radius", "mass-acc": "mass accuracy"}[m]
                        for m in measure)


def _names(keys):
    """What a run did not log, without articles: `radius` -> 'radius'."""
    names = {"radius": "radius", "ms2_ppm": "'Optimised mass accuracy' (MS2)",
             "ms1_ppm": "'Recommended MS1 mass accuracy setting' (MS1)"}
    return " or ".join(names[k] for k in keys)


def _values(p, measure):
    """What one probe logged, as the job log shows it: 'radius 7, MS2 14 ppm, MS1 4.1 ppm'."""
    out = []
    if "window" in measure and p.get("radius") is not None:
        out.append(f"radius {p['radius']}")
    if "mass-acc" in measure:
        for key, level in (("ms2_ppm", "MS2"), ("ms1_ppm", "MS1")):
            if p.get(key) is not None:
                out.append(f"{level} {ppm_text(p[key])} ppm")
    return ", ".join(out)


if __name__ == "__main__":
    main()
