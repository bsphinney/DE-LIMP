# QC runs around a project (`qc_bracket.py`)

The Core injects a QC standard on each instrument between projects, and STAN
(https://ucd.stan-proteomics.org) searches and scores every one. A project's results are only as
good as the instrument was while its samples ran. `qc_bracket.py` finds, for each instrument the
project used, the QC runs around the samples, grades them, and gives one verdict.

**The verdict is staff-only** (Brett, 2026-10-01). It and everything behind it go to staff:

- the session's `logs/` folder, which is Core-internal (never delivered, kept out of the
  session zip);
- the run registry, which only the Core group reads.

No client deliverable carries a verdict word. A check or concern holds the client delivery
until a staff member records an acknowledgement (see "The delivery gate").

## When to run it

- **Only on a Core staff route:** a CoreOmics submission (step 1c, fetched with the Core's
  CoreOmics key) or a Core HT submission (step 1a), the routes the other staff-only steps
  key on. **Never for anyone else**, including a user who says the data came off the Core's
  instruments: this skill is public and clients run it. Skip the step for them, and never
  show them a verdict.
- **After step 2, before the delivery (step 12b).** Running it before step 9 gives staff the
  verdict while there is still time to re-inject samples; the client report does not change
  either way.
- **Not for data from elsewhere.** STAN records the instrument model, not its serial number. A
  timsTOF HT at another lab would be matched to the Core's timsTOF HT QC, and that answer
  would be meaningless.

```
python3 scripts/qc_bracket.py --session <session>
python3 scripts/qc_bracket.py --files /data/*.d /data/*.raw --acquisition-json acq.json --out qc.json
```

With `--session`, it reads:

- the session's raw-file list: `input/raw_files.txt`, or else the search's own record of the
  files it read;
- `input/acquisition.json`, which is step 2's `detect_acquisition.py` output (save it there at
  step 3b).

It writes `logs/qc_bracket.json` (the record, `schema_version` 1 — the Core's board reads it;
a field is never renamed or removed without bumping it) and `logs/qc_bracket.md` (the staff
page), both mode 0640 like `logs/decisions.md`. The record also goes to stdout, and the
plain-English summary to stderr.

**Where to run it:** anywhere that can reach `https://ucd.stan-proteomics.org`, such as a laptop
or a HIVE login node. It makes a few small anonymous HTTP requests. Reading the files is the only
heavy part, and only for a `.raw` that is not in `acquisition.json`: that takes one
ThermoRawFileParser call, about 3-7 s. On a cluster login node it refuses to read more than 5
of them. In that case, save step 2's JSON or run the script under `srun`.

## Exit codes, and what to do

| exit | meaning | do |
|---|---|---|
| 0 | checked; `verdict` is set | tell the staff member the summary |
| 2 | usage error — **no raw-file list** in the session (no `input/raw_files.txt` and no search record of its files), no files — or refused (too many `.raw` to read on a login node) | write the raw-file list; pass step 2's JSON, or use `srun` |
| 3 | **STAN could not be reached**, or did not answer with its run list | **no verdict.** Never say the QC was fine. The delivery is held. Re-run when STAN answers |
| 4 | no file gave both an instrument and an acquisition time | read `not_checked`; fix the reason (often a `.raw` that ThermoRawFileParser could not read) |

A record is written for exits 0, 3 and 4 (`status`: `ok` / `stan_unreachable` /
`nothing_checked`). **Exit 2 writes none**, and with no record the delivery is held (below).
`qc_bracket.py gate` and `ack` (below) exit 0 when the client report may go out, and 2 when it
is held or the input was bad.

## What it checks

For each instrument it checks:

- **the nearest QC run before** the first sample, and how long before;
- **every QC run during** the project, meaning between the first and last sample;
- **the nearest QC run after** the last sample, or "none yet".

Each QC run is compared with the same instrument's QC runs over the **30 days before the
project**. The comparison uses runs of the same acquisition mode and throughput (samples per
day). A 35-minute and a 120-minute gradient differ 1.7-fold in IDs on the same instrument. If
there are fewer than 3 such runs, it compares against every run of that mode and says so. With
fewer than 3 runs of that mode at all (a sparse schedule, or the month after a downtime), there
is nothing to compare with, and the run is graded **check** ("too few baseline runs"), never
normal.

| grade | when |
|---|---|
| concern | no identifications; or IDs under 50% of the 30-day median |
| check | IDs under 80% of the median; MS1 or MS2 mass error over 1.5× the median and at least 1 ppm above it; or FWHM peak width over 1.5× the median |

**IPS is shown, not graded.** Each run's IPS appears with STAN's dashboard colour: green ≥ 80,
amber ≥ 60, red below that. The instrument's 30-day median IPS is shown as well. IPS compares a
run with STAN's April 2026 reference cohort, not with the instrument's own recent runs. On
2026-10-01, 62% of STAN's 4,683 QC runs scored red. The Exploris 480's QC at 38 samples/day had
a median IPS of 28 over the past year. If IPS were graded, most projects on that instrument would
be called a concern while its QC ran exactly as usual. Whether a red IPS that has lasted for
months should raise a `check` is Brett's decision; it is not built in.

**`gate_result` is never used.** STAN records "pass" for every run, including runs that identified
nothing (STAN `docs/qc_gating_and_slack_summary.md`).

Some STAN rows are left out, and each is listed in `ignored_rows`:

- runs named as blanks (`blank...HELA`: a blank before a QC), even when STAN's QC-name rule
  matched them;
- runs with impossible dates (STAN holds a Fusion Lumos run dated 1980-01-02);
- runs an operator hid in STAN.

## The verdict

The QC runs that **speak for the project** are:

- the runs during it;
- the nearest run before and the nearest run after, if each is within 7 days.

A flagged run is **isolated** when the QC runs just before and just after it in STAN are both
normal and both within 7 days of it (`ISOLATION_DAYS`). That usually means a bad injection,
with the instrument working on both sides of it. An isolated run counts one level lower. A
normal run weeks away clears nothing.

| verdict | when |
|---|---|
| **concern** | a run that speaks for the project is graded concern |
| **check** | one is graded check, or an isolated concern; or no QC run within 7 days of the project (a "QC desert": the instrument's state while the samples ran was not measured) |
| **no QC on record** | STAN has no QC run for that instrument (for example an instrument STAN does not track) |
| **good** | otherwise |

The overall verdict is the worst across instruments, in this order: concern > no QC on record >
check > good.

**Samples between normal QC runs.** For every flagged run that speaks for the project, the record
lists `samples_between_normal_qc`. These are the project's files acquired after the last normal
QC run before it and before the first normal one after it. They are the samples a staff member
should look at first.

**How often each verdict occurs.** These rules were replayed over STAN's last year of QC and its
maintenance log on 2026-10-01, using 1,000 random projects of 2 hours to 3 days per instrument.

| instrument | good | check | concern | note |
|---|---|---|---|---|
| timsTOF HT | 43% | 29% | 28% | 13% of its QC runs identified under half its usual, often several in a row |
| Fusion Lumos | 40% | 35% | 24% | |
| Exploris 480 | 41% | 46% | 13% | 16% of projects had no QC run within 7 days |

A `check` is common, and it is meant for staff to look at. It is not a finding about the client's
data.

## The maintenance log

Long QC gaps are often the instrument being down, not QC that was missed. STAN keeps a
maintenance log: `GET /api/instruments/<model>/events`, which is anonymous. It records column
changes and clogs, source cleaning, calibration, LC service and PM, and downtime with an end
date. For each instrument, `qc_bracket.py` reads the log and does four things.

1. **Lists every event within 30 days either side of the project.** For example:
   - `column change 2026-07-30 (date only)`
   - `instrument downtime 2026-06-25 to 2026-07-21`
2. **Explains every gap of more than 7 days** between the QC runs around the project, using the
   events that overlap it. The explanation quotes the log and never guesses a cause. When no
   event overlaps the gap, it says so plainly: "No QC for 28 days (2026-06-24 to 2026-07-22);
   STAN has no maintenance record for this period." For an instrument with an empty log, it
   adds "none is logged for this instrument at all".
3. **Marks the first QC run after each event as the post-maintenance check**, wherever that run
   is shown. If project files ran after the event with no other QC run between them and it,
   then:
   - that QC run speaks for the project however far away it is;
   - a flag on it is never discounted as isolated, because its neighbour before it is the
     instrument before the work.

   A change long before the project, with normal QC runs since, does not count this way.
4. **Raises the verdict to a check** for files acquired after an event and before any QC run,
   or acquired inside a logged downtime.

**Coverage on 2026-10-01:**

| instrument | events logged |
|---|---|
| timsTOF HT | 11 |
| Fusion Lumos | 4, none from May to September |
| Exploris 480 | 0 |

The Orbitrap downtimes are not logged yet, so "no record" does not mean "no downtime".

**Dates with no time.** A bare date, or the noon-UTC time that STAN's older form stores a date
as, means a whole Pacific day. Files acquired on that day are "order unknown". The first QC run
"on or after the day" may have run before the work, and the summary says so.

**`first_run` pins the change.** When an event's `first_run` names one of the project's files or
one of STAN's QC runs, the change is pinned to that run.

**Privacy.** Notes, operator and creator are never copied out, because they can name people.
`first_run` is used only for matching.

**When the log cannot be read**, the verdict rests on the QC runs alone, and every gap says it is
unexplained.

## The delivery gate

`core_submission.py deliver` (step 12b) asks `qc_bracket.delivery_gate()` before an analysis
delivery, dry run included.

**What refuses the delivery outright:** a record for a different file list from the session's.
The record stores the sha256 of the sorted list of files it checked (`files.list_sha256`), and
the gate compares it with the session's raw-file list now. If files were added or removed since
the check (a late batch, re-injections kept, a re-search), it says "re-run step 8e". No
acknowledgement releases this, and `ack` refuses to record one.

**What holds the delivery until acknowledged:**
- any instrument's verdict is **check** or **concern**;
- any file that could not be placed in time (`not_checked`): the QC around those samples was
  never looked at;
- the check could not run (STAN unreachable, nothing placeable);
- the record is unreadable;
- there is no record at all.

**What goes through:** **good** and **no QC on record**. The instruments with no QC on record
are in the gate deliver prints to the staff member and in `logs/qc_gate.json`.

To release a held delivery, a staff member reads `logs/qc_bracket.md`, then records:
```
python3 scripts/qc_bracket.py ack --session <S> --by <HIVE login> --note "<one line>"
python3 scripts/qc_bracket.py ack --session <S> --by <HIVE login> --note "<one line>" \
    --no-record-reason "<why the check could not run>"     # only when there is no record
python3 scripts/qc_bracket.py gate --session <S>     # exit 0 = may go out, 2 = held
```

**Who may acknowledge.** `--by` must be a HIVE login on the Core staff list:
- the list is `$CORE_STAFF_FILE`, by default `/quobyte/proteomics-grp/.config/core_staff.txt`,
  with one login per line and `#` comments;
- it is trusted only when it is one plain file (not a link, not hard-linked) that a Core admin
  (`scripts/core_admins.txt`) owns and nobody else can write. This is `notes.py`'s rule for
  `SENDERS`.

With no trusted list, `ack` refuses names that look like an agent's (`claude`, `assistant`,
`agent`, `bot`, `gpt`…), warns, and marks the acknowledgement `NOT VERIFIED -- the Core staff
list is not set up`. **This stops an agent acknowledging by accident. It is not proof that a
person did it**: that comes later, with the board's Duo approval. Record an acknowledgement
only on a staff member's own words, never in an agent's name.

The acknowledgement:
- goes into `logs/qc_bracket_ack.json`, append-only, carrying the staff-only marker;
- records the login, the computer account, the staff-list state, when (UTC), the note, the
  verdict, the reason it held and any `--no-record-reason`;
- is tied to the record's sha256, so a re-run of `qc_bracket.py` needs a new acknowledgement.
  One made with no record (`--no-record-reason`, required then) never releases a record
  written later;
- needs a one-line note of 1-300 characters, saying what was reviewed and decided.

`ack` refuses a record for another file list (re-run step 8e instead). A damaged
acknowledgement file never counts.

**What deliver writes.** The full gate it applied (verdict, reason, who acknowledged, the note)
goes to `logs/qc_gate.json`, mode 0640. `delivery.json` sits in the service tree, which other
HIVE users can read, so it holds only `qc_gate: {proceed, record_sha256, ack_at}`; it is 0640
and never goes into the session zip.

## Telling people

**Core staff** get the whole record:

- the staff page `logs/qc_bracket.md` (and its copy in the run registry, `qc/qc_bracket.md`):
  the summary, the table, the maintenance log;
- the QC run names in `qc_bracket.json`;
- the samples between normal QC runs;
- the 30-day context.

For a **concern**, a staff member decides whether to re-inject the samples acquired around the
failed QC before anything goes to the client.

**The client** sees nothing from this check by default:

- the report has no QC section;
- the Methods say only: "The Core monitors instrument performance with routine HeLa digest
  quality-control standards; QC records are kept by the UC Davis Proteomics Core." This is a
  statement of the Core's practice, true whatever this project's QC looked like. The earlier
  "acquired regularly on each instrument used" was false for an instrument with a 27-day QC
  gap, or with no QC on record.

The sentence is the same whatever the verdict, and it is printed whenever the project has a QC
record path, even when STAN was unreachable or the record is unreadable: a sentence that came
and went with the check's outcome would tell the client something the verdict does not.

Telling a client about a check or concern is a staff decision. When staff decide to, use plain
words: leave out IPS, FWHM, ppm and QC run names, and do not make it sound like an alarm.

- **check:** "One of the instrument's routine QC checks around your samples was looked at by
  Core staff," then what they found.
- **concern:** say what was found ("a QC run on the day your samples ran identified far fewer
  peptides than usual"), which samples were affected, and what the Core is doing about it
  (re-run, or why the data stand).
- **QC desert:** "The instrument's QC standard was not run within a week of your samples, so we
  cannot show its performance on those days." If the maintenance log explains the gap (for
  example the instrument was down and repaired), say what the log says and that the first QC
  run after the repair is the check. Never offer a cause the log does not record.
- **not checked** (STAN unreachable): never say "it was fine".

## What it cannot know, and says

- **Instrument match is by model.** See "When to run it".
- **STAN's time for a QC run may be the end of the run.** STAN does not record whether a run's
  time came from the file's header or, when it could not read that, from the file's modification
  time (`stan/pipeline/hive_process.py`). For the Core's Orbitrap QC runs it is often the
  modification time, which is the end of the run, up to one gradient length after its start. Two
  runs checked on 2026-10-01 were like this. The project's own files are timed from the start of
  acquisition.
- **A Thermo `.raw` records no time zone** (see below).
- **Files with no readable time or instrument** are listed in `not_checked`, never guessed.

## Time zones (`acq_time.py`)

All comparisons are made in UTC. Times are shown in the Core's local time (Pacific), which is the
instruments' clock and the one STAN's dashboard shows.

- **Bruker `.d`:** `analysis.tdf` GlobalMetadata `AcquisitionDateTime` is ISO 8601 with the
  instrument PC's UTC offset, so it is unambiguous.
- **Thermo `.raw`:** ThermoRawFileParser's metadata `NCIT:C69199` "Content Creation Date"
  (`FileHeader.CreationDate.ToString()`) has no time zone. It is printed in the parser's
  culture: `09/30/2026 03:06:18` on HIVE (invariant culture), or `9/30/2026 3:06:18 AM` on a
  US-English Windows.
  - It was checked on HIVE (srun job 24264231). An Exploris 480 QC run on a 120-minute gradient
    printed 03:06:18, and its file was last written at 05:06:20 PDT. A Fusion Lumos run printed
    13:39:43 and was last written at 15:39:46 PDT.
  - The value was the same with `TZ=UTC` and `TZ=Asia/Tokyo`. So it is the instrument PC's own
    wall clock at the **start** of acquisition, and it is read as `America/Los_Angeles`. Every
    record that carries such a time says so.
  - A day-first culture (en-GB: `30/09/2026`) prints the same shape and cannot be told apart.
    `check_against_mtime` flags a reading that falls after the file was last written.
- **No time-zone database:** Windows Pythons often have none (it needs the `tzdata` package).
  Then the US Pacific DST rule, in force since 2007, stands in for it, and the note says so.

`detect_acquisition.py` (step 2) records these for every file as `acquired_at`,
`acquired_at_source` and `acquired_at_note`.

## Where it goes

| Where | Who reads it | What it holds |
|---|---|---|
| `<session>/logs/qc_bracket.json`, `qc_bracket.md`, `qc_bracket_ack.json`, `qc_gate.json` (deliver's snapshot) | staff, the Core's board | everything |
| run registry (`record_run.py`): `qc/` copies, `run_record.json` `analysis.instrument_qc`, `SEARCH_LOG.md` "Instrument QC around this project (staff only -- never delivered)" | staff, the Core's board | the verdict, the summary and the delivery gate (held or not, who acknowledged) |
| `delivery.json` (`deliver`, in the session or work dir, 0640, never zipped) | staff | `qc_gate`: `proceed`, `record_sha256`, `ack_at` only |
| Methods (`make_methods.py --qc-bracket`; steps 1c.7 and 9d pass it, and so does `finalize`) | the client | ONE fixed sentence about the Core's practice, printed whenever a QC record path is given: whatever the verdict, also when STAN was unreachable or the record unreadable |
| `reproducibility/run_manifest.json` (`provenance.py --qc-bracket`) | the client | the record's path, sha256 and `checked_at` only |
| `Analysis_Report.html` / `.md` / `.pdf`, README, AGENTS.md | the client | nothing |

**Every staff-only output carries a marker** (`qc_bracket:staff-only`): it is the record's and
the acknowledgements' first JSON key, the staff page's first line, and in deliver's gate
snapshot. The session zip and `deliver` leave out **any file carrying it, whatever it is called
and wherever it is**. `--out qc.json` can put a record anywhere, so names alone are not enough.
`deliver` lists such a file as `[SKIPPED] ... staff-only`.
