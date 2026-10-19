# The audio discussion of the results ("Signal to Noise")

A two-host, roughly 20-minute audio episode in which a cell biologist and a statistician talk
through ONE finished analysis. It sits beside the report, and a "Listen" card near the top of
`Analysis_Report.html` links to it. **You write the script** from the report you have already
read. `scripts/make_podcast.py` then checks it against the report, turns it into audio, and
links it in. It never writes a word of the conversation. A CoreOmics submission gets one by
default; internal Core work and anyone else are offered one ("When it is made").

**Who it is for: the collaborator who submitted the samples** -- picture a biologist who knows
their own system well but not mass spectrometry or statistics. The episode has three jobs:

1. **What their data say**: the findings, how confident to be in each, and the caveats.
2. **How their samples were measured and analysed**: the instrument, the acquisition, the
   search, the database and the DE model, in their own study's terms.
3. **How proteomics works, taught with their data as the examples**, so they can read their
   own report. A biologist who has never run a mass spectrometer should finish the episode
   knowing what a precursor, a protein group and a 1% FDR are.

**Define every term the first time it is said**, in a clause, in the host's own words: FDR,
log2 fold change (log2FC), adjusted p (adj.P), the IgG control (or whatever the control is),
protein group, precursor, DIA, and the per-group detection counts (`Detected_<group>`: detected in
k of n runs). A term used before it is defined loses this listener.

**Discuss only what the lab was given.** Every fact comes from the files in the session's
`output/` folder, which the lab receives: the report, `methods.md`, `SAMPLE_QUALITY.md`,
`AUDIT.md`. The Listen card says the episode "discusses this report and methods; items beyond
them are listed in the transcript", so that has to be true. A Core note, a draft or a scratch
file the lab cannot open is an **extra source**. Pass it as `--extra-source <file> --label
"<what it is>"`, and list anything taken only from it under Claims beyond the report, naming
it (e.g. "41 entries removed: from the v2 notes (Core, not delivered)"). `check` fails a
number, symbol or phrase found only in an extra source until such a bullet names it, and the
Listen card adds "It also drew on: <label> (not in this folder)".

Why it is built this way: when a model wrote the script from the report, it passed a number
check and still made things up ("limpa is our Core's custom extension to limma"; "Kcnb2, the
Kv2.1 beta subunit"). It also ran 30 minutes and read like a lecture. You have read the report,
the tables and the figures, so you write the script. The tool then checks the script's tokens
against the report -- numbers, symbol-like words, capitalised names -- which catches invented
ones, not wrong sentences (see "What check cannot catch" below).

## When it is made

- **Only after the report is final**: `Analysis_Report.html` is written, and any expert-review
  or data-quality changes are in. A podcast of a draft goes stale.
- **A CoreOmics submission: by default.** Core staff analysing a CoreOmics submission (the
  session answers one, `submission_report.core_run`). Brett (Core director) made the podcast
  the default on 2026-10-15. Once the report is final, make it without asking, and tell the
  staff member in one line that you are. If they say no, skip it. The consent is his standing
  approval, which `--core-run` records ("Privacy and consent").
- **Internal Core work: offer it.** An analysis of the Core's own samples with no submission (a
  method test, QC; `session.py init --internal`) is offered in one line when its results are
  written up, never made by default. On a yes, `--core-run internal` records the same standing
  approval.
- **Anyone else: offer it when the results go to a collaborator or the lab.** A user outside the
  Core, with their own data, is offered it in one line then. For example: "Want a ~20-min audio
  discussion of these results for the lab (AI-generated, synthetic voices), linked at the top of
  the report? It also explains how the measurement and the statistics work, using their data."
  Make it only on a yes. Their data, their consent: cloud voices need it too ("Privacy and
  consent").
- **It never holds up the report or the delivery.** No key, Gemini unreachable or out of quota,
  a render or a job that fails: say so in one line ("No podcast this time: Gemini's daily quota
  is used up.") and deliver without it. On HIVE, submit the job as soon as `job` writes it. If
  it has not ended when everything else is ready to deliver, do not wait: cancel it (`scancel
  <jobid>`; the chunks it made are cached), say so in one line, and finalize and deliver
  without it. A job left running would change the report while `deliver` copies it.
- **Adding it later.** Make it (submitting the job again resumes from its cache), run
  `finalize` again, then deliver again with a new `--label` (SKILL.md step 12b): the new delivery
  goes beside the first, which is never changed.
- If the analysis is preliminary (a v1 with a v2 planned), say so in the disclosure turn.

## Where it is made

Where the report was finished, so the podcast's checks read the report it links into.

- **On HIVE** when the report was built there: this computer cannot run the skill's Python
  (Windows staff: `check_access.sh` → `local_python3.usable: false`), or any session of record
  on HIVE whose report was made there. You write the script here and put it beside the report
  in `$S/output/podcast/`. `check` runs there through `hive_exec.sh`, and a small SLURM job
  (`make_podcast.py job`) renders it in Gemini voices, verifies it, links it and shares it
  ("Commands", "On HIVE"). A Core run uses the Core's shared key. Anyone else needs their own
  key on HIVE, `~/.config/ucdavis-proteomics/gemini_key` with mode 600, saved in their own
  terminal and never in this chat; `say` does not exist there.
- **Here** when the report was built here: `local` mode, or a copy pulled from a HIVE session of
  record (SKILL.md step 1c.7 on a Mac, or any hive_remote user who built the report here). Make
  the podcast before the finished files go back. Then put `podcast/`, the linked
  `Analysis_Report.html` and `.md`, and `Analysis_Report_with_audio.html` into the session of
  record (step 1c.7's push lines), since finalize and `deliver` read them there. As in
  "Commands", "Here": for a Core run use **`--tts gemini --core-run`**. The Gemini voices are
  far better than `say`, the consent is standing, and `verify` (the only check of the audio)
  uses Gemini anyway. It needs a key on this computer; `make_podcast.py key --core-run` says
  whether there is one. With none, use `--tts say`: it is offline and sends nothing, but with
  no key `verify` cannot run either, so say that someone must listen to it.

## Privacy and consent (read before any cloud voice)

Two things can leave the machine, and nothing else: the **transcript** (render: the turns,
after pronunciation substitutions, plus one line of voice direction) and the **rendered
audio** (verify: downsampled, for transcription). The report, tables and figures are never
sent. Even so, both describe unpublished client data, so both need the consent below.

What Google's Gemini API terms say (checked 2026-09-25 at
<https://ai.google.dev/gemini-api/terms>; "Effective March 23, 2026", page last updated
2026-04-28):

- **Unpaid Services** (the free tier, "unpaid quota in Gemini API"): "Google uses the content you
  submit to the Services and any generated responses to provide, improve, and develop Google
  products and services and machine learning technologies". Also: "To help with quality and
  improve our products, human reviewers may read, annotate, and process your API input and
  output." And: "Do not submit sensitive, confidential, or personal information to the Unpaid
  Services."
- **Paid Services** (the API "through a Cloud Project associated with an active billing
  account"): "Google doesn't use your prompts … or responses to improve our products". Also:
  "Google logs prompts and responses for a limited period of time, solely for detecting and
  preventing violations of the Prohibited Use Policy … and any required legal or regulatory
  disclosures."

So:

- **A Core run: the consent is standing.** Brett (Core director) ruled on 2026-09-28 that the
  Core director's approval is enough to send a Core run's transcript and audio to the TTS and
  ASR service; the collaborator's lab is not asked separately. On 2026-10-15 he made the
  podcast the default for CoreOmics submissions. Pass **`--core-run`** to `render`, `verify` and
  `job`. It records that approval in `podcast.json` (`cloud_tts_consent`, with `core_run` saying
  which kind) in the one wording kept in `make_podcast.py` (`CORE_CLOUD_OK`); never type it
  yourself. `--core-run` refuses a session that answers no CoreOmics submission
  (`submission_report.core_run`). For internal Core work that the staff member said yes to,
  pass `--core-run internal`; it refuses a session whose `session.json` does not say internal
  work (`analysis_for`, from `session.py init --internal`), and a session for a CoreOmics
  submission (attached, or named in `analysis_for`): a client's data is never labelled the
  Core's own. A submission named in `session.json` but not attached is a `--core-run` refusal
  that says to run `submission_report.py attach` first. Never use either for anyone else's
  data.
- **Only Core staff can claim a Core run.** Every member of proteomics-grp can read the Core's
  key, but the approval is for Core staff. So `--core-run` also needs the person running it on
  the Core staff list (`$CORE_STAFF_FILE`, default `/quobyte/proteomics-grp/.config/
  core_staff.txt`, trusted only when a Core admin owns it, `staff.load_staff`). On HIVE that is
  this account, by its uid, never `$USER`. On a laptop HIVE is asked, through `hive_exec.sh`,
  about the account the saved HIVE login reaches. No usable list is a no. `podcast.json` says
  `"by": "Core staff"` and how it was checked; the login goes to the session's staff-only
  `logs/podcast.staff.json` (`staff.py`: names stay with the staff, never delivered or zipped).
- The terms above still apply to whatever tier the key is on: on the free tier Google may use
  the text and audio, and people may read or hear it. Brett accepted that for Core runs.
- **Anyone else's data: ask them.** It is theirs. Say what is sent and what the key's tier
  means, and record their answer: `--cloud-ok "<who agreed>, <date>, <paid|free> key"`, written
  to `podcast.json` as `cloud_tts_consent`. For unpublished data suggest a **paid-tier key**
  (billing enabled on the key's Cloud project), or **`--tts say`**, which runs offline on macOS
  and sends nothing. A free-tier key is for published or synthetic data only.
- Without consent the gemini backend sends nothing and exits, explaining why. `--cloud-ok no`
  (or `false`, `0`, `none`, `declined`) is a refusal. `--core-run` and `--cloud-ok` cannot be
  given together.
- **Keys**, the first found: `GEMINI_API_KEY`; `~/.config/ucdavis-proteomics/gemini_key`
  (`chmod 600`), beside the skill's other per-user settings. On HIVE, where every home folder is
  open to every account, a personal key others can read is **refused**, not used
  (`core_submission.key_mode_problem`, the skill's one rule for key files): chmod 600 it, and if
  it was readable for a while, make a new key. On a laptop the tool only warns. And, **for a
  Core run only**, the Core's shared key on HIVE,
  `/quobyte/proteomics-grp/.config/gemini_key`. Brett installs that one: owner brettsp, group
  proteomics-grp, mode 640, in the Core's `.config` folder (mode 2750). It is used only while a
  Core admin (`scripts/core_admins.txt`) owns the file and its folder and nobody else can write
  either (the Core's rule for its files on HIVE, `slack_collab.trust_problem`), and only while
  accounts outside the group cannot read it; otherwise render names the reason and sends
  nothing. Group-readable is how the Core's members use it, so that is no warning.
  `make_podcast.py key [--core-run]` says where the key would come from, never the key.
- Never read, print, copy or paste a key, and never ask for one in the chat. The tool never
  prints the key and scrubs it from every error with the skill's one list of secret patterns
  (`notify_slack.redact`).

## The hosts

- **Maya, cell biologist.** Vivid, reaches for metaphors, and drives the biology ("so what does
  this mean in a neuron?").
- **Leo, statistician.** Dry and skeptical. He loves QC, design and variance, and he deflates
  over-reading ("what's the denominator?").

The tension carries the episode. Maya gets excited, Leo checks the evidence, and they converge
on what the data support. Give each their own voice and running bits: recurring metaphors and
callbacks are good. The default show name is **Signal to Noise**. Keep the split roughly even:
`check` warns when either host has under 35% of the words.

**The hosts must not claim a real research specialty.** "I study membrane contact sites in
neurons" gives a synthetic voice an authority it does not have. They introduce themselves by
role: "I'm Maya, the biologist of the pair", "I'm Leo, the statistician". `check` fails "I
study…", "I research…", "in my lab", "as a neuroscientist" and the like.

**Leo's standard deflations** -- the over-readings every episode should head off where they
apply:
- enriched is not direct binding (a pulldown, and even more a cross-linked one, finds
  neighbours);
- not significant is not absent, and not unchanged (it may be underpowered);
- the quantities are relative between samples, not absolute amounts;
- with n = 3 per group only large, consistent changes can reach significance (low power), and a
  huge fold change against a near-empty control is a presence call, not a precise ratio.

## Structure

Use 6–10 segments separated by `---`, and aim for **2,500–3,200 words (~17–21 min)**. `check`
warns outside 2,200–3,600.

1. **Cold open.** A hook from the data, in two or three turns.
2. **AI disclosure in the first 3 turns**, spoken plainly: "this is an AI-generated discussion
   …, our voices are synthetic, and the written report is the record" (plus "preliminary" if it
   is). `check` fails without the words "AI-generated" in one of the first 3 turns.
3. **The design**: samples, groups, controls, replicates, and why that design.
4. **How proteomics works, with your data** (required: one full segment, early). Walk their
   samples through the pipeline, using this run's own numbers (runs, precursors, protein
   groups, FDR), taken from the report and the methods:
   - what LC-MS/MS does: digestion to peptides, chromatography, the mass spectrometer;
   - DIA and, on a timsTOF, dia-PASEF, versus DDA;
   - precursors versus protein groups, and why proteins come as groups;
   - what a 1% FDR means, and the target-decoy idea;
   - library-free search and match-between-runs;
   - detected versus inferred values (the detection-probability model), and the per-group
     detection counts that say how often a protein was actually seen in each group;
   - empirical Bayes (why n = 3 is workable);
   - multiple testing (adjusted p).

   Teach the ones that apply to this analysis and skip the rest (for example, MBR for a DDA
   run without it). The general explanations are claims beyond the report, so list them there.
   `check` warns about any of these topics the transcript never mentions.
5. **Did the experiment work?** Whatever the design offers as evidence: identifications per run
   and how even they are, CVs within groups, whether the PCA separates the groups (and not the
   batches), known markers or positive controls behaving as expected (for a pulldown: the bait
   tops its own list and known partners appear), and anything the QC flagged.
6. **The findings, in the report's own order**, with the report's numbers and how confident to
   be in each.
7. **The main caveat or confound** the report raises.
8. **A data-quality "detective story"**, if the report has one (a contaminant, a thin run, a
   batch).
9. **A stats "nerd moment"**: one idea from segment 4 taken deeper on a real hit, for example
   an inferred value behind a big fold change, empirical Bayes with n = 3, or blocking.
10. **Lightning round**: one follow-up experiment each.
11. **What to do with this, then the sign-off** (short):
    - which files to open: `README.html` (the folder's front page), then the report
      (`Analysis_Report.html`);
    - how to tier hits: the per-group detection columns in the `DE_*.csv` tables
      (`Detected_<group>`: detected in k of n runs; `Evidence`), or the tier file if the
      analysis wrote one. A hit seen in every run of its group is solid; one carried by
      inferred values is a follow-up;
    - what to validate first;
    - **Core runs only** (the session answers a CoreOmics submission, so the report ends with
      the Core's "How did we do?" survey): the sign-off may say, in words, "the report ends
      with a short survey; tell us what you thought of this episode". Never speak the address.
      Outside the Core, never mention a survey. `check` never asks for it; it FAILs an
      invitation to a survey ("survey" near "tell us", "feedback", "how did we do", "5-minute"
      or "questionnaire") unless the session is a Core run (`submission_report.core_run`, as
      link and share ask it). "MS1 survey scan" is the DDA term and always fine.

    **Never tier by PropObs.** It is the observed fraction of a protein's precursors over ALL
    runs, so a protein found only in its own group -- exactly the hit a pulldown is looking for
    -- scores low. `check` fails advice to tier, rank, sort, filter or prioritise by it,
    however it is written or spoken ("PropObs", "prop obs", "proportion observed"), within a
    clause; saying not to ("Never tier by PropObs", "rank by fold change, not PropObs") is fine.
    `check` warns when the last two segments never point at README.html, the report, the
    detection columns or the tier file, or validation.

Merge neighbours to stay within 6–10 segments. For example, fold the lightning round into the
close, or the detective story into the caveat.

## Fidelity rules

Every number and data claim comes from the report or the other sources you pass to `check`.
General-knowledge colour is allowed only when it is listed under `## Claims beyond the report`:
background biology, textbook mechanisms, arithmetic glosses ("about 1,000-fold", "about
half"), metaphors, and each host's speculative follow-up. **Never state a protein's function
from memory as if the data showed it.** Describe a figure only from what the report says about
it, or what you can see in its PNG.

What `check` tests, token by token (a token not in the sources passes only if it appears in the
Claims section):

- **Numbers**, kept as digits in the transcript (215, 6,112, 5e-15, −1.3):
  - An integer must match exactly: "29%" does not match a report's "29.4%". Only an unsigned
    count of 10 or less ("4 baits") is exempt, and never one followed by %, percent, fold or ×.
    "3k", "6K", "2M" and "1B" are 3,000, 6,000, 2,000,000 and 1,000,000,000, and are checked.
  - A decimal may be rounded to fewer places: 10.6 matches 10.62.
  - A sign must match: "−2.68" does not match "+2.68". An unsigned number ("fell by 1.3")
    matches either.
  - A p-value keeps its power of ten: say "5.65 times 10 to the minus 10", not "5.65". A bare
    power of ten ("10 to the minus 15" for 5.2e-15) is accepted as an order of magnitude.
- **Quantities in words** -- numbers over ten and their plurals and -fold forms ("thousands",
  "hundredfold"), "tenfold", "a dozen", "twice", "half" -- pass only as the report's own
  phrase. The phrase is the word, any little words after it, and the next content word, within
  one sentence, and it is matched after normalising, within one sentence or table cell of the
  sources. "Half of the runs" passes when the report says "fewer than
  half of the runs"; "half the proteins were inferred" fails when it does not. Otherwise use
  the report's digits, or list your phrase as a claim.
- **Symbols** -- anything with digits, three or more capitals, inner capitals or a Greek letter
  (Jph3, FKBP12.6, SERCA, IgG, TNF-α), and a capitalised name joined to a number by a hyphen
  (IL-6, COVID-19) -- must be spelled as the report spells them (case does not matter). A symbol is a whole word: SOD does not match "sodium". Do not write ordinary
  words in ALL CAPS for emphasis; use *italics*.
- **Capitalised words mid-sentence** (Gapdh, Actb, a person, a place) must be in the sources,
  the claims, a host's name or the show's name. A capitalised word that starts a sentence
  cannot be told from an ordinary one ("Gapdh went up." / "The gel ran."), so it is not
  failed; `check.txt` lists the ones not in the sources as INFO -- read that line.
- **Names**: no person's name except the hosts; the lab name ("the Dickson lab") is the most
  you may use. Pass people's full and first names to `check --forbid-name` (the PI, the
  submitter, lab members, Core staff). They are searched in the transcript, the header lines,
  the claims, the styles and the Pronunciation table.
- **Pronunciation rows** go to the voices without being checked against the report, so a row
  may only re-spell its written form. Its digits -- numerals and number words, read the way
  people say them ("sep two fifty" is 250) -- must be the written form's digits in the same
  order, so `2.68 → two point eight six` fails. A capitalised word or symbol that is not part
  of the written form (`Ryr2 → Gapdh`) fails too. Every row is listed in `check.txt`. Write phonetic respellings in lowercase ("teck R").
- **The disclosure** ("AI-generated") is in one of the first 3 turns, and **no host claims a
  real specialty**.
- **No draft history.** The listener never saw a draft of the report, so no host says what an
  earlier draft said, got wrong or fell for, or that something "was corrected". A staff report:
  "the first draft of this report fell for it" passed `check` and reached a podcast script.
  The hosts talk about what the report finds and what was checked, in the present tense: "the
  review checked whether a loading difference could explain this, and it does not", never "the
  first draft fell for it, and we corrected that". `client_text.py` holds the phrases, in two
  tiers. `check` FAILS on the unambiguous ones ("first draft", "an earlier version of this
  report", "fell for it", "we got this wrong", "the report originally said", "we've updated the
  wording", "... was incorrect, and it has been changed", "honest correction") and WARNS on
  words that are draft history only in context ("corrected", "fixed", "updated" in a sentence
  that also speaks of a review, a mistake, a draft or the report; "the draft"; "an earlier
  pass"; "initially reported as"). Science wording trips neither ("the p-values are corrected",
  "batch effects were corrected", "fixed modification", "a draft genome", "fell for this group
  of mice", "previously reported"). What neither tier catches is the meaning check's (f).

### The meaning check: a separate reader, before render

`check` matches tokens. A staff episode (2026-10-16) passed it with 0 problems and 0 warnings,
and a separate read of every line against the report then found: a hit from one comparison told
as part of another comparison's story and called "nearly absent" (it was detected in 2 of 3
cultures); proteins
put right after "It survives everything" whose robustness table says `robust` = no or "not
testable"; the report's "not a main driver" spoken as "aren't driving"; a tool presented as a
database the client already uses; a number given for "the shifted cultures" that belongs to one
of them. So **after `check` passes, an agent (or person) who did not write the script reads every
sentence against the report** before anything is rendered, and checks it for:

- **(a)** the right comparison, group and direction for every number and named protein;
- **(b)** every hedge the report gives, kept ("not a main driver" is not "isn't driving");
- **(c)** the robust status of every named protein: one spoken next to "robust", "survives" or
  "holds up" must be robust (`robust` = yes) in that comparison's robustness table;
- **(d)** no superlative or totality the report does not make ("the only", "all", "every",
  "the strongest");
- **(e)** no claim outside the report and the Claims list;
- **(f)** no draft history: nothing an earlier draft or version said, got wrong or had changed.

The reviewer's brief, with every sentence numbered by its line and the robustness tables found
in `output/tables/`: `make_podcast.py meaning <script> --brief` (also written as
`podcast/meaning_brief.md`). Spawn the reviewer with it (one `Agent` call; not the agent that
wrote the script). The writer fixes what it finds, then `check` again. The reviewer records the
verdict; the tool takes the script's sha256 itself:

```bash
python3 scripts/make_podcast.py meaning $P/podcast_script.md --verdict PASS \
    --reviewer "<who: an agent or person that did not write the script>" --findings <their file>
```

That writes `podcast/meaning_check.md`: the script's sha256, `verdict: PASS` or `FAIL`, the
reviewer, the sha256 of the brief, of each delivered source `check.txt` names and of the
robustness tables, the checks, the mechanical robustness lines and the findings. It is recorded
only after `check` passes for this script and its sources and after `--brief` was run for this
script, and it refuses a reviewer whose name is the script's `Author:` line and an empty
findings file. **`render`, `link`, `share` and `job` refuse without a PASS that still holds**:
for the script as it is now, the report and other sources as they were read, and the same
robustness tables. An edited script, a rewritten report or new tables need `check` and the
meaning check again. `--no-meaning-check "<why>"` (render, job) goes on without it and is
recorded in `podcast.json` `meaning_check` (`--unchecked` overrides only the token check at
render; at link and share it also carries a meaning check that went stale, as a warning). An
episode rendered before this check existed (no `meaning_check` in its `podcast.json`) is linked
with a warning. A job or requeue written by 2.13.1 runs the new render and refuses without a
meaning check: record one, then submit again.

What `check` now catches of it mechanically: a protein named **in the same clause** as a robust
claim ("robust", "survives", "holds up") that every row of the robustness tables calls not
robust is a FAIL; robust in some comparisons only, or a `robust` cell it does not understand, is
a WARN. A clause is cut at `.`, `;`, `!`, `?` and at "but", "while", "whereas", "unlike",
"except" ("Jph3 survives every check; Ryr2 doesn't" names Ryr2 in another clause), and at a
comma that opens a second clause (", and", ", so", ", which": "Gapdh survives every check, and
Ryr2 goes up threefold" says nothing about Ryr2), but never at a colon ("It survives everything:
Ryr2, Jph3, and Actb" is about all three), in a list of names before the claim ("Ryr2, Jph3, and
Actb all survive") or at ", and so does Ryr2". A negation before the claim ("Neither Ryr2 nor
Actb is robust", "Not one of them ... is robust") or beside a name ("Ryr2 and Actb don't", "drop
out", "unlike Ryr2") is not a claim; one after the claim word is ("Ryr2 holds up, which doesn't
happen often"). A gene symbol is matched as written, so SET, WAS or MAX in the table never match
"set", "was" or "max"; one with a digit, or written in capitals, matches in either case (Jph3 /
JPH3, ACTB for a mouse table's Actb). A common English word that starts a sentence ("Was that
robust?", "Set against that, ...") is read as English, not as the mouse gene; in capitals (SET)
it is the gene. **The tables** are every CSV under `output/tables/` (three folders
deep) with a `robust` column, whatever made them: `yes`/`no`, and `TRUE`/`FALSE`, `1`/`0` or
"not testable ..." as R and pandas write them, matched on `Genes`, `Protein.Names` and
`Protein.Group`, read with or without a byte-order mark; `check.txt` lists them. A superlative
("the only", "the strongest", "the top", "none of" ...) whose words are not in the sources is a
WARN for the meaning check to confirm. The mechanical lines are a floor, not the check: the
meaning reviewer reads robust status across sentences too.

### What check cannot catch

Read the script against the report for these; a PASS says nothing about them:

- small integers: an unsigned count of 10 or less ("4 baits") is not checked;
- context: a real number or protein attached to the wrong protein, contrast, group or figure;
- relational words: "higher", "more than", "most", "only", "the top hit";
- false statements built from true numbers and true names;
- lowercase symbols ("gapdh") and lowercase respellings in the Pronunciation table;
- a capitalised word that starts a sentence ("Gapdh went up."): the INFO line lists them;
- numbers with a unit other than %, fold, ×, k, M and B ("3 kDa", "2 µg");
- meaning: a caveat the report makes can be dropped, and speculation can be spoken as fact;
- draft history in words the phrase list does not have ("we had that backwards before").

## The script file (`<session>/output/podcast/podcast_script.md`)

```markdown
# Signal to Noise — <episode title>

- **Show:** Signal to Noise
- **Title:** <episode title, shown on the Listen card>
- **Hosts:** Maya (cell biologist), Leo (statistician)
- **Voices:** Maya=Kore, Leo=Charon          (optional; Gemini prebuilt voices, these are the defaults)
- **Say voices:** Maya=Samantha, Leo=Daniel  (optional; macOS `say`, these are the defaults)
- **Styles:** Maya: curious and energetic; Leo: calm, precise and dryly funny  (optional; how each host sounds on Gemini)
- **Written by:** <model>, <date>, from <the sources>

## Pronunciation

| Written | Spoken |
|---|---|
| Jph3 | J P H 3 |
| Tecr | teck R |

## Claims beyond the report

- <each statement that is not in the report: background biology, metaphors, glosses, speculation>

## Transcript

<!-- TRANSCRIPT START -->
**MAYA:** One turn per line, the whole turn on that ONE line.

**LEO:** Blank lines between turns are fine.

---

**MAYA:** A line with only --- starts a new segment.
<!-- TRANSCRIPT END -->
```

- The header lines may be plain (`Show: …`) or bulleted and bold as above.
  `Hosts: MAYA (cell biologist, voice Kore) · LEO (statistician, voice Charon)` also works.
- There must be exactly two hosts, and each turn must start with `**NAME:**`. A wrapped turn or
  an unknown speaker is an error, never silently dropped.
- `## Pronunciation` is required, but its table may be empty. `## Claims beyond the report` is
  required: write bullets, or `None`.
- **Pronunciation** applies only to the text sent to the voices; the transcript keeps the
  written form. Entries match whole tokens (no letter or digit on either side), longest first,
  in one pass. Your table adds to the built-in defaults and wins on the same spelling. The
  built-in defaults are:
  - `Kv2.1` → "K V two point one" (also `Kv2.2`, `Kv2`)
  - `IgG` → "I G G"
  - `DIA-NN` → "D I A N N"
  - `dia-PASEF` → "dia pasef"
  - `timsTOF` → "tims toff"
  - `log2FC` → "log two fold change"; `log2` → "log two"
  - `adj.P` / `adj.P.Val` → "adjusted p"
  - `m/z` → "m over z"
  - `1/K0` → "one over K zero"
  - `ER` → "E R" (whole word only)
  - `PropObs` → "prop obs"
  - `LC-MS/MS`, `MS/MS`
  - `e.g.`, `i.e.`, `vs`

  Numbers are made speakable automatically: `5e-15` and `5×10⁻¹⁵` → "5 times 10 to the minus
  15"; `−1.3` → "minus 1.3"; `≥`, `≤`, `±` and `~` become words. Add only the study's own
  terms. `render --dry-run` prints exactly what the voices will read.

## Commands

### Here (a Mac, or `local` mode)

```bash
S=<session>; P=$S/output/podcast
# 1. write $P/podcast_script.md (above), then check it until it passes
python3 scripts/make_podcast.py check $P/podcast_script.md \
    --source $S/output/AI_Analysis_Report.md $S/output/methods.md $S/output/AUDIT.md \
             $S/output/SAMPLE_QUALITY.md \
    --forbid-name "<PI full name>" "<submitter full name>" "<first names>"
#    -> $P/check.txt. Fix every FAIL line, re-run until PASS. Read the WARN lines.
# 1b. the meaning check: a reader who did not write the script (an Agent call) reads every
#     sentence against the report from the brief; fix, check again; they record the verdict
python3 scripts/make_podcast.py meaning $P/podcast_script.md --brief --source $S/output/AI_Analysis_Report.md
python3 scripts/make_podcast.py meaning $P/podcast_script.md --verdict PASS \
    --reviewer "<who>" --findings <their findings file>
# 2. preview what will be spoken (sends nothing)
python3 scripts/make_podcast.py render $P/podcast_script.md --tts say --dry-run | less
# 3. render. A Core run: Gemini voices under the standing approval (verify runs at the end)
python3 scripts/make_podcast.py render $P/podcast_script.md --tts gemini --core-run
#    ... anyone else, once they agreed:
python3 scripts/make_podcast.py render $P/podcast_script.md --tts gemini \
    --cloud-ok "<who agreed>, <date>, <paid|free> key"
#    ... or offline (macOS; no key needed, nothing sent):
python3 scripts/make_podcast.py render $P/podcast_script.md --tts say
# 4. verify: the ASR round trip (automatic after a consented gemini render; by hand after say)
python3 scripts/make_podcast.py verify $P/podcast_script.md --core-run   # or --cloud-ok "<who, when>"
#    -> $P/verify.txt. Listen to every segment it names; fix, re-render, verify again.
# 5. link it into the report, README and AGENTS.md (safe to re-run)
python3 scripts/make_podcast.py link $S/output
# 6. the one file to send: the report with the audio and transcript built in (finalize also
#    makes it; safe to re-run)
python3 scripts/make_podcast.py share $S/output
```

A session of record on HIVE whose report was built here (SKILL.md step 1c.7 on a Mac, or a
hive_remote user): make the podcast on the local copy before the finished files go back. The
push (step 1c.7's lines) puts `podcast/`, the linked `Analysis_Report.html` and `.md`, and
`Analysis_Report_with_audio.html` into the session of record, where finalize and `deliver`
read them.

### On HIVE (the report was built there; Windows staff)

Every step runs in the session of record `$S` on HIVE, once the finished report is there
(`Analysis_Report.html`, `AI_Analysis_Report.md`, `methods.md`, `AUDIT.md`,
`SAMPLE_QUALITY.md` in `$S/output/`). You write the script on this computer; nothing else
runs here.

```bash
# 1. write podcast_script.md here, put it beside the report, check it there until it passes
bash scripts/hive_exec.sh "mkdir -p $S/output/podcast"
bash scripts/hive_exec.sh --put <here>/podcast_script.md "$S/output/podcast/"
bash scripts/hive_exec.sh "python3 ~/proteomics-pipeline/scripts/make_podcast.py check \
    $S/output/podcast/podcast_script.md --source $S/output/AI_Analysis_Report.md \
    $S/output/methods.md $S/output/AUDIT.md $S/output/SAMPLE_QUALITY.md \
    --forbid-name '<PI full name>' '<submitter full name>' '<first names>'"
bash scripts/hive_exec.sh "cat $S/output/podcast/check.txt"   # fix here, --put again, check again
# 1b. the meaning check (a reader who did not write the script; see above), recorded on HIVE
bash scripts/hive_exec.sh "python3 ~/proteomics-pipeline/scripts/make_podcast.py meaning \
    $S/output/podcast/podcast_script.md --brief --source $S/output/AI_Analysis_Report.md"
bash scripts/hive_exec.sh --put <here>/findings.md "$S/output/podcast/"
bash scripts/hive_exec.sh "python3 ~/proteomics-pipeline/scripts/make_podcast.py meaning \
    $S/output/podcast/podcast_script.md --verdict PASS --reviewer '<who>' \
    --findings $S/output/podcast/findings.md"
# 2. preview what will be spoken (sends nothing; light, so no job)
bash scripts/hive_exec.sh "python3 ~/proteomics-pipeline/scripts/make_podcast.py render \
    $S/output/podcast/podcast_script.md --tts gemini --dry-run"
# 3. the job: render in Gemini voices (verify follows), link, share -- submit it at once
bash scripts/hive_exec.sh "python3 ~/proteomics-pipeline/scripts/make_podcast.py job \
    $S/output/podcast/podcast_script.md --core-run"
#    internal Core work they said yes to: --core-run internal
#    anyone else: --cloud-ok "<who agreed>, <date>, <paid|free> key" (their own key on HIVE)
bash scripts/hive_exec.sh '<the "submit" line it printed>'
# 4. watch it like any job, with the "log" it printed; when it ends, read verify.txt
bash scripts/watch_run.sh --hive --slurm <jobid> --log "<the log it printed, with the job id>"
bash scripts/hive_exec.sh "tail -n 30 <that log>; cat $S/output/podcast/verify.txt"
```

- **`job` checks first and writes nothing on a problem**, each a one-line refusal: the script
  must be `output/podcast/podcast_script.md` in a session's output folder, its check must hold,
  the consent must be there (`--core-run` with the person on the Core staff list, or
  `--cloud-ok` for anyone else), a python must have
  google-genai >= 2.25.0 (the pipeline env's, from `setup.json`; `bash
  ~/proteomics-pipeline/scripts/setup.sh` installs it, and `setup.sh --check` reports it as
  `podcast`), and there must be a key (for a Core run, the Core's shared key).
- **The job** is `~/.proteomics-pipeline/jobs/podcast_<PROT or session>_<hash>.sh`, its log
  `..._<jobid>.log` beside it; `job` prints both. That folder is yours alone (mode 700). Core
  session folders are group-writable with no sticky bit, so a script there could be swapped for
  another member's before `sbatch` reads it, and the job would run as you. The script is written
  without following a link, and a link where it goes is refused. The job asks for 1 CPU, 2 GB
  and 2 hours, on the queue `run_search.slurm_queue` picks for this account (as every job the
  skill writes; `--requeue` on a preemptible one, with `--open-mode=append` so a requeued job
  adds to its log). Its `--output` path is in double quotes, so a space in it is fine (single
  quotes leave no log in a real run). It is a login shell (`#!/bin/bash -l`), so it can
  `module load ffmpeg/7.0.2` for the AAC; when no ffmpeg can be found, its log says so in a
  `[WARN]` (the episode then stays a ~58 MB `podcast.wav`, and `share` cannot build). Then `render --tts gemini --core-run` (verify follows by itself), `link` and `share`;
  a failed step stops it there. Submitting it again resumes: finished chunks are cached. An
  episode takes minutes; rate limits can stretch that to an hour.
- **Never run render, verify or share on the login node**: they are minutes of network calls.
  HIVE's compute nodes reach Google's API (checked 2026-10-15). `check`, the preview, `link`
  and `key` are light, so they may run through `hive_exec.sh` directly.
- **After the job, the HIVE report carries the Listen card.** Never put a local copy of
  `Analysis_Report.html` over it. To change the report later, rebuild it
  (`make_analysis_html.py` keeps the card), then `check` and the job again: only changed chunks
  are rendered again.
- Then finalize (step 12) and deliver (step 12b), which find `podcast/` and
  `Analysis_Report_with_audio.html` in `$S/output`. **They never wait for it.** A job that
  failed: one line to the staff member, then finalize and deliver without it. A job still
  running when everything else is ready: `scancel` it, one line, then finalize and deliver
  without it, and add the audio later ("When it is made").

- **Sources.** Pass the delivered report text you wrote from, all from the session's `output/`
  folder: `AI_Analysis_Report.md` (or the report's `Analysis_Report.md` twin, or
  `Analysis_Report.html`, whose embedded images and scripts are ignored), `methods.md` (the
  teaching segment's instrument, search and FDR numbers come from it), `AUDIT.md` and
  `SAMPLE_QUALITY.md`. A `--source` outside `output/` (or inside `output/podcast/`) fails;
  give it as `--extra-source <file> --label "<what it is>"` instead (see "Discuss only what the
  lab was given"). `check` finds `output/` as the folder above the script's `podcast/` folder,
  or from `--output-dir`, and only accepts a folder that holds the report
  (`Analysis_Report.html`) or the DE record (`tables/de_provenance.json`): pointing it at the
  session or its parent would let a draft beside `output/` count as delivered. Do **not** pass the DE tables: with thousands of numbers in the
  sources, almost any number would match something, and `check` warns when that happens.
- **No absolute paths are delivered.** `check.txt` records the script by name and each source
  by its path relative to `check.txt` (an extra source by its file name and label), each with
  its sha256; forbidden names are counted, not listed. `podcast.json` copies those. A render
  into another folder rewrites the relative paths so the check still holds there.
- **Render refuses** unless `check.txt` says PASS for this exact script and every source is
  unchanged since (each is hashed without any podcast block, so link's own Listen line does
  not count). `--unchecked` overrides this and is recorded in `podcast.json`; use it only when
  the user asks. `link` refuses the same way (the card vouches for the episode), with its own
  `--unchecked`.
- **Gemini** needs `google-genai`: `setup.sh` installs it into the pipeline env (conda-forge,
  >= 2.25.0), or `python3 -m pip install google-genai`; the tool says so if it is missing.
  - **Model order.** Names are read from the API (`models.list`), never assumed: a pro TTS
    model, then `gemini-3.8-flash-tts`, then `gemini-2.5-flash-preview-tts`. Any other TTS
    model (a preview or lite one) is used only with `--model`.
  - **Two APIs, chosen from the model name**
    (<https://ai.google.dev/gemini-api/docs/speech-generation>):
    - 2.x models take one text prompt through `generate_content`.
    - Gemini 3 and later take `interactions.create`: one text part per turn, annotated with
      its speaker and that host's style, and `speech_config` `{"mode": "conversational",
      "speakers": [...]}`. The audio comes back as base64 WAV. `generate_content` fails on
      these models with a 400 ("must specify speaker names for each part").
    - Each host's style comes from an optional `Styles:` header line, for example
      `Styles: Maya: curious and energetic; Leo: calm, precise and dryly funny` (these are the
      defaults).
  - **Handing over.** On a fresh render, a model that is missing, has no quota, or rejects the
    request before it has made its first chunk hands over to the next.
  - **One model per episode.** Once a model has made a chunk, the episode stays on it, because
    a second model's voices sound different. A resume (chunks already in the cache) and
    `--redo` are pinned to the model that made the cache: they never hand over. If that model
    is out of quota, render stops, names the chunk and keeps the cache. Re-run later to
    resume, or pass `--model <other>` to re-render every chunk with that model. Only an
    explicit `--model` render prunes another model's cached chunks.
  - **A model that stops answering is left, and every chunk is made again.** A staff render
    (2026-10-16) sat 28 minutes on a stalled `gemini-3.8-flash-tts` while
    `gemini-2.5-flash-preview-tts` answered a two-line test in 2.9 s; its first chunks were
    cached from the stalled model, so it could never hand over. Now:
    - **The two-line test (preflight).** Before the first chunk a render must make, each
      candidate model in turn (the cache's model first) gets one tiny two-speaker request with
      a 60 s deadline; the first that answers is used. A rate limit counts as an answer; a
      missing model, no quota, a rejected request or no answer moves on; a key problem stops.
      If the model that made the cached chunks does not answer, the next one makes **every
      chunk again**, so the voices never mix. Nothing to make, no test. Each test is one request
      against the key's daily quota (a free-tier key's TTS quota is small), and a resume walks
      the models until one answers (up to 60 s each). `--no-preflight` skips it;
      `podcast.json` `tts.preflight` records which.
    - **Mid-episode.** When a chunk's request goes unanswered twice (`--request-timeout`), the
      render moves to the next candidate model and makes **every chunk again** with it, and
      says so in the log; with no other model it stops with the resume message, as before.
      `tts.switched` and the warnings record each hand-over. The next model then hands over
      like any first model: one that is missing, out of quota or rejects the request is passed
      over (never "made this episode's cached chunks", which it did not).
    - **A stall is remembered.** A model that stopped answering mid-episode (twice on a chunk)
      is written to `podcast/.cache/stalled.json`; for 12 hours a resume does not pin it, even
      when it made the most cached chunks, and the model used makes every chunk again. So a
      resume cannot repeat the same timeouts. A miss at the two-line test is not remembered
      (one short outage would otherwise shut a model out for 12 hours). The remembered models
      are still tried, last, before a render gives up -- on a free-tier key the pro models
      never answer, so they may be all there is -- and one that answers again is forgotten.
      When nothing answers, the message names them; `--model <name>` tries a model whatever
      happened before.
    - **No silent retries in the SDK.** The client is built with a per-request HTTP timeout
      (`--request-timeout`) and `HttpRetryOptions(attempts=1)` (`tts.http`). Read from
      google-genai 2.31.0: `generate_content` then never retries; the `interactions` client
      (Gemini 3) takes `attempts` as its retry count, so it keeps one retry (without the option
      it made three, obeying a server's Retry-After with nothing logged). The render's own
      deadline bounds every request either way.
  - **The key's tier.** When a 429 names a free-tier quota (`generate_content_free_tier_requests`,
    `…-FreeTier`), render says the key's Google project is on the FREE tier and records
    `key_tier: "free"` in `podcast.json` (`null` when no answer named a tier: never guessed).
    A hand-typed `--cloud-ok` note that says "paid" then gets one plain warning that the note
    does not match the key, in the log and the warnings. It never stops a render: the Core's
    standing approval covers free-tier Core podcasts (2026-09-28).
  - Gemini 3 TTS needs a google-genai with `client.interactions` (2.25.0 has it). With an older
    one the 3.x models are left out, and render says to run `python3 -m pip install -U
    google-genai`.
  - **Unreadable audio is not retried.** An answer is parsed once; audio that cannot be read
    stops the chunk, because asking again would be billed again. The RIFF header is parsed by
    hand, so WAVE_FORMAT_EXTENSIBLE and streaming (size 0) headers work on Python 3.9 too. An
    empty answer is asked once more.
  - **Rate limits.** A 429 on a per-minute quota, or one that carries a `retryDelay`, is waited
    out, not fatal. Render reads the delay and the quota id from the error details, prints
    `[render] rate-limited, waiting Ns`, sleeps the delay + 5 s, and retries the same chunk,
    up to 10 times. All the waits in one render are capped by `--rate-wait-budget` (default 30
    min), after which it stops with the resume message, as it does at once for a daily quota
    (`…PerDay…`).
  - **Pacing.** Requests to one model are spaced at least 20 s apart (`--min-interval`).
  - **A request that does not answer.** Each request logs `[render] chunk k/N ...: sent to
    <model>` and, while it waits, `still waiting for <model> (Ns)` every minute. With no answer
    in `--request-timeout` seconds (default 600; a chunk takes about 3 minutes) it is asked once
    more, then the next model takes over (above); with no other model render stops with the
    resume message, so the job fails where its watcher sees it instead of hanging until its
    time limit. The time lost to unanswered requests in one render
    is capped as well (`--request-timeout-budget`, default 30 min). Re-running resumes from the
    cached chunks. (An abandoned request may still finish and be billed: the retry can cost that
    chunk twice.)
  - An episode is one request per segment (a segment over 500 words is split at turn
    boundaries), so about 8–12 requests, and at least 3–4 minutes with pacing. Google does not
    publish the free-tier TTS limits. They are per Cloud project and shown at
    <https://aistudio.google.com/rate-limit>.
- **Verify: the audio heard back.** You cannot listen, so `verify` has a Gemini text model
  transcribe the rendered episode and compares that with the spoken script (the text after
  pronunciation).
  - The audio goes to Google downsampled to 16 kHz mono AAC at 32 kbps (about 6 MB for 24
    min), through a 16 kHz WAV (afconvert needs one; ffmpeg works too). The transcriber is
    `gemini-2.5-flash`, falling back to the newest 3.x flash text model, with the prompt
    "Transcribe verbatim. Write numbers as digits."
  - Both texts are normalised to words: spelled-out numbers become digits ("five thousand and
    twenty-four" -> 5024), and spelled-out letters are joined ("I G G" -> igg). They are then
    compared with difflib (`autojunk=False`).
  - It reports the word match ratio; every script span of 6+ words not heard (dropped or
    garbled audio); and every number in the script that is not in the transcript, with ±60
    characters of the transcript where it should be. That context tells a voice's misread
    from the ASR's mishearing. Counts of 10 or less are not listed.
  - It writes `verify.txt`, `verify_transcript.txt` and a `verify` block in `podcast.json`
    (ratio, gaps, numbers_not_heard, segments_to_check, model). It WARNs when the ratio is
    below 0.93, or when there is a gap or a number not heard, and names the segments to
    listen to.
  - Every render removes the last render's `verify.txt` and `verify_transcript.txt`, so an
    old OK never stands beside new audio. When this render's verify does not run (it failed,
    for example on a 429; `--no-verify`; an offline `say` render), `verify.txt` says `NOT RUN`
    with the reason, and so does the `verify` block: nothing has heard the audio back.
  - It runs at the end of every `render --tts gemini` given consent (`--core-run` or
    `--cloud-ok`; `--no-verify` skips it). It never fails the render and never touches the
    report. It **sends the audio** to Google, so the same consent rules apply: nothing without
    that consent. Run by hand, `verify` needs the same `--core-run` or `--cloud-ok`.
  - Calibration, PROT_0756: a first render scored 0.969 with no gaps and found one real
    misread (the voice said "5,244" for "5,024"), which a Pronunciation row fixed. The clean
    26-minute episode on gemini-3.8-flash-tts then scored ratio 0.965, coverage 0.957, 0 gaps,
    no numbers flagged. A perfectly heard episode scores about 0.96-0.97, not 1: a respelled
    symbol ("rye R 2") is not written the way the ASR writes it ("Ryr2").
  - A number glued to letters is still heard: both texts are split at letter/digit
    boundaries ("LRS124" and "LRS-124" are both "lrs 124"), except inside 5e-15 and 3k.
  - **Fixing what it finds.** A misread number gets a Pronunciation row that spells it out,
    for example `| 5,024 | five thousand and twenty-four |`. Then re-run `check` and
    `render`: only the chunk with that line is re-made. Dropped or garbled audio with the
    text unchanged: `render --redo <segment> ...` re-makes just those segments. Then `verify`
    again. The ASR can mishear too, so a flagged line is a place to listen, not proof of a
    fault, and a clean result is not a listen.
- **say** (macOS) uses one voice per host (Samantha and Daniel by default), runs offline and
  needs no key. A 22-minute episode rendered in about 2 minutes.
- **Re-runs are cheap.** Every chunk is cached as `podcast/.cache/<sha256>.wav`, keyed on its
  spoken text, voices and model, so editing one line re-makes one chunk. After an edit,
  re-run `check`, then `render`. A chunk's warnings are kept beside it (`<sha256>.json`), so a
  resumed render still reports them, and a successful render removes cached chunks the script
  no longer uses.
- **Audio.** A two-note chime opens and closes the episode, with 0.7 s between segments. Each
  chunk is levelled to about −20 dBFS, peaking at −1 dBFS or lower. The file is `podcast.m4a`
  (AAC, via afconvert or ffmpeg), or `podcast.wav` when neither encoder is installed.
- **Link** puts the Listen card at the top of the report's reading column (or at an explicit
  `<!-- podcast:slot -->`). It shows a player, the title, the length, the disclosure and the
  transcript link. When the report is printed, the player is hidden and the file names are
  printed instead. Link also adds a line under the Markdown report's title, an entry in
  `README.html` / `README.md` and a section in `AGENTS.md`.
  - Everything link adds sits between `<!-- podcast:start -->` and `<!-- podcast:end -->` and
    is replaced on a re-run.
  - A Core run's report ends with the Core's "How did we do?" survey line. It is written when
    the report is made, usually before the podcast. Link rewords it to ask about the podcast
    too (`core_submission.refresh_feedback`, the one wording) in the HTML and the `.md`,
    before the PDF reprint. It never adds one to a report that has none.
  - `make_analysis_html.py` and `session.py finalize` keep the card and the entries on their
    own whenever `podcast/podcast.json` exists. The podcast can never stop the report of
    record: a `podcast.json` that is unreadable or has a wrong field is reported once as a
    `[WARN]` naming the file and the reason, and the report is made without the card; `link`
    says "podcast.json exists but is invalid: …".
  - If an `Analysis_Report.pdf` is older than the HTML, link reprints it with
    `html_to_pdf.print_report`, so the PDF shows the card's print text (the audio's file
    name). When it cannot reprint (no browser, say), the old PDF is renamed
    `Analysis_Report.stale.pdf` -- kept, so it never passes for current -- and link reports it
    as `[SKIPPED]` with the reason. Without `html_to_pdf.py` at all, an `[INFO]` line says how
    to reprint it by hand.

## Sending the report with the audio (`share`)

The Listen card's links are relative (`podcast/podcast.m4a`, `podcast/transcript.html`), so
`Analysis_Report.html` sent by itself has no audio. `share` writes
`output/Analysis_Report_with_audio.html`: the same report, with the card built by the same code
in its share mode. The audio is inside the file as a `data:audio/mp4` URI, the transcript is
inside it as a "Read the transcript" panel (an `iframe` whose `srcdoc` holds `transcript.html`
without its player or links), and nothing points outside the file. This is the file to send
when someone asks for "the report with the audio".

- **Audio.** Re-encoded for speech to AAC, 24 kHz mono, `--kbps` 32 to 64 (default 48), from a
  24 kHz WAV, by afconvert and otherwise ffmpeg; with neither, `podcast.m4a` goes in as rendered
  (larger). A 21-minute episode at 48 kbps is about 7.5 MB of audio, 10 MB as base64: the
  PROT_0756 v2 report came to 20.0 MB.
- **Checked before it is kept.** The written file is read back. The base64 must decode to
  exactly the audio that was put in, and that audio must be an MP4 with an AAC (`mp4a`) track of
  the episode's length (`podcast.json` `duration_s`, within 1 s or 1%). afconvert or ffmpeg must
  decode it back to that length, and it must have sound in it. There must also be no reference
  to `podcast/` anywhere. Then, when a Chrome, Chromium or Edge is found (the browser
  `html_to_pdf.py` prints with), a headless one opens the file and plays the audio for 1.5 s,
  muted, over the DevTools pipe with a 60 s deadline. It also counts the transcript's turns and
  any file the page fetched. (A `--virtual-time-budget` probe hangs on a `data:` URI audio, so
  it is not used.) Any failure removes the file and exits 1. `--no-browser` skips the playback.
- **Size.** An email attachment travels as base64, about 4/3 larger, so it prints both sizes.
  Above ~18 MB of file (~25 MB attached) it warns that many email systems will refuse it and
  suggests Bioshare (`core_submission.py deliver` puts it in the share), or `--kbps 32` with
  the sizes that would give.
- **Its survey line.** On a Core run the report's "How did we do?" line comes along, with
  `src=share` (the copy made to be forwarded). README and AGENTS.md name the file only while
  it is current: a stale copy is never advertised.
- **Refuses** while the podcast's check does not hold, like `link` (`--unchecked` overrides and
  is recorded in the file). Other links in the report that point at files beside it are listed
  as a `[WARN]`: they will not open when the file travels alone.
- **Current or not.** The file carries `<!-- podcast:share report=… podcast=… audio=… -->`, the
  sha256 of `Analysis_Report.html`, `podcast.json` and the audio it was built from, so a copy is
  checked by content, not by file times. `session.py finalize` leaves a current one alone and
  builds a missing or out-of-date one (`[OK]`/`[SKIPPED]` in MANIFEST.txt). If it cannot, it
  renames the old one `Analysis_Report_with_audio.stale.html`, which is scratch and never
  zipped or sent. `deliver` copies the file only while it is current.
- **Print.** The printed card hides the player and the transcript panel, and names the file the
  audio is in. Print `Analysis_Report.html` for paper.

## What ends up in `output/podcast/`

`podcast.m4a`, `podcast_script.md` (the script with its claims ledger), `transcript.html`
(disclosure, player, transcript, claims, how it was made -- the pronunciation table is in
`check.txt`, not here), `podcast.json` (show, title, hosts and voices, TTS backend and exact
model, consent and, for a Core run, `core_run` (`submission` with its PROT number, or
`internal`), script sha256, sources as relative paths or labelled file names with their
sha256, words, duration, `ai_generated: true`, and the `verify` block), `check.txt`, and after
verify `verify.txt` and `verify_transcript.txt`. For a Core run, who made it (their login) is in
the session's staff-only `logs/podcast.staff.json`, never delivered or zipped. A podcast made on
HIVE leaves its job script and log in `~/.proteomics-pipeline/jobs/` of whoever submitted it.

`.cache/` is scratch. So is any `*.part` file, and `podcast.wav` when `podcast.m4a` sits beside it
(`render --keep-wav`). One rule (`scripts/scratch_files.py`) keeps all three out of the session
zip and `OUTPUT_FILES.md`; never deposit them. The Core run-registry record copies what the
Listen card points at: the audio, the transcript, the script, `check.txt` and `podcast.json`.

## Known limitations

- **You cannot hear the audio.** `verify` hears it back through an ASR and catches dropped
  audio and misread numbers. It does not catch pronunciation, pacing or tone, and the ASR can
  mishear, so ask for a person to listen before it is shared, and say so when you hand it
  over. A Gemini render also warns when a chunk's audio is far shorter or longer than its text
  implies.
- **Synthetic voices.** There are no real interruptions or laughter. Gemini voices a whole
  segment at once, so it flows better than one call per turn.
- **Figures only through the report**, or through the PNGs you looked at.
- **The check is mechanical.** It catches invented tokens -- numbers, symbols, capitalised names,
  quantities in words -- plus a missing disclosure, a claimed specialty, forbidden names and
  pronunciation rows that add content. It cannot catch the list under "What check cannot
  catch", so your own read against the report is the other half.

`make_podcast.py verify --model <model>` (repeatable; `render --verify-model` for the verify a
render runs) picks the Gemini text model(s) that transcribe the audio (the default is in
`make_podcast.py`).
