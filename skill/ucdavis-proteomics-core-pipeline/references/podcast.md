# Optional: an audio discussion of the results ("Signal to Noise")

A two-host, roughly 20-minute audio episode in which a cell biologist and a statistician talk
through ONE finished analysis. It sits beside the report, and a "Listen" card near the top of
`Analysis_Report.html` links to it. **You write the script** from the report you have already
read. `scripts/make_podcast.py` then checks it against the report, turns it into audio, and
links it in. It never writes a word of the conversation.

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

## When to offer it

- **Only after the report is final**: `Analysis_Report.html` is written, and any expert-review
  or data-quality changes are in. A podcast of a draft goes stale.
- **When you deliver the results to a collaborator**, offer it in one line. Never make one by
  default. For example: "Want a ~20-min audio discussion of these results for the lab
  (AI-generated, synthetic voices), linked at the top of the report? It also explains how the
  measurement and the statistics work, using their data." Make it only on a yes. Cloud voices
  need consent too (next section).
- If the analysis is preliminary (a v1 with a v2 planned), say so in the disclosure turn.

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

- For **unpublished client data**, use a **paid-tier key** (billing enabled on the key's
  Cloud project), or **`--tts say`**, which runs offline on macOS and sends nothing.
- A free-tier key is for published or synthetic data only.
- **Who approves: the Core director.** Brett ruled on 2026-09-28 that the Core director's
  approval is enough to send a collaborator's transcript and audio to the TTS and ASR service;
  the collaborator's lab is not asked separately. Ask the director, say which tier the key is
  on, and record the answer: `--cloud-ok "Brett (Core director), 2026-09-28, free-tier key"`.
  It is written to `podcast.json` as `cloud_tts_consent`. Without `--cloud-ok`, the gemini
  backend sends nothing and exits, explaining why. The free-tier terms above still apply:
  Google may use the text and audio, and people may read or hear it.
- Key: `GEMINI_API_KEY`, or `~/.config/ucdavis-proteomics/gemini_key` (`chmod 600`; the tool
  warns if others can read it), beside the skill's other per-user settings. Never paste a key
  into the chat. The tool never prints the key and scrubs it from every error with the skill's
  one list of secret patterns (`notify_slack.redact`).
- `--cloud-ok` records consent; `--cloud-ok no` (or `false`, `0`, `none`, `declined`) is a
  refusal, and nothing is sent.

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
    - what to validate first.

    **Never tier by PropObs.** It is the observed fraction of a protein's precursors over ALL
    runs, so a protein found only in its own group -- exactly the hit a pulldown is looking for
    -- scores low. `check` fails "tier/rank/sort by PropObs" advice.
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

### What check cannot catch

Read the script against the report for these; a PASS says nothing about them:

- small integers: an unsigned count of 10 or less ("4 baits") is not checked;
- context: a real number or protein attached to the wrong protein, contrast, group or figure;
- relational words: "higher", "more than", "most", "only", "the top hit";
- false statements built from true numbers and true names;
- lowercase symbols ("gapdh") and lowercase respellings in the Pronunciation table;
- a capitalised word that starts a sentence ("Gapdh went up."): the INFO line lists them;
- numbers with a unit other than %, fold, ×, k, M and B ("3 kDa", "2 µg");
- meaning: a caveat the report makes can be dropped, and speculation can be spoken as fact.

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

```bash
S=<session>; P=$S/output/podcast
# 1. write $P/podcast_script.md (above), then check it until it passes
python3 scripts/make_podcast.py check $P/podcast_script.md \
    --source $S/output/AI_Analysis_Report.md $S/output/methods.md $S/output/AUDIT.md \
             $S/output/SAMPLE_QUALITY.md \
    --forbid-name "<PI full name>" "<submitter full name>" "<first names>"
#    -> $P/check.txt. Fix every FAIL line, re-run until PASS. Read the WARN lines.
# 2. preview what will be spoken (sends nothing)
python3 scripts/make_podcast.py render $P/podcast_script.md --tts say --dry-run | less
# 3. render: offline ...
python3 scripts/make_podcast.py render $P/podcast_script.md --tts say
#    ... or Gemini voices, with the user's consent recorded (verify runs at the end)
python3 scripts/make_podcast.py render $P/podcast_script.md --tts gemini \
    --cloud-ok "<who agreed>, <date>, <paid|free> key"
# 4. verify: the ASR round trip (automatic after a consented gemini render; by hand after say)
python3 scripts/make_podcast.py verify $P/podcast_script.md --cloud-ok "<who agreed>, <date>"
#    -> $P/verify.txt. Listen to every segment it names; fix, re-render, verify again.
# 5. link it into the report, README and AGENTS.md (safe to re-run)
python3 scripts/make_podcast.py link $S/output
```

- **Sources.** Pass the delivered report text you wrote from, all from the session's `output/`
  folder: `AI_Analysis_Report.md` (or the report's `Analysis_Report.md` twin, or
  `Analysis_Report.html`, whose embedded images and scripts are ignored), `methods.md` (the
  teaching segment's instrument, search and FDR numbers come from it), `AUDIT.md` and
  `SAMPLE_QUALITY.md`. A `--source` outside `output/` (or inside `output/podcast/`) fails;
  give it as `--extra-source <file> --label "<what it is>"` instead (see "Discuss only what the
  lab was given"). `check` finds `output/` as the folder above the script's `podcast/` folder,
  or from `--output-dir`. Do **not** pass the DE tables: with thousands of numbers in the
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
- **Gemini** needs `google-genai` (`python3 -m pip install google-genai`); the tool says so if
  it is missing.
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
  - It runs at the end of every `render --tts gemini` given `--cloud-ok` (`--no-verify`
    skips it). It never fails the render and never touches the report. It **sends the
    audio** to Google, so the same consent rules apply: nothing without `--cloud-ok`.
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
  - `make_analysis_html.py` and `session.py finalize` keep the card and the entries on their
    own whenever `podcast/podcast.json` exists. The podcast can never stop the report of
    record: a `podcast.json` that is unreadable or has a wrong field is reported once as a
    `[WARN]` naming the file and the reason, and the report is made without the card; `link`
    says "podcast.json exists but is invalid: …".
  - If an `Analysis_Report.pdf` is older than the HTML, link reprints it with `html_to_pdf.py`
    where the skill has it, so the PDF shows the card's print text (the audio's file name).
    Otherwise it prints an `[INFO]` line saying how to reprint it by hand.

## What ends up in `output/podcast/`

`podcast.m4a`, `podcast_script.md` (the script with its claims ledger), `transcript.html`
(disclosure, player, transcript, claims, how it was made -- the pronunciation table is in
`check.txt`, not here), `podcast.json` (show, title, hosts and voices, TTS backend and exact
model, consent, script sha256, sources as relative paths or labelled file names with their
sha256, words, duration, `ai_generated: true`, and the `verify` block), `check.txt`, and after
verify `verify.txt` and `verify_transcript.txt`.

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
