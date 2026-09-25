# Optional: an audio discussion of the results ("Signal to Noise")

A two-host, roughly 20-minute audio episode in which a cell biologist and a statistician talk
through ONE finished analysis. It sits beside the report, and a "Listen" card near the top of
`Analysis_Report.html` links to it. **You write the script** from the report you have already
read. `scripts/make_podcast.py` then checks it against the report, turns it into audio, and
links it in. It never writes a word of the conversation.

**Who it is for: the collaborator who submitted the samples.** The episode has three jobs:

1. **What their data say**: the findings, how confident to be in each, and the caveats.
2. **How their samples were measured and analysed**: the instrument, the acquisition, the
   search, the database and the DE model, in their own study's terms.
3. **How proteomics works, taught with their data as the examples**, so they can read their
   own report. A biologist who has never run a mass spectrometer should finish the episode
   knowing what a precursor, a protein group and a 1% FDR are.

Why it is built this way: when a model wrote the script from the report, it passed a number
check and still made things up ("limpa is our Core's custom extension to limma"; "Kcnb2, the
Kv2.1 beta subunit"). It also ran 30 minutes and read like a lecture. You have read the report,
the tables and the figures, so you write the script. The tool then checks every number and
every gene-like symbol against the report.

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

Only the **transcript** ever leaves the machine: the turns, after pronunciation substitutions,
plus one line of voice direction. The report, tables and figures are never sent. Even so, the
transcript describes unpublished client data.

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
- **Ask the user**, and tell them which tier the key is on. Record the answer:
  `--cloud-ok "Brett, 2026-09-25, paid key"`. It is written to `podcast.json` as
  `cloud_tts_consent`. Without `--cloud-ok`, the gemini backend sends nothing and exits,
  explaining why.
- Key: `GEMINI_API_KEY`, or `~/.config/podcast/gemini_key` (`chmod 600`; the tool warns if
  others can read it). Never paste a key into the chat. The tool never prints the key and scrubs
  it from every error.

## The hosts

- **Maya, cell biologist.** Vivid, reaches for metaphors, and drives the biology ("so what does
  this mean in a neuron?").
- **Leo, statistician.** Dry and skeptical. He loves QC, design and variance, and he deflates
  over-reading ("what's the denominator?").

The tension carries the episode. Maya gets excited, Leo checks the evidence, and they converge
on what the data support. Give each their own voice and running bits: recurring metaphors and
callbacks are good. The default show name is **Signal to Noise**. Keep the split roughly even:
`check` warns when either host has under 35% of the words.

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
   - detected versus inferred values (the detection-probability model, PropObs);
   - empirical Bayes (why n = 3 is workable);
   - multiple testing (adjusted p).

   Teach the ones that apply to this analysis and skip the rest (for example, MBR for a DDA
   run without it). The general explanations are claims beyond the report, so list them there.
   `check` warns about any of these topics the transcript never mentions.
5. **Did the experiment work?** The positive controls: the bait tops its own pulldown, known
   partners appear, and the QC shows the runs are sound.
6. **The findings, in the report's own order**, with the report's numbers and how confident to
   be in each.
7. **The main caveat or confound** the report raises.
8. **A data-quality "detective story"**, if the report has one (a contaminant, a thin run, a
   batch).
9. **A stats "nerd moment"**: one idea from segment 4 taken deeper on a real hit, for example
   an inferred value behind a big fold change, empirical Bayes with n = 3, or blocking.
10. **Lightning round**: one follow-up experiment each.
11. **What to do with this, then the sign-off** (short):
    - which files to open: `Analysis_Report.html` first, then the `DE_*.csv` tables;
    - how to tier hits by PropObs (measured in most samples = solid; mostly inferred = follow
      up);
    - what to validate first.

    `check` warns when the last two segments never point at a file, PropObs or validation.

Merge neighbours to stay within 6–10 segments. For example, fold the lightning round into the
close, or the detective story into the caveat.

## Fidelity rules (the check enforces the mechanical half)

- **Every number and data claim comes from the report or the other sources you pass to
  `check`.** Quote the report's precision:
  - An integer must match exactly. "29%" does not match a report's "29.4%".
  - A decimal may be rounded to fewer places: 10.6 matches 10.62.
  - A p-value keeps its power of ten. Say "5.65 times 10 to the minus 10", not "5.65". A bare
    power of ten, such as "10 to the minus 15" for 5.2e-15, is accepted as an order of
    magnitude.
- **Numbers stay as digits in the transcript** (215, 6,112, 5e-15, −1.3). The pronunciation
  step decides how they are spoken. A spelled-out number over ten ("twenty", "a hundred",
  "thousand") fails the check because it cannot be verified.
- **Never state a protein's function from memory as if the data showed it.** General-knowledge
  colour is allowed only when it is listed under `## Claims beyond the report`. That covers
  background biology, textbook mechanisms, arithmetic glosses ("about 1,000-fold"),
  metaphors, and each host's speculative follow-up. `check` accepts a symbol or number that is
  not in the sources only if it appears there.
- **Gene and protein symbols** (Jph3, Kv2.1, FKBP12.6, C1qa, IgG …) must be spelled as the
  report spells them (case does not matter) or be listed as claims. Do not write ordinary words
  in ALL CAPS for emphasis: `check` reads them as symbols. Use *italics*.
- **Figures**: describe only what the report says about a figure, or what you can see in its
  PNG. Never invent a feature of a plot.
- **Names**: no person's name except the hosts. The lab name ("the Dickson lab") is the most you
  may use. Pass the people's full and first names to `check --forbid-name`: the PI, the
  submitter, lab members and Core staff.
- `check` compares tokens, not meaning. A sentence built from real numbers can still say
  something the report does not. Read your transcript against the report once more before
  rendering.

## The script file (`<session>/output/podcast/podcast_script.md`)

```markdown
# Signal to Noise — <episode title>

- **Show:** Signal to Noise
- **Title:** <episode title, shown on the Listen card>
- **Hosts:** Maya (cell biologist), Leo (statistician)
- **Voices:** Maya=Kore, Leo=Charon          (optional; Gemini prebuilt voices, these are the defaults)
- **Say voices:** Maya=Samantha, Leo=Daniel  (optional; macOS `say`, these are the defaults)
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
#    ... or Gemini voices, with the user's consent recorded
python3 scripts/make_podcast.py render $P/podcast_script.md --tts gemini \
    --cloud-ok "<who agreed>, <date>, <paid|free> key"
# 4. link it into the report, README and AGENTS.md (safe to re-run)
python3 scripts/make_podcast.py link $S/output
```

- **Sources.** Pass the report text you wrote from: `AI_Analysis_Report.md` (or the report's
  `Analysis_Report.md` twin, or `Analysis_Report.html`, whose embedded images and scripts are
  ignored), `methods.md` (the teaching segment's instrument, search and FDR numbers come from
  it), `AUDIT.md`, `SAMPLE_QUALITY.md`, and any other document you drew on. Do **not**
  pass the DE tables. With thousands of numbers in the sources, almost any number would match
  something; `check` warns when that happens.
- **Render refuses** unless `check.txt` says PASS for this exact script (its sha256).
  `--unchecked` overrides this and is recorded in `podcast.json`; use it only when the user
  asks.
- **Gemini** needs `google-genai` (`python3 -m pip install google-genai`); the tool says so if
  it is missing.
  - Model names are read from the API: a pro TTS model first, then
    `gemini-2.5-flash-preview-tts`, then any other TTS model.
  - A model that is missing, out of quota, or throttled twice before its first chunk hands
    over to the next.
  - Once a model has made a chunk, the episode stays on it, because a second model's voices
    sound different. If that model runs out of quota partway through, render stops, names the
    chunk and keeps the cache. Re-run later to resume, or pass `--model <other>` to re-render
    every chunk with that model.
  - An episode is one request per segment (a segment over 500 words is split at turn
    boundaries), so about 8–12 requests. Google does not publish the free-tier TTS limits.
    They are per Cloud project and shown at <https://aistudio.google.com/rate-limit>. A small
    daily quota can stop a render partway; it resumes on the next run.
- **say** (macOS) uses one voice per host (Samantha and Daniel by default), runs offline and
  needs no key. A 22-minute episode rendered in about 2 minutes.
- **Re-runs are cheap.** Every chunk is cached as `podcast/.cache/<sha256>.wav`, keyed on its
  spoken text, voices and model, so editing one line re-makes one chunk. After an edit,
  re-run `check`, then `render`.
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
    own whenever `podcast/podcast.json` exists.
  - If an `Analysis_Report.pdf` is older than the HTML, link reprints it with `html_to_pdf.py`
    where the skill has it, so the PDF shows the card's print text (the audio's file name).
    Otherwise it prints an `[INFO]` line saying how to reprint it by hand.

## What ends up in `output/podcast/`

`podcast.m4a`, `podcast_script.md` (the script with its claims ledger), `transcript.html`
(disclosure, player, transcript, claims, how it was made), `podcast.json` (show, title, hosts
and voices, TTS backend and exact model, consent, script sha256, sources with their sha256,
words, duration, `ai_generated: true`) and `check.txt`.

`.cache/` is scratch. The session zip leaves it out; never deposit it.

## Known limitations

- **You cannot hear the audio.** Pronunciation, pacing and any clipped turn are unverified until
  a person listens. Say so when you hand it over, and ask for a listen before it is shared.
  A Gemini render warns when a chunk's audio is far shorter or longer than its text implies.
- **Synthetic voices.** There are no real interruptions or laughter. Gemini voices a whole
  segment at once, so it flows better than one call per turn.
- **Figures only through the report**, or through the PNGs you looked at.
- **The check is mechanical.** It catches invented numbers and symbols, spelled-out numbers, a
  missing disclosure and forbidden names. It does not catch a wrong sentence built from real
  numbers, so your own read against the report is the other half.
