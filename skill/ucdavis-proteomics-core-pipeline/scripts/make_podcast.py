#!/usr/bin/env python3
"""
make_podcast.py -- an audio discussion of a finished report: two synthetic hosts, a cell
biologist and a statistician, talk through the results ("Signal to Noise"). Made by default for
a Core run (--core-run), on request for anyone else (references/podcast.md).

The agent running the skill WRITES the script itself, from the finished report it has read,
following references/podcast.md. This script never writes a word of the conversation: it
CHECKS the script, RENDERS it to audio and LINKS the audio into the outputs.

  check  SCRIPT.md --source FILE [FILE ...] [--forbid-name NAME ...]
         -> check.txt beside the script; exit 1 on any FAIL. Token by token, not meaning:
            numbers (only unsigned counts of 10 or less are exempt; sign kept), symbol-like
            words and capitalised words mid-sentence must be in the sources or listed under
            "Claims beyond the report"; quantities in words fail; Pronunciation rows may only
            re-spell; the AI disclosure is in the first 3 turns; no host claims a specialty; no
            forbidden name appears anywhere shown or sent. check.txt lists what it cannot catch.
  meaning SCRIPT.md --brief | --verdict PASS|FAIL --reviewer WHO --findings FILE
         -> the meaning check, by a reader who did not write the script: every sentence against
            the report for (a) the right comparison, group and direction, (b) the report's hedges
            kept, (c) robust status of every named protein, (d) no superlative the report does not
            make, (e) no claim outside the report, (f) no draft history. --brief first (the
            reviewer works from it); --verdict, after a passing check, writes meaning_check.md
            with THIS script's sha256 and its sources' and robustness tables'; render, link,
            share and job refuse without its PASS while all of them hold.
  render SCRIPT.md [--out DIR] --tts gemini|say [--cloud-ok [NOTE] | --core-run [internal]]
         [--model M] [--keep-wav]
         -> DIR/podcast.m4a (podcast.wav when neither afconvert nor ffmpeg is present),
            transcript.html, podcast.json, podcast_script.md. Refuses unless check.txt says PASS
            for this exact script and unchanged sources (--unchecked overrides; podcast.json
            then says so). Every chunk is cached as DIR/.cache/<sha256>.wav, so a re-run resumes
            and editing one line re-synthesizes only its chunk; unused chunks are pruned after
            a successful render. --dry-run prints what would be spoken and exits.
  verify SCRIPT.md (--cloud-ok NOTE | --core-run [internal]) [--out DIR] [--model M]
         -> the ASR round trip: Gemini transcribes the rendered audio (16 kHz AAC 32 kbps) and
            difflib compares it with the spoken script. verify.txt, verify_transcript.txt and a
            `verify` block in podcast.json: word match ratio, spans of 6+ words not heard,
            numbers not heard (with transcript context), segments to listen to. Runs by itself
            after a consented --tts gemini render (--no-verify skips); never fails the render.
  link   OUTDIR
         -> a "Listen" card near the top of OUTDIR/Analysis_Report.html, a line near the top of
            the Markdown report, an entry in README.html / README.md and AGENTS.md. Idempotent:
            the card sits between <!-- podcast:start --> and <!-- podcast:end --> and is replaced
            on a re-run. A file that is not there is skipped with [INFO]. Refuses while the
            check does not hold (--unchecked overrides); rewords a Core run's "How did we do?"
            survey line to ask about the podcast too (never adds one); reprints an older
            Analysis_Report.pdf.
            The hooks the report calls never raise: a bad podcast.json is one [WARN].
  share  OUTDIR [--kbps 48] [--unchecked] [--no-browser]
         -> OUTDIR/Analysis_Report_with_audio.html: the report with the audio (AAC, 32-64 kbps
            mono, as a data: URI) and the transcript (an iframe srcdoc) built in -- one file to
            send by itself; nothing in it points outside it. The same Listen card code as link.
            Checked before it is kept: the embedded bytes decode back from the written file to
            an MP4 audio track of the episode's length (its boxes, then afconvert/ffmpeg), and a
            headless Chromium plays it when one is here. Warns above ~18 MB (~25 MB attached:
            email sends base64) with Bioshare as the way to send it. Refuses while the check
            does not hold, like link.
  job    SCRIPT.md (--core-run [internal] | --cloud-ok NOTE) [--model M] [--python PY]
         -> ~/.proteomics-pipeline/jobs/podcast_<name>_<hash>.sh, a SLURM job (1 CPU, 2 GB, 2 h)
            that runs render --tts gemini (verify follows), link and share -- the route on HIVE,
            where minutes of network calls must not run on the login node. Written only when
            what the job needs is here: a passing check, consent, a python with google-genai >=
            GENAI_MIN (setup.sh installs it) and a key; in a folder only you can change (a
            group-writable session folder would let another member swap the script). Loads
            HIVE's ffmpeg module for the AAC, and says so in its log when it cannot.
  key    [--core-run]
         -> where the Gemini key would come from, and whether the Core's shared key can be
            used; never the key. Exit 0 when there is one.

Privacy: only the final transcript (render: the turns, after pronunciation substitutions) and
the rendered audio (verify, downsampled) ever leave the machine -- never the report -- and only
with consent, recorded in podcast.json: --cloud-ok NOTE (the person whose data it is agreed), or
--core-run for a Core run, which records CORE_CLOUD_OK (the Core director's standing approval;
the one place that text lives). --core-run needs a CoreOmics submission (or, with `internal`,
session.json saying internal Core work) and a person on the Core staff list; podcast.json says
"Core staff", and their login goes to the staff-only logs/podcast.staff.json. --tts say (macOS) is offline. The Gemini key is read from
GEMINI_API_KEY, else ~/.config/ucdavis-proteomics/gemini_key (on HIVE refused when others can
read it), else -- for a Core run only -- the
Core's shared key on HIVE (CORE_KEY_FILE, used only while the Core's trust rule holds and other
accounts cannot read it). It is never printed or logged and is scrubbed from every error message
(notify_slack.redact: the skill's one list of secret patterns).

Stdlib only at import time: google-genai is imported by the gemini backend alone. No numpy and
no ffmpeg requirement -- audio is assembled with array/wave; afconvert (macOS) or ffmpeg, when
present, encodes the AAC.
"""
import argparse
import array
import base64
import binascii
import csv
import datetime
import decimal
import difflib
import hashlib
import html
import io
import json
import math
import os
import re
import select
import shlex
import shutil
import stat
import struct
import subprocess
import sys
import tempfile
import threading
import time
import traceback
import urllib.parse
import wave

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
from notify_slack import REDACTED, redact as _redact    # noqa: E402  the ONE secret list
from report_files import SHARE_NAME, SHARE_STALE    # noqa: E402  names defined once
from client_text import draft_history, draft_history_loose  # noqa: E402  the ONE list

SHOW = "Signal to Noise"
HOSTS = (("Maya", "cell biologist"), ("Leo", "statistician"))
GEMINI_VOICES = ("Kore", "Charon")                 # first host, second host
SAY_VOICES = ("Samantha", "Daniel")                # installed on macOS by default (say -v '?')
FALLBACK_TTS_MODEL = "gemini-2.5-flash-preview-tts"
KNOWN_TTS_MODELS = ("gemini-2.5-pro-preview-tts", "gemini-3.8-flash-tts", FALLBACK_TTS_MODEL)
HOST_STYLES = ("curious and energetic", "calm, precise and dryly funny")   # first, second host
MIN_INTERVAL = 20.0                                 # s between requests to one Gemini model
RATE_WAITS = 10                                     # per-minute 429 waits allowed per chunk
RATE_WAIT_BUDGET_MIN = 30.0                         # rate-limit waits allowed per render, minutes
# One TTS request's answer is waited for at most this long, then asked once more, then the
# render stops with the resume message: a request once hung for 20+ minutes with no log line,
# holding the job until a person cancelled it (2026-10; a chunk takes about 3 minutes).
REQUEST_TIMEOUT = 600.0                             # s
REQUEST_TIMEOUTS = 2                                # unanswered requests allowed per chunk
HEARTBEAT = 60.0                                    # s between "still waiting" log lines
# ...and all the time lost to unanswered requests in one render is capped too, like the rate-limit
# waits: eleven chunks each hanging once would otherwise reach the job's time limit (review of
# 2.13.1), where SLURM's TIMEOUT says less than the resume message does
REQUEST_TIMEOUT_BUDGET_MIN = 30.0                   # minutes
# A model that stops answering: after REQUEST_TIMEOUTS unanswered requests on one chunk the render
# moves to the next candidate model and makes EVERY chunk again with it (one episode never mixes
# two models' voices). And before the first chunk it must make, each candidate is asked a two-line
# test first and the first that answers is used. A staff render (2026-10-16) sat 28 minutes on a
# stalled gemini-3.8-flash-tts while gemini-2.5-flash-preview-tts answered the same test in 2.9 s;
# the episode's first chunks were cached from the stalled model, so it could never hand over.
PREFLIGHT_TIMEOUT = 60.0                            # s for the two-line test
# A model that stopped answering is remembered in the cache (STALL_FILE) for STALL_MEMORY_H
# hours: a resume then never pins it again, so it cannot repeat the same timeouts (review of
# 2.13.2: every resume went back to the stalled model, waited, and stopped).
STALL_FILE = "stalled.json"
STALL_MEMORY_H = 12.0
# What Google's answers said about the key's billing tier (GeminiTTS.note_tier), kept in the
# cache for TIER_MEMORY_DAYS under a fingerprint of the key (never the key): a 402 stops the
# render, so the render that records podcast.json is a later one (staff report, skill 2.11.2).
TIER_FILE = "key_tier.json"
TIER_MEMORY_DAYS = 7.0
PREFLIGHT_LINES = ("This is a short sound check.", "And this is the second voice.")
TERMS_URL = "https://ai.google.dev/gemini-api/terms"
KEY_FILE = os.path.join("~", ".config", "ucdavis-proteomics", "gemini_key")
# The Core's shared key on HIVE (Brett installs it: owner brettsp, group proteomics-grp, mode
# 640). Read last, for a Core run only (--core-run), and only while slack_collab.trust_problem
# holds (a Core admin owns it and its folder; nobody else can write either) and other accounts
# cannot read it. Tests replace this constant; nothing in the environment can.
CORE_KEY_FILE = "/quobyte/proteomics-grp/.config/gemini_key"
# THE consent of record for a Core run (--core-run): Brett ruled on 2026-09-28 that the Core
# director's approval is enough to send a Core run's transcript and audio to Gemini, and made
# the podcast the default for Core runs on 2026-10-15. Docs point here; never retype it.
CORE_CLOUD_OK = ("Brett Phinney (Core director): standing approval for Core runs (2026-09-28; "
                 "podcast by default since 2026-10-15)")
CORE_RUN_KINDS = ("submission", "internal")
VERIFY_FILES = ("verify.txt", "verify_transcript.txt")      # what verify writes beside the audio
# Gemini 3 TTS needs client.interactions (2.25.0 has it). setup.sh installs at least this into
# the pipeline env (its GENAI_MIN; tests/test_make_podcast_hive.py keeps the two equal).
GENAI_MIN = "2.25.0"
# The podcast job (job): one network-bound process; an episode renders in minutes, a rate-limited
# one in up to an hour. HIVE's ffmpeg (verified 2026-10-15, AAC encoder present) is a module, so
# the job is a login shell (#!/bin/bash -l), where `module` exists.
JOB_HOURS, JOB_MEM = 2, "2G"
HIVE_FFMPEG_MODULE = "ffmpeg/7.0.2"

RATE = 24000                                        # Hz, 16-bit mono throughout
WPM = 150                                           # spoken words per minute, for estimates
WORDS_LO, WORDS_HI = 2200, 3600
SEGMENTS_LO, SEGMENTS_HI = 6, 10
CHUNK_WORDS = 500                                   # a segment longer than this is split
GAP_SEGMENT, GAP_CHUNK, GAP_TURN = 0.7, 0.3, 0.25   # seconds of silence

START, END = "<!-- podcast:start -->", "<!-- podcast:end -->"
T_START, T_END = "<!-- TRANSCRIPT START -->", "<!-- TRANSCRIPT END -->"
DISCLOSURE = ("AI-generated, in synthetic voices. It discusses this report and methods; items "
              "beyond them are listed in the transcript. The report is the record — check "
              "anything important against it.")


def extras_note(man, where="folder"):
    """" It also drew on: <labels> (not in this folder)." when the check had extra sources."""
    labs = [x.get("label") or x.get("file") for x in (man.get("sources") or [])
            if isinstance(x, dict) and x.get("delivered") is False]
    labs = [str(x) for x in labs if x]
    return f" It also drew on: {'; '.join(labs)} (not in this {where})." if labs else ""
DISCLOSURE_RE = re.compile(r"\bAI[-\s]generated\b|\bgenerated\s+(?:by|with)\s+(?:an?\s+)?"
                           r"(?:AI|artificial\s+intelligence)\b", re.I)

# Written -> spoken, applied to the TTS text only (the transcript keeps the written form).
# Whole tokens only (no letter or digit either side), longest first; a script's own
# "## Pronunciation" table adds to these and wins on the same Written form.
PRONOUNCE = [
    ("Kv2.1", "K V two point one"), ("Kv2.2", "K V two point two"), ("Kv2", "K V two"),
    ("IgG", "I G G"), ("DIA-NN", "D I A N N"), ("dia-PASEF", "dia pasef"),
    ("diaPASEF", "dia pasef"), ("timsTOF", "tims toff"), ("log2FC", "log two fold change"),
    ("log2", "log two"), ("log10", "log ten"), ("adj.P.Val", "adjusted p value"),
    ("adj.P", "adjusted p"), ("m/z", "m over z"), ("1/K0", "one over K zero"), ("ER", "E R"),
    ("PropObs", "prop obs"), ("LC-MS/MS", "L C M S M S"), ("LC-MS", "L C M S"),
    ("MS/MS", "M S M S"), ("e.g.", "for example"), ("i.e.", "that is"), ("vs.", "versus"),
    ("vs", "versus"),
]

# Acronyms a conversation uses that are not claims about the data. Anything else that looks
# like a symbol must be in the sources or under "Claims beyond the report".
GENERIC_TOKENS = {"ai", "ok", "dna", "rna", "mrna", "pcr", "pdf", "html", "csv", "png", "tv",
                  "phd", "fyi", "usa"}

# The listener is the collaborator who submitted the samples; one job of the episode is to teach
# how proteomics works with their own data (references/podcast.md). check warns -- never fails --
# on a topic the transcript never touches: not every analysis has MBR or dia-PASEF.
TEACHING = [
    ("what LC-MS/MS does", r"\bLC-?MS|mass spec|chromatograph|\bpeptides?\b"),
    ("DIA / dia-PASEF", r"\bDIA\b|dia-?PASEF|data[- ]independent"),
    ("precursors vs protein groups", r"precursor|protein groups?\b"),
    ("what 1% FDR means", r"\bFDR\b|false discovery"),
    ("library-free search / match-between-runs",
     r"library[- ]free|match(?:ing)?[- ]between[- ]runs|\bMBR\b|spectral librar"),
    ("detected vs inferred values", r"\binferred\b|detection[- ]probability|PropObs"),
    ("empirical Bayes", r"empirical Bayes|borrow\w* (?:strength|information)|moderated t"),
    ("multiple testing", r"multiple[- ]testing|Benjamini|adjusted p|\bFDR\b"),
]
# The close tells them what to do next: which file to open, how to tier hits, what to validate.
NEXT_STEPS = (r"README|Analysis[_ ]Report|Detected_|detection (?:count|column)|detected in \d+ "
              r"of \d+|\bEvidence\b|tier file|\bvalidat")
# Advice the brief forbids: PropObs is the observed fraction of a protein's precursors over ALL
# runs, so a hit found only in its bait's pulldowns scores low. Tier with the per-group
# detection columns (Detected_<group> k/n, Evidence) or the tier file instead.
# PropObs however it is written or spoken: PropObs, prop obs, prop-obs, "proportion observed".
_PROPOBS = r"(?:\bprop\.?[\s_-]*obs\b|\bproportion(?:\s+of)?\s+observ\w*)"
_ADVISE = r"\b(?:tier|rank|sort|filter|triage|prioriti[sz]e)\w*"
# within one clause (no . ! ? ; :), in either order
PROPOBS_TIER = re.compile(rf"{_ADVISE}[^.!?;:]{{0,80}}{_PROPOBS}|{_PROPOBS}[^.!?;:]{{0,80}}{_ADVISE}",
                          re.I)
# "Never tier by PropObs", "rank by fold change, not PropObs": advice against it is fine.
_NEGATION = re.compile(r"\b(?:never|not|no|nor|don[’']?t|do not|doesn[’']?t|won[’']?t|"
                       r"shouldn[’']?t|avoid\w*|instead of|rather than)\b", re.I)


def propobs_advice(text):
    """Each clause that advises tiering / ranking / sorting hits by PropObs, unless a negation
    comes before it in its clause or inside it."""
    for m in PROPOBS_TIER.finditer(text):
        start = max(text.rfind(c, 0, m.start()) for c in ".!?;:") + 1
        if not _NEGATION.search(text[start:m.end()]):
            yield m.group(0)

# A quantity said in words cannot be checked against the sources: numbers over ten, their
# plurals and -fold forms, N-fold, dozen(s), twice, half. A phrase listed under "Claims beyond the
# report" is allowed (an arithmetic gloss, say).
SPELLED_NUMBER = re.compile(
    r"\b(?:(?:eleven|twelve|thirteen|fourteen|fifteen|sixteen|seventeen|eighteen|nineteen|"
    r"twenty|thirty|forty|fifty|sixty|seventy|eighty|ninety|hundred|thousand|million|billion)"
    r"(?:s|[- ]?fold)?"
    r"|(?:one|two|three|four|five|six|seven|eight|nine|ten)[- ]?fold"
    r"|(?:a\s+)?dozens?|twice|half|halves)\b", re.I)
# Little words between a quantity word and the word it counts: "half OF THE runs".
_QUANT_FILLER = {"of", "the", "a", "an", "as", "than", "more", "less", "many", "much", "to",
                 "in", "on", "at", "by", "for", "all", "our", "their", "its", "these", "those",
                 "this", "that"}


# Where a phrase ends: the end of a sentence or clause, a table cell, a line.
_BOUNDARY = re.compile(r"[.!?;:](?=\s|$)|\||\n")


def _words_norm(text):
    """Lowercase words, anything else a single space, padded, with " # " at every sentence,
    clause, table-cell or line boundary: a phrase is searched within one stretch of prose."""
    parts = (" ".join(re.findall(r"[a-z0-9]+", p)) for p in _BOUNDARY.split((text or "").lower()))
    return " " + " # ".join(p for p in parts if p) + " "


def quantity_phrase(text, m):
    """The phrase a quantity word stands in: the word, any little words after it, and the next
    content word -- "hundreds of times", "half of the runs", "thousands of proteins" -- within
    its sentence. None when the sentence ends first (nothing to anchor it)."""
    rest = re.findall(r"[a-z0-9]+", _BOUNDARY.split(text[m.end():].lower())[0])
    tail = []
    for w in rest[:6]:
        tail.append(w)
        if w not in _QUANT_FILLER:
            return " ".join(re.findall(r"[a-z0-9]+", m.group(0).lower()) + tail)
    return None


# A host must not claim a real research specialty ("I study membrane contact sites in neurons"):
# it lends the synthetic voice an authority it does not have. "I'm the biologist of the pair".
SPECIALTY = re.compile(r"\bI(?:'m| am)?\s+(?:study|studied|research|specialise|specialize|"
                       r"work on|run a lab)\b|\b(?:in\s+)?my\s+(?:lab|research|thesis|postdoc|"
                       r"PhD)\b|\bas an?\s+(?:neuroscientist|biochemist|cell biologist|"
                       r"statistician|professor|expert)\b", re.I)
# Said in check.txt, the brief and SKILL.md: the check matches tokens, not meaning.
CANNOT_CATCH = (
    "small integers: an unsigned count of 10 or less (\"4 baits\") is not checked",
    "context: a real number or protein attached to the wrong protein, contrast, group or "
    "figure passes",
    "relational words: \"higher\", \"more than\", \"most\", \"only\", \"the top hit\" are "
    "not checked",
    "false statements built from true numbers and true names",
    "lowercase symbols (\"gapdh\") and lowercase respellings in the Pronunciation table",
    "a capitalised word that starts a sentence (\"Gapdh went up.\"): check.txt lists the ones "
    "not in the sources as INFO, to read",
    "numbers with a unit or suffix other than %, fold, x, k, M and B (\"3 kDa\", \"2 µg\") are "
    "matched as plain numbers, and a count of 10 or less with one is not checked",
    "meaning: a caveat the report makes can be dropped, and speculation can be spoken as fact",
    "draft history in words client_text.py's phrase list does not have (\"we had that "
    "backwards before\")",
)
# --cloud-ok with one of these means no.
NO_CONSENT = {"", "false", "no", "n", "0", "none", "null", "off", "declined", "decline",
              "refused", "refuse"}
# Words a spoken form may use to say a digit (Pronunciation rows): each maps to the digits it
# stands for, which the written form must contain.
NUMBER_WORDS = {"zero": "0", "oh": "0", "one": "1", "two": "2", "three": "3", "four": "4",
                "five": "5", "six": "6", "seven": "7", "eight": "8", "nine": "9", "ten": "10",
                "eleven": "11", "twelve": "12", "thirteen": "13", "fourteen": "14",
                "fifteen": "15", "sixteen": "16", "seventeen": "17", "eighteen": "18",
                "nineteen": "19", "twenty": "2", "thirty": "3", "forty": "4", "fifty": "5",
                "sixty": "6", "seventy": "7", "eighty": "8", "ninety": "9", "hundred": "",
                "thousand": "", "million": ""}

_SECRETS = []
_sleep = time.sleep                                 # tests replace these two
_now = time.monotonic


# ----------------------------------------------------------------------------- small helpers
def log(msg):
    print(scrub(msg), file=sys.stderr, flush=True)


def _never_raise(default):
    """For the hooks the report of record calls (make_analysis_html.py, session_docs.py): any
    error is reported as a [WARN] and the report goes on without the podcast."""
    def wrap(fn):
        def inner(*args, **kw):
            try:
                return fn(*args, **kw)
            except Exception as e:
                log(f"[WARN] podcast left out of the report ({fn.__name__}: {type(e).__name__}: {e})")
                return default(*args) if callable(default) else default
        inner.__name__, inner.__doc__ = fn.__name__, fn.__doc__
        return inner
    return wrap


def scrub(text):
    """Every secret out of a string: the key read at run time, and every secret-shaped substring
    in notify_slack's list (Google API keys, key=..., tokens, webhooks, private keys)."""
    s = str(text)
    for k in _SECRETS:
        if k:
            s = s.replace(k, REDACTED)
    return _redact(s)


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def now_iso():
    return datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat()


def minutes_label(seconds):
    if not seconds:
        return "length not recorded"
    m = seconds / 60.0
    return "under a minute" if m < 1 else f"about {m:.0f} min"


def rel(path, base):
    return os.path.relpath(path, base).replace(os.sep, "/")


def href(path):
    return urllib.parse.quote(path, safe="/")


def esc(s):
    return html.escape(str(s if s is not None else ""))


# ----------------------------------------------------------------------------- the script
class Turn(object):
    def __init__(self, speaker, text, line, seg):
        self.speaker, self.text, self.line, self.seg = speaker, text, line, seg
        self.words = len(text.split())


class Script(object):
    """A parsed podcast script. `problems` are format errors (line, message)."""

    def turns(self):
        return [t for seg in self.segments for t in seg]


HEADER_RE = re.compile(r"^([A-Za-z][A-Za-z \-]{0,30}?)\s*:\s*(.*?)\s*$")
TURN_RE = re.compile(r"^\*\*\s*([A-Za-z][A-Za-z .'\-]{0,40}?)\s*(?::\s*\*\*|\*\*\s*:)\s*(.*?)\s*$")


def _pairs(value):
    return [(m.group(1).strip(), m.group(2).strip()) for m in
            re.finditer(r"([A-Za-z][A-Za-z .'\-]*?)\s*[=:]\s*([^,;]+)", value or "")]


def _cells(row):
    parts = re.split(r"(?<!\\)\|", row.strip().strip("|"))
    return [p.replace("\\|", "|").strip().strip("`").strip() for p in parts]


@_never_raise(lambda text, *rest: text)
def strip_block(text):
    """The text without any <!-- podcast:start --> ... <!-- podcast:end --> block (and the line
    breaks around it), so a report that carries a Listen line reads as it did before."""
    return re.sub(r"\n?[ \t]*" + re.escape(START) + r".*?" + re.escape(END) + r"[ \t]*\n?", "\n",
                  text or "", flags=re.S)


def parse_script(path):
    with open(path, "rb") as fh:
        raw = fh.read()
    text = raw.decode("utf-8", errors="replace")
    s = Script()
    s.path, s.text, s.sha256, s.problems = os.path.abspath(path), text, sha256_bytes(raw), []
    lines = text.splitlines()

    head, h1 = {}, None
    for ln in lines:
        if ln.startswith("## ") or ln.strip() == T_START:
            break
        if ln.startswith("# ") and h1 is None:
            h1 = ln[2:].strip()
            continue
        m = HEADER_RE.match(re.sub(r"^\s*[-*]\s+", "", ln).replace("**", ""))
        if m and m.group(2):
            head.setdefault(m.group(1).strip().lower(), m.group(2))
    s.header = head
    s.show = head.get("show") or SHOW
    s.title = head.get("title") or h1
    s.author = head.get("author") or head.get("script author") or head.get("written by")
    s.style = head.get("style")
    # "Maya (cell biologist), Leo (statistician)" or "MAYA (cell biologist, voice Kore) · LEO
    # (statistician, voice Charon)": a voice named in the parentheses is the Gemini voice
    hosts, host_voice = [], {}
    for part in re.split(r"[,;·|](?![^()]*\))", head.get("hosts", "")):
        m = re.match(r"\s*([A-Za-z][A-Za-z .'\-]*?)\s*(?:\((.*?)\))?\s*$", part)
        if not (m and m.group(1)):
            continue
        name = m.group(1).strip()
        name = name.title() if name.isupper() else name
        bits = [b.strip() for b in (m.group(2) or "").split(",") if b.strip()]
        for b in bits:
            v = re.match(r"(?i)voice\s*[:=]?\s*(\S+)$", b)
            if v:
                host_voice[name] = v.group(1)
        hosts.append((name, ", ".join(b for b in bits if not re.match(r"(?i)voice\b", b))))
    s.hosts = hosts if head.get("hosts") else list(HOSTS)
    if len(s.hosts) != 2:
        s.problems.append((0, f"Hosts: needs exactly two hosts, got {len(s.hosts)} "
                              f"({head.get('hosts')!r})"))
    names = {h.lower(): h for h, _ in s.hosts}

    def voices(value, defaults, named=None):
        v = {h: defaults[i] if i < len(defaults) else defaults[-1]
             for i, (h, _) in enumerate(s.hosts)}
        v.update(named or {})
        for name, voice in _pairs(value):
            if name.lower() in names:
                v[names[name.lower()]] = voice
        return v
    s.gemini_voices = voices(head.get("voices") or head.get("gemini voices"), GEMINI_VOICES,
                             host_voice)
    # "Styles: Maya: curious and energetic; Leo: calm, precise and dryly funny" -- how each host
    # sounds (Gemini 3 TTS takes it per turn; the 2.x prompt line uses it too)
    s.styles = {}
    for part in (head.get("styles") or "").split(";"):
        m = re.match(r"\s*([A-Za-z][A-Za-z .'\-]*?)\s*[=:]\s*(.+?)\s*$", part)
        if m and m.group(1).lower() in names:
            s.styles[names[m.group(1).lower()]] = m.group(2)
    s.say_voices = voices(head.get("say voices"), SAY_VOICES)

    secs = {}
    parts = re.split(r"(?m)^##\s+(.+?)\s*$", text)
    for i in range(1, len(parts), 2):
        secs[parts[i].strip().lower()] = parts[i + 1]
    s.sections = secs

    s.pronunciation, s.has_pronunciation = [], "pronunciation" in secs
    for row in (secs.get("pronunciation") or "").splitlines():
        if not row.strip().startswith("|"):
            continue
        c = _cells(row)
        if len(c) < 2 or not c[0] or c[0].lower() == "written" or re.match(r"^:?-{2,}:?$", c[0]):
            continue
        s.pronunciation.append((c[0], c[1]))

    ckey = next((k for k in secs if k.startswith("claims beyond the report")), None)
    s.has_claims = ckey is not None
    body = (secs.get(ckey) or "").split(T_START)[0] if ckey else ""   # before the transcript
    body = re.sub(r"(?s)<!--.*?-->", "", body)
    s.claims = [re.sub(r"^\s*[-*]\s+", "", ln).strip() for ln in body.splitlines()
                if re.match(r"^\s*[-*]\s+\S", ln)]
    s.claims = [c for c in s.claims if c.rstrip(".").lower() != "none"]
    if ckey and not s.claims and body.strip().strip("-* ").rstrip(".").lower() != "none":
        s.problems.append((0, "'## Claims beyond the report' is empty: list each claim as a "
                              "bullet, or write 'None'"))
    s.claims_text = "\n".join(s.claims)

    s.segments = []
    st = [i for i, ln in enumerate(lines) if ln.strip() == T_START]
    en = [i for i, ln in enumerate(lines) if ln.strip() == T_END]
    if not st or not [e for e in en if e > st[-1]]:
        s.problems.append((0, f"no transcript: put the turns between a line '{T_START}' and a "
                              f"line '{T_END}'"))
        return s
    a = st[-1]
    b = [e for e in en if e > a][0]
    cur = []
    for i in range(a + 1, b):
        ln = lines[i].strip()
        if not ln or re.match(r"^<!--.*-->$", ln):
            continue
        if ln == "---":
            if cur:
                s.segments.append(cur)
                cur = []
            continue
        m = TURN_RE.match(ln)
        if not m or not m.group(2):
            s.problems.append((i + 1, "not a turn: each turn is ONE line starting "
                                      + " or ".join(f"**{h.upper()}:**" for h, _ in s.hosts)
                                      + f" -- {ln[:70]!r}"))
            continue
        who = names.get(m.group(1).strip().lower())
        if not who:
            s.problems.append((i + 1, f"speaker {m.group(1)!r} is not one of the hosts "
                                      f"({', '.join(h for h, _ in s.hosts)})"))
            continue
        cur.append(Turn(who, m.group(2), i + 1, len(s.segments)))
    if cur:
        s.segments.append(cur)
    if not s.segments:
        s.problems.append((a + 1, "the transcript has no turns"))
    return s


# ----------------------------------------------------------------------------- pronunciation
def pronouncer(script_table=()):
    """-> speak(text): the TTS text. Built-in defaults plus the script's table (which wins),
    applied in ONE pass, longest Written form first, whole tokens only."""
    table = dict(PRONOUNCE)
    table.update(dict(script_table))
    keys = sorted((k for k in table if k), key=len, reverse=True)
    rx = re.compile(r"(?<![A-Za-z0-9])(?:" + "|".join(re.escape(k) for k in keys)
                    + r")(?![A-Za-z0-9])") if keys else None

    def speak(text):
        t = re.sub(r"`([^`]*)`", r"\1", text)
        t = re.sub(r"\*\*?([^*]+)\*\*?", r"\1", t)        # emphasis marks are not spoken
        if rx:
            t = rx.sub(lambda m: table[m.group(0)], t)
        return speak_numbers(t)
    return speak


def _sci_words(mant, sign, exp):
    lead = "" if mant in ("1", "1.0") else f"{mant} times "
    return f"{lead}10 to the {'minus ' if sign in ('-', '−') else ''}{exp}"


_SUPER = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺", "0123456789-+")


def speak_numbers(t):
    """Symbols a voice reads badly: scientific notation, minus signs, ≥ ≤ ± ~ ×."""
    t = re.sub(r"[⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+", lambda m: "^" + m.group(0).translate(_SUPER), t)
    t = re.sub(r"(\d+(?:\.\d+)?)\s*[×x*]\s*10\s*\^\s*([+\-−]?)(\d+)",
               lambda m: _sci_words(m.group(1), m.group(2), m.group(3)), t)
    t = re.sub(r"(?<![\w.])(\d+(?:\.\d+)?)[eE]([+\-−]?)(\d+)(?!\w|\.\d)",
               lambda m: _sci_words(m.group(1), m.group(2), m.group(3)), t)
    t = re.sub(r"(?<![\w.])10\s*\^\s*([+\-−]?)(\d+)",
               lambda m: _sci_words("1", m.group(1), m.group(2)), t)
    t = re.sub(r"(^|[\s(=])[−\-](?=\d)", r"\1minus ", t)
    t = t.replace("−", "-")
    t = re.sub(r"\s*≥\s*", " at least ", t)
    t = re.sub(r"\s*≤\s*", " at most ", t)
    t = re.sub(r"\s*±\s*", " plus or minus ", t)
    t = re.sub(r"(^|\s)~\s*(?=\d)", r"\1about ", t)
    t = re.sub(r"(\d)\s*×", r"\1 times", t)
    t = re.sub(r"\s<\s", " less than ", t)
    t = re.sub(r"\s>\s", " greater than ", t)
    return re.sub(r"[ \t]{2,}", " ", t).strip()


# ----------------------------------------------------------------------------- the number check
D = decimal.Decimal
_NUM = re.compile(
    r"(?<![A-Za-z0-9_.])(?:"
    r"(?P<m1>\d+(?:\.\d+)?)[eE](?P<e1>[+-]?\d+)(?!\d|\.\d)"
    r"|(?P<m2>\d+(?:\.\d+)?)\s*(?:[x×*·]|times)\s*10\s*(?:\^\s*|\*\*\s*|\s+to\s+the\s+"
    r"(?:power\s+of\s+)?)(?P<n2>minus\s+|negative\s+)?(?P<e2>[+-]?\d+)"
    r"|10\s*(?:\^\s*|\*\*\s*|\s+to\s+the\s+(?:power\s+of\s+)?)(?P<n3>minus\s+|negative\s+)?"
    r"(?P<e3>[+-]?\d+)"
    r"|(?P<p>\d+(?:\.\d+)?|\.\d+)"
    r")", re.I)


def normalize_numbers(t):
    """Unicode minus, superscript exponents and thousands separators to plain ASCII, so
    "6,112", "−1.3" and "5×10⁻¹⁵" compare as 6112, -1.3 and 5×10^-15."""
    t = (t or "").replace("−", "-").replace("‒", "-")
    t = re.sub(r"[–—]", " ", t)                    # en/em dash: a range, not a minus
    t = re.sub(r"[⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+", lambda m: "^" + m.group(0).translate(_SUPER), t)
    return re.sub(r"(?<![\d.,])\d{1,3}(?:,\d{3})+(?!\d)", lambda m: m.group(0).replace(",", ""), t)


class Num(object):
    """kind: plain | sci | pow10. value: the absolute value. sign: '-', '+' or None as written
    (a plain number only). dec: decimals written (the mantissa's, for sci). quant: followed by
    %, percent, fold or x, so a measurement even when small."""

    def __init__(self, text, value, kind, dec, sign=None, quant=False):
        self.text, self.value, self.kind, self.dec = text, value, kind, dec
        self.sign, self.quant = sign, quant

    def cls(self):
        return "-" if self.sign == "-" else "+"


def _decimals(s):
    return len(s.split(".", 1)[1]) if "." in s else 0


# "-1.3", "(−0.9)", "minus 1.3" are negative; "LRS-124" and "10-20" are not (a letter or a digit
# before the hyphen).
_SIGN_BEFORE = re.compile(r"(?:(?<![A-Za-z0-9_.])([+-])|\b(minus|negative|plus)\s+)$", re.I)
_QUANT_AFTER = re.compile(r"\s*(?:%|percent\b|-?\s*fold\b|×|x\b)", re.I)
# 3k, 6K, 2M, 1B: the suffix multiplies, and the number is a measurement, never a count.
_SUFFIX = re.compile(r"([kKMB])(?![A-Za-z0-9])")
_SUFFIX_X = {"k": 1000, "K": 1000, "M": 1000000, "B": 1000000000}


def numbers_in(text):
    out = []
    t = normalize_numbers(text)
    for m in _NUM.finditer(t):
        g = m.groupdict()
        try:
            if g["p"] is not None:
                p = g["p"]
                sm = _SIGN_BEFORE.search(t[max(0, m.start() - 16):m.start()])
                sign = None
                if sm:
                    sign = sm.group(1) or ("+" if sm.group(2).lower() == "plus" else "-")
                sx = _SUFFIX.match(t, m.end())
                if sx:
                    v = abs(D(p)) * _SUFFIX_X[sx.group(1)]
                    out.append(Num(p + sx.group(1), v, "plain", max(0, -v.normalize().as_tuple().exponent),
                                   sign, True))
                    continue
                out.append(Num((sign or "") + p if sign == "-" else p, abs(D(p)), "plain",
                               _decimals(p), sign, bool(_QUANT_AFTER.match(t, m.end()))))
                continue
            mant = g["m1"] or g["m2"] or "1"
            e = g["e1"] or g["e2"] or g["e3"]
            exp = -abs(int(e)) if (g["n2"] or g["n3"]) else int(e)
            v = abs(D(mant).scaleb(exp))
            kind = "pow10" if D(mant) == 1 else "sci"
            out.append(Num(m.group(0).strip(), v, kind, _decimals(mant)))
        except (decimal.InvalidOperation, ValueError):
            continue
    return out


def _sci(v):
    """(exponent, mantissa in [1, 10)) of a positive Decimal."""
    k = v.adjusted()
    return k, v.scaleb(-k)


def _rounded(v, d):
    q = D(1).scaleb(-d)
    return {v.quantize(q, rounding=decimal.ROUND_HALF_UP), v.quantize(q, rounding=decimal.ROUND_HALF_EVEN)}


class NumberBook(object):
    """Every number in the sources, and what it may legitimately be spoken as: itself, rounded
    to fewer decimals (never to an integer), or -- for a p-value -- its order of magnitude. A
    number spoken with a sign must match the sign in the sources (-2.68 is not +2.68); one
    spoken without a sign ("fell by 1.3") matches either."""

    def __init__(self):
        self.abs, self.signed, self.sci_rounded, self.exponents = set(), set(), {}, set()
        self.rounded, self.rounded_signed = {}, {}
        self.mantissas = {}                     # 5.65 -> "5.65e-10": for a helpful message

    def add_text(self, text):
        for n in numbers_in(text):
            self.add(n)

    def add(self, n):
        v = n.value
        self.abs.add(v)
        self.signed.add((n.cls(), v))
        if not v:
            return
        if n.kind != "plain" or v < D("0.01"):
            k, mant = _sci(v)
            self.exponents.add(k)
            self.mantissas.setdefault(mant, n.text)
            md = max(0, -mant.normalize().as_tuple().exponent)
            for d in range(0, md):
                for r in _rounded(mant, d):
                    self.sci_rounded.setdefault(d, set()).add((k, r))
        if n.kind == "plain":
            for d in range(1, n.dec):
                for r in _rounded(v, d):
                    self.rounded.setdefault(d, set()).add(r)
                    self.rounded_signed.setdefault(d, set()).add((n.cls(), r))

    def verdict(self, n):
        """-> 'trivial' | 'exact' | 'rounded' | 'magnitude' | None (not in the sources). Only an
        unsigned integer of 10 or less that is not a percentage, a fold or a multiple is
        trivial (a count in speech: "4 baits")."""
        v = n.value
        if n.kind == "plain" and n.dec == 0 and v <= 10 and n.sign is None and not n.quant:
            return "trivial"
        if (v in self.abs) if n.sign is None else ((n.cls(), v) in self.signed):
            return "exact"
        if n.kind == "plain" and n.dec >= 1:
            if ((v in self.rounded.get(n.dec, ())) if n.sign is None else
                    ((n.cls(), v) in self.rounded_signed.get(n.dec, ()))):
                return "rounded"
        if n.kind in ("sci", "pow10") and v:
            k, mant = _sci(v)
            md = max(0, -mant.normalize().as_tuple().exponent)
            if (k, mant) in self.sci_rounded.get(md, ()):
                return "rounded"
            if n.kind == "pow10" and k in self.exponents:
                return "magnitude"
        return None


# ----------------------------------------------------------------------------- symbol check
_SYM = re.compile(r"(?<![A-Za-z0-9α-ωΑ-Ω])[A-Za-zα-ωΑ-Ω][A-Za-z0-9α-ωΑ-Ω]*(?:\.\d+)*")
_TITLE = re.compile(r"^[A-Z][a-z]+$")
# A capitalised name, a hyphen and a number: IL-6, COVID-19. "IL" alone is too short to be a
# symbol and "6" alone is a trivial count, so the pair is checked whole.
_HYPHEN_SYM = re.compile(r"(?<![A-Za-z0-9\-])[A-Z][A-Za-z]{0,7}-\d+[A-Za-zα-ω]?(?![A-Za-z0-9])")


def is_symbol(tok):
    if tok.lower() in GENERIC_TOKENS:
        return False
    if re.search(r"[α-ωΑ-Ω]", tok):
        return True                                   # TNF-α, PKCδ, Gβγ
    if any(c.isdigit() for c in tok):
        return True                                   # Jph3, Kv2.1, FKBP12.6, C1qa, RYR2
    if len(tok) >= 3 and tok.isupper():
        return True                                   # SERCA, VAPA, BSA
    return bool(re.search(r"[a-z][A-Z]", tok))        # IgG, timsTOF, PropObs


def symbol_found(tok, hay):
    """Case-insensitive, as a whole token. One ending in a digit must not run on into more
    digits (Kcnb2 is not Kcnb20); one ending in a letter must not run on into more letters (SOD
    is not sodium, APP is not applied) but may take a plural or digits (RyR -> RyRs, RyR2)."""
    t = tok.lower()
    tail = r"(?![0-9])" if t[-1].isdigit() else r"(?:e?s)?(?![a-zα-ω])"
    if re.search(r"(?<![a-z0-9α-ω])" + re.escape(t) + tail, hay):
        return True
    return len(t) > 3 and t.endswith("s") and symbol_found(tok[:-1], hay)


def title_tokens(text, initial=False):
    """Capitalised words that do not start a sentence: Gapdh, Western, a person, a place. With
    initial=True, the ones that do (after . ! ? or at the start of the turn), which check lists
    but cannot judge: "Gapdh went up." and "The gel ran." look the same."""
    for m in re.finditer(r"(?<![A-Za-z0-9'’_α-ωΑ-Ω])([A-Z][a-z]+)(?![A-Za-z0-9α-ωΑ-Ω])", text):
        pre = re.sub(r"[\s\"“”‘’'(\[*_]+$", "", text[:m.start()])
        starts = not pre or pre[-1] in ".!?…"
        if starts == initial:
            yield m.group(1)


def spoken_digits(spoken):
    """The digits a spoken form says, in order: numerals as written, number words read the way
    verify reads them ("two point eight six" -> 286, "sep two fifty" -> 250, "zero seven five
    six" -> 0756, "oh" as 0)."""
    words = [[("zero" if w.lower() == "oh" else w.lower()), 0, 0]
             for w in re.findall(r"[A-Za-z]+|\d+", spoken or "")]
    return "".join(c for t in _spoken_numbers(words) for c in t[0] if c.isdigit())


def pronunciation_problems(written, spoken):
    """What a Pronunciation row must not do: the spoken form is sent to the voices and never
    checked against the report, so it may only re-spell the written form. Its digits -- numerals
    and number words, read in order -- must be the written form's digits in the same order
    (2.68 is not "two point eight six"), and a symbol or capitalised word in it must be part of
    the written form. A lowercase respelling ("teck R" for Tecr) cannot be judged."""
    probs = []
    wdig = re.sub(r"\D", "", written)
    sdig = spoken_digits(spoken)
    if sdig != wdig:
        probs.append(f"adds the number {sdig}, which the written form does not have" if not wdig
                     else f"drops the digits {wdig}" if not sdig
                     else f"says the digits {sdig}, not {wdig} in order")
    for m in _SYM.finditer(spoken):
        tok = m.group(0)
        if not (is_symbol(tok) or _TITLE.match(tok) or (len(tok) == 1 and tok.isupper())):
            continue
        if tok.lower() not in written.lower():
            probs.append(f"'{tok}' is not a spelling of {written!r}")
    return probs


# ----------------------------------------------------------------------------- sources
def source_text(path):
    """A source as check hashes it: its text without any podcast block, so `link` adding its
    Listen line to the report does not make the check stale -- and without the Core's survey
    line, which link rewords to ask about the podcast (core_submission.without_feedback)."""
    with open(path, "rb") as fh:
        text = strip_block(fh.read().decode("utf-8", errors="replace"))
    try:
        import core_submission
    except ImportError:
        return text
    return core_submission.without_feedback(text)


def source_sha(text):
    """The hash check records for a source: its source_text with whitespace runs collapsed, so
    the blank lines around a removed Listen block do not count as an edit."""
    return sha256_bytes(re.sub(r"\s+", " ", text).strip().encode("utf-8"))


def load_source(path):
    """-> (source_sha, text) with the parts that are not prose removed: embedded images (a
    base64 blob is full of digit runs), script/style, tags, and any podcast block."""
    t = source_text(path)
    sha = source_sha(t)
    t = re.sub(r"data:[\w/+.\-]+;base64,[A-Za-z0-9+/=]+", " ", t)
    if os.path.splitext(path)[1].lower() in (".html", ".htm"):
        t = re.sub(r"(?is)<(script|style)\b.*?</\1\s*>", " ", t)
        t = re.sub(r"(?s)<!--.*?-->", " ", t)
        t = re.sub(r"<[^>]+>", " ", t)
        t = html.unescape(t)
    return sha, t


# ----------------------------------------------------------------------------- check
# "survey" as an invitation to the Core's feedback survey: near a feedback word, and never
# "survey scan(s)" -- the DDA term (an MS1 survey scan) the teaching segment may well use
SURVEY = re.compile(r"\bsurveys?\b(?![\s-]+scans?\b)", re.I)
_FEEDBACK_WORDS = re.compile(r"how did we do|tell us|feedback|(?:5|five)[\s-]*minutes?\b|"
                             r"questionnaire", re.I)


def survey_invites(text):
    """The places `text` invites someone to a survey: "survey" (not "survey scan") within 120
    characters of a feedback word. -> [snippet]."""
    out = []
    for m in SURVEY.finditer(text):
        a, b = max(0, m.start() - 120), m.end() + 120
        if _FEEDBACK_WORDS.search(text[a:b]):
            out.append(" ".join(text[max(0, m.start() - 40):m.end() + 40].split()))
    return out


# ---------------------------------------------------------------- meaning (beyond tokens)
# check matches tokens. A staff episode (2026-10-16) passed it with a true number attached to the
# wrong comparison, a hedge dropped ("not a main driver" spoken as "aren't driving"), proteins put
# right after "It survives everything" that the session's robustness tables call not robust, and
# a tool presented as one the client uses. Two things follow: what can be caught mechanically is
# (a robust claim next to a protein whose table says it is not: FAIL; a superlative the sources do
# not make: WARN), and a separate reviewer reads every sentence against the report before render
# (meaning_check.md, below): render, link, share and job refuse without its PASS for this script.
ROBUST_CLAIM = re.compile(r"\b(?:robust(?:ly)?|surviv(?:e|es|ed|ing)|holds?\s+up|held\s+up|"
                          r"stands?\s+up|stood\s+up)\b", re.I)
# A claim covers its own clause: a sentence, cut at ; and at the words that turn a contrast
# ("Jph3 survives every check; Ryr2 doesn't", "Jph3 is robust, unlike Ryr2"). A colon does not
# cut: "It survives everything: Ryr2 and Jph3" names what the claim is about.
_CLAUSE_CUT = re.compile(r"[.;!?](?=\s|$)|\b(?:but|while|whereas|although|though|unlike|except|"
                         r"however|yet)\b", re.I)
# A negation before the claim in its clause turns it around ("Neither Ryr2 nor Actb is robust",
# "Not one of them ... is robust"); so does one beside a name ("Ryr2 and Actb don't").
_CLAUSE_NEG = re.compile(r"\b(?:not|no|never|neither|nor|none|without|n[’']t|\w+n[’']t|fails?\s+to|"
                         r"failed\s+to)\b", re.I)
_NAME_NEG_AFTER = re.compile(r"^[^.;!?]{0,40}?\b(?:\w+n[’']t|do(?:es)?\s+not|did\s+not|is\s+not|"
                             r"are\s+not|was\s+not|were\s+not|drops?\s+out|dropped\s+out|"
                             r"falls?\s+away|fell\s+away|fails?|failed)\b", re.I)
_NAME_NEG_BEFORE = re.compile(r"\b(?:not|neither|nor|unlike|except|but\s+not)\W+(?:\w+\W+){0,2}$",
                              re.I)
# What a `robust` cell says. R's write.csv writes TRUE/FALSE, pandas True/False; a cell that
# says neither is not understood (a WARN, never a FAIL).
ROBUST_YES = {"yes", "y", "true", "t", "1", "robust"}
ROBUST_NO = {"no", "n", "false", "f", "0", "not robust"}
_NAME_TOK = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.\-]*[A-Za-z0-9]|[A-Za-z0-9]")


def robust_value(v):
    """'yes', 'no' or None (not understood) for one `robust` cell."""
    v = " ".join(str(v or "").split()).lower()
    if v in ROBUST_YES:
        return "yes"
    if v in ROBUST_NO or v.startswith("not "):
        return "no"
    return None
SUPERLATIVE = re.compile(r"\b(?:the\s+only|only\s+one|the\s+(?:strongest|biggest|largest|"
                         r"highest|lowest|best|most|top)|every\s+single|none\s+of)\b", re.I)
MEANING_FILE = "meaning_check.md"
# What the reviewer checks in every sentence (the staff report's (a)-(e), and (f) from the 2.13.2
# review: draft history no phrase list catches). Said in the brief,
# in meaning_check.md and in references/podcast.md.
MEANING_CHECKS = (
    ("a", "the right comparison, group and direction for every number and named protein: the "
          "report's sentence or table it comes from says so"),
    ("b", "every hedge the report gives is kept: \"not a main driver\" never becomes \"isn't "
          "driving\", \"may\" never becomes \"does\""),
    ("c", "the robust status of every named protein: one spoken next to \"robust\", \"survives\" "
          "or \"holds up\" is robust in that comparison's robustness table"),
    ("d", "no superlative or totality the report does not make: \"the only\", \"all\", "
          "\"every\", \"the strongest\""),
    ("e", "no claim outside the report and the Claims list: a tool or database presented as one "
          "the client uses, a host answering their own question as if it were a finding"),
    ("f", "no draft history: nothing an earlier draft or version said, got wrong or had changed "
          "(the listener never saw one); say what the report finds and what was checked"),
)


def robustness_files(out, depth=3):
    """Every CSV under <out>/tables (to `depth` folders) whose header has a `robust` column --
    the session's own robustness tables, whatever made them. Sorted paths; [] when none."""
    root = os.path.join(out or "", "tables")
    if not out or not os.path.isdir(root):
        return []
    found, base = [], root.rstrip(os.sep).count(os.sep)
    for d, dirs, files in os.walk(root):
        if d.count(os.sep) - base >= depth:
            dirs[:] = []
        for f in sorted(files):
            path = os.path.join(d, f)
            if not f.lower().endswith(".csv"):
                continue
            try:
                with open(path, newline="", encoding="utf-8-sig", errors="replace") as fh:
                    head = next(csv.reader(fh), [])
            except (OSError, csv.Error):
                continue
            if "robust" in {h.strip().lower() for h in head}:
                found.append(path)
    return sorted(found)


def robustness_digest(out):
    """One sha256 of the robustness tables (their paths under <out> and contents), or "none":
    check.txt and meaning_check.md record it, so new or changed tables make both stale."""
    files = robustness_files(out)
    if not files:
        return "none"
    h = hashlib.sha256()
    for path in files:
        h.update(os.path.relpath(path, out).replace(os.sep, "/").encode("utf-8") + b"\0")
        h.update(sha256_file(path).encode("ascii") + b"\n")
    return h.hexdigest()


def robustness_tables(out, depth=3):
    """{protein name as written: [(table, robust cell, name)]} from robustness_files: names from
    Genes, Protein.Names and Protein.Group (";"-separated). Read as utf-8-sig (Excel and readr
    write a byte-order mark). {} when there are none."""
    found = {}
    for path in robustness_files(out, depth):
        try:
            with open(path, newline="", encoding="utf-8-sig", errors="replace") as fh:
                rows = csv.DictReader(fh)
                col = {h.strip().lower(): h for h in (rows.fieldnames or [])}
                for row in rows:
                    v = (row.get(col["robust"]) or "").strip()
                    for k in ("genes", "protein.names", "protein.group"):
                        for n in (row.get(col[k]) or "").split(";") if k in col else ():
                            n = n.strip()
                            if n:
                                found.setdefault(n, []).append(
                                    (os.path.relpath(path, out).replace(os.sep, "/"), v, n))
        except (OSError, csv.Error):
            continue
    return found


# A comma does not end a clause, except before a second, independent one (", and Ryr2 goes up",
# ", while ...", ", which ..."): never ", and so does Ryr2", which carries the claim on, and never
# inside a list a colon opened ("It survives everything: Ryr2, Jph3, and Actb").
_COMMA_CUT = re.compile(r",\s+(?:and|so|while|which|whereas)\b(?!\s+so\s+(?:does|did|do|is|are|was|"
                        r"were)\b)", re.I)
# English words that are also (mouse) gene symbols, read as English at the start of a sentence:
# "Was that robust?", "Set against that, ..." (review of 2.13.2, round 2).
SENTENCE_START_WORDS = {"was", "set", "max", "son", "rest", "impact", "cat", "arc", "star", "can",
                        "far", "ran", "gas", "spin", "pole", "tan", "lot", "not", "per", "top",
                        "act", "aim", "bad", "cap", "end", "fat", "fit", "ill", "key", "kit",
                        "mad", "met", "net", "pan", "ret", "tor", "via", "wars", "wise", "yes"}


def _robust_rows(robust, tok, sentence_start=False):
    """The table rows for a token: as written (a symbol that is also an English word -- SET, WAS,
    MAX -- is never matched by "set" or "was"); a token with a digit, or in capitals (ACTB for a
    mouse table's Actb), in either case. The first word of a sentence that is a common English
    word ("Was", "Set"; not in capitals: SET is the gene) is read as English."""
    if sentence_start and tok.lower() in SENTENCE_START_WORDS and not tok.isupper():
        return None
    rows = robust.get(tok)
    if rows or not (re.search(r"\d", tok) or (tok.isupper() and len(tok) >= 3)):
        return rows
    low = tok.lower()
    return [r for n, rs in robust.items() if n.lower() == low for r in rs] or None


def robust_claims(text, robust):
    """-> [(level, protein, claim, detail)]: a protein named in the same clause as a robust claim
    ("robust", "survives", "holds up") that no negation turns around, whose robustness tables
    say not robust in every row (FAIL), robust in some only (WARN), or something not understood
    (WARN). The meaning check reads the rest."""
    out, seen = [], set()
    if not robust:
        return out
    starts = {0} | {m.end() for m in re.finditer(r"[.!?]\s+", text)}     # sentence starts
    cuts = [0] + [m.end() for m in _CLAUSE_CUT.finditer(text)] + [len(text)]
    for a, b in zip(cuts, cuts[1:]):
        clause = text[a:b]
        for m in ROBUST_CLAIM.finditer(clause):
            if _CLAUSE_NEG.search(clause[:m.start()]):
                continue
            # the claim's reach: its clause, less a second one a comma opens -- after the claim,
            # unless a colon opened a list; before it, unless what precedes is a list of names
            # ("Ryr2, Jph3, and Actb all survive")
            begin, end = 0, len(clause)
            c = _COMMA_CUT.search(clause, m.end())
            if c and ":" not in clause[m.end():c.start()]:
                end = c.start()
            for c in _COMMA_CUT.finditer(clause, 0, m.start()):
                words = [w for w in _NAME_TOK.findall(clause[begin:c.start()])
                         if w.lower() not in ("and", "or")]
                if not words or not all(_robust_rows(robust, w) for w in words):
                    begin = c.end()
            for t in _NAME_TOK.finditer(clause, begin, end):
                rows = _robust_rows(robust, t.group(0), sentence_start=(a + t.start()) in starts)
                if not rows or rows[0][2] in seen:
                    continue
                # a negation beside the name -- after it, but never past the claim word ("Ryr2
                # holds up, which doesn't happen often" is a claim)
                after = clause[t.end():m.start()] if t.start() < m.start() else clause[t.end():end]
                if _NAME_NEG_AFTER.search(after) or _NAME_NEG_BEFORE.search(clause[:t.start()]):
                    continue
                seen.add(rows[0][2])
                said = [(tbl, v, robust_value(v)) for tbl, v, _ in rows]
                yes = [x for x in said if x[2] == "yes"]
                no = [x for x in said if x[2] == "no"]
                odd = [x for x in said if x[2] is None]
                where = "; ".join(f"{tbl}: robust = {v or 'blank'}" for tbl, v, _ in (no + odd)[:3])
                if odd and not yes:
                    out.append(("WARN", rows[0][2], m.group(0),
                                f"robust value not understood -- {where}; the meaning check "
                                "decides"))
                elif not yes:
                    out.append(("FAIL", rows[0][2], m.group(0), where))
                elif no or odd:
                    out.append(("WARN", rows[0][2], m.group(0),
                                f"robust in {', '.join(t for t, _, _ in yes[:3])}; not in {where}"))
    return out


def meaning_path(script_path):
    return os.path.join(os.path.dirname(os.path.abspath(script_path)), MEANING_FILE)


def read_meaning(s):
    """-> (status, problem) from the meaning_check.md beside the script: ("PASS", None) only for
    a PASS recorded for THIS script's sha256. Else the status ("missing", "FAIL", "STALE", ...)
    and why render, link, share and job refuse."""
    path = meaning_path(s.path)
    try:
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
    except OSError:
        return "missing", (f"no {MEANING_FILE} beside {os.path.basename(s.path)}: the meaning "
                           "check has not been recorded")
    verdict = re.search(r"(?m)^verdict: (\w+)", text)
    sha = re.search(r"(?m)^script_sha256: (\w+)", text)
    if not verdict or verdict.group(1) != "PASS":
        return (verdict.group(1) if verdict else "unreadable"), f"{MEANING_FILE} does not say PASS"
    if not sha or sha.group(1) != s.sha256:
        return "STALE", f"the script changed after {MEANING_FILE} was recorded"
    # ...and what it was read against: each delivered source and the robustness tables (a
    # meaning check from before these lines were written has none to compare)
    for m in re.finditer(r"(?m)^source: (.+) sha256=(\w+)$", text):
        try:
            now = source_sha(source_text(os.path.join(os.path.dirname(path), m.group(1))))
        except OSError:
            return "STALE", f"{os.path.basename(m.group(1))} is missing"
        if now != m.group(2):
            return "STALE", (f"{os.path.basename(m.group(1))} changed after {MEANING_FILE} was "
                             "recorded")
    tables = re.search(r"(?m)^robust_tables_sha256: (\w+)$", text)
    out = output_dir_of(s.path)
    if tables and out and robustness_digest(out) != tables.group(1):
        return "STALE", f"the robustness tables changed after {MEANING_FILE} was recorded"
    return "PASS", None


def meaning_brief(s, robust, sources):
    """The reviewer's brief: what to check in every sentence, the sources, the robustness tables,
    and the transcript, sentence by sentence with its line."""
    L = [f"# Meaning check brief: {os.path.basename(s.path)}", f"script_sha256: {s.sha256}", "",
         "You did NOT write this script. Read every sentence against the report and the files "
         "below, and check each for:", ""]
    L += [f"- ({k}) {t}" for k, t in MEANING_CHECKS]
    L += ["", "Sources: " + (", ".join(sources) if sources else "the files check.txt lists"),
          "Robustness tables (a `robust` column): " + (", ".join(sorted({r[0] for v in robust.values()
                                                                        for r in v}))
                                                      if robust else "none in output/tables"),
          "", "For each problem write one line: `line N: <what is wrong> -> <the fix>` (mark it "
          "`resolved` once the script is fixed). Then record the verdict, the script's sha256 "
          "is taken by the tool:", "",
          f"    make_podcast.py meaning {os.path.basename(s.path)} --verdict PASS|FAIL "
          "--reviewer \"<who, not the writer>\" --findings <your findings file>", "",
          "## The transcript, sentence by sentence", ""]
    for t in s.turns():
        for sent in re.split(r"(?<=[.!?])\s+", t.text.strip()):
            if sent:
                L.append(f"- line {t.line} {t.speaker}: {sent}")
    return "\n".join(L) + "\n"


def cmd_meaning(a):
    """--brief: print (and write beside the script) what the reviewer checks. --verdict: record
    the reviewer's verdict for THIS script (its sha256 is computed here, never typed)."""
    s = parse_script(a.script)
    out = output_dir_of(s.path, getattr(a, "output_dir", None))
    robust = robustness_tables(out) if out else {}
    if a.brief:
        text = meaning_brief(s, robust, [rel(x, os.path.dirname(s.path)) for x in a.source or []])
        with open(os.path.join(os.path.dirname(s.path), "meaning_brief.md"), "w",
                  encoding="utf-8") as fh:
            fh.write(text)
        if hasattr(sys.stdout, "reconfigure"):           # a cp1252 console: "α", "→" print as ?
            try:
                sys.stdout.reconfigure(errors="replace")
            except (ValueError, OSError):
                pass
        print(text)
        return 0
    if not a.verdict:
        log("[meaning] give --brief (what to check) or --verdict PASS|FAIL with --reviewer and "
            "--findings")
        return 2
    who = " ".join((a.reviewer or "").split())
    if not who:
        log("[meaning] --reviewer is required: who read the script (an agent or a person that "
            "did not write it)")
        return 2
    if s.author and who.lower() == " ".join(str(s.author).split()).lower():
        log(f"[meaning] refusing: the reviewer ({who}) is the script's author: the meaning check "
            "is a separate reader's")
        return 2
    try:
        with open(a.findings, encoding="utf-8") as fh:
            findings = fh.read().strip()
    except (OSError, TypeError) as e:
        log(f"[meaning] --findings: {e}")
        return 2
    if not findings:
        log("[meaning] --findings is empty: list each problem found (or write `None`)")
        return 2
    # read against what passed the check: its sources (and the tables) are recorded below, so a
    # report rewritten after this verdict makes it stale
    status, srcs, why = read_check(s)
    if why:
        log(f"[meaning] refusing: {why}. check.txt must PASS for this script and its sources "
            "first: the meaning check reads the same report")
        return 2
    brief = os.path.join(os.path.dirname(s.path), "meaning_brief.md")
    try:
        with open(brief, encoding="utf-8") as fh:
            brief_ok = f"script_sha256: {s.sha256}" in fh.read()
    except OSError:
        brief_ok = False
    if not brief_ok:
        log(f"[meaning] refusing: no meaning_brief.md for this script: run `make_podcast.py "
            f"meaning {a.script} --brief` first -- the reviewer works from it")
        return 2
    mech = []
    for t in s.turns():
        for level, name, claim, detail in robust_claims(t.text, robust):
            mech.append(f"- {level}: line {t.line}: {name} next to '{claim}' ({detail})")
    if a.verdict == "PASS" and any(x.startswith("- FAIL") for x in mech):
        log("[meaning] refusing PASS: " + "; ".join(x[2:] for x in mech if x.startswith("- FAIL")))
        return 2
    L = [f"# Meaning check: {os.path.basename(s.path)}",
         f"script_sha256: {s.sha256}", f"verdict: {a.verdict}", f"reviewer: {who}",
         f"checked: {now_iso()}", f"brief_sha256: {sha256_file(brief)}"]
    L += [f"source: {x['file']} sha256={x['sha256']}" for x in srcs if x.get("delivered")]
    L += [f"robust_tables_sha256: {robustness_digest(out) if out else 'none'}", "",
          "Every sentence was read against the report for:"]
    L += [f"- ({k}) {t}" for k, t in MEANING_CHECKS]
    L += ["", "## Mechanical (robustness tables)", ""] + (mech or ["- none"])
    L += ["", "## Findings", "", findings, ""]
    path = meaning_path(s.path)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    print(json.dumps({"meaning_check": path, "verdict": a.verdict, "script_sha256": s.sha256},
                     indent=2))
    return 0 if a.verdict == "PASS" else 1


def _core_run(session):
    """submission_report.core_run: the one test of a UC Davis Core run. False when it cannot be
    asked -- never invite a user outside the Core to the Core's survey."""
    try:
        import submission_report
    except ImportError:
        return False
    return bool(submission_report.core_run(session))


def _as_source(x):
    """A source as check uses it: a dict with path, sha, text, kind ('delivered' | 'extra'),
    label and ref (what check.txt records). A (path, sha, text) tuple is a delivered source."""
    if isinstance(x, dict):
        return x
    path, sha, text = x
    return {"path": path, "sha": sha, "text": text, "kind": "delivered", "label": None,
            "ref": os.path.basename(path)}


def check(s, sources, forbid=(), core=False, robust=None):
    """-> dict with fails / warns / infos (lists of str) and stats. `sources`: dicts from
    cmd_check (or (path, sha, text) tuples, taken as delivered). A token found only in an
    extra source (one the lab was not given) passes only when a Claims bullet names it and that
    source. `core`: a UC Davis Core run (submission_report.core_run of the session, as link and
    share ask it), whose report ends with the Core's survey -- the only kind of episode that may
    invite listeners to it."""
    fails, warns, infos = [], [], []
    for line, msg in s.problems:
        fails.append((f"line {line}: " if line else "") + msg)
    if not s.title:
        warns.append("no 'Title:' in the header (and no '# ' heading): the card will say "
                     "'Untitled'")
    if not s.has_pronunciation:
        fails.append("missing section '## Pronunciation' (a table 'Written | Spoken'; it may "
                     "have no rows)")
    if not s.has_claims:
        fails.append("missing section '## Claims beyond the report' (bullets, or 'None')")
    srcs = [_as_source(x) for x in sources]
    if not [x for x in srcs if x["kind"] == "delivered"]:
        fails.append("no delivered --source given: nothing to check the numbers and symbols "
                     "against")

    def pool(kind):
        book, hay = NumberBook(), []
        for x in srcs:
            if x["kind"] == kind:
                book.add_text(x["text"])
                hay.append(normalize_numbers(x["text"]).lower())
        hay = "\n".join(hay)
        return book, hay, _words_norm(hay)
    book, hay, phrases = pool("delivered")
    xs = []                                  # each extra source on its own, to name the right one
    for x in srcs:
        if x["kind"] == "extra":
            b = NumberBook()
            b.add_text(x["text"])
            h = normalize_numbers(x["text"]).lower()
            xs.append({"label": x.get("label") or os.path.basename(x["path"]), "book": b, "hay": h,
                       "phrases": _words_norm(h), "names": [n.lower() for n in (
                           x.get("label"), os.path.basename(x["path"])) if n]})

    def extra_hits(test):
        return [x for x in xs if test(x)]

    def extra_only(where, what, pred, hits):
        """`what` is only in the extra sources `hits`: pass when a Claims bullet says it (pred)
        and names one of them."""
        labels = ", ".join(repr(x["label"]) for x in hits)
        if any(pred(c) and any(n in c.lower() for x in hits for n in x["names"]) for c in s.claims):
            infos.append(f"{where}: {what} is only in {labels} (not delivered); it is listed under "
                         "Claims beyond the report with its source")
        else:
            fails.append(f"{where}: {what} is only in {labels}, not in the files the lab has -- "
                         "list it under Claims beyond the report naming that source, or drop it")
    claims_hay = normalize_numbers(s.claims_text).lower()
    disclosed = NumberBook()
    disclosed.add_text(s.claims_text)
    if len(book.abs) > 20000:
        warns.append(f"the sources hold {len(book.abs):,} distinct numbers: a large numeric "
                     "table makes the number check weak -- pass the report's text (e.g. "
                     "AI_Analysis_Report.md, AUDIT.md), not the DE tables")

    turns = s.turns()
    for i, t in enumerate(turns, 1):
        where = f"turn {i} (line {t.line}, {t.speaker.upper()})"
        for n in numbers_in(t.text):
            v = book.verdict(n)
            if v in ("exact", "trivial"):
                continue
            if v == "rounded":
                infos.append(f"{where}: {n.text} matches a source value rounded to its precision")
            elif v == "magnitude":
                infos.append(f"{where}: {n.text} is spoken as an order of magnitude of a source "
                             "value (same power of ten)")
            elif xs and extra_hits(lambda x, n=n: x["book"].verdict(n) in (
                    "exact", "rounded", "magnitude")):
                extra_only(where, f"number {n.text}", lambda c, v=n.value: any(
                    b.value == v for b in numbers_in(c)),
                    extra_hits(lambda x, n=n: x["book"].verdict(n) in ("exact", "rounded", "magnitude")))
            elif disclosed.verdict(n) in ("exact", "trivial"):
                infos.append(f"{where}: {n.text} is not in the sources; it is listed under "
                             "Claims beyond the report")
            else:
                hint = book.mantissas.get(n.value) if n.kind == "plain" else None
                fails.append(f"{where}: number {n.text} is not in the sources"
                             + (f" (they have {hint}: say the power of ten too)" if hint else "")
                             + f" -- {_ctx(t.text, n.text)}")
        for m in SPELLED_NUMBER.finditer(t.text):
            # The report's own phrase is fine: "fewer than half of the runs" in the sources lets
            # "half of the runs" through, not "half the proteins were inferred".
            ph = quantity_phrase(t.text, m)
            if ph and f" {ph} " in phrases:
                continue
            hits = extra_hits(lambda x, ph=ph: f" {ph} " in x["phrases"]) if ph else []
            if hits:
                extra_only(where, f"'{ph}'", lambda c, ph=ph: f" {ph} " in _words_norm(c), hits)
                continue
            if re.search(r"\b" + re.escape(m.group(0).lower()) + r"\b", s.claims_text.lower()):
                infos.append(f"{where}: '{m.group(0)}' is listed under Claims beyond the report")
                continue
            fails.append(f"{where}: quantity in words '{ph or m.group(0)}' is not in the sources "
                         "-- use the report's own phrase or its digits (the pronunciation step "
                         "handles speech), or list it under Claims beyond the report if it is "
                         "your own gloss")
        for adv in propobs_advice(t.text):
            fails.append(f"{where}: '{adv}' -- do not tier or rank hits by PropObs: it is "
                         "the observed fraction over ALL runs, so a hit found only in its own "
                         "group scores low. Point to the per-group detection columns "
                         "(Detected_<group>, Evidence) or the tier file instead")
        for inv in ([] if core else survey_invites(t.text)):
            fails.append(f"{where}: '{inv}' -- invites listeners to a survey, but this is not a "
                         "UC Davis Core run: only a Core run's report ends with the Core's "
                         "feedback survey. Take it out")
        for m in SPECIALTY.finditer(t.text):
            fails.append(f"{where}: '{m.group(0)}' -- a host must not claim a real research "
                         "specialty (it lends a synthetic voice false authority); say \"I'm the "
                         "biologist of the pair\" or \"I'm the statistician\"")
        for level, name, claim, detail in robust_claims(t.text, robust):
            (fails if level == "FAIL" else warns).append(
                f"{where}: {name} is named next to '{claim}', but the robustness tables say "
                f"otherwise ({detail}) -- name only the robust proteins there, or say which "
                "comparison it holds in")
        for m in SUPERLATIVE.finditer(t.text):
            ph = quantity_phrase(t.text, m) or m.group(0).lower()
            if f" {ph} " not in phrases:
                warns.append(f"{where}: '{ph}' -- a superlative the sources do not say in these "
                             "words: the meaning check must confirm the report makes it")
        for snip in draft_history(t.text):
            fails.append(f"{where}: draft history '{snip}' -- the listener never saw a draft: "
                         "say what the report finds and what was checked, in the present tense, "
                         "never what an earlier draft said, got wrong or corrected "
                         "(references/podcast.md)")
        for snip in draft_history_loose(t.text):
            warns.append(f"{where}: possible draft history '{snip}' -- if it tells of an earlier "
                         "draft or a fix, say what the report finds instead; the meaning check "
                         "(f) decides")

    seen = {}
    for i, t in enumerate(turns, 1):
        for m in _SYM.finditer(t.text):
            tok = m.group(0)
            if is_symbol(tok) and tok not in seen:
                seen[tok] = (i, t)
        for m in _HYPHEN_SYM.finditer(t.text):             # IL-6, COVID-19, LRS-124
            if m.group(0) not in seen:
                seen[m.group(0)] = (i, t)
    for tok, (i, t) in seen.items():
        where = f"turn {i} (line {t.line}, {t.speaker.upper()})"
        if symbol_found(tok, hay):
            continue
        hits = extra_hits(lambda x, tok=tok: symbol_found(tok, x["hay"]))
        if hits:
            extra_only(where, tok, lambda c, tok=tok: symbol_found(tok, normalize_numbers(c).lower()),
                       hits)
            continue
        if symbol_found(tok, claims_hay):
            infos.append(f"{where}: {tok} is not in the sources; it is listed under Claims beyond "
                         "the report")
            continue
        fails.append(f"{where}: {tok} looks like a gene/protein symbol or acronym and is not in "
                     "the sources -- use the report's own spelling, list it under Claims beyond "
                     "the report, or (if it is emphasis) write it in lowercase or *italics*")

    # Capitalised words mid-sentence (Gapdh, Actb, Western, a person): symbol-like in speech but
    # not caught above. Allowed: in the sources, in the claims, a host's name, the show's name.
    allowed = {h.lower() for h, _ in s.hosts} | {w.lower() for w in re.findall(r"[A-Za-z]+",
                                                                                s.show or "")}
    caps = {}
    for i, t in enumerate(turns, 1):
        for tok in title_tokens(t.text):
            if tok.lower() not in allowed and tok.lower() not in GENERIC_TOKENS and tok not in caps:
                caps[tok] = (i, t)
    initial = {}
    for i, t in enumerate(turns, 1):
        for tok in title_tokens(t.text, initial=True):
            if (tok.lower() not in allowed and tok not in initial and tok not in caps
                    and not symbol_found(tok, hay) and not symbol_found(tok, claims_hay)):
                initial[tok] = i
    if initial:
        infos.append(f"capitalised at the start of a sentence, so not checked, and not in the "
                     f"sources ({len(initial)}): " + ", ".join(f"{k} (turn {v})" for k, v in
                                                              initial.items()))
    for tok, (i, t) in caps.items():
        where = f"turn {i} (line {t.line}, {t.speaker.upper()})"
        if symbol_found(tok, hay):
            continue
        hits = extra_hits(lambda x, tok=tok: symbol_found(tok, x["hay"]))
        if hits:
            extra_only(where, tok, lambda c, tok=tok: symbol_found(tok, normalize_numbers(c).lower()),
                       hits)
            continue
        if symbol_found(tok, claims_hay):
            infos.append(f"{where}: {tok} is not in the sources; it is listed under Claims beyond "
                         "the report")
            continue
        fails.append(f"{where}: {tok} is capitalised mid-sentence and is not in the sources, the "
                     "claims, a host's name or the show's name -- a gene, a person or a place the "
                     "report does not mention? Use the report's spelling, list it under Claims "
                     "beyond the report, or lowercase it if it is an ordinary word")

    # Pronunciation rows: the spoken form goes to the voices unchecked, so it may only re-spell.
    pron = []
    for w, sp in s.pronunciation:
        probs = pronunciation_problems(w, sp)
        pron.append((w, sp, probs))
        for pr in probs:
            fails.append(f"pronunciation {w!r} -> {sp!r}: {pr}")

    if turns and not any(DISCLOSURE_RE.search(t.text) for t in turns[:3]):
        fails.append("no AI disclosure in the first 3 turns: one host must say, early and plainly, "
                     "that this is an AI-generated discussion (e.g. 'AI-generated')")

    for name in forbid:
        name = name.strip()
        if not name:
            continue
        rx = re.compile(r"(?<![A-Za-z])" + re.escape(name) + r"(?![A-Za-z])", re.I)
        hits = [f"header '{k}'" for k, v in s.header.items() if rx.search(v)]
        if s.title and rx.search(s.title) and not any(h == "header 'title'" for h in hits):
            hits.insert(0, "the title")
        hits += [f"turn {i} (line {t.line})" for i, t in enumerate(turns, 1) if rx.search(t.text)]
        if rx.search(s.claims_text):
            hits.append("Claims beyond the report")
        hits += [f"pronunciation {w!r} -> {sp!r}" for w, sp in s.pronunciation
                 if rx.search(w) or rx.search(sp)]
        hits += [f"the style for {h}" for h, st in s.styles.items() if rx.search(st)]
        if hits:
            fails.append(f"forbidden name {name!r} appears in: {', '.join(hits)}")

    said = "\n".join(t.text for t in turns)
    missing = [name for name, rx in TEACHING if not re.search(rx, said, re.I)]
    if turns and missing:
        warns.append("the collaborator is never taught: " + "; ".join(missing) + " (the 'How "
                     "proteomics works, with your data' segment; skip what does not apply here)")
    tail = "\n".join(t.text for seg in s.segments[-2:] for t in seg)
    if turns and not re.search(NEXT_STEPS, tail, re.I):
        warns.append("the last two segments never say what to do with this: README.html, the "
                     "report, the per-group detection columns (Detected_<group>, Evidence) or the "
                     "tier file, and what to validate first")

    words = sum(t.words for t in turns)
    by = {}
    for t in turns:
        by[t.speaker] = by.get(t.speaker, 0) + t.words
    stats = {"words": words, "turns": len(turns), "segments": len(s.segments),
             "minutes": round(words / float(WPM), 1),
             "share": {h: round(100.0 * by.get(h, 0) / words) if words else 0 for h, _ in s.hosts}}
    if turns and not WORDS_LO <= words <= WORDS_HI:
        warns.append(f"{words:,} words (about {words / float(WPM):.0f} min): the brief asks for "
                     f"{WORDS_LO:,}-{WORDS_HI:,} (target 2,500-3,200)")
    if s.segments and not SEGMENTS_LO <= len(s.segments) <= SEGMENTS_HI:
        warns.append(f"{len(s.segments)} segments: the brief asks for {SEGMENTS_LO}-{SEGMENTS_HI} "
                     "('---' lines)")
    for h, _ in s.hosts:
        if turns and not by.get(h):
            fails.append(f"host {h} never speaks")
        elif words and by.get(h, 0) < 0.35 * words:
            warns.append(f"{h} speaks only {100.0 * by[h] / words:.0f}% of the words")
    return {"status": "FAIL" if fails else "PASS", "fails": fails, "warns": warns,
            "infos": infos, "stats": stats, "pron": pron,
            "robust_tables": sorted({r[0] for v in (robust or {}).values() for r in v})}


def _ctx(text, token):
    t = normalize_numbers(text)
    i = t.find(token)
    if i < 0:
        return repr(text[:80])
    a, b = max(0, i - 45), min(len(t), i + len(token) + 45)
    return repr(("…" if a else "") + t[a:b] + ("…" if b < len(t) else ""))


def check_report(s, res, sources, forbid):
    st = res["stats"]
    L = [f"Podcast script check: {res['status']} ({len(res['fails'])} problem(s), "
         f"{len(res['warns'])} warning(s))",
         f"status: {res['status']}",
         f"script: {os.path.basename(s.path)}",
         f"script_sha256: {s.sha256}",
         f"checked: {now_iso()}"]
    for x in (_as_source(y) for y in sources):          # relative paths / basenames only
        if x["kind"] == "delivered":
            L.append(f"source: {x['ref']} sha256={x['sha']}")
        else:
            L.append(f"extra-source: {x['ref']} sha256={x['sha']} label=\"{x.get('label') or ''}\"")
    names = [n for n in forbid if n.strip()]
    L.append(f"forbidden names: {len(names)} given (not listed here)" if names else
             "forbidden names: none given")
    rt = res.get("robust_tables") or []
    L.append("robustness tables: " + (", ".join(rt) if rt else "none (no CSV under output/tables "
                                                               "has a `robust` column)"))
    if res.get("robust_digest"):
        L.append(f"robust_tables_sha256: {res['robust_digest']}")
    L += ["", f"words: {st['words']:,} · turns: {st['turns']} · segments: {st['segments']} · "
              f"about {st['minutes']} min at {WPM} wpm",
          "speaking share: " + " · ".join(f"{h} {v}%" for h, v in st["share"].items())]
    for title, items in (("FAIL", res["fails"]), ("WARN", res["warns"]), ("INFO", res["infos"])):
        if items:
            L += ["", title] + [f"- {x}" for x in items]
    rows = res.get("pron") or []
    L += ["", f"PRONUNCIATION ({len(rows)} row(s) from the script, applied to speech only; the "
              "built-in defaults are not listed)"]
    L += [f"- {w} -> {sp}" + (f"   [FAIL: {'; '.join(pr)}]" if pr else "") for w, sp, pr in rows]
    if res["status"] == "PASS":
        L += ["", "PASS means: every number (except unsigned counts of 10 or less), every "
                  "symbol-like token and every capitalised word mid-sentence in the transcript was "
                  "found in the sources or is listed under Claims beyond the report; the AI "
                  "disclosure is spoken early; no forbidden name appears; every pronunciation row "
                  "only re-spells its written form; no turn narrates an earlier draft "
                  "(client_text.draft_history); no protein the robustness tables call not "
                  "robust is named next to a robust claim. It does NOT mean the script is right. "
                  f"The meaning check ({MEANING_FILE}) reads it for meaning before render."]
    L += ["", "What check cannot catch -- read the script against the report for these:"]
    L += [f"- {x}" for x in CANNOT_CATCH]
    return "\n".join(L) + "\n"


def output_dir_of(script_path, given=None):
    """The session's output folder: --output-dir, else the folder above a podcast/ folder."""
    if given:
        return os.path.realpath(given)
    d = os.path.dirname(os.path.abspath(script_path))
    return os.path.realpath(os.path.dirname(d)) if os.path.basename(d) == "podcast" else None


def output_dir_problem(d):
    """None when `d` is a session's output folder -- it holds the report
    (Analysis_Report.html) or the DE record (tables/de_provenance.json) -- else why not. A
    folder above it (the session, its parent) would let a draft count as delivered."""
    if not os.path.isdir(d):
        return f"{d} is not a folder"
    if (os.path.isfile(os.path.join(d, "Analysis_Report.html"))
            or os.path.isfile(os.path.join(d, "tables", "de_provenance.json"))):
        return None
    return (f"{os.path.basename(d) or d} is not a session's output folder: it holds neither "
            "Analysis_Report.html nor tables/de_provenance.json")


def cmd_check(a):
    s = parse_script(a.script)
    sdir = os.path.dirname(s.path)
    out = output_dir_of(s.path, a.output_dir)
    pod = os.path.realpath(sdir)
    sources, problems = [], []
    if not out:
        problems.append("cannot tell the session's output folder: keep the script in "
                        "<session>/output/podcast/ or pass --output-dir")
    elif output_dir_problem(out):
        problems.append(("--output-dir" if a.output_dir else "the folder above podcast/") + ": "
                        + output_dir_problem(out) + ". Point it at <session>/output, the folder the "
                        "lab receives; nothing is counted as delivered until then")
        out = None
    for p in a.source:
        if not out:
            continue                              # the one problem above says why
        if not os.path.isfile(p):
            problems.append(f"source not found: {p}")
            continue
        rp = os.path.realpath(p)
        if not out or not rp.startswith(out + os.sep) or rp.startswith(pod + os.sep):
            problems.append(f"source {os.path.basename(p)} is not a delivered file: a --source must "
                            f"be in the session's output folder the lab receives (the report, "
                            f"methods.md, SAMPLE_QUALITY.md, AUDIT.md), not the podcast folder. "
                            f"Give anything else as --extra-source {os.path.basename(p)} "
                            f"--label \"what it is\"")
            continue
        sha, text = load_source(p)
        sources.append({"path": rp, "sha": sha, "text": text, "kind": "delivered", "label": None,
                        "ref": os.path.relpath(rp, pod).replace(os.sep, "/")})
    extras, labels = list(a.extra_source or []), list(a.label or [])
    if len(extras) != len(labels):
        problems.append(f"{len(extras)} --extra-source but {len(labels)} --label: give each extra "
                        "source its own --label \"what it is\"")
    for p, lab in zip(extras, labels):
        if not os.path.isfile(p):
            problems.append(f"extra source not found: {p}")
            continue
        sha, text = load_source(p)
        sources.append({"path": os.path.realpath(p), "sha": sha, "text": text, "kind": "extra",
                        "label": lab.strip() or os.path.basename(p), "ref": os.path.basename(p)})
    res = check(s, sources, a.forbid_name or [], core=bool(out and _core_run(_session_of(out))),
                robust=robustness_tables(out) if out else None)
    res["robust_digest"] = robustness_digest(out) if out else None
    for pr in reversed(problems):
        res["fails"].insert(0, pr)
    if problems:
        res["status"] = "FAIL"
    path = os.path.join(sdir, "check.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(check_report(s, res, sources, a.forbid_name or []))
    print(json.dumps({"status": res["status"], "fails": len(res["fails"]),
                      "warnings": len(res["warns"]), "check": path, **res["stats"]}, indent=2))
    for f in res["fails"][:40]:
        log(f"[FAIL] {f}")
    for w in res["warns"]:
        log(f"[WARN] {w}")
    return 0 if res["status"] == "PASS" else 1


def read_check(s):
    """-> (status, sources, problem) from the check.txt beside the script. The check holds only
    while the script AND every source are unchanged: a source is re-hashed as check hashed it
    (source_text: its text without a podcast block), so link's Listen line does not count as a
    change and a regenerated or edited report does."""
    path = os.path.join(os.path.dirname(s.path), "check.txt")
    try:
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
    except OSError:
        return None, [], f"no check.txt beside {os.path.basename(s.path)}"
    status = re.search(r"(?m)^status: (\w+)", text)
    sha = re.search(r"(?m)^script_sha256: (\w+)", text)
    srcs = [{"file": m.group(1), "sha256": m.group(2), "delivered": True}
            for m in re.finditer(r"(?m)^source: (.+) sha256=(\w+)$", text)]
    extra = [{"file": m.group(1), "sha256": m.group(2), "delivered": False, "label": m.group(3)}
             for m in re.finditer(r'(?m)^extra-source: (.+) sha256=(\w+) label="(.*)"$', text)]
    if not status or status.group(1) != "PASS":
        return (status.group(1) if status else None), srcs, "check.txt does not say PASS"
    if not sha or sha.group(1) != s.sha256:
        return "STALE", srcs, "the script changed after check.txt was written"
    if not srcs:
        return "STALE", srcs + extra, "check.txt names no delivered source"
    tables = re.search(r"(?m)^robust_tables_sha256: (\w+)$", text)
    out = output_dir_of(s.path)
    if tables and out and robustness_digest(out) != tables.group(1):
        return "STALE", srcs + extra, "the robustness tables changed after check.txt was written"
    for src in srcs:                     # relative to check.txt; extra sources are not re-read
        try:
            now = source_sha(source_text(os.path.join(os.path.dirname(path), src["file"])))
        except OSError:
            return "STALE", srcs + extra, f"source {src['file']} is missing"
        if now != src["sha256"]:
            return "STALE", srcs + extra, (f"source {os.path.basename(src['file'])} changed after "
                                           "check.txt was written")
    return "PASS", srcs + extra, None


# ----------------------------------------------------------------------------- audio
def _arr(data):
    """Native-order 16-bit samples: what the wave module reads and writes (it converts to and
    from a WAV file's little-endian itself, on a big-endian host too)."""
    a = array.array("h")
    a.frombytes(data[: len(data) // 2 * 2])
    return a


def _arr_le(data):
    """Raw little-endian 16-bit PCM (Gemini's audio/L16) as native samples: the one place that
    swaps, and only on a big-endian host."""
    a = _arr(data)
    if sys.byteorder == "big":
        a.byteswap()
    return a


def silence(sec):
    return array.array("h", [0]) * int(RATE * sec)


def resample(a, src, dst=RATE):
    if src == dst or not a:
        return a
    n = int(len(a) * dst / float(src))
    step = src / float(dst)
    out = array.array("h", [0]) * n
    last = len(a) - 1
    for i in range(n):
        x = i * step
        j = int(x)
        f = x - j
        out[i] = int(a[j] * (1 - f) + a[min(j + 1, last)] * f)
    return out


def parse_wav_bytes(data):
    """A RIFF/WAVE byte string -> array('h') at RATE, parsed by hand rather than with the wave
    module. It takes PCM and WAVE_FORMAT_EXTENSIBLE with a PCM sub-format, which Python 3.9's
    wave module refuses. A data chunk whose size is 0 or 0xFFFFFFFF (a streaming header), or
    runs past the end, is read to the end. 16-bit only; the first channel of several. Raises
    ValueError for anything else."""
    if len(data) < 12 or data[:4] != b"RIFF" or data[8:12] != b"WAVE":
        raise ValueError("not a RIFF/WAVE file")
    pos, fmt = 12, None
    while pos + 8 <= len(data):
        cid, size = data[pos:pos + 4], struct.unpack("<I", data[pos + 4:pos + 8])[0]
        body = pos + 8
        if cid == b"fmt ":
            if body + 16 > len(data):
                raise ValueError("truncated fmt chunk")
            tag, ch, rate, _, _, bits = struct.unpack("<HHIIHH", data[body:body + 16])
            if tag == 0xFFFE and size >= 40 and body + 26 <= len(data):
                tag = struct.unpack("<H", data[body + 24:body + 26])[0]   # the sub-format
            fmt = (tag, ch, rate, bits)
        elif cid == b"data":
            if fmt is None:
                raise ValueError("data chunk before fmt")
            tag, ch, rate, bits = fmt
            if tag != 1 or bits != 16 or not ch or not rate:
                raise ValueError(f"not 16-bit PCM (format {tag:#x}, {bits}-bit, {ch} channel(s))")
            end = len(data) if size in (0, 0xFFFFFFFF) or body + size > len(data) else body + size
            pcm = data[body:end]
            a = _arr_le(pcm[: len(pcm) // (2 * ch) * (2 * ch)])
            return resample(a[::ch] if ch > 1 else a, rate)
        if size in (0xFFFFFFFF,):
            break
        pos = body + size + (size & 1)                       # chunks are word-aligned
    raise ValueError("no data chunk")


def decode_audio(data, mime=""):
    """Audio bytes -> array('h') at RATE: a RIFF/WAVE file (parse_wav_bytes), or raw 16-bit
    little-endian PCM when the mime type says so (Gemini 2.x: "audio/L16;codec=pcm;rate=24000").
    Anything else -- audio/wav without a RIFF header included -- is a ValueError, never read as
    PCM."""
    if data[:4] == b"RIFF":
        return parse_wav_bytes(data)
    mime = (mime or "").lower()
    if "wav" in mime or not re.search(r"l16|pcm", mime):
        raise ValueError(f"not a RIFF file and not raw PCM (mime type {mime or 'not given'})")
    m = re.search(r"rate=(\d+)", mime)
    return resample(_arr_le(data), int(m.group(1)) if m else RATE)


def _from_wave(w):
    ch, sw, rate, n = w.getnchannels(), w.getsampwidth(), w.getframerate(), w.getnframes()
    if sw != 2:
        raise ValueError(f"expected 16-bit audio, got {8 * sw}-bit")
    a = _arr(w.readframes(n))
    if ch > 1:
        a = a[::ch]
    return resample(a, rate)


def read_wav(path):
    with wave.open(path, "rb") as w:
        return _from_wave(w)


def write_wav(path, a, rate=RATE):
    tmp = path + ".part"
    with wave.open(tmp, "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(rate)
        w.writeframes(a.tobytes())
    os.replace(tmp, path)


def chime():
    """Two soft notes (C5 then G5, exponential decay), for the start and the end."""
    n, off = int(RATE * 1.4), int(RATE * 0.22)
    out = array.array("h", [0]) * n
    for i in range(n):
        t = i / float(RATE)
        v = 0.18 * math.exp(-3 * t) * math.sin(2 * math.pi * 523.25 * t)
        if i >= off:
            u = (i - off) / float(RATE)
            v += 0.16 * math.exp(-3 * u) * math.sin(2 * math.pi * 783.99 * u)
        out[i] = int(32767 * v * min(1.0, i / (RATE * 0.005)))
    return out


def normalize(a, target_rms=3300.0, peak_cap=29000.0, max_gain=8.0):
    """Speech to about -20 dBFS RMS (measured on the voiced samples), never past -1 dBFS peak,
    so chunks from separate TTS calls play at one loudness."""
    if not a:
        return a
    voiced = [x for x in a[::4] if x > 400 or x < -400]
    if not voiced:
        return a
    rms = math.sqrt(sum(x * x for x in voiced) / float(len(voiced)))
    peak = float(max(max(a), -min(a))) or 1.0
    g = min(target_rms / rms, peak_cap / peak, max_gain)
    if abs(g - 1.0) < 0.03:
        return a
    out = array.array("h")
    for i in range(0, len(a), 65536):
        out.extend(array.array("h", [int(x * g) for x in a[i:i + 65536]]))
    return out


def encode_aac(wav_path, m4a_path):
    """-> ("<tool> AAC <n> kbps", None) or (None, why). afconvert (macOS) first, then ffmpeg;
    96 kbps first, then 64: afconvert refuses anything above 64 kbps for 24 kHz mono AAC
    ("Couldn't set audio converter property", measured on macOS 26)."""
    tried = []
    for tool in ("afconvert", "ffmpeg"):
        if not shutil.which(tool):
            continue
        for kbps in (96, 64):
            cmd = (["afconvert", "-f", "m4af", "-d", "aac", "-b", str(kbps * 1000), wav_path,
                    m4a_path] if tool == "afconvert" else
                   ["ffmpeg", "-y", "-loglevel", "error", "-i", wav_path, "-c:a", "aac", "-b:a",
                    f"{kbps}k", m4a_path])
            try:
                r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
            except (OSError, subprocess.SubprocessError) as e:
                tried.append(f"{tool} {kbps} kbps: {e}")
                continue
            if r.returncode == 0 and os.path.isfile(m4a_path) and os.path.getsize(m4a_path) > 1000:
                return f"{tool} AAC {kbps} kbps", None
            tried.append(f"{tool} {kbps} kbps: exit {r.returncode} "
                         f"{(r.stderr or r.stdout).strip()[:160]}")
    return None, ("; ".join(tried) if tried else "neither afconvert nor ffmpeg is installed")


# ----------------------------------------------------------------------------- chunks
class Chunk(object):
    def __init__(self, seg, turns):
        self.seg, self.turns = seg, turns
        self.words = sum(t.words for t in turns)


def make_chunks(segments, limit=CHUNK_WORDS):
    """One chunk per segment; a segment over `limit` words is split at turn boundaries into
    near-equal parts (a single turn is never split)."""
    out = []
    for si, seg in enumerate(segments):
        total = sum(t.words for t in seg)
        n = max(1, -(-total // limit))
        target = total / float(n)
        cur, cw, made = [], 0, 0
        for t in seg:
            if cur and made < n - 1 and cw + t.words / 2.0 > target:
                out.append(Chunk(si, cur))
                made += 1
                cur, cw = [], 0
            cur.append(t)
            cw += t.words
        if cur:
            out.append(Chunk(si, cur))
    return out


# ----------------------------------------------------------------------------- TTS backends
class RenderError(Exception):
    pass


class SwitchModel(Exception):
    pass


class EmptyAudio(RenderError):
    """A response with no audio in it: retried once, unlike an unreadable one."""


class RequestTimeout(Exception):
    """A TTS request with no answer within --request-timeout seconds."""


class RestartRender(Exception):
    """The model stopped answering and the render moved to another: every chunk is made again
    with the new model, so the episode never mixes two models' voices."""


def gemini_client(genai, types, key, timeout_s):
    """genai.Client with a per-request HTTP timeout and no silent SDK retries -> (client, record).
    Read from google-genai 2.31.0's source (HIVE, 2026-10-16): HttpOptions.timeout is in
    milliseconds; models.generate_content retries only when retry_options is set; the interactions
    client (Gemini 3 TTS) retries by itself -- 3 retries when retry_options is unset, with a
    server's Retry-After obeyed while nothing is logged -- and takes HttpRetryOptions.attempts as
    its retry count, so attempts=1 leaves it one retry, the fewest this SDK allows (attempts=0 is
    turned into 1 when the client is built). timed_request's deadline bounds every request either
    way. An SDK without these options keeps its own defaults, and the record says so."""
    want = {"timeout_ms": int(timeout_s * 1000), "retry_attempts": 1}
    try:
        opts = types.HttpOptions(timeout=want["timeout_ms"],
                                 retry_options=types.HttpRetryOptions(attempts=1))
    except (AttributeError, TypeError, ValueError) as e:
        return genai.Client(api_key=key), {
            "set": False, "note": f"this google-genai has no HttpOptions / HttpRetryOptions "
                                  f"({type(e).__name__}): its own timeout and retries apply"}
    try:
        return genai.Client(api_key=key, http_options=opts), dict(want, set=True)
    except TypeError as e:
        return genai.Client(api_key=key), {
            "set": False, "note": f"this google-genai's Client takes no http_options ({e}): "
                                  "its own timeout and retries apply"}


class GeminiTTS(object):
    """Gemini multi-speaker TTS, one request per chunk. Model names are read from the API
    (models.list), never assumed. Default order: a pro TTS model, then gemini-3.8-flash-tts,
    then gemini-2.5-flash-preview-tts. A model missing (404), without quota, or rejecting the
    request before it has made a chunk hands over to the next. Once a model has made a chunk
    the episode stays on it, because a second model's voices sound different -- unless it stops
    answering (preflight, stall_switch): then the next model makes every chunk again.

    Two APIs, chosen from the model name (https://ai.google.dev/gemini-api/docs/speech-generation,
    read 2026-09-25):
      2.x   models.generate_content: one text prompt ("TTS the following conversation ...:" +
            "Maya: ..." lines) and a MultiSpeakerVoiceConfig. Returns raw 24 kHz PCM.
      3.x+  interactions.create: one text part per turn, each annotated with speech_metadata
            {speaker, style}; speech_config {"mode": "conversational", "speakers": [...]}.
            Returns base64 audio/wav. generate_content fails on these models with 400
            "Multi-speaker generation requests must specify speaker names for each part".

    Rate limits: a 429 on a per-minute quota (no "PerDay" quotaId), or one carrying a
    retryDelay, waits that long + 5 s and retries the same chunk, up to RATE_WAITS times. A
    daily or zero quota stops the render with the resume message. Requests to one model are
    at least --min-interval seconds apart (default 20)."""
    name, cloud = "gemini", True

    def __init__(self, s, a):
        self.hosts = [h for h, _ in s.hosts]
        self.voices = dict(s.gemini_voices)
        self.styles = {h: s.styles.get(h) or HOST_STYLES[min(i, len(HOST_STYLES) - 1)]
                       for i, h in enumerate(self.hosts)}
        # The 2.x prompt's opening line. Unchanged from the first release, so a resumed render
        # finds its cached chunks (the cache key hashes this text).
        self.style = s.style or (
            f"TTS the following conversation between {self.hosts[0]} and {self.hosts[1]}, two "
            f"hosts of a lively, natural science podcast. {self.hosts[0]} sounds "
            f"{self.styles[self.hosts[0]]}; {self.hosts[1]} sounds {self.styles[self.hosts[1]]}")
        try:
            from google import genai
            from google.genai import types
        except ImportError:
            raise RenderError("[INFO] the gemini backend needs the google-genai package: "
                              "python3 -m pip install google-genai (or use --tts say)")
        self.key = read_key(core=bool(getattr(a, "core_run", None)))
        # which key the cache's TIER_FILE speaks for: a fingerprint, never the key
        self.key_id = hashlib.sha256(self.key.encode("utf-8")).hexdigest()[:16]
        self.types = types
        self.request_timeout = float(getattr(a, "request_timeout", None) or REQUEST_TIMEOUT)
        self.client, self.http = gemini_client(genai, types, self.key, self.request_timeout)
        if not self.http.get("set"):
            log(f"[INFO] {self.http['note']}")
        self.user_models = list(a.model or [])
        self.candidates = self.user_models or self.discover()
        self.model = self.candidates[0]
        self.min_interval = float(getattr(a, "min_interval", MIN_INTERVAL))
        # what the key is (note_tier): "prepaid" once a 402 shows prepay billing, "free" once a
        # 429 names a free-tier quota (never assumed), with the evidence; and whether a
        # hand-typed --cloud-ok note says otherwise ("paid") -- recorded, never a stop
        self.key_tier, self.key_tier_basis, self.tier_warning = None, None, None
        self.tier_evidence, self.prior_evidence, self.rate_logged = [], [], False
        note = getattr(a, "cloud_ok", None)
        self.paid_note = (note if isinstance(note, str) and re.search(r"\bpaid\b", note, re.I)
                          and not re.search(r"\b(?:un|not\s+|non-?)paid\b", note, re.I) else None)
        self.stalled, self.switches = [], []     # models that stopped answering; the hand-overs
        self.unusable = []                       # missing, no quota or rejected in this render
        self.remembered = {}                     # stalled in an earlier render (STALL_FILE)
        self.cache_dir = None
        self.preflight_record = None
        self.timeout_budget = 60.0 * float(getattr(a, "request_timeout_budget", None)
                                           or REQUEST_TIMEOUT_BUDGET_MIN)
        self.timed_out = 0.0                 # s lost to unanswered requests in this render
        self.wait_budget = 60.0 * float(getattr(a, "rate_wait_budget", RATE_WAIT_BUDGET_MIN))
        self.waited = 0.0
        self.pinned = False                  # set by prepare() / --redo: no model hand-over
        self.last_call = {}
        if not hasattr(self.client, "interactions"):
            three = [m for m in self.candidates if self.api(m) == "interactions"]
            if three and len(three) == len(self.candidates):
                raise RenderError(f"{', '.join(three)} need client.interactions, which this "
                                  "google-genai does not have: python3 -m pip install -U "
                                  "google-genai")
            if three:
                self.candidates = [m for m in self.candidates if m not in three]
                self.model = self.candidates[0]
                log(f"[WARN] this google-genai has no client.interactions, which Gemini 3 TTS "
                    f"needs: leaving out {', '.join(three)}. To use them: python3 -m pip install "
                    "-U google-genai")
        self.produced, self.models_used, self.warnings = 0, [], []

    @staticmethod
    def api(model):
        """'interactions' for Gemini 3 and later TTS models, 'generate_content' for 2.x."""
        m = re.match(r"(?:models/)?gemini-(\d+)", model or "")
        return "interactions" if m and int(m.group(1)) >= 3 else "generate_content"

    def discover(self):
        try:
            names = []
            for m in self.client.models.list():
                n = (getattr(m, "name", "") or "").split("/")[-1]
                acts = getattr(m, "supported_actions", None) or []
                if "tts" in n.lower() and (not acts or "generateContent" in acts):
                    names.append(n)
        except Exception as e:
            log(f"[WARN] could not list the Gemini models ({scrub(e)[:200]}); trying "
                f"{', '.join(KNOWN_TTS_MODELS)}")
            names = list(KNOWN_TTS_MODELS)
        if not names:
            names = list(KNOWN_TTS_MODELS)

        def ver(n):
            m = re.search(r"gemini-(\d+(?:\.\d+)?)", n)
            return float(m.group(1)) if m else 0.0
        pro = sorted((n for n in names if "-pro" in n), key=ver, reverse=True)
        # the 3.x flash model the docs show for multi-speaker: not a preview, not lite
        flash3 = sorted((n for n in names if self.api(n) == "interactions" and "-flash" in n
                         and "lite" not in n and "preview" not in n), key=ver, reverse=True)
        out = pro + flash3[:1] + [FALLBACK_TTS_MODEL]
        other = [n for n in names if n not in out]
        if other:
            log(f"[render] other TTS models, used only with --model: {', '.join(other)}")
        return out

    def remembered_stalls(self):
        """{model: record} of the models that stopped answering in this cache within
        STALL_MEMORY_H hours."""
        try:
            with open(os.path.join(self.cache_dir, STALL_FILE), encoding="utf-8") as fh:
                got = json.load(fh)
        except (OSError, ValueError, TypeError):
            return {}
        now = datetime.datetime.now(datetime.timezone.utc)
        out = {}
        for m, rec in (got.items() if isinstance(got, dict) else ()):
            try:
                at = datetime.datetime.fromisoformat(str(rec.get("at")))
                if at.tzinfo is None:
                    at = at.replace(tzinfo=datetime.timezone.utc)
            except (ValueError, AttributeError):
                continue
            if (now - at).total_seconds() < STALL_MEMORY_H * 3600:
                out[m] = rec
        return out

    def forget_stall(self, model):
        """`model` answers again: take it out of STALL_FILE (never fatal)."""
        self.remembered.pop(model, None)
        path = os.path.join(self.cache_dir or "", STALL_FILE)
        try:
            with open(path, encoding="utf-8") as fh:
                got = json.load(fh)
            if isinstance(got, dict) and got.pop(model, None) is not None:
                with open(path + ".part", "w", encoding="utf-8") as fh:
                    json.dump(got, fh, indent=2)
                os.replace(path + ".part", path)
        except (OSError, ValueError):
            pass

    def remember_stall(self, model, where, why):
        """Record in the cache that `model` stopped answering (never fatal)."""
        if not self.cache_dir:
            return
        path = os.path.join(self.cache_dir, STALL_FILE)
        try:
            with open(path, encoding="utf-8") as fh:
                got = json.load(fh)
            got = got if isinstance(got, dict) else {}
        except (OSError, ValueError):
            got = {}
        got[model] = {"at": now_iso(), "where": where, "why": why}
        try:
            with open(path + ".part", "w", encoding="utf-8") as fh:
                json.dump(got, fh, indent=2)
            os.replace(path + ".part", path)
        except OSError as e:
            log(f"[WARN] could not record that {model} stopped answering ({e})")

    def prepare(self, chunks, speak, cache_dir):
        """Stay on the model the cache was made with: resuming must not mix voices -- unless that
        model stopped answering within STALL_MEMORY_H hours: then it is not resumed with, and the
        model that is used makes every chunk again."""
        self.cache_dir = cache_dir
        self.prior_evidence = self._load_tier()
        self._settle_tier()
        if self.user_models:
            return
        recent = self.remembered_stalls()
        if all(m in recent for m in self.candidates):
            recent = {}                         # every model stalled: each gets another chance
        for m in recent:
            if m in self.candidates and m not in self.stalled:
                self.stalled.append(m)
        self.remembered = {m: r for m, r in recent.items() if m in self.candidates}
        have = {m: sum(os.path.isfile(os.path.join(cache_dir, cache_key(self, c, speak, m) +
                                                   ".wav")) for c in chunks)
                for m in self.candidates}
        for m in [m for m in self.candidates if m in recent and have[m]]:
            log(f"[render] not resuming with {m}: it stopped answering at {recent[m].get('at')} "
                f"({recent[m].get('where')}); the model used makes every chunk again, so the "
                "episode has one model's voices")
        usable = [m for m in self.candidates if m not in recent]
        self.candidates = usable + [m for m in self.candidates if m not in usable]
        self.model = self.candidates[0]
        best = max(usable, key=lambda m: have[m])
        if have[best]:
            self.candidates.remove(best)
            self.candidates.insert(0, best)
            self.model = best
            self.pinned = True                   # a resume never hands over to other voices
            log(f"[render] resuming with {best}: {have[best]} of {len(chunks)} chunk(s) cached")

    def prompt(self, chunk, speak):
        return self.style + ":\n\n" + "\n".join(f"{t.speaker}: {speak(t.text)}" for t in chunk.turns)

    def parts(self, chunk, speak):
        """The 3.x input: one text part per turn, annotated with its speaker and style."""
        return [{"type": "text", "text": speak(t.text),
                 "annotations": [{"type": "speech_metadata", "speaker": t.speaker,
                                  "style": self.styles[t.speaker]}]} for t in chunk.turns]

    def speech_config(self):
        return {"mode": "conversational",
                "speakers": [{"speaker": h, "voice": self.voices[h]} for h in self.hosts]}

    def payload(self, chunk, speak, model=None):
        model = model or self.model
        if self.api(model) == "interactions":
            return {"backend": "gemini", "api": "interactions", "model": model, "rate": RATE,
                    "speech_config": self.speech_config(), "parts": self.parts(chunk, speak)}
        return {"backend": "gemini", "model": model, "rate": RATE,
                "voices": [[h, self.voices[h]] for h in self.hosts],
                "text": self.prompt(chunk, speak)}

    def config(self):
        T = self.types
        return T.GenerateContentConfig(
            response_modalities=["AUDIO"],
            speech_config=T.SpeechConfig(multi_speaker_voice_config=T.MultiSpeakerVoiceConfig(
                speaker_voice_configs=[T.SpeakerVoiceConfig(
                    speaker=h, voice_config=T.VoiceConfig(
                        prebuilt_voice_config=T.PrebuiltVoiceConfig(voice_name=self.voices[h])))
                    for h in self.hosts])))

    @staticmethod
    def audio(r):
        for cand in (getattr(r, "candidates", None) or []):
            for p in (getattr(getattr(cand, "content", None), "parts", None) or []):
                d = getattr(p, "inline_data", None)
                if d is not None and getattr(d, "data", None):
                    data = d.data
                    if isinstance(data, str):
                        data = base64.b64decode(data)
                    return decode_audio(data, getattr(d, "mime_type", "") or "")
        why = [str(getattr(c, "finish_reason", "")) for c in (getattr(r, "candidates", None) or [])]
        raise EmptyAudio(f"the response held no audio (finish reason: {', '.join(why) or 'none'})")

    @staticmethod
    def interaction_audio(r):
        """interaction.output_audio.data: base64 audio/wav (24 kHz mono 16-bit) by default;
        the RIFF header is read, not spliced, so the rate comes from the file itself."""
        out = getattr(r, "output_audio", None)
        data = getattr(out, "data", None) if out is not None else None
        if not data:
            raise EmptyAudio(f"the interaction held no audio (status: "
                             f"{getattr(r, 'status', None) or 'not given'})")
        if isinstance(data, str):
            data = data.encode("ascii", "replace")
        data = bytes(data)
        if data[:4] != b"RIFF" and re.fullmatch(rb"[A-Za-z0-9+/=\s]+", data):
            try:                                         # base64, line-wrapped or not
                data = base64.b64decode(data, validate=False)
            except (binascii.Error, ValueError) as e:
                raise ValueError(f"the audio is not valid base64 ({e})")
        mime = getattr(out, "mime_type", None) or ""
        rate = getattr(out, "sample_rate", None)
        if rate and "rate=" not in str(mime):
            mime = f"{mime};rate={rate}"
        return decode_audio(bytes(data), str(mime))

    def request(self, chunk, speak):
        """One billed request; the answer is parsed by parse(), outside the retry loop."""
        if self.api(self.model) == "interactions":
            return self.client.interactions.create(
                model=self.model, input=[{"type": "user_input", "content": self.parts(chunk, speak)}],
                response_format={"type": "audio"},
                generation_config={"speech_config": self.speech_config()})
        return self.client.models.generate_content(model=self.model,
                                                   contents=self.prompt(chunk, speak),
                                                   config=self.config())

    def parse(self, r):
        return (self.interaction_audio(r) if self.api(self.model) == "interactions" else
                self.audio(r))

    def timed_request(self, chunk, speak, label, timeout=None):
        """request() in a worker thread, waited for at most `timeout` (request_timeout) seconds,
        with a log line when it is sent and every HEARTBEAT seconds while it is out: a hung request
        shows in the log and ends as RequestTimeout, never an endless wait. The abandoned thread
        is a daemon, so it cannot keep the process alive."""
        timeout = self.request_timeout if timeout is None else timeout
        box = {}

        def work():
            try:
                box["r"] = self.request(chunk, speak)
            except BaseException as e:            # handed to the caller as it was raised
                box["e"] = e
        t = threading.Thread(target=work, name="tts-request", daemon=True)
        start = _now()
        log(f"[render] {label}: sent to {self.model} (an answer is waited for up to "
            f"{timeout:.0f} s)")
        t.start()
        while True:
            left = timeout - (_now() - start)
            t.join(max(0.0, min(HEARTBEAT, left)))
            if not t.is_alive():
                break
            waited = _now() - start
            if waited >= timeout:
                raise RequestTimeout(f"no answer from {self.model} in {waited:.0f} s")
            log(f"[render] {label}: still waiting for {self.model} ({waited:.0f} s)")
        if "e" in box:
            raise box["e"]
        return box["r"]

    def _has_next(self):
        return (self.model in self.candidates
                and self.candidates.index(self.model) + 1 < len(self.candidates))

    def _next_model(self):
        """The next candidate after this model that has not stopped answering and was not found
        missing, out of quota or rejected in this render, or None."""
        if self.model not in self.candidates:
            return None
        return next((m for m in self.candidates[self.candidates.index(self.model) + 1:]
                     if m not in self.stalled and m not in self.unusable), None)

    def note_tier(self, err):
        """What an error says about the key's billing tier, as evidence: a 402 (prepay credits
        used up) shows PREPAID billing; a 429 that names a free-tier quota, the FREE tier. Said
        in the log when the tier it settles on changes, recorded (podcast.json key_tier,
        key_tier_evidence) and kept for this key in the cache (TIER_FILE). A --cloud-ok note that
        says "paid" gets one plain warning when the tier is free. Never a stop by itself: the
        Core's standing approval covers free-tier Core podcasts (2026-09-28)."""
        if err.get("kind") == "billing":
            tier = "prepaid"
            basis = (f"a 402 from {self.model}: the project's prepay credits are depleted, so it "
                     "is on prepaid (paid) billing")
        elif err.get("free_tier"):
            tier = "free"
            basis = f"a 429 from {self.model} named the free-tier quota {err['free_tier']}"
        else:
            return
        if any((e["tier"], e["model"]) == (tier, self.model) for e in self.tier_evidence):
            return
        was = self.key_tier
        self.tier_evidence.append({"tier": tier, "model": self.model, "at": now_iso(),
                                   "basis": basis})
        self._settle_tier()
        self._save_tier()
        if self.key_tier != was:
            log(f"[render] the Gemini key's Google project is "
                + ("on PREPAID billing" if self.key_tier == "prepaid" else "on the FREE tier")
                + f" ({self.key_tier_basis})")
            if self.tier_warning:
                log(f"[WARN] {self.tier_warning}")
        elif tier == "free" and self.key_tier == "prepaid":
            log(f"[render] {basis}; not taken as the key's tier: billing outvotes it "
                f"({self.key_tier_basis})")

    def _settle_tier(self):
        """key_tier from the evidence, this render's and this key's earlier renders' (TIER_FILE):
        "prepaid" when any shows prepay billing, else "free" when any names a free-tier quota,
        else None; the latest of that kind is the basis. A free-tier quota never outvotes
        billing: the staff render's project answered a 429 naming a free-tier quota (limit 0) and,
        later the same day, a 402 that showed prepaid billing."""
        every = self.prior_evidence + self.tier_evidence
        for tier in ("prepaid", "free"):
            got = [e for e in every if e["tier"] == tier]
            if got:
                e = got[-1]
                self.key_tier = tier
                self.key_tier_basis = e["basis"] + (
                    "" if e in self.tier_evidence else f" (an earlier render, {e['at']})")
                break
        else:
            self.key_tier, self.key_tier_basis = None, None
        self.tier_warning = None
        if self.key_tier == "free" and self.paid_note:
            self.tier_warning = (f"the --cloud-ok note says the key is paid ({self.paid_note!r}), but "
                                 f"its Google project is on the free tier ({self.key_tier_basis}): "
                                 "the note does not match the key; correct it for the next render")

    def key_tier_evidence(self):
        """Every answer the tier rests on, earlier renders' for this key first."""
        return self.prior_evidence + self.tier_evidence

    def _load_tier(self):
        """This key's evidence from TIER_FILE, the last TIER_MEMORY_DAYS of it ([] for another
        key or none)."""
        try:
            with open(os.path.join(self.cache_dir, TIER_FILE), encoding="utf-8") as fh:
                got = json.load(fh)
        except (OSError, ValueError, TypeError):
            return []
        if not isinstance(got, dict) or got.get("key_id") != self.key_id:
            return []
        now = datetime.datetime.now(datetime.timezone.utc)
        out = []
        for e in got.get("evidence") or []:
            try:
                at = datetime.datetime.fromisoformat(str(e["at"]))
                if at.tzinfo is None:
                    at = at.replace(tzinfo=datetime.timezone.utc)
                if (now - at).total_seconds() < TIER_MEMORY_DAYS * 86400 \
                        and e["tier"] in ("prepaid", "free"):
                    out.append({k: e[k] for k in ("tier", "model", "at", "basis")})
            except (KeyError, TypeError, ValueError):
                continue
        return out

    def _save_tier(self):
        """Keep this key's evidence in the cache (never fatal)."""
        if not self.cache_dir:
            return
        path = os.path.join(self.cache_dir, TIER_FILE)
        rec = {"key_id": self.key_id, "tier": self.key_tier, "basis": self.key_tier_basis,
               "evidence": self.key_tier_evidence()[-20:]}
        try:
            with open(path + ".part", "w", encoding="utf-8") as fh:
                json.dump(rec, fh, indent=2)
            os.replace(path + ".part", path)
        except OSError as e:
            log(f"[WARN] could not record the key's tier in the cache ({e})")

    def stall_switch(self, label, why):
        """This model stopped answering: hand over to the next candidate, which makes EVERY chunk
        again (RestartRender), or None when there is no other."""
        was = self.model
        if was not in self.stalled:
            self.stalled.append(was)
        self.remember_stall(was, label, why)
        nxt = self._next_model()
        if not nxt:
            return None
        note = (f"{was} stopped answering ({label}: {why}); switched to {nxt}, which made every "
                "chunk again, so the episode has one model's voices")
        self.switches.append({"from": was, "to": nxt, "at": label, "why": why, "note": note})
        log(f"[render] {was} stopped answering ({label}: {why}). Switching to {nxt} and making "
            f"EVERY chunk again with it, so the episode never mixes two models' voices")
        self.model, self.produced, self.models_used = nxt, 0, []
        self.pinned = False        # the next model made nothing yet: a model it cannot use hands over
        return RestartRender(note)

    def switch_notes(self):
        """One line per hand-over in this render, for podcast.json's warnings."""
        return [x["note"] for x in self.switches]

    def preflight(self, cached_from=None):
        """Before the first chunk this render must make: a two-line test to each candidate in
        turn, from the current one, on a PREFLIGHT_TIMEOUT deadline; the first that answers is
        used. A rate limit counts as an answer (the normal loop waits it out); a missing model,
        no quota, a rejected request, an error or no answer moves on; a key problem stops. A
        model other than the one the cache was made with (`cached_from`) makes every chunk
        again. -> the model chosen; RenderError when none answers."""
        probe = Chunk(0, [Turn(self.hosts[0], PREFLIGHT_LINES[0], 0, 0),
                          Turn(self.hosts[1], PREFLIGHT_LINES[1], 0, 0)])
        start = self.model if self.model in self.candidates else self.candidates[0]
        order = [m for m in self.candidates[self.candidates.index(start):]
                 if m not in self.stalled and m not in self.unusable]
        # then, before giving up, the models an earlier render remembered as stalled: an outage
        # may be over (round 2: on a free-tier key the pro models never answer, so passing over
        # the remembered flash models left nothing)
        later = [m for m in self.candidates if m in self.remembered and m not in order
                 and m not in self.unusable]
        tried = []
        for m in order + later:
            if m in later:
                self.stalled = [x for x in self.stalled if x != m]
                log(f"[render] preflight: {m} stopped answering earlier "
                    f"({self.remembered[m].get('at')}); trying it once more")
            self.model = m
            self._pace()
            t0 = _now()
            try:
                self.parse(self.timed_request(probe, lambda x: x, "preflight",
                                              timeout=PREFLIGHT_TIMEOUT))
            except RequestTimeout:
                # not remembered: one 60 s miss in a short outage must not keep the model out
                # of every render for STALL_MEMORY_H (review of 2.13.2, round 2); only a stall
                # mid-episode is (stall_switch)
                self.stalled.append(m)
                tried.append(f"{m}: no answer in {PREFLIGHT_TIMEOUT:.0f} s")
                continue
            except EmptyAudio as e:
                tried.append(f"{m}: {e}")
                continue
            except Exception as e:
                err = error_info(e)
                self.note_tier(err)
                if err["kind"] == "billing":            # credits are per project: every model
                    raise RenderError(billing_stop("preflight", m, err))
                if err["kind"] == "fatal":
                    raise RenderError(f"preflight: {m}: {why_line(err, 300)}")
                if err["kind"] != "rate":
                    self.unusable.append(m)
                    tried.append(f"{m}: {err['kind']} ({why_line(err, 120)})")
                    continue
            secs = round(_now() - t0, 1)
            self.preflight_record = {"model": m, "seconds": secs, "skipped": tried}
            if m in later:
                log(f"[render] preflight: {m} stopped answering earlier; it answers now")
                self.forget_stall(m)
            if tried:
                log(f"[render] preflight: {'; '.join(tried)}; using {m}, which answered the "
                    f"two-line test in {secs:g} s")
            else:
                log(f"[render] preflight: {m} answered the two-line test in {secs:g} s")
            if cached_from and m != cached_from:
                note = (f"{cached_from} made this episode's cached chunks but did not answer the "
                        f"preflight; {m} makes every chunk again, so the episode has one model's "
                        "voices")
                self.switches.append({"from": cached_from, "to": m, "at": "preflight",
                                      "why": "; ".join(tried), "note": note})
                self.pinned = False              # none of its chunks are made yet
                log(f"[render] {note}")
            return m
        before = [m for m in self.candidates if m in self.remembered]
        raise RenderError(
            "no Gemini TTS model answered a two-line test (" + "; ".join(tried) + ")"
            + (f"; {', '.join(before)} also stopped answering in an earlier render "
               f"(podcast/.cache/{STALL_FILE})" if before else "")
            + ". The chunks already made are cached: re-run later"
            + (" (this key's project is on the free tier: its quotas are small)"
               if self.key_tier == "free" else "")
            + ", name a model with --model <name> (it is tried whatever happened before), or "
            "--tts say")

    def _advance(self):
        nxt = None if self.produced else self._next_model()
        if not nxt:
            return False
        self.model = nxt
        return True

    def _pace(self):
        """Keep requests to one model at least min_interval seconds apart."""
        last = self.last_call.get(self.model)
        if last is not None:
            gap = self.min_interval - (_now() - last)
            if gap > 0.5:
                log(f"[render] pacing: {gap:.0f}s before the next request to {self.model}")
                _sleep(gap)
        self.last_call[self.model] = _now()

    def _stop(self, label, err, why):
        return RenderError(
            f"{label}: {self.model} {why} ({why_line(err)}). The chunks already made are "
            f"cached: re-run the same command later to resume with {self.model}, or pass "
            "--model <another> to render EVERY chunk with that model (voices differ between "
            "models, so one episode never mixes them).")

    def synth(self, chunk, speak, label):
        errors, waits, empty, retried_len, timeouts = 0, 0, 0, False, 0
        while True:
            self._pace()
            try:
                r = self.timed_request(chunk, speak, label)
            except RequestTimeout as e:
                timeouts += 1
                self.timed_out += self.request_timeout
                if self.timed_out >= self.timeout_budget:
                    raise RenderError(
                        f"{label}: {e}; unanswered requests have cost {self.timed_out / 60:.0f} "
                        f"min in this render (--request-timeout-budget "
                        f"{self.timeout_budget / 60:g} min). The chunks already made are cached: "
                        "re-run the same command to resume")
                if timeouts >= REQUEST_TIMEOUTS:
                    switch = self.stall_switch(label, f"{timeouts} requests of "
                                                      f"{self.request_timeout:.0f} s unanswered")
                    if switch:
                        raise switch
                    raise RenderError(
                        f"{label}: {self.model} did not answer within {self.request_timeout:.0f} s, "
                        f"{timeouts} times ({e}). The chunks already made are cached: re-run the "
                        "same command to resume (--request-timeout to wait longer)")
                log(f"[render] {label}: {e}; asking once more")
                continue
            except Exception as e:
                err = error_info(e)
                self.note_tier(err)
                kind = err["kind"]
                if kind == "billing":
                    raise RenderError(billing_stop(label, self.model, err))
                if kind in ("missing", "exhausted", "rejected") and self.model not in self.unusable:
                    self.unusable.append(self.model)
                if kind in ("missing", "exhausted", "rejected") and not self.produced \
                        and not self.pinned:
                    was = self.model
                    if self._advance():
                        why = {"missing": "not available", "exhausted": "no quota",
                               "rejected": "request rejected"}[kind]
                        log(f"[render] {was}: {why} ({why_line(err)}); trying {self.model}")
                        raise SwitchModel()
                    raise RenderError(
                        f"{label}: no Gemini TTS model could be used (tried "
                        f"{', '.join(self.candidates)}; last: {why_line(err)}). Quotas are per Cloud "
                        "project (see https://aistudio.google.com/rate-limit): re-run later "
                        "(finished chunks are cached), use a paid-tier key, or --tts say.")
                if kind in ("missing", "exhausted", "rejected") and self.pinned and not self.produced:
                    raise self._stop(label, err, "made this episode's cached chunks, so the "
                                                 "render stays on it, but it failed")
                if kind == "exhausted":
                    raise self._stop(label, err, "is out of its daily quota")
                if kind == "rate":
                    waits += 1
                    wait = (err["delay"] if err["delay"] is not None else 60.0) + 5.0
                    if waits > RATE_WAITS:
                        raise self._stop(label, err, f"is still rate-limited after {RATE_WAITS} "
                                                     "waits for this chunk")
                    if self.waited + wait > self.wait_budget:
                        raise self._stop(label, err, f"has cost {self.waited / 60:.0f} min of "
                                                     f"rate-limit waits, and the next "
                                                     f"{wait:.0f} s would pass --rate-wait-budget "
                                                     f"({self.wait_budget / 60:g} min)")
                    self.waited += wait
                    if not self.rate_logged:
                        # what Google said, once: a misread error shows here, not after 11 min
                        self.rate_logged = True
                        log(f"[render] {label}: Google answered {err['code'] or '?'} "
                            f"{err['status'] or '?'}: {(err['message'] or err['msg'])[:300]}")
                    log(f"[render] rate-limited, waiting {wait:.0f}s ({label}, {self.model}"
                        + (f", {err['quota']}" if err["quota"] else "") + f"; wait {waits} of "
                        f"{RATE_WAITS})")
                    _sleep(wait)
                    continue
                errors += 1
                if kind in ("fatal", "rejected") or errors >= 4:
                    raise RenderError(f"{label} failed after {errors} error(s) on "
                                      f"{self.model}: {why_line(err, 300)}. The chunks already made are "
                                      "cached; re-run the same command to resume.")
                wait = err["delay"] if err["delay"] is not None else min(60.0, 5.0 * 2 ** (errors - 1))
                log(f"[render] {label}: {kind} error on {self.model} ({why_line(err, 120)}); retrying in "
                    f"{wait:.0f} s")
                _sleep(wait)
                continue
            # Parsed outside the retry: the same answer to the same request would be billed again.
            try:
                pcm = self.parse(r)
            except EmptyAudio as e:
                empty += 1
                if empty > 1:
                    raise RenderError(f"{label}: {self.model} returned no audio twice ({e}); the "
                                      "chunks already made are cached")
                log(f"[WARN] {label}: {e}; asking once more")
                continue
            except (ValueError, struct.error, binascii.Error, EOFError) as e:
                raise RenderError(f"{label}: {self.model} answered with audio that could not be "
                                  f"read ({type(e).__name__}: {e}); not retried, since the same "
                                  "request would be billed again")
            dur, expect = len(pcm) / float(RATE), chunk.words * 60.0 / WPM
            if expect > 5 and not 0.45 <= dur / expect <= 2.5:
                note = (f"{label}: {dur:.0f} s of audio for about {expect:.0f} s of text on "
                        f"{self.model}")
                if not retried_len:
                    retried_len = True
                    log(f"[WARN] {note}; asking once more")
                    continue
                self.warnings.append(note + " (kept; listen to this part)")
                log(f"[WARN] {note}; kept")
            if self.model not in self.models_used:
                self.models_used.append(self.model)
            return pcm


def why_line(err, n=200):
    """An error for a log line with the useful part first: the quota id and the retry delay,
    which the SDK's message buries after a long preamble, then the message itself."""
    parts = [f"quota {err['quota']}"] if err.get("quota") else []
    if err.get("delay") is not None:
        parts.append(f"retry in {err['delay']:g}s")
    return "; ".join(parts + [err["msg"][:n]])


def billing_stop(where, model, err):
    """The message a 402 stops a render with: what it is, what clears it, what is kept --
    with Google's own words, which carry its link."""
    said = (f"{err.get('code') or '?'} {err.get('status') or '?'}: {err['message'][:300]}"
            if err.get("message") else why_line(err, 300))
    return (f"{where}: Gemini prepay credits are depleted ({model}: {said}). "
            "Waiting cannot clear this, and every model draws on the same project: add prepay "
            "credits to the key's Google project in AI Studio (its Billing page), then re-run the "
            "same command. The chunks already made are cached.")


def _error_details(e):
    """The structured part of an SDK error: google-genai's APIError keeps the response JSON in
    .details; the interactions client's errors keep it in .body."""
    out = []
    for attr in ("details", "body", "response_json"):
        v = getattr(e, attr, None)
        if isinstance(v, (bytes, bytearray)):
            v = v.decode("utf-8", "replace")
        if isinstance(v, str):
            try:
                v = json.loads(v)
            except ValueError:
                pass
        if v is not None and not callable(v):
            out.append(v)
    return out


def _seconds(v):
    m = re.match(r"\s*(\d+(?:\.\d+)?)\s*s?\s*$", str(v))
    return float(m.group(1)) if m else None


# A 402 / prepay credits used up. Never the bare word "billing": Google's ordinary 429 says
# "please check your plan and billing details".
_PREPAY = re.compile(r"\bprepay(?:ment)?\b|\bcredits?\s+(?:are|is|have\s+been)\s+(?:depleted|"
                     r"exhausted|used\s+up)\b", re.I)
_RPC_STATUS = re.compile(r"\b(RESOURCE_EXHAUSTED|PERMISSION_DENIED|INVALID_ARGUMENT|NOT_FOUND|"
                         r"UNAVAILABLE|INTERNAL|DEADLINE_EXCEEDED|UNAUTHENTICATED|"
                         r"FAILED_PRECONDITION)\b")


def error_info(e):
    """-> {kind, msg, delay, quota, free_tier, code, status, message}. kind: missing | billing (402:
    the project's prepay credits are used up -- no wait clears it) | exhausted (daily or zero
    quota) | rate (per-minute, or any 429 that gives a retryDelay) | rejected (400: the request
    does not fit this model) | fatal (key, permission) | transient. The retry delay and the
    quota come from google.rpc.RetryInfo / QuotaFailure in the error details when present, else
    from the text."""
    msg = scrub(f"{type(e).__name__}: {e}")
    code = next((c for c in (getattr(e, "code", None), getattr(e, "status_code", None))
                 if isinstance(c, int)), None)
    found = {"ids": [], "values": [], "delay": None, "metrics": [], "status": None,
             "message": None}

    def walk(x):
        if isinstance(x, dict):
            for k, v in x.items():
                if k == "quotaId" and isinstance(v, str):
                    found["ids"].append(v)
                elif k == "quotaMetric" and isinstance(v, str):
                    found["metrics"].append(v.rsplit("/", 1)[-1])
                elif k == "quotaValue":
                    found["values"].append(str(v))
                elif k == "retryDelay" and found["delay"] is None:
                    found["delay"] = _seconds(v)
                elif k == "status" and isinstance(v, str):
                    found["status"] = found["status"] or v
                elif k == "message" and isinstance(v, str):
                    found["message"] = found["message"] or v
                else:
                    walk(v)
        elif isinstance(x, (list, tuple)):
            for v in x:
                walk(v)
    details = _error_details(e)
    walk(details)
    text = msg + " " + scrub(json.dumps(details, default=str))[:6000]
    if code is None:
        m = re.search(r"\b(4\d\d|5\d\d)\b", text)
        code = int(m.group(1)) if m else None
    if found["delay"] is None:
        m = re.search(r"retry(?:Delay)?['\"]?\s*(?:in|:)?\s*['\"]?(\d+(?:\.\d+)?)\s*s", text, re.I)
        found["delay"] = float(m.group(1)) if m else None
    if found["delay"] is not None:
        found["delay"] = min(max(found["delay"], 1.0), 300.0)
    quota = ", ".join(dict.fromkeys(found["ids"]))
    # The key's project is on Google's FREE tier when a quota it hit says so: quotaMetric
    # generate_content_free_tier_requests, quotaId ...-FreeTier (render records key_tier)
    tier = [q for q in found["metrics"] + found["ids"] if re.search(r"free[_-]?tier", q, re.I)]
    if not tier and re.search(r"\w*free[_-]?tier\w*", text, re.I):
        tier = [re.search(r"\w*free[_-]?tier\w*", text, re.I).group(0)]
    status = found["status"] or next(iter(_RPC_STATUS.findall(text)), None)
    info = {"msg": msg, "delay": found["delay"], "quota": quota,
            "free_tier": ", ".join(dict.fromkeys(tier)) or None, "code": code, "status": status,
            "message": scrub(found["message"]) if found["message"] else None}
    if code == 404 or "NOT_FOUND" in text:
        return dict(info, kind="missing")
    # Before RESOURCE_EXHAUSTED, which a 402 also says: read as a rate limit, it waited 10 x 65 s
    # per chunk on an error no wait can clear (staff report, skill 2.11.2)
    if code == 402 or _PREPAY.search(text):
        return dict(info, kind="billing")
    if code == 429 or "RESOURCE_EXHAUSTED" in text:
        zero = "0" in found["values"] or re.search(r"limit:\s*0\b", text)
        daily = (any("perday" in q.lower() for q in found["ids"]) if found["ids"] else
                 bool(re.search(r"per ?day|daily", text, re.I)) and found["delay"] is None)
        return dict(info, kind="exhausted" if zero or daily else "rate")
    if code in (401, 403) or re.search(r"PERMISSION_DENIED|API key not valid|UNAUTHENTICATED",
                                       text):
        return dict(info, kind="fatal")
    if code == 400 or "INVALID_ARGUMENT" in text:
        return dict(info, kind="rejected")
    return dict(info, kind="transient")


def core_key():
    """The Core's shared key (CORE_KEY_FILE) -> (key, None); (None, None) when there is no such
    file here (a laptop, or not installed); (None, why it is not used). Used only while the
    Core's rule for its files on HIVE holds (slack_collab.trust_problem: /quobyte/proteomics-grp
    is group-writable and not sticky, so any member could put a folder with their own key in
    its place, and every Core episode would then go to their Gemini project) and no account
    outside the group can read it. Group-readable is the point of it, so that is no warning."""
    path = CORE_KEY_FILE
    if not os.path.lexists(path):
        return None, None
    try:
        import slack_collab as sc                  # the Core's trust rule and its admins
    except ImportError as e:
        return None, f"{path} is not used: the Core's trust rule cannot be checked ({e})"
    info = sc._look(path, True)
    if info.get("error"):
        return None, (f"{path} is not readable ({info['error']}): it is for members of "
                      "proteomics-grp")
    why = sc.trust_problem(info, sc._folder(path), sc.core_admins())
    if why:
        return None, f"{path} is not used: {why}"
    if int(info.get("mode") or 0) & 0o007:
        return None, f"{path} is not used: accounts outside its group can read it (chmod 640)"
    key = (info.get("text") or "").strip()
    return (key, None) if key else (None, f"{path} is empty")


def _key_mode_problem(path):
    """core_submission.key_mode_problem, the skill's one rule for a key file: on HIVE every home
    folder is open to every account, so a key file others can read is refused, not used.
    Elsewhere None (a laptop's home is its user's)."""
    try:
        from core_submission import key_mode_problem
    except ImportError as e:
        return f"its mode cannot be checked ({e})"
    return key_mode_problem(path)


def find_key(core=False):
    """-> (key, where it came from, why the Core's shared key was not used or None). The order:
    GEMINI_API_KEY; ~/.config/ucdavis-proteomics/gemini_key (on HIVE refused, not used, when
    other accounts can read it); for a Core run only, the Core's shared key (core_key). Never
    printed."""
    key = (os.environ.get("GEMINI_API_KEY") or "").strip()
    if key:
        return key, "GEMINI_API_KEY", None
    path = os.path.expanduser(KEY_FILE)
    if os.path.isfile(path):
        bad = _key_mode_problem(path)
        if bad:
            return None, None, (f"{KEY_FILE} is not used: {bad}. Make it yours alone: chmod 600 "
                                f"{KEY_FILE}; if others could read it for a while, make a new "
                                "key in Google AI Studio and save that one")
        mode = os.stat(path).st_mode & 0o777
        if mode & 0o077:
            log(f"[WARN] {KEY_FILE} is readable by other accounts (mode {oct(mode)[2:]}): "
                f"chmod 600 {KEY_FILE}")
        with open(path, encoding="utf-8") as fh:
            key = fh.read().strip()
        if key:
            return key, KEY_FILE, None
    if not core:
        return None, None, (f"{CORE_KEY_FILE} is used only for a Core run (--core-run)"
                            if os.path.lexists(CORE_KEY_FILE) else None)
    key, why = core_key()
    return (key, CORE_KEY_FILE, None) if key else (None, None, why)


def read_key(core=False):
    """The Gemini key (find_key), remembered for scrub(). `core`: a Core run, which may use the
    Core's shared key. Never printed."""
    key, _src, why = find_key(core)
    if not key:
        raise RenderError(f"no Gemini API key: set GEMINI_API_KEY or put it in {KEY_FILE} "
                          "(chmod 600)" + (f"; {why}" if why else ""))
    _SECRETS.append(key)
    return key


def say_voices():
    try:
        out = subprocess.run(["say", "-v", "?"], capture_output=True, text=True, timeout=60).stdout
    except (OSError, subprocess.SubprocessError):
        return []
    return [(m.group(1).strip(), m.group(2)) for m in
            re.finditer(r"(?m)^(.+?)\s+([a-z]{2}_[A-Z]{2,3})\s+#", out)]


class SayTTS(object):
    """macOS `say`, one call per turn, offline: nothing leaves the machine."""
    name, cloud = "say", False

    def __init__(self, s, a):
        if not shutil.which("say"):
            raise RenderError("--tts say needs macOS `say` (not found); use --tts gemini")
        installed = say_voices()
        names = [n for n, _ in installed]
        english = [n for n, loc in installed if loc.startswith("en_")]
        self.hosts = [h for h, _ in s.hosts]
        self.voices = {}
        for i, h in enumerate(self.hosts):
            v = s.say_voices.get(h)
            if names and v not in names:
                alt = [n for n in english if n not in self.voices.values()]
                log(f"[WARN] say voice {v!r} for {h} is not installed; using "
                    f"{alt[0] if alt else names[0]!r}")
                v = alt[0] if alt else names[0]
            self.voices[h] = v
        self.model = "macOS say"
        self.produced, self.models_used, self.warnings = 0, [self.model], []

    def prepare(self, chunks, speak, cache_dir):
        return

    def prompt(self, chunk, speak):
        return "\n".join(f"{t.speaker}: {speak(t.text)}" for t in chunk.turns)

    def payload(self, chunk, speak, model=None):
        return {"backend": "say", "rate": RATE, "gap": GAP_TURN,
                "voices": [[h, self.voices[h]] for h in self.hosts],
                "turns": [[t.speaker, speak(t.text)] for t in chunk.turns]}

    def synth(self, chunk, speak, label):
        out = array.array("h")
        with tempfile.TemporaryDirectory() as td:
            for k, t in enumerate(chunk.turns):
                txt, wav = os.path.join(td, "t.txt"), os.path.join(td, "t.wav")
                with open(txt, "w", encoding="utf-8") as fh:
                    fh.write(speak(t.text))
                # A turn takes ~0.6 s; now and then `say` stalls with no CPU (seen once in a
                # 135-turn render, not reproducible), so a stalled call is retried, not waited on.
                r = None
                for attempt in range(3):
                    try:
                        r = subprocess.run(["say", "-v", self.voices[t.speaker],
                                            "--file-format=WAVE", "--data-format=LEI16@24000",
                                            "-o", wav, "-f", txt],
                                           capture_output=True, text=True, timeout=60)
                        break
                    except subprocess.TimeoutExpired:
                        log(f"[WARN] {label}: say stalled on the turn at line {t.line}; retrying")
                if r is None:
                    raise RenderError(f"{label}: say stalled 3 times on the turn at line {t.line}; "
                                      "re-run to resume (finished chunks are cached)")
                if r.returncode or not os.path.isfile(wav):
                    raise RenderError(f"{label}: say failed on the turn at line {t.line}: "
                                      f"{r.stderr.strip()[:200]}")
                if k:
                    out.extend(silence(GAP_TURN))
                out.extend(read_wav(wav))
        return out


BACKENDS = {"gemini": GeminiTTS, "say": SayTTS}


def cache_key(backend, chunk, speak, model=None):
    blob = json.dumps(backend.payload(chunk, speak, model), sort_keys=True, ensure_ascii=False)
    return sha256_bytes(blob.encode("utf-8"))


# ----------------------------------------------------------------------------- consent
def _on_hive():
    """core_submission.on_hive, the skill's test (its CORE_ON_HIVE switch is for tests)."""
    try:
        from core_submission import on_hive
    except ImportError:
        return os.path.isdir("/quobyte")
    return on_hive()


def _uid_login():
    """This account's login, by its uid: never $USER or $LOGNAME, which anyone can set."""
    try:
        import pwd
        return pwd.getpwuid(os.getuid()).pw_name
    except (ImportError, KeyError, AttributeError):
        return None


def _staff_check_on_hive():
    """From a laptop: ask HIVE, through hive_exec.sh, whether the account this computer's HIVE
    login reaches is on the Core staff list -- staff.py's own check, of `id -un` there (what the
    ssh key proves), never a name this computer claims. -> (login, how, None) or (None, None,
    why)."""
    cmd = 'python3 ~/proteomics-pipeline/scripts/staff.py check --by "$(id -un)"'
    try:
        r = subprocess.run(["bash", os.path.join(HERE, "hive_exec.sh"), cmd],
                           capture_output=True, text=True, timeout=120)
        got = json.loads(r.stdout[r.stdout.find("{"):])      # past anything a login shell prints
    except (OSError, subprocess.SubprocessError, ValueError) as e:
        return None, None, ("the Core staff list is on HIVE and could not be asked from here "
                            f"({type(e).__name__}): this computer needs its HIVE login "
                            "(check_access.sh saves it) and the skill on HIVE "
                            "(hive_exec.sh --put-skill)")
    if not isinstance(got, dict):
        return None, None, "HIVE's answer about the Core staff list was not readable"
    checked = str(got.get("checked") or "")
    if got.get("ok") and checked.startswith("on the Core staff list"):
        return got.get("by"), checked + ", asked on HIVE", None
    return None, None, str(got.get("error") or checked or "HIVE did not confirm it")


def core_staff():
    """Is the person running this on the Core staff list? -> (login, how it was checked, None)
    or (None, None, why not). The Core's Gemini key is readable by every member of
    proteomics-grp, but the standing approval and the key are for Core staff: on HIVE this
    account (by uid) must be on the list staff.load_staff trusts (CORE_STAFF_FILE, owned by a
    Core admin); on a laptop HIVE is asked about the account its key logs in as. Fails closed:
    no usable list is no."""
    if not _on_hive():
        return _staff_check_on_hive()
    try:
        import staff
    except ImportError as e:
        return None, None, f"the Core staff list cannot be checked ({e})"
    login = _uid_login()
    logins, why = staff.load_staff()
    if logins is None:
        return None, None, f"the Core staff list cannot be used ({why})"
    if not login or login not in logins:
        return None, None, (f"{login or 'this account'} is not on the Core staff list "
                            f"({staff.staff_file()}); a Core admin adds staff to it")
    return login, f"on the Core staff list ({staff.staff_file()})", None


def core_run_consent(script_path, kind):
    """--core-run: the consent of record for a Core run -> (CORE_CLOUD_OK, the Core run, None),
    or (None, None, why not). `submission` (the default): the session the script sits in answers
    a CoreOmics submission (submission_report.core_run, the one test); `internal`: its
    session.json says it is internal Core work (analysis_for, from session.py init --internal).
    Either way the person running it must be Core staff (core_staff). The Core run names the
    ROLE; its "login" is for the staff-only record (record_core_run), never podcast.json.
    Anyone else's data is never a Core run: they are asked, and --cloud-ok records the answer."""
    out = output_dir_of(script_path)
    session = _session_of(out) if out else None
    work, got = None, None
    if session:
        try:
            import session as session_mod
            work = session_mod.analysis_for(session)
        except ImportError:
            work = None
        try:
            import submission_report
            got = submission_report.core_run(session)
        except ImportError:
            got = None
    named = [str(x) for x in ((work or {}).get("submissions") or [])] \
        if isinstance(work, dict) else []
    if kind == "internal":
        if got or named:              # a client's data is never labelled as the Core's own
            which = (got or {}).get("prot") or ", ".join(named) or "a CoreOmics submission"
            return None, None, (f"--core-run internal: this session is for {which}, a CoreOmics "
                                "submission, so its podcast is that submission's: pass "
                                "--core-run.")
        if not (isinstance(work, dict) and work.get("internal")):
            return None, None, ("--core-run internal: this session's session.json does not say "
                                "it is internal Core work (analysis_for, written by session.py "
                                "init --internal \"<its board title>\"). For a CoreOmics "
                                "submission pass --core-run; for anyone else's data, ask them "
                                "and pass --cloud-ok \"<who agreed, when>\".")
        rec = {"kind": "internal", "title": work["internal"]}
    else:
        if not got and named:
            return None, None, (f"--core-run: session.json says this session is for "
                                f"{', '.join(named)}, but no submission record is attached: run "
                                f"submission_report.py attach --session {session} --record "
                                f"<the submission's folder> (SKILL.md step 1c.7), then --core-run "
                                "again.")
        if not got:
            return None, None, ("--core-run: this session answers no CoreOmics submission "
                                "(submission_report.core_run, asked of the session above the "
                                "script's podcast/ folder). For internal Core work pass "
                                "--core-run internal; for anyone else's data, ask them and "
                                "pass --cloud-ok \"<who agreed, when>\".")
        rec = {"kind": "submission", "prot": got.get("prot")}
    login, how, why = core_staff()
    if not login:
        return None, None, f"--core-run is for Core staff: {why}."
    try:
        from staff import ROLE
    except ImportError:
        ROLE = "Core staff"
    rec.update(by=ROLE, staff_check=how, login=login)
    return CORE_CLOUD_OK, rec, None


def record_core_run(script_path, rec, event):
    """Take the login out of a Core run (`rec`, as core_run_consent made it) and append who it
    was to the session's staff-only record, logs/podcast.staff.json: names and logins stay with
    the staff (staff.py), and podcast.json says "Core staff". -> `rec` without the login."""
    rec = dict(rec or {})
    login = rec.pop("login", None)
    if not login:
        return rec or None
    out = output_dir_of(script_path)
    where = os.path.join(_session_of(out) if out else os.path.dirname(os.path.abspath(
        script_path)), "logs", "podcast.staff.json")
    try:
        import staff
        staff.append(where, dict(event=event, by=login, staff_check=rec.get("staff_check"),
                                 at=now_iso(), **{k: v for k, v in rec.items()
                                                  if k in ("kind", "prot", "title")}))
    except Exception as e:                      # the podcast never stops for its record
        log(f"[WARN] the staff-only record of who made this Core run's podcast was not written "
            f"({type(e).__name__}: {e})")
    return rec


def consent_of(a):
    """The consent a render, verify or job sends under -> (consent, the Core run or None, why
    there is none or None): --core-run records CORE_CLOUD_OK; --cloud-ok a note (or True); a
    refusal word ("no", "declined" ...) is no consent."""
    kind = getattr(a, "core_run", None)
    if kind:
        return core_run_consent(a.script, kind)
    consent = getattr(a, "cloud_ok", False)
    if isinstance(consent, str) and consent.strip().lower() in NO_CONSENT:
        return False, None, f"--cloud-ok {consent!r} is a refusal."
    return consent, None, None


# ----------------------------------------------------------------------------- render
def render_chunks(backend, chunks, speak, cache):
    """Every chunk, cached or made now -> (pieces, n_cached, n_new, used cache keys, warnings).
    RestartRender (the model stopped answering and the next one takes over) propagates: the
    caller starts again from the first chunk, with the new model."""
    pieces, n_cached, n_new, used, warnings = [], 0, 0, set(), []
    for i, c in enumerate(chunks, 1):
        label = f"chunk {i}/{len(chunks)} (segment {c.seg + 1})"
        while True:
            key = cache_key(backend, c, speak)
            path = os.path.join(cache, key + ".wav")
            side = os.path.join(cache, key + ".json")
            if os.path.isfile(path):
                try:
                    pcm = read_wav(path)
                    n_cached += 1
                    backend.produced += 1
                    if backend.model not in backend.models_used:
                        backend.models_used.append(backend.model)
                    try:
                        with open(side, encoding="utf-8") as fh:
                            warnings += [str(w) for w in (json.load(fh).get("warnings") or [])]
                    except (OSError, ValueError, AttributeError):
                        pass                             # a chunk cached with no warnings
                    log(f"[render] {label}: cached")
                    break
                except (wave.Error, EOFError, OSError, ValueError):
                    os.remove(path)
            before = len(backend.warnings)
            try:
                pcm = backend.synth(c, speak, label)
            except SwitchModel:
                continue
            new_w = backend.warnings[before:]
            write_wav(path, pcm)
            with open(side, "w", encoding="utf-8") as fh:
                json.dump({"label": label, "model": backend.model, "warnings": new_w}, fh)
            warnings += new_w
            n_new += 1
            backend.produced += 1
            log(f"[render] {label}: {len(pcm) / float(RATE):.0f} s from {backend.model}")
            break
        used.add(key)
        pieces.append((c.seg, pcm))
    return pieces, n_cached, n_new, used, warnings


def cmd_render(a):
    s = parse_script(a.script)
    if s.problems:
        for line, msg in s.problems:
            log(f"[FAIL] {'line %d: ' % line if line else ''}{msg}")
        return 2
    backend_cls = BACKENDS[a.tts]
    speak = pronouncer(s.pronunciation)
    chunks = make_chunks(s.segments)
    words = sum(c.words for c in chunks)
    if a.dry_run:                                  # a preview sends nothing: no gate needed
        for i, c in enumerate(chunks, 1):
            print(f"===== chunk {i}/{len(chunks)} (segment {c.seg + 1}, {c.words} words)")
            print("\n".join(f"{t.speaker}: {speak(t.text)}" for t in c.turns))
        return 0
    status, sources, why = read_check(s)
    if why and not a.unchecked:
        log(f"[render] refusing: {why}. Run `make_podcast.py check {a.script} --source <the "
            "report files>` and fix every flagged line first (or pass --unchecked, which "
            "podcast.json then records).")
        return 2
    consent, core_rec, refusal = consent_of(a) if backend_cls.cloud else (False, None, None)
    if backend_cls.cloud and not consent:
        if refusal:
            log(f"[render] not sending anything: {refusal}")
        if not getattr(a, "core_run", None):
            log(f"[render] not sending anything. --tts {a.tts} sends the transcript -- "
                f"{words:,} words of the turns, after pronunciation substitutions; never the "
                f"report -- to Google's Gemini API, and verify then sends the rendered AUDIO "
                f"(downsampled) for transcription. On a free-tier key Google may use both to "
                f"improve its products and human reviewers may read it ({TERMS_URL}). Once the "
                f"user has agreed, re-run with --cloud-ok \"who agreed, when\" (a Core run: "
                f"--core-run, the Core director's standing approval), or use --tts say "
                f"(offline, macOS).")
        return 2
    m_status, m_why = read_meaning(s)
    m_override = " ".join((getattr(a, "no_meaning_check", None) or "").split())
    if m_why and not m_override:
        log(f"[render] refusing: {m_why}. A reader who did not write the script checks every "
            f"sentence against the report first: `make_podcast.py meaning {a.script} --brief`, "
            f"then `make_podcast.py meaning {a.script} --verdict PASS --reviewer \"<who>\" "
            "--findings <file>` (or pass --no-meaning-check \"<why>\", which podcast.json "
            "records).")
        return 2
    out = os.path.abspath(a.out or os.path.dirname(s.path))
    cache = os.path.join(out, ".cache")
    os.makedirs(cache, exist_ok=True)
    backend = backend_cls(s, a)
    backend.prepare(chunks, speak, cache)
    redo = set(getattr(a, "redo", None) or [])
    for c in chunks:                              # re-make these segments, text unchanged
        if c.seg + 1 in redo:
            for ext in (".wav", ".json"):
                f = os.path.join(cache, cache_key(backend, c, speak) + ext)
                if os.path.exists(f):
                    os.remove(f)
    if redo:
        log(f"[render] re-making segment(s) {', '.join(str(x) for x in sorted(redo))}")
        if hasattr(backend, "pinned") and not getattr(a, "model", None):
            backend.pinned = True                 # the re-made chunks join the others' voices
    log(f"[render] {len(chunks)} chunk(s), {words:,} words -> about {words / float(WPM):.0f} min "
        f"with {backend.name}")
    # A two-line test before the first chunk this render must make (PREFLIGHT_TIMEOUT): a model
    # that does not answer is passed over before the episode depends on it. Nothing to make, no
    # test; --no-preflight skips it (and podcast.json says so).
    preflight = {"run": False, "why": ("every chunk is cached" if hasattr(backend, "preflight")
                                       else f"no two-line test for {backend.name} (offline)")}
    if hasattr(backend, "preflight"):
        todo = [c for c in chunks if not os.path.isfile(os.path.join(
            cache, cache_key(backend, c, speak) + ".wav"))]
        if getattr(a, "no_preflight", False):
            preflight = {"run": False, "why": "--no-preflight"}
        elif todo:
            cached_from = backend.model if len(todo) < len(chunks) else None
            backend.preflight(cached_from=cached_from)
            preflight = dict(run=True, **(backend.preflight_record or {}))

    # Each chunk: <key>.wav, and <key>.json with the warnings made when it was synthesized, so
    # a resumed render still reports them. A model that stops answering mid-episode hands over
    # (RestartRender), and every chunk is made again with the next one.
    while True:
        try:
            pieces, n_cached, n_new, used, warnings = render_chunks(backend, chunks, speak, cache)
            break
        except RestartRender:
            continue
    warnings = list(backend.switch_notes()) + warnings if hasattr(backend, "switch_notes") \
        else warnings
    if getattr(backend, "tier_warning", None):
        warnings.append(backend.tier_warning)

    audio = array.array("h")
    audio.extend(chime())
    audio.extend(silence(0.5))
    prev = None
    for seg, pcm in pieces:
        if prev is not None:
            audio.extend(silence(GAP_SEGMENT if seg != prev else GAP_CHUNK))
        audio.extend(normalize(pcm))
        prev = seg
    audio.extend(silence(0.6))
    audio.extend(chime())
    duration = round(len(audio) / float(RATE), 1)

    wav = os.path.join(out, "podcast.wav")
    m4a = os.path.join(out, "podcast.m4a")
    write_wav(wav, audio)
    tool, why_not = encode_aac(wav, m4a)
    if tool:
        audio_name = "podcast.m4a"
        if not a.keep_wav:
            os.remove(wav)
    else:
        audio_name = "podcast.wav"
        if os.path.exists(m4a):
            os.remove(m4a)                             # an older render's; not this audio
        log(f"[INFO] no AAC encoder ({why_not}); keeping {audio_name}")

    script_copy = os.path.join(out, "podcast_script.md")
    if os.path.abspath(s.path) != script_copy:
        shutil.copyfile(s.path, script_copy)
        for name in ("check.txt", MEANING_FILE):      # both name the script's sha and sources
            src = os.path.join(os.path.dirname(s.path), name)
            if not os.path.isfile(src) or os.path.abspath(src) == os.path.join(out, name):
                continue
            with open(src, encoding="utf-8") as fh:              # sources stay relative, rebased
                text = fh.read()
            text = re.sub(r"(?m)^source: (.+) (sha256=\w+)$", lambda m: "source: " + (
                m.group(1) if os.path.isabs(m.group(1)) else rel(os.path.normpath(os.path.join(
                    os.path.dirname(s.path), m.group(1))), out)) + " " + m.group(2), text)
            with open(os.path.join(out, name), "w", encoding="utf-8") as fh:
                fh.write(text)
    man = {
        "show": s.show, "title": s.title,
        "hosts": [{"name": h, "role": r or None, "voice": backend.voices.get(h)} for h, r in s.hosts],
        "tts": {"backend": backend.name, "model": (backend.models_used or [backend.model])[0],
                "models_used": backend.models_used, "sample_rate": RATE,
                "sent_to_cloud": bool(backend.cloud),
                # the SDK's per-request timeout and retries as set (gemini_client), the two-line
                # test before the first chunk, and any hand-over from a model that stopped
                # answering (every chunk made again by the next one)
                "http": getattr(backend, "http", None), "preflight": preflight,
                "switched": list(getattr(backend, "switches", None) or []),
                "what_was_sent": ("the transcript turns after pronunciation substitutions, plus "
                                  "one style line; then, for verify (unless --no-verify), the "
                                  "rendered audio downsampled to 16 kHz; never the report"
                                  if backend.cloud else "nothing (offline)")},
        "cloud_tts_consent": (consent if backend.cloud else False),
        # "prepaid" once a 402 showed prepay billing, "free" once a 429 named a free-tier quota
        # (this render, or this key's earlier ones); null when no answer said which tier
        "key_tier": getattr(backend, "key_tier", None),
        "key_tier_basis": getattr(backend, "key_tier_basis", None) or (
            "no answer from Google named the key's tier in this render" if backend.cloud else None),
        "key_tier_evidence": (backend.key_tier_evidence()
                              if hasattr(backend, "key_tier_evidence") else None),
        "core_run": record_core_run(a.script, core_rec, "render") if core_rec else None,
        "script": "podcast_script.md", "script_sha256": s.sha256,
        "script_author": s.author,
        "check": {"status": status or "not run", "file": "check.txt",
                  "overridden": bool(why and a.unchecked), "reason": why},
        # the separate reader's sentence-by-sentence check (meaning_check.md, this script's sha)
        "meaning_check": {"status": m_status, "file": MEANING_FILE,
                          "overridden": bool(m_why and m_override), "reason": m_why,
                          "override_reason": m_override if m_why else None},
        "sources": sources,
        "words": words, "turns": len(s.turns()), "segments": len(s.segments),
        "chunks": len(chunks), "chunks_cached": n_cached, "chunks_synthesized": n_new,
        "duration_s": duration, "created": now_iso(), "ai_generated": True,
        "audio": audio_name, "encoder": tool, "transcript": "transcript.html",
        "warnings": warnings,
        "made_by": "make_podcast.py (UC Davis Proteomics Core pipeline skill)",
    }
    for name in VERIFY_FILES:                     # the last render's: never this audio's check
        if os.path.exists(os.path.join(out, name)):
            os.remove(os.path.join(out, name))
    with open(os.path.join(out, "podcast.json"), "w", encoding="utf-8") as fh:
        json.dump(man, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    with open(os.path.join(out, "transcript.html"), "w", encoding="utf-8") as fh:
        fh.write(transcript_html(s, man, out))
    print(json.dumps({"audio": os.path.join(out, audio_name), "duration_s": duration,
                      "minutes": round(duration / 60.0, 1), "backend": backend.name,
                      "model": man["tts"]["model"], "chunks": len(chunks),
                      "chunks_cached": n_cached, "chunks_synthesized": n_new,
                      "transcript": os.path.join(out, "transcript.html"),
                      "manifest": os.path.join(out, "podcast.json"),
                      "cache_pruned": prune_cache(cache, used, backend.model,
                                                  explicit=bool(getattr(a, "model", None))),
                      "warnings": warnings}, indent=2))
    # The ASR round trip. It sends the audio, so only under this render's cloud consent. It
    # never fails the render and never touches the report.
    retry = (f"`make_podcast.py verify {script_copy} "
             f"{'--core-run' if core_rec else '--cloud-ok ...'}`")
    if backend.cloud and consent and not getattr(a, "no_verify", False):
        try:
            run_verify(s.path if os.path.dirname(s.path) == out else script_copy, out,
                       getattr(a, "verify_model", None), consent, core=bool(core_rec))
        except Exception as e:
            why = scrub(e)
            log(f"[WARN] verify skipped: {why} (the render is fine; run {retry} to retry)")
            verify_not_run(out, f"verify failed: {why}", retry)
    elif backend.cloud:
        verify_not_run(out, "--no-verify was given", retry)
    else:
        verify_not_run(out, f"an offline render ({backend.name}) sends nothing",
                       "`make_podcast.py verify` with consent (it sends the audio to Google)")
    return 0


def verify_not_run(out, why, how):
    """This render's audio has not been heard back: verify.txt says so and why, and podcast.json
    gets a `verify` block with status NOT RUN -- never the last render's OK left in place
    (record_run copies verify.txt). Never fails the render."""
    try:
        with open(os.path.join(out, "verify.txt"), "w", encoding="utf-8") as fh:
            fh.write(f"Podcast audio check (ASR round trip): NOT RUN\nreason: {why}\n"
                     f"checked: {now_iso()}\nNothing has heard this audio back. Listen to it "
                     f"before it is shared, or run {how}.\n")
        mpath = os.path.join(out, "podcast.json")
        with open(mpath, encoding="utf-8") as fh:
            man = json.load(fh)
        man["verify"] = {"status": "NOT RUN", "reason": why, "created": now_iso(),
                         "report": "verify.txt"}
        with open(mpath + ".part", "w", encoding="utf-8") as fh:
            json.dump(man, fh, indent=2, ensure_ascii=False)
            fh.write("\n")
        os.replace(mpath + ".part", mpath)
    except (OSError, ValueError, AttributeError, TypeError) as e:
        log(f"[WARN] could not record that verify did not run ({type(e).__name__}: {e})")


def prune_cache(cache, used, model, explicit=False):
    """After a successful render: remove cached chunks this script no longer uses (an edited
    line's old chunk) and stray *.part files. A chunk made by another model -- or one with no
    record of its model -- is kept unless this render ran on an explicit --model: those chunks
    may still be the rest of an episode that is being resumed. -> number removed."""
    n, files, made_by = 0, sorted(os.listdir(cache)), {}
    for f in files:                                   # every chunk's model, before any delete
        if f.endswith(".json"):
            try:
                with open(os.path.join(cache, f), encoding="utf-8") as fh:
                    made_by[f.split(".", 1)[0]] = json.load(fh).get("model")
            except (OSError, ValueError, AttributeError):
                pass
    for f in files:
        stem = f.split(".", 1)[0]
        if f in (STALL_FILE, TIER_FILE):
            continue                             # the memories are not a chunk's sidecar
        if f.endswith(".part"):
            pass
        elif stem in used:
            continue
        elif not explicit and made_by.get(stem) != model:
            continue
        try:
            os.remove(os.path.join(cache, f))
            n += 1
        except OSError as e:
            log(f"[WARN] could not remove {f} from the cache: {e}")
    if n:
        log(f"[render] removed {n} cached file(s) this script no longer uses")
    return n


# ----------------------------------------------------------------------------- verify
VERIFY_PROMPT = "Transcribe verbatim. Write numbers as digits."
VERIFY_TEXT_MODEL = "gemini-2.5-flash"            # then the newest 3.x flash text model
VERIFY_RATIO = 0.93                               # WARN below this word match ratio
VERIFY_GAP = 6                                    # a script span this long not heard is a gap
VERIFY_INLINE_MAX = 18 * 1024 * 1024              # Gemini takes ~20 MB of inline data
_VTOK = re.compile(r"\d[\d,]*(?:\.\d+)?|[A-Za-z]+(?:'[A-Za-z]+)?")
_UNITS = {w: str(i) for i, w in enumerate(
    "zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen "
    "fifteen sixteen seventeen eighteen nineteen".split())}
_TENS = {w: str(10 * i) for i, w in enumerate(
    "_ _ twenty thirty forty fifty sixty seventy eighty ninety".split()) if w != "_"}


class VerifyError(Exception):
    pass


# Letter/digit boundaries, except inside 5e-15 and a 3k / 2M suffix.
_LETTER_DIGIT = re.compile(r"(?<=[A-Za-z])(?=\d)|(?<=\d)(?=[A-Za-z])(?![eE][+\-−]?\d)"
                           r"(?![kKMB](?![A-Za-z]))")


def split_letters_digits(text):
    """"LRS124" -> "LRS 124", "KCNB1" -> "KCNB 1"; 5e-15 and 3k stay whole."""
    return _LETTER_DIGIT.sub(" ", text or "")


def verify_tokens(text):
    """-> [(word, start, end)]: text as comparable words, split at letter/digit boundaries
    ("LRS124" and "LRS-124" are both lrs 124, "KCNB1" is kcnb 1). Lowercase; numbers without thousands
    separators; number words as digits ("forty six" -> 46, "two point one" -> 2.1); "minus"
    and "plus" dropped (the transcript writes -1.3); runs of capital single letters with only
    spaces between them joined ("K V" -> kv, "J P H 3" -> jph 3), so a spelled-out symbol and
    the transcript's IgG meet."""
    words = []
    for m in _VTOK.finditer(text or ""):
        w = m.group(0)
        words.append([w.replace(",", "") if w[0].isdigit() else w.lower(), m.start(), m.end()])
    out = _spoken_numbers(words)
    merged = []
    for t in out:
        if (len(merged) >= 2 and merged[-1][0] == "point" and merged[-2][0].isdigit()
                and t[0].isdigit()):
            merged[-2:] = [[merged[-2][0] + "." + t[0], merged[-2][1], t[2]]]
        else:
            merged.append(t)
    merged = [t for t in merged if t[0] not in ("minus", "plus")]
    joined, run = [], []

    def flush():
        if len(run) >= 2:
            joined.append(["".join(r[0] for r in run), run[0][1], run[-1][2]])
        else:
            joined.extend(run)
        run[:] = []
    for t in merged:
        letter = len(t[0]) == 1 and t[0].isalpha() and text[t[1]].isupper()
        if letter and (not run or not text[run[-1][2]:t[1]].strip()):
            run.append(t)
            continue
        flush()
        if letter:
            run.append(t)
        else:
            joined.append(t)
    flush()
    return [tuple(t) for t in joined]


_MULT = {"hundred": 100, "thousand": 1000, "million": 1000000}


def _spoken_numbers(words):
    """Runs of number words as one number, the way a transcript that "writes numbers as digits"
    has them. Read the way people say numbers: "five thousand and twenty four" -> 5024, "forty
    six" -> 46, "two fifty" -> 250, "nineteen ninety" -> 1990, and digit by digit, "zero seven
    five six" -> 0756. A tens word takes a following unit; any other pair of number words
    without a multiplier between them starts a new group of digits."""
    out, i = [], 0
    isnum = lambda w: w in _UNITS or w in _TENS or w in _MULT        # noqa: E731
    while i < len(words):
        if not isnum(words[i][0]):
            out.append(words[i])
            i += 1
            continue
        j, run = i, []
        while j < len(words) and (isnum(words[j][0]) or (      # "hundred and five"
                words[j][0] == "and" and run and run[-1] in _MULT and j + 1 < len(words)
                and isnum(words[j + 1][0]))):
            if words[j][0] != "and":
                run.append(words[j][0])
            j += 1
        groups, total, cur, last = [], 0, None, None
        for w in run:
            if w in _MULT:
                cur = (cur or 1) * _MULT[w]
                if _MULT[w] >= 1000:
                    total, cur = total + cur, 0
                last = "mult"
                continue
            v = int(_TENS.get(w) or _UNITS[w])
            if last == "tens" and v < 10:
                cur, last = cur + v, "unit"                   # forty six
                continue
            if last in ("unit", "teen", "tens"):              # two | fifty, seven | five
                groups.append(str(total + (cur or 0)))
                total, cur = 0, None
            cur = (cur or 0) + v
            last = "tens" if w in _TENS else ("teen" if v >= 10 else "unit")
        groups.append(str(total + (cur or 0)))
        out.append(["".join(groups), words[i][1], words[j - 1][2]])
        i = j
    return out


def compare_audio_text(s, speak, transcript):
    """The spoken script (after pronunciation) against what the ASR heard.
    -> dict: ratio, gaps (script spans of >= VERIFY_GAP words not heard), numbers_not_heard
    (each with +-60 characters of the transcript where it should have been), segments."""
    stoks, smeta = [], []
    for si, seg in enumerate(s.segments):
        for t in seg:
            for w, _, _ in verify_tokens(speak(t.text)):
                stoks.append(w)
                smeta.append((si + 1, t.line))
    tt = verify_tokens(transcript)
    ttoks = [w for w, _, _ in tt]
    sm = difflib.SequenceMatcher(None, stoks, ttoks, autojunk=False)
    ops = sm.get_opcodes()
    matched = sum(i2 - i1 for tag, i1, i2, _, _ in ops if tag == "equal")

    def heard_at(i):
        """The transcript character offset aligned with script word i."""
        for tag, i1, i2, j1, j2 in ops:
            if i1 <= i < i2 or (i1 == i2 == i):
                j = j1 + (i - i1) if tag == "equal" else j1
                j = min(j, len(tt) - 1)
                return tt[j][1] if tt else 0
        return len(transcript)

    def around(off):
        a, b = max(0, off - 60), min(len(transcript), off + 60)
        return ("…" if a else "") + transcript[a:b].replace("\n", " ") + ("…" if b < len(transcript) else "")

    gaps = []
    for tag, i1, i2, j1, j2 in ops:
        if tag in ("delete", "replace") and i2 - i1 >= VERIFY_GAP:
            gaps.append({"segment": smeta[i1][0], "line": smeta[i1][1], "words": i2 - i1,
                         "script": " ".join(stoks[i1:i2])[:300],
                         "heard": " ".join(ttoks[j1:j2])[:300]})

    # An ASR writes "LRS124" where the script says "LRS-124": a number glued to letters is
    # still heard. Both texts' words were split at letter/digit boundaries by verify_tokens;
    # the heard numbers are read from the transcript split the same way, and unsplit.
    heard = ({n.value for n in numbers_in(transcript)} |
             {n.value for n in numbers_in(split_letters_digits(transcript))})
    missing, seen = [], set()
    first_tok = {}
    for k, meta in enumerate(smeta):
        first_tok.setdefault(meta[1], k)
    for si, seg in enumerate(s.segments):
        for t in seg:
            for n in numbers_in(t.text):
                if (n.kind == "plain" and n.dec == 0 and n.value <= 10 and n.sign is None
                        and not n.quant) or n.value in heard or (t.line, n.value) in seen:
                    continue
                seen.add((t.line, n.value))
                start = first_tok.get(t.line, 0)
                digits = str(n.value)
                k = next((k for k in range(start, len(stoks)) if smeta[k][1] == t.line and
                          stoks[k] == digits), start)
                missing.append({"number": n.text, "segment": si + 1, "line": t.line,
                                "context": around(heard_at(k))})
    segs = sorted({g["segment"] for g in gaps} | {m["segment"] for m in missing})
    ratio = round(sm.ratio(), 3)
    return {"ratio": ratio, "coverage": round(matched / float(len(stoks)), 3) if stoks else 0.0,
            "words_script": len(stoks), "words_heard": len(ttoks), "gaps": gaps,
            "numbers_not_heard": missing, "segments_to_check": segs,
            "status": "WARN" if (ratio < VERIFY_RATIO or gaps or missing) else "OK"}


def _verify_audio(audio, tmp):
    """The episode as Gemini will take it inline: 16 kHz mono, then AAC at 32 kbps (ADTS,
    audio/aac) through afconvert or ffmpeg -- via a 16 kHz WAV, which afconvert needs. Without
    an encoder, the 16 kHz WAV itself when it fits. -> (bytes, mime, description)."""
    wav16 = os.path.join(tmp, "verify16k.wav")
    if audio.lower().endswith(".wav"):
        write_wav(wav16, resample(read_wav(audio), RATE, 16000), rate=16000)
    else:
        for cmd in (["afconvert", "-f", "WAVE", "-d", "LEI16@16000", "-c", "1", audio, wav16],
                    ["ffmpeg", "-y", "-loglevel", "error", "-i", audio, "-ar", "16000", "-ac",
                     "1", wav16]):
            if shutil.which(cmd[0]) and subprocess.run(cmd, capture_output=True,
                                                       timeout=900).returncode == 0:
                break
        else:
            raise VerifyError(f"cannot decode {os.path.basename(audio)}: neither afconvert nor "
                              "ffmpeg worked")
    aac = os.path.join(tmp, "verify.aac")
    for cmd in (["afconvert", "-f", "adts", "-d", "aac", "-b", "32000", wav16, aac],
                ["ffmpeg", "-y", "-loglevel", "error", "-i", wav16, "-c:a", "aac", "-b:a", "32k",
                 "-f", "adts", aac]):
        if (shutil.which(cmd[0]) and subprocess.run(cmd, capture_output=True,
                                                    timeout=900).returncode == 0
                and os.path.isfile(aac) and os.path.getsize(aac) > 1000):
            with open(aac, "rb") as fh:
                return fh.read(), "audio/aac", "16 kHz mono AAC 32 kbps"
    if os.path.getsize(wav16) <= VERIFY_INLINE_MAX:
        with open(wav16, "rb") as fh:
            return fh.read(), "audio/wav", "16 kHz mono WAV (no AAC encoder)"
    raise VerifyError("no AAC encoder, and the 16 kHz WAV is too large to send inline")


class Transcriber(object):
    """A Gemini TEXT model transcribes the audio: gemini-2.5-flash, then the newest 3.x flash
    text model models.list offers. Same key, scrubbing and rate-limit handling as rendering."""

    def __init__(self, models=None, core=False):
        try:
            from google import genai
            from google.genai import types
        except ImportError:
            raise VerifyError("the verify step needs the google-genai package: python3 -m pip "
                              "install google-genai")
        self.types = types
        self.client = genai.Client(api_key=read_key(core))
        self.candidates = list(models or []) or self.discover()

    def discover(self):
        names = []
        try:
            for m in self.client.models.list():
                n = (getattr(m, "name", "") or "").split("/")[-1]
                acts = getattr(m, "supported_actions", None) or []
                if (re.match(r"^gemini-3(?:\.\d+)?-flash$", n) and
                        (not acts or "generateContent" in acts)):
                    names.append(n)
        except Exception as e:
            log(f"[WARN] could not list the Gemini models ({scrub(e)[:160]})")
        names.sort(key=lambda n: float(re.search(r"gemini-(\d+(?:\.\d+)?)", n).group(1)),
                   reverse=True)
        return [VERIFY_TEXT_MODEL] + names[:1]

    def transcribe(self, data, mime):
        last = None
        for model in self.candidates:
            waits = errors = 0
            while True:
                try:
                    r = self.client.models.generate_content(
                        model=model, contents=[self.types.Part.from_bytes(data=data, mime_type=mime),
                                               VERIFY_PROMPT],
                        config=self.types.GenerateContentConfig(temperature=0))
                    text = getattr(r, "text", None)
                    if not text:
                        raise VerifyError("the response held no text")
                    return text, model
                except Exception as e:
                    err = error_info(e)
                    last = f"{model}: {why_line(err)}"
                    if err["kind"] == "rate" and waits < 5:
                        waits += 1
                        wait = (err["delay"] if err["delay"] is not None else 60.0) + 5.0
                        log(f"[verify] rate-limited, waiting {wait:.0f}s ({model})")
                        _sleep(wait)
                        continue
                    if err["kind"] == "transient" and errors < 2:
                        errors += 1
                        _sleep(5.0 * errors)
                        continue
                    if err["kind"] in ("fatal", "billing"):
                        raise VerifyError(last)
                    log(f"[verify] {last}; trying the next model")
                    break
        raise VerifyError(f"no model could transcribe the audio (last: {last})")


def verify_report(res, audio, how, size, model):
    L = [f"Podcast audio check (ASR round trip): {res['status']}",
         f"audio: {audio} (sent to Google as {how}, {size / 1e6:.1f} MB)",
         f"transcribed by: {model} (\"{VERIFY_PROMPT}\")",
         f"checked: {now_iso()}",
         f"word match ratio: {res['ratio']} (difflib, normalised words; WARN below {VERIFY_RATIO})"
         f" · script words: {res['words_script']:,} · heard: {res['words_heard']:,} · "
         f"script words matched: {res['coverage']:.1%}"]
    if res["gaps"]:
        L += ["", f"GAPS: script spans of {VERIFY_GAP}+ words not heard (dropped or garbled "
                  "audio, or an ASR slip)"]
        L += [f"- segment {g['segment']}, line {g['line']} ({g['words']} words): "
              f"\"{g['script']}\" -- heard: \"{g['heard'] or '(nothing)'}\"" for g in res["gaps"]]
    if res["numbers_not_heard"]:
        L += ["", "NUMBERS NOT HEARD (a count of 10 or less is not listed) -- the transcript "
                  "where each should be, to tell a TTS misread from an ASR mishearing"]
        L += [f"- {m['number']} (segment {m['segment']}, line {m['line']}): \"{m['context']}\""
              for m in res["numbers_not_heard"]]
    if res["segments_to_check"]:
        segs = " ".join(str(x) for x in res["segments_to_check"])
        L += ["", f"SEGMENTS TO LISTEN TO: {segs}",
              "A number the voice misread: add a Pronunciation row for it (e.g. | 5,024 | five "
              "thousand and twenty-four |), re-run check, then render -- only that chunk is "
              "re-made. Dropped or garbled audio with the text unchanged: render --redo "
              f"{segs}. Then verify again."]
    L += ["", "The ASR can mishear too: a flagged line is a place to listen, not proof of an "
              "audio fault. A clean result is not a listen either."]
    return "\n".join(L) + "\n"


def run_verify(script_path, out=None, models=None, consent=True, core=False):
    """Transcribe the rendered audio and compare it with the spoken script. Writes
    verify.txt, verify_transcript.txt and a `verify` block in podcast.json. -> the result dict.
    `core`: a Core run, which may use the Core's shared key. Raises VerifyError when it cannot
    run; never touches the report."""
    s = parse_script(script_path)
    out = os.path.abspath(out or os.path.dirname(s.path))
    mpath = os.path.join(out, "podcast.json")
    try:
        with open(mpath, encoding="utf-8") as fh:
            man = json.load(fh)
        audio = man["audio"]
        if not (isinstance(audio, str) and re.match(_FILE_RX["audio"], audio)):
            raise ValueError(f"'audio' is {audio!r}")
    except (OSError, ValueError, KeyError, TypeError) as e:
        raise VerifyError(f"no usable podcast.json in {out} ({type(e).__name__}: {e}): render "
                          "first")
    apath = os.path.join(out, audio)
    if not os.path.isfile(apath):
        raise VerifyError(f"{audio} is not in {out}")
    with tempfile.TemporaryDirectory() as td:
        data, mime, how = _verify_audio(apath, td)
    if len(data) > VERIFY_INLINE_MAX:
        raise VerifyError(f"the downsampled audio is {len(data) / 1e6:.0f} MB, over the "
                          "inline limit")
    log(f"[verify] sending {len(data) / 1e6:.1f} MB ({how}) for transcription")
    text, model = Transcriber(models, core).transcribe(data, mime)
    res = compare_audio_text(s, pronouncer(s.pronunciation), text)
    res.update(model=model, audio_sent=how, created=now_iso(), cloud_consent=consent,
               report="verify.txt", transcript="verify_transcript.txt")
    with open(os.path.join(out, "verify_transcript.txt"), "w", encoding="utf-8") as fh:
        fh.write(text.rstrip() + "\n")
    with open(os.path.join(out, "verify.txt"), "w", encoding="utf-8") as fh:
        fh.write(verify_report(res, audio, how, len(data), model))
    man["verify"] = res
    tmp = mpath + ".part"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(man, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    os.replace(tmp, mpath)
    if res["status"] == "WARN":
        log(f"[WARN] verify: word match ratio {res['ratio']}, {len(res['gaps'])} gap(s), "
            f"{len(res['numbers_not_heard'])} number(s) not heard -- listen to segment(s) "
            f"{', '.join(str(x) for x in res['segments_to_check']) or '(none named)'}; see "
            f"{os.path.join(out, 'verify.txt')}")
    else:
        log(f"[verify] word match ratio {res['ratio']}, no gaps, every number heard")
    return res


def cmd_verify(a):
    consent, core_rec, refusal = consent_of(a)
    if not consent:
        if refusal:
            log(f"[verify] not sending anything: {refusal}")
        if not getattr(a, "core_run", None):
            log(f"[verify] not sending anything. verify sends the rendered AUDIO (downsampled) "
                f"to Google's Gemini API for transcription. On a free-tier key Google may use it "
                f"to improve its products and human reviewers may hear it ({TERMS_URL}). Once "
                f"the user has agreed, re-run with --cloud-ok \"who agreed, when\" (a Core run: "
                f"--core-run).")
        return 2
    try:
        res = run_verify(a.script, a.out, a.model, consent, core=bool(core_rec))
    except (VerifyError, RenderError) as e:
        log(f"[verify] could not verify: {scrub(e)}")
        return 1
    if core_rec:
        record_core_run(a.script, core_rec, "verify")
    print(json.dumps({k: res[k] for k in ("status", "ratio", "coverage", "words_script",
                                          "words_heard", "segments_to_check", "model")},
                     indent=2))
    return 0


# ----------------------------------------------------------------------------- pages
TRANSCRIPT_CSS = """
.pc-turn{margin:.55rem 0;max-width:var(--measure,72ch)}
.pc-turn b{display:inline-block;min-width:3.6rem;font-size:.78rem;text-transform:uppercase;letter-spacing:.06em;color:var(--accent,#2463c9)}
.pc-turn.h2 b{color:var(--warning,#a8590a)}
.pc-seg{border-top:1px solid var(--line,#dfe3e8);padding-top:.4rem;margin-top:.9rem}
.pc-seg:first-of-type{border-top:0}
main.doc audio{display:block;width:100%;max-width:40rem;margin:.8rem 0 1rem}
@media print{main.doc audio{display:none}}
"""


def _inline(t):
    t = esc(t)
    return re.sub(r"(?<![*\w])\*([^*]+)\*(?!\w)", r"<em>\1</em>",
                  re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", t))


def transcript_html(s, man, out, embedded=False):
    """transcript.html for the podcast folder `out`. `embedded`: the copy built into the
    shareable report (share) -- no player (the card's is right above it) and nothing that
    points out of the file."""
    sys.path.insert(0, HERE)
    import report_style as rs
    title = f"{s.show} — {s.title or 'Untitled'}"
    report = os.path.join(os.path.dirname(out), "Analysis_Report.html")
    what = ("this report" if embedded else
            '<a href="../Analysis_Report.html">the analysis report</a>' if os.path.isfile(report)
            else "the analysis report")
    where = " (in the podcast folder of the full results)" if embedded else ""
    second = s.hosts[1][0] if len(s.hosts) > 1 else None
    segs = []
    for seg in s.segments:
        segs.append('<div class="pc-seg">' + "".join(
            f'<p class="pc-turn{" h2" if t.speaker == second else ""}"><b>{esc(t.speaker)}</b> '
            f"{_inline(t.text)}</p>" for t in seg) + "</div>")
    claims = ("<ul>" + "".join(f"<li>{_inline(c)}</li>" for c in s.claims) + "</ul>" if s.claims
              else "<p>None listed.</p>")
    tts = man["tts"]
    rows = [["Hosts", esc(", ".join(f"{h['name']} ({h['role']})" if h.get("role") else h["name"]
                                    for h in man["hosts"]))],
            ["Voices", esc(", ".join(f"{h['name']}: {h['voice']}" for h in man["hosts"]))],
            ["Speech", esc(f"{tts['backend']} ({', '.join(tts['models_used']) or tts['model']})"
                           + ("; the transcript text was sent to this cloud service"
                              if tts["sent_to_cloud"] else "; made offline"))],
            ["Length", esc(f"{minutes_label(man['duration_s'])} ({man['words']:,} words, "
                           f"{man['turns']} turns)")],
            ["Script check", esc(f"{man['check']['status']}"
                                 + (" (overridden with --unchecked)" if man["check"]["overridden"]
                                    else "") + " — check.txt" + where)],
            ["Made", esc(man["created"])],
            ["Script sha256", f"<code>{esc(s.sha256[:16])}…</code>"]]
    body = (f"<header><h1>{esc(title)}</h1><p>An AI-generated audio discussion of {what}: "
            f"two synthetic hosts talk through the results.</p></header>\n"
            + rs.callout("info", f"<p>{esc(DISCLOSURE + extras_note(man, 'file' if embedded else 'folder'))} The hosts are fictional, their voices are "
                                 "synthetic, and the script was written by an AI from the report; "
                                 "anything said beyond the report is listed below.</p>",
                         title="AI-generated")
            + ("" if embedded else
               f'<audio controls preload="none" src="{href(man["audio"])}"></audio>')
            + f"<style>{TRANSCRIPT_CSS}</style>"
            + '<section><h2 id="transcript">Transcript</h2>' + "".join(segs) + "</section>"
            + '<section><h2 id="claims">Claims beyond the report</h2><p>Statements the script '
              "makes that do not come from the report (general knowledge, metaphors, "
              "speculation). Check these before relying on them.</p>" + claims + "</section>"
            + '<section><h2 id="made">How this was made</h2>' + rs.table(["", ""], rows)
            + "<p>Pronunciation substitutions change only what the voices say; the transcript "
              f"above keeps the written form. They are listed in check.txt{where}.</p></section>")
    if hasattr(rs, "document"):
        return rs.document(title, body)
    return rs.page(title, body)


_HEADPHONES = ('<svg aria-hidden="true" width="18" height="18" viewBox="0 0 24 24" fill="none" '
               'stroke="currentColor" stroke-width="2" stroke-linecap="round" '
               'stroke-linejoin="round"><path d="M3 18v-6a9 9 0 0 1 18 0v6"/><path d="M21 19a2 2 0 '
               '0 1-2 2h-1a2 2 0 0 1-2-2v-3a2 2 0 0 1 2-2h3zM3 19a2 2 0 0 0 2 2h1a2 2 0 0 0 2-2v-3a2 '
               '2 0 0 0-2-2H3z"/></svg>')

CARD_CSS = (
    ".pc-card{--pc:var(--accent,#2463c9);background:var(--surface,#fff);color:var(--fg,#17191c);"
    "border:1px solid var(--line,#dfe3e8);border-left:5px solid var(--pc);border-radius:12px;"
    "padding:.85rem 1.1rem;margin:0 0 1.4rem;max-width:calc(var(--measure,72ch) + 4rem);"
    "box-shadow:var(--shadow,none)}"
    ".pc-card .pc-h{display:flex;flex-wrap:wrap;align-items:center;gap:.3rem .55rem;"
    "margin:0 0 .45rem;font-weight:650;color:var(--fg,#17191c)}"
    ".pc-card .pc-h svg{color:var(--pc)}"
    ".pc-card .pc-k{font-size:.72rem;text-transform:uppercase;letter-spacing:.07em;"
    "border:1px solid var(--pc);color:var(--pc);border-radius:999px;padding:.02rem .5rem}"
    ".pc-card .pc-d{color:var(--muted,#5d6570);font-weight:400;font-size:.9rem}"
    ".pc-card audio{display:block;width:100%;margin:.2rem 0 .45rem}"
    ".pc-card p{margin:.2rem 0;font-size:.9rem;color:var(--muted,#5d6570);max-width:none}"
    ".pc-card .pc-print{display:none}"
    ".pc-card details{margin:.4rem 0 0}.pc-card summary{cursor:pointer;font-size:.9rem;"
    "color:var(--pc)}.pc-card iframe{display:block;width:100%;height:32rem;margin-top:.4rem;"
    "border:1px solid var(--line,#dfe3e8);border-radius:8px;background:var(--surface,#fff)}"
    "@media print{.pc-card audio{display:none}.pc-card .pc-print{display:block}"
    ".pc-card details{display:none}"
    ".pc-card{box-shadow:none;break-inside:avoid;-webkit-print-color-adjust:exact;"
    "print-color-adjust:exact}}")


_FILE_RX = {"audio": r"^[A-Za-z0-9_.\-]+\.(?:m4a|wav)$", "transcript": r"^[A-Za-z0-9_.\-]+\.html$",
            "script": r"^[A-Za-z0-9_.\-]+\.md$"}
_WARNED = set()


def manifest_problem(outdir):
    """-> (manifest, None) for a usable podcast/podcast.json; (None, None) when there is none;
    (None, why) when it exists but cannot be used -- unreadable, truncated, or a field of the
    wrong type (a Listen card built from it would crash or point anywhere)."""
    path = os.path.join(outdir, "podcast", "podcast.json")
    if not os.path.exists(path):
        return None, None
    try:
        with open(path, encoding="utf-8") as fh:
            man = json.load(fh)
    except (OSError, ValueError) as e:
        return None, f"{path}: not readable JSON ({type(e).__name__}: {e})"
    if not isinstance(man, dict):
        return None, f"{path}: not a JSON object"
    for k, rx in _FILE_RX.items():
        v = man.get(k)
        if k == "audio" and v is None:
            return None, f"{path}: no 'audio' file name"
        if v is not None and not (isinstance(v, str) and re.match(rx, v)):
            return None, f"{path}: '{k}' must be a plain file name in podcast/, got {v!r}"
    d = man.get("duration_s")
    if d is not None and (isinstance(d, bool) or not isinstance(d, (int, float)) or d < 0):
        return None, f"{path}: 'duration_s' must be a number of seconds, got {d!r}"
    for k in ("show", "title"):
        if man.get(k) is not None and not isinstance(man[k], str):
            return None, f"{path}: '{k}' must be text, got {man[k]!r}"
    if man.get("check") is not None and not isinstance(man["check"], dict):
        return None, f"{path}: 'check' must be an object"
    return man, None


def load_manifest(outdir):
    """The manifest, or None. An invalid one is reported once (a [WARN] naming the file and the
    reason) and left out: the optional podcast never stops the report of record."""
    man, why = manifest_problem(outdir)
    if why and why not in _WARNED:
        _WARNED.add(why)
        log(f"[WARN] podcast.json exists but is invalid, so the podcast is left out: {why}")
    return man


def has_podcast(outdir):
    """Does this run have a podcast -- a usable podcast/podcast.json under `outdir`? The one test
    (the Listen card, the survey line's "and the podcast", README, the delivery's wording): an
    invalid podcast.json is no podcast, as it gets no Listen card."""
    return bool(load_manifest(outdir))


def _paths(outdir, base, man):
    pdir = os.path.join(outdir, "podcast")
    return (rel(os.path.join(pdir, man["audio"]), base),
            rel(os.path.join(pdir, man.get("transcript") or "transcript.html"), base))


def _name(man):
    return f"{man.get('show') or SHOW} — {man.get('title') or 'Untitled'}"


@_never_raise("")
def listen_card_html(outdir, base_dir=None, share=None):
    """The "Listen" card (with its markers) for a page in `base_dir` (default: `outdir`, where
    Analysis_Report.html sits). '' when `outdir` has no podcast/podcast.json. `share` (from
    build_share: the audio as a data: URI, the transcript page, the provenance mark) makes the
    shareable report's card: the audio and the transcript built in, no link out of the file."""
    man = load_manifest(outdir)
    if not man:
        return ""
    head = (f"{START}<style>{CARD_CSS}</style>"
            f'<aside class="pc-card" aria-label="Audio discussion of these results">'
            f'<div class="pc-h">{_HEADPHONES}<span class="pc-k">Listen</span>'
            f"<span>{esc(_name(man))}</span>"
            f'<span class="pc-d">{esc(minutes_label(man.get("duration_s")))}</span></div>')
    if share:
        return (head + share["mark"]
                + f'<audio controls preload="none" src="{share["audio"]}">Your browser cannot '
                  "play the built-in audio.</audio>"
                + f"<p>{esc(DISCLOSURE + extras_note(man, 'file'))} The transcript, with the "
                  "claims that go beyond the report, is below. The audio and the transcript are "
                  "built into this file, so it can be sent on its own.</p>"
                + f'<p class="pc-print">The audio and its transcript are built into '
                  f"{esc(SHARE_NAME)}; open it in a browser to listen.</p>"
                + '<details><summary>Read the transcript</summary><iframe title="Podcast '
                  f'transcript" srcdoc="{esc(share["transcript"])}"></iframe></details>'
                + f"</aside>{END}")
    audio, tr = _paths(outdir, base_dir or outdir, man)
    return (head
            + f'<audio controls preload="none" src="{href(audio)}">Open '
            f'<a href="{href(audio)}">{esc(audio)}</a> to listen.</audio>'
            f"<p>{esc(DISCLOSURE + extras_note(man))} <a href=\"{href(tr)}\">Read the transcript</a> "
            f"(with the claims that go beyond the report). The audio is a separate file in the podcast folder "
            f"beside this report.</p>"
            f'<p class="pc-print">Audio: {esc(audio)} · transcript: {esc(tr)} (beside this '
            f"report)</p></aside>{END}")


@_never_raise("")
def listen_md(outdir, base_dir=None):
    man = load_manifest(outdir)
    if not man:
        return ""
    audio, tr = _paths(outdir, base_dir or outdir, man)
    return "\n".join([START, f"> **Listen:** *{man.get('show') or SHOW}* — "
                             f"{man.get('title') or 'Untitled'} ({minutes_label(man.get('duration_s'))}): "
                             f"[{audio}]({href(audio)}). {DISCLOSURE}{extras_note(man)} Transcript: [{tr}]({href(tr)}).",
                      END])


@_never_raise("")
def readme_item_md(outdir, base_dir):
    """The README "Start here" bullet, paths relative to `base_dir`. '' without a podcast.
    session_docs.py and `link` both use this, so the wording lives here only."""
    man = load_manifest(outdir)
    if not man:
        return ""
    audio, tr = _paths(outdir, base_dir, man)
    return (f"- [Audio discussion of these results]({href(audio)}) — AI-generated, two synthetic "
            f"voices, {minutes_label(man.get('duration_s'))}; the report is the record "
            f"([transcript]({href(tr)}))")


@_never_raise("")
def share_item_md(outdir, base_dir):
    """The README "Start here" bullet for the shareable report (share), after the podcast's.
    '' when it is not there."""
    path = os.path.join(outdir, SHARE_NAME)
    if not (load_manifest(outdir) and _share_current(outdir)):
        return ""
    return (f"- [{SHARE_NAME}]({href(rel(path, base_dir))}) — **to send the report with its "
            f"audio, send this one file**: the same report with the audio discussion and its "
            f"transcript built in ({_mb(os.path.getsize(path))})")


@_never_raise([])
def agents_md_lines(outdir, base_dir):
    """AGENTS.md's section on the podcast. [] without a podcast."""
    man = load_manifest(outdir)
    if not man:
        return []
    pdir = rel(os.path.join(outdir, "podcast"), base_dir)
    share = os.path.join(outdir, SHARE_NAME)
    return ["## Audio discussion (podcast): a derivative, not a record", "",
            f"`{pdir}/` holds an AI-generated audio discussion of this analysis "
            f"(*{man.get('show') or SHOW}*: two synthetic hosts, a script written by an AI from "
            "the report). It is a derivative of the report and NOT authoritative: never cite it or "
            "take a number from it -- use the report and the tables. Its transcript is "
            f"`{pdir}/transcript.html`; the script, with its \"Claims beyond the report\" ledger "
            f"(statements not in the report), is `{pdir}/podcast_script.md`; the automated "
            f"fidelity check is `{pdir}/check.txt`; how it was made (TTS service, model, voices, "
            f"consent) is `{pdir}/podcast.json`."] + ([
            "", f"`{rel(share, base_dir)}` is the same report with this audio and its transcript "
            "built in, one file for a person to send on; it is mostly the audio as base64. Read "
            f"`{rel(os.path.join(outdir, 'Analysis_Report.html'), base_dir)}` instead, and never "
            "treat the copy as a second source."] if _share_current(outdir) else [])


def _md_inline_html(t):
    t = esc(t)
    t = re.sub(r"\[([^\]]+)\]\(([^)\s]+)\)", r'<a href="\2">\1</a>', t)
    return re.sub(r"(?<![*\w])\*([^*]+)\*(?!\w)", r"<em>\1</em>", t)


def _upsert(text, block, pos, sep=""):
    """Replace an existing marker block, else insert `block` + `sep` at `pos` (a callable ->
    index, or None when there is no place). A re-run replaces exactly the block, so the file
    ends up the same however often it runs."""
    i = text.find(START)
    j = text.find(END, i + 1) if i >= 0 else -1
    if i >= 0 and j >= 0:
        return text[:i] + block + text[j + len(END):], "replaced"
    k = pos(text)
    if k is None:
        return text, None
    return text[:k] + block + sep + text[k:], "added"


def _html_slot(doc):
    """Near the top: an explicit <!-- podcast:slot -->; else the top of <main> (after a header
    that opens it); else after the first </header>; else after <body>."""
    m = re.search(r"<!--\s*podcast:slot\s*-->", doc)
    if m:
        return m.end()
    m = re.search(r"<main\b[^>]*>", doc, re.I)
    if m:
        h = re.match(r"\s*<header\b.*?</header>", doc[m.end():], re.I | re.S)
        return m.end() + (h.end() if h else 0)
    for rx in (r"</header>", r"<body\b[^>]*>"):
        m = re.search(rx, doc, re.I)
        if m:
            return m.end()
    return 0


@_never_raise(lambda text, *rest: text)
def add_listen_card(doc, outdir, share=None):
    """`doc` with the Listen card near the top; unchanged when `outdir` has no podcast. The
    hook make_analysis_html.py calls, so regenerating the report keeps the card. `share`: the
    shareable report's card in its place (build_share)."""
    card = listen_card_html(outdir, share=share)
    return _upsert(doc, card, _html_slot)[0] if card else doc


@_never_raise(lambda text, *rest: text)
def add_listen_md(md, outdir):
    """The Markdown twin of add_listen_card: `md` with the Listen line under its title;
    unchanged when `outdir` has no podcast. For a writer of Analysis_Report.md."""
    line = listen_md(outdir)
    return _upsert(md, line, _md_slot, "\n\n")[0] if line else md


def _md_slot(text):
    m = re.search(r"(?m)^# .*\n", text)
    if not m:
        return 0
    k = m.end()
    b = re.match(r"[ \t]*\n", text[k:])
    return k + (b.end() if b else 0)


def _mentions(text):
    """True when the podcast is already listed outside a marker block (session_docs.py)."""
    return bool(re.search(r"podcast/(?:transcript\.html|podcast\.json|podcast\.(?:m4a|wav))",
                          strip_block(text)))


def _readme_html_slot(doc):
    m = re.search(r"<li>(?:(?!</li>).)*Analysis_Report\.html(?:(?!</li>).)*</li>", doc, re.S)
    if m:
        return m.end()
    m = re.search(r"<ul>", doc)
    return m.end() if m else None


def _readme_md_slot(text):
    m = re.search(r"(?m)^- .*Analysis_Report\.html.*\n", text)
    if m:
        return m.end()
    m = re.search(r"(?m)^## Start here[ \t]*\n(?:[ \t]*\n)?", text)
    return m.end() if m else None


def _agents_slot(text):
    for h in ("## Do not", "## Where this lives on HIVE"):
        m = re.search(r"(?m)^" + re.escape(h), text)
        if m:
            return m.start()
    return len(text)


def vouch(out, man, unchecked=False):
    """The Listen card vouches for the episode, so link and share go on only while its check
    still holds. -> (script, refuse, warn): the parsed script (None when it is missing); why to
    refuse (None: go on); what to warn about when going on although the check does not hold
    (--unchecked now, or the episode was rendered with it)."""
    path = os.path.join(out, "podcast", man.get("script") or "podcast_script.md")
    s = parse_script(path) if os.path.isfile(path) else None
    stale = (read_check(s)[2] if s else f"no {os.path.basename(path)} beside podcast.json")
    overridden = bool((man.get("check") or {}).get("overridden"))
    if stale and not (unchecked or overridden):
        return s, (f"refusing: {stale}. Re-run `make_podcast.py check` (and `render`, if the "
                   "script changed) so the episode matches the report, or pass --unchecked."), None
    mc = man.get("meaning_check")
    m_warn = None
    if s and isinstance(mc, dict):                    # rendered since the meaning check (2.13.2)
        m_why = read_meaning(s)[1]
        if m_why and not (mc.get("overridden") or unchecked):
            return s, (f"refusing: {m_why}. Record the meaning check for this script "
                       "(`make_podcast.py meaning <script> --verdict ...`), and render again if "
                       "the script changed, or pass --unchecked."), None
        m_warn = m_why and f"meaning check: {m_why}"
    elif s and ("meaning", out) not in _WARNED:
        _WARNED.add(("meaning", out))
        log("[WARN] this episode was rendered before the meaning check existed (skill 2.13.2): "
            "nothing on record says a separate reader checked its sentences against the report")
    warn = stale and stale + (" (rendered with --unchecked)" if overridden else " (--unchecked)")
    if m_warn:
        warn = (warn + "; " if warn else "") + m_warn + (
            " (rendered with --no-meaning-check)" if mc.get("overridden") else " (--unchecked)")
    return s, None, warn


def cmd_link(a):
    out = os.path.abspath(a.outdir)
    man, why = manifest_problem(out)
    if why:
        log(f"[link] podcast.json exists but is invalid: {why}. Re-render, or fix the file.")
        return 2
    if not man:
        log(f"[link] no podcast/podcast.json under {out}: render first, and point link at the "
            "folder that holds Analysis_Report.html and podcast/")
        return 2
    refuse, warn = vouch(out, man, a.unchecked)[1:]
    if refuse:
        log(f"[link] {refuse}")
        return 2
    if warn:
        log(f"[WARN] linking an episode whose check does not hold: {warn}")
    done = []

    def edit(path, fn, what):
        if not os.path.isfile(path):
            done.append(("INFO", path, "not there; skipped"))
            return
        with open(path, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
        new, how = fn(text)
        if how is None:
            done.append(("INFO", path, f"no place for {what}; skipped"))
            return
        if new != text:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(new)
        done.append(("OK", path, f"{what} {how}"))

    edit(os.path.join(out, "Analysis_Report.html"),
         lambda t: _upsert(t, listen_card_html(out), _html_slot), "Listen card")
    for name in ("Analysis_Report.md", "AI_Analysis_Report.md"):
        edit(os.path.join(out, name), lambda t: _upsert(t, listen_md(out), _md_slot, "\n\n"),
             "Listen line")                   # add_listen_md, with the added/replaced verdict
    # README / AGENTS.md sit at the session root: the folder above output/
    roots = [os.path.dirname(out) if os.path.basename(out) == "output" else out]
    for root in roots:
        item = readme_item_md(out, root)

        def readme_html(t, item=item):
            if _mentions(t):
                return t, "already listed (session_docs.py)"
            block = f"{START}<li>{_md_inline_html(item[2:])}</li>{END}"
            return _upsert(t, block, _readme_html_slot)

        def readme_md(t, item=item):
            if _mentions(t):
                return t, "already listed (session_docs.py)"
            return _upsert(t, f"{START}\n{item}\n{END}", _readme_md_slot, "\n")

        def agents(t, root=root):
            if _mentions(t):
                return t, "already described (session_docs.py)"
            block = START + "\n" + "\n".join(agents_md_lines(out, root)) + "\n" + END
            return _upsert(t, block, _agents_slot, "\n\n")
        edit(os.path.join(root, "README.html"), readme_html, "file-list entry")
        edit(os.path.join(root, "README.md"), readme_md, "file-list entry")
        edit(os.path.join(root, "AGENTS.md"), agents, "podcast section")
    done += refresh_feedback(out)                     # before the PDF is reprinted
    done.append(refresh_pdf(os.path.join(out, "Analysis_Report.html")))
    for level, path, note in done:
        print(f"[{level}] {rel(path, os.path.dirname(out))}: {note}")
    return 0


def _session_of(out):
    """The session folder of an output folder (README and AGENTS.md sit there too)."""
    return os.path.dirname(out) if os.path.basename(out) == "output" else out


def refresh_feedback(out):
    """A Core run's report ends with the Core's "How did we do?" survey line, written when the
    report was made -- usually before the podcast, so it asks about the report alone. Reword it
    to ask about the podcast too (core_submission.refresh_feedback / feedback_line: the one
    wording) in Analysis_Report.html and its .md twin. Only a line that is there: a report
    outside the Core has none, and link never adds one. -> [(level, path, note)]."""
    try:
        import core_submission
        import submission_report
    except ImportError as e:
        return [("INFO", os.path.join(out, "Analysis_Report.html"),
                 f"feedback line not checked: {e}")]
    if not submission_report.core_run(_session_of(out)):   # the one test of a Core run
        return []                                           # (link: never add a line)
    done = []
    for name, fmt in (("Analysis_Report.html", "html"), ("Analysis_Report.md", "md")):
        path = os.path.join(out, name)
        if not os.path.isfile(path):
            continue
        with open(path, encoding="utf-8", errors="replace") as fh:
            text = fh.read()
        new, n = core_submission.refresh_feedback(text, fmt, podcast=has_podcast(out))
        if n:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(new)
            done.append(("OK", path, "the feedback survey line now asks about the podcast too"))
    return done


def refresh_pdf(html_path):
    """Analysis_Report.pdf is printed from the HTML (html_to_pdf.py, where the skill has it):
    once link has put the card in the HTML, an older PDF lacks it. Reprint it with
    html_to_pdf.print_report; the card's print style shows the audio's file name in place of
    the player. When it cannot reprint, print_report renames the old PDF
    Analysis_Report.stale.pdf (kept, so it never passes for current) and says SKIPPED.
    -> (level, path, note)"""
    pdf = os.path.splitext(html_path)[0] + ".pdf"
    if not os.path.isfile(pdf):
        return "INFO", pdf, "no PDF beside the report; nothing to reprint"
    if not os.path.isfile(html_path) or os.path.getmtime(pdf) >= os.path.getmtime(html_path):
        return "OK", pdf, "up to date with the HTML"
    how = "open Analysis_Report.html in a browser, Print, Save as PDF"
    try:
        import html_to_pdf
    except ImportError:
        return "INFO", pdf, f"older than the HTML, so it has no Listen card; to reprint it, {how}"
    status, note = html_to_pdf.print_report(html_path, pdf)
    if status == "OK":
        return "OK", pdf, f"reprinted with the Listen card ({note})"
    if status == "SKIPPED":
        return "SKIPPED", pdf, f"NOT reprinted: {note}"
    return "INFO", pdf, f"older than the HTML and NOT reprinted: {note}"


# ----------------------------------------------------------------------------- share
SHARE_KBPS = 48                         # AAC, 24 kHz mono: speech stays clear from 32 to 64 kbps
SHARE_KBPS_RANGE = (32, 64)             # afconvert takes at most 64 kbps for 24 kHz mono AAC
EMAIL_LIMIT_MB = 25.0                   # many email systems refuse a message larger than ~25 MB
SHARE_WARN_MB = 18.0                    # ... and an attachment travels as base64, ~4/3 larger
SHARE_MARK = re.compile(r"<!-- podcast:share ([^>]*?) -->")
_DATA_AUDIO = re.compile(r'<audio\b[^>]*?\bsrc="data:(audio/[\w.+-]+);base64,([A-Za-z0-9+/=]*)"')
_REF = re.compile(r"""\s(?:src|href|poster|srcset|action)\s*=\s*(?:"([^"]*)"|'([^']*)'|([^\s>"']+))""",
                  re.I)
_CSS_URL = re.compile(r"""url\(\s*['"]?([^'")]+)""")
_SAFE_REF = re.compile(r"^(?:#|data:|https?:|mailto:|javascript:|about:)", re.I)
_EMBED_SAFE = re.compile(r"^(?:#|data:|https?:)", re.I)     # inside an iframe srcdoc
_SRCDOC = re.compile(r"""\ssrcdoc\s*=\s*(?:"([^"]*)"|'([^']*)')""", re.I)
_POD_REF = re.compile(r"(?:^|/)podcast/|\.(?:m4a|wav|mp3|aac)(?:[?#]|$)", re.I)


class ShareError(Exception):
    """share could not make a file it can vouch for. `code`: 2 a refusal (nothing to share, or
    the check does not hold), 1 a failure."""

    def __init__(self, msg, code=1):
        super().__init__(msg)
        self.code = code


def _mb(n):
    return f"{n / 1e6:.1f} MB"


def _attached(n):
    """The size of an n-byte file as an email attachment: base64, 4 bytes for every 3."""
    return (n + 2) // 3 * 4


def _share_current(outdir):
    """share_state's verdict, False when it cannot be read: README and AGENTS.md name the file
    only while it is the current report and episode, never a stale copy."""
    try:
        return share_state(outdir)[0] is True
    except OSError:
        return False


def local_refs(doc, safe=_SAFE_REF):
    """What in an HTML page points at another file -- a relative path, file:, a CSS url() --
    and so breaks when the page is sent by itself. data:, #anchors, http(s) and mailto do not.
    The page itself only: an iframe's srcdoc is a page of its own (embedded_refs)."""
    doc = _SRCDOC.sub(' srcdoc=""', doc)
    out = []
    for m in _REF.finditer(doc):
        v = html.unescape(next(g for g in m.groups() if g is not None)).strip()
        if v and not safe.match(v):
            out.append(v)
    out += [m.group(1).strip() for m in _CSS_URL.finditer(doc)
            if not safe.match(m.group(1).strip())]
    return out


def embedded_refs(doc):
    """The same inside every iframe srcdoc of `doc`, nested ones too, where only #anchors,
    data: and http(s) may appear. A srcdoc resolves a relative src against the page that holds
    it, so transcript.html's own player (src="podcast.m4a") showed a second player reading
    "Error" in the hand-made prototype (Brett, 2026-09-28), and its ../Analysis_Report.html link
    went nowhere."""
    out = []
    for m in _SRCDOC.finditer(doc):
        inner = html.unescape(m.group(1) if m.group(1) is not None else m.group(2))
        out += local_refs(inner, _EMBED_SAFE) + embedded_refs(inner)
    return out


def _boxes(data, start=0, end=None):
    """(type, body start, end) of each ISO base-media box in data[start:end]."""
    end = len(data) if end is None else end
    i = start
    while i + 8 <= end:
        size, kind = struct.unpack(">I4s", data[i:i + 8])
        head = 8
        if size == 1:
            size, head = struct.unpack(">Q", data[i + 8:i + 16])[0], 16
        elif size == 0:
            size = end - i
        if size < head or i + size > end:
            raise ValueError(f"box {kind.decode('latin-1')!r} runs past the end of the file")
        yield kind, i + head, i + size
        i += size


def _box(data, span, kind):
    for k, a, b in _boxes(data, *span):
        if k == kind:
            return a, b
    raise ValueError(f"no {kind.decode()} box")


def mp4_info(data):
    """What a browser needs to play an .m4a, read from its boxes (ISO/IEC 14496-12): an ftyp
    first, a moov holding an audio track (handler 'soun') whose sample entry is mp4a (AAC), and
    the samples (mdat). -> {"brand", "codec", "rate", "duration_s", "mdat_bytes"}; ValueError
    for anything else. (The sample entry's channel count is a fixed 2 whatever the audio.)"""
    try:
        if data[4:8] != b"ftyp":
            raise ValueError("not an MP4 file (it does not start with an ftyp box)")
        top = {}
        for k, a, b in _boxes(data):
            top.setdefault(k, (a, b))
        if b"moov" not in top:
            raise ValueError("no moov box (the index a player reads first)")
        mdat = top.get(b"mdat")
        if not mdat or mdat[1] - mdat[0] < 1000:
            raise ValueError("no audio samples (the mdat box is missing or empty)")
        for k, a, b in _boxes(data, *top[b"moov"]):
            if k != b"trak":
                continue
            mdia = _box(data, (a, b), b"mdia")
            h = _box(data, mdia, b"hdlr")[0]
            if data[h + 8:h + 12] != b"soun":
                continue
            m = _box(data, mdia, b"mdhd")[0]
            scale, dur = (struct.unpack(">IQ", data[m + 20:m + 32]) if data[m] == 1 else
                          struct.unpack(">II", data[m + 12:m + 20]))
            stsd = _box(data, _box(data, _box(data, mdia, b"minf"), b"stbl"), b"stsd")
            entry = next(_boxes(data, stsd[0] + 8, stsd[1]), None)
            if not entry or not scale:
                raise ValueError("the audio track has no sample description or timescale")
            info = {"brand": data[top[b"ftyp"][0]:top[b"ftyp"][0] + 4].decode("latin-1").strip(),
                    "codec": entry[0].decode("latin-1"),
                    "rate": struct.unpack(">I", data[entry[1] + 24:entry[1] + 28])[0] >> 16,
                    "duration_s": dur / float(scale), "mdat_bytes": mdat[1] - mdat[0]}
            if info["codec"] != "mp4a":
                raise ValueError(f"the audio is {info['codec']!r}, not AAC (mp4a)")
            return info
        raise ValueError("no audio track in the moov box")
    except struct.error as e:
        raise ValueError(f"a truncated MP4 box ({e})")


def _tool(cmd, tried):
    """Run an audio tool when it is installed; True when it exited 0. A failure is added to
    `tried` (for the message that names what was attempted)."""
    if not shutil.which(cmd[0]):
        return False
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    except (OSError, subprocess.SubprocessError) as e:
        tried.append(f"{cmd[0]}: {e}")
        return False
    if r.returncode:
        tried.append(f"{cmd[0]}: exit {r.returncode} {(r.stderr or r.stdout).strip()[:160]}")
    return r.returncode == 0


def share_audio(audio, kbps, tmp):
    """The episode for the shareable report: AAC at `kbps`, 24 kHz mono, in an .m4a -- through a
    24 kHz WAV (afconvert encodes AAC from PCM only), with afconvert, else ffmpeg. With neither,
    podcast.m4a as rendered (larger). -> (bytes, how)."""
    wav, m4a, tried = os.path.join(tmp, "speech.wav"), os.path.join(tmp, "speech.m4a"), []
    if audio.lower().endswith(".wav"):
        wav = audio                                     # render writes 24 kHz mono 16-bit
    elif not any(_tool(cmd, tried) and os.path.isfile(wav) for cmd in (
            ["afconvert", "-f", "WAVE", "-d", f"LEI16@{RATE}", "-c", "1", audio, wav],
            ["ffmpeg", "-y", "-loglevel", "error", "-i", audio, "-ar", str(RATE), "-ac", "1",
             wav])):
        wav = None
    for cmd in ([] if wav is None else [
            ["afconvert", "-f", "m4af", "-d", "aac", "-b", str(kbps * 1000), wav, m4a],
            ["ffmpeg", "-y", "-loglevel", "error", "-i", wav, "-c:a", "aac", "-b:a", f"{kbps}k",
             "-movflags", "+faststart", m4a]]):
        if _tool(cmd, tried) and os.path.isfile(m4a) and os.path.getsize(m4a) > 1000:
            with open(m4a, "rb") as fh:
                return fh.read(), f"AAC {kbps} kbps, {RATE // 1000} kHz mono, by {cmd[0]}"
    why = "; ".join(tried) or "neither afconvert nor ffmpeg is installed"
    if audio.lower().endswith(".m4a"):
        with open(audio, "rb") as fh:
            return fh.read(), f"{os.path.basename(audio)} as rendered, not re-encoded ({why})"
    raise ShareError(f"cannot make AAC from {os.path.basename(audio)} ({why}): install ffmpeg")


def decode_back(data, tmp):
    """Decode `data` -- the audio taken back out of the written page -- with afconvert or ffmpeg,
    to 8 kHz mono. -> (seconds, peak sample), or None when neither tool is here."""
    src, back, tried = os.path.join(tmp, "embedded.m4a"), os.path.join(tmp, "decoded.wav"), []
    with open(src, "wb") as fh:
        fh.write(data)
    for cmd in (["afconvert", "-f", "WAVE", "-d", "LEI16@8000", "-c", "1", src, back],
                ["ffmpeg", "-y", "-loglevel", "error", "-i", src, "-ar", "8000", "-ac", "1", back]):
        if _tool(cmd, tried) and os.path.isfile(back):
            peak = 0
            with wave.open(back, "rb") as w:
                n, rate = w.getnframes(), w.getframerate()
                while True:
                    a = _arr(w.readframes(1 << 16))
                    if not a:
                        break
                    peak = max(peak, max(a), -min(a))
            return n / float(rate), peak
    if tried:
        raise ShareError("the built-in audio does not decode: " + "; ".join(tried))
    return None


def _close(a, b):
    return abs(a - b) <= max(1.0, 0.01 * b)


PLAY_JS = r"""(async () => { try {
  const a = document.querySelector('.pc-card audio');
  if (!a) return {error: 'no audio player in the Listen card'};
  if (!a.canPlayType('audio/mp4; codecs="mp4a.40.2"')) return {skip: 'this browser has no AAC decoder'};
  const wait = (ev, ms) => new Promise((ok, no) => {
    const t = setTimeout(() => no({timeout: 'no ' + ev + ' within ' + ms / 1000 + ' s'}), ms);
    a.addEventListener(ev, () => { clearTimeout(t); ok(); }, {once: true});
    a.addEventListener('error', () => { clearTimeout(t);
      no(new Error('media error ' + (a.error ? a.error.code + ' ' + a.error.message : ''))); }, {once: true});
  });
  a.preload = 'auto'; const meta = wait('loadedmetadata', 20000); a.load(); await meta;
  a.muted = true; const t0 = a.currentTime; let played = 0, playErr = null;
  try { await a.play(); await new Promise(r => setTimeout(r, 1500)); played = a.currentTime - t0; a.pause(); }
  catch (e) { playErr = String(e); }
  const f = document.querySelector('.pc-card iframe');
  const d = f && f.contentDocument;
  return {duration: a.duration, played: played, playErr: playErr,
          turns: d ? d.querySelectorAll('.pc-turn').length : null,
          players: d ? d.querySelectorAll('audio, video').length : null,
          fetched: performance.getEntriesByType('resource').map(e => e.name.slice(0, 120))};
} catch (e) { return e && e.timeout ? e : {error: String(e && e.message || e)}; } })()"""


def browser_play(path, timeout=60.0):
    """Open the page in a headless Chromium (html_to_pdf.find_browser, the browser the PDF is
    printed with) and play the built-in audio for 1.5 s, muted, over the DevTools pipe
    (--remote-debugging-pipe: JSON messages ending in NUL on fds 3 and 4). A --virtual-time-budget
    probe hangs on a data: URI audio (2026-09-28); this has a real deadline. -> (result, None),
    or (None, why the test could not run)."""
    if os.name != "posix":
        return None, "the browser test runs on macOS and Linux only"
    try:
        import html_to_pdf
    except ImportError:
        return None, "html_to_pdf.py, which finds the browser, is not beside this script"
    exe = html_to_pdf.find_browser()
    if not exe:
        return None, "no Chrome, Chromium or Edge found"
    profile = tempfile.mkdtemp(prefix="podcast-play-")
    to_b, from_me = os.pipe()                   # we write commands -> the browser's fd 3
    from_b, to_me = os.pipe()                   # the browser writes replies -> its fd 4

    def fds():
        a, b = os.dup(to_b), os.dup(to_me)      # dup first: a pipe may already sit on 3 or 4
        os.dup2(a, 3)
        os.dup2(b, 4)
        os.set_inheritable(3, True)             # dup2 onto itself keeps close-on-exec
        os.set_inheritable(4, True)
    try:
        proc = subprocess.Popen(
            [exe, "--headless=new", "--remote-debugging-pipe", f"--user-data-dir={profile}",
             "--no-first-run", "--no-default-browser-check", "--mute-audio",
             "--autoplay-policy=no-user-gesture-required", "about:blank"],
            preexec_fn=fds, close_fds=False, start_new_session=True,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except OSError as e:
        for fd in (to_b, from_me, from_b, to_me):
            os.close(fd)
        html_to_pdf._remove_profile(profile)
        return None, f"the browser did not start ({e})"
    os.close(to_b)
    os.close(to_me)
    deadline, buf, ids = time.monotonic() + timeout, b"", [0]

    def call(method, params=None, session=None):
        nonlocal buf
        ids[0] += 1
        msg = {"id": ids[0], "method": method, "params": params or {}}
        if session:
            msg["sessionId"] = session
        os.write(from_me, json.dumps(msg).encode() + b"\0")
        while True:
            while b"\0" in buf:
                raw, buf = buf.split(b"\0", 1)
                reply = json.loads(raw.decode("utf-8", "replace"))
                if reply.get("id") == ids[0]:
                    if "error" in reply:
                        raise RuntimeError(f"{method}: {reply['error'].get('message')}")
                    return reply.get("result") or {}
            left = deadline - time.monotonic()
            if left <= 0 or not select.select([from_b], [], [], left)[0]:
                raise TimeoutError(f"no answer from the browser within {timeout:.0f} s")
            chunk = os.read(from_b, 1 << 16)
            if not chunk:
                raise RuntimeError("the browser closed the DevTools pipe")
            buf += chunk
    try:
        target = call("Target.createTarget", {"url": "file://" + urllib.parse.quote(
            os.path.abspath(path))})["targetId"]
        session = call("Target.attachToTarget", {"targetId": target, "flatten": True})["sessionId"]
        while call("Runtime.evaluate", {"expression": "location.protocol + document.readyState",
                                        "returnByValue": True}, session
                   )["result"].get("value") != "file:complete":
            if time.monotonic() > deadline:
                raise TimeoutError(f"the page did not load within {timeout:.0f} s")
            time.sleep(0.2)
        r = call("Runtime.evaluate", {"expression": PLAY_JS, "awaitPromise": True,
                                      "returnByValue": True}, session)
        return r["result"].get("value") or {"error": "the page returned nothing"}, None
    except (OSError, ValueError, KeyError, RuntimeError, TimeoutError) as e:
        return None, f"the browser test did not finish: {e}"
    finally:
        for fd in (from_me, from_b):
            os.close(fd)
        html_to_pdf._stop(proc)
        html_to_pdf._remove_profile(profile)


def share_keys(out, man):
    """What the shareable report is built from: Analysis_Report.html, podcast.json and the audio
    (sha256, 16 hex). Its mark records them, so a copy on HIVE or in a delivery can be checked
    without trusting file times."""
    pdir = os.path.join(out, "podcast")
    return (f"report={sha256_file(os.path.join(out, 'Analysis_Report.html'))[:16]} "
            f"podcast={sha256_file(os.path.join(pdir, 'podcast.json'))[:16]} "
            f"audio={sha256_file(os.path.join(pdir, man['audio']))[:16]}")


def share_state(outdir):
    """Is OUTDIR/Analysis_Report_with_audio.html the current report and episode? -> (None, why)
    with no podcast/podcast.json (nothing to share); (False, why) when it is missing, built from
    another report, podcast.json or audio, or the podcast's check no longer holds; (True, what)."""
    out = os.path.abspath(outdir)
    man, why = manifest_problem(out)
    if not man:                             # the reason without the session's path: it is delivered
        return ((False, "podcast/podcast.json is not usable: " + why.split(": ", 1)[-1]) if why
                else (None, "no podcast was made"))
    path = os.path.join(out, SHARE_NAME)
    if not os.path.isfile(path):
        return False, "not built yet"
    with open(path, encoding="utf-8", errors="replace") as fh:
        m = SHARE_MARK.search(fh.read())
    try:
        now = share_keys(out, man)
    except OSError as e:
        return False, f"cannot read what it is built from ({e.strerror or type(e).__name__})"
    if not m or not m.group(1).startswith(now):
        return False, "built from another version of the report or the episode"
    if "unchecked=1" not in m.group(1):
        refuse = vouch(out, man)[1]
        if refuse:
            return False, "the podcast's check no longer holds"
    return True, "current: built from this report and episode"


def build_share(outdir, kbps=SHARE_KBPS, unchecked=False, play=True):
    """Write OUTDIR/Analysis_Report_with_audio.html: Analysis_Report.html with the Listen card
    (listen_card_html, share mode) holding the audio as a data: URI and the transcript as an
    iframe srcdoc. Kept only once the file as written holds up: the audio taken back out of it
    is the audio put in, an MP4 AAC track of the episode's length that afconvert/ffmpeg decode
    to that length, nothing in it points at podcast/ and nothing in the transcript's srcdoc
    points outside the file, and (play) a headless Chromium plays it.
    Written via a .part; intermediates live in a temporary folder outside the session.
    -> {"path", "bytes", "lines": [(level, text)], "note"}. ShareError otherwise."""
    out = os.path.abspath(outdir)
    lo, hi = SHARE_KBPS_RANGE
    if not lo <= kbps <= hi:
        raise ShareError(f"--kbps {kbps}: use {lo} to {hi} (speech, 24 kHz mono AAC)", 2)
    man, why = manifest_problem(out)
    if why:
        raise ShareError(f"podcast.json exists but is invalid: {why}. Re-render, or fix the file.", 2)
    if not man:
        raise ShareError(f"no podcast/podcast.json under {out}: render first, and point share at "
                         "the folder that holds Analysis_Report.html and podcast/", 2)
    report = os.path.join(out, "Analysis_Report.html")
    audio = os.path.join(out, "podcast", man["audio"])
    for f in (report, audio):
        if not os.path.isfile(f):
            raise ShareError(f"no {rel(f, os.path.dirname(out))}: nothing to build it from", 2)
    s, refuse, warn = vouch(out, man, unchecked)
    if refuse:
        raise ShareError(refuse, 2)
    if not s:
        raise ShareError("no podcast_script.md beside podcast.json: the transcript is built from it", 2)
    lines = [("WARN", f"the podcast's check does not hold: {warn}")] if warn else []
    mark = f"<!-- podcast:share {share_keys(out, man)} kbps={kbps}{' unchecked=1' if warn else ''} -->"
    path, want = os.path.join(out, SHARE_NAME), man.get("duration_s")
    part = path + ".part"
    try:
        with tempfile.TemporaryDirectory(prefix="podcast-share-") as tmp:
            data, how = share_audio(audio, kbps, tmp)
            try:
                info = mp4_info(data)
            except ValueError as e:
                raise ShareError(f"the audio to build in is not playable MP4/AAC ({how}): {e}")
            if want and not _close(info["duration_s"], want):
                raise ShareError(f"the audio to build in is {info['duration_s']:.1f} s long; "
                                 f"podcast.json says {want:.1f} s")
            page = transcript_html(s, man, os.path.dirname(audio), embedded=True)
            with open(report, encoding="utf-8", errors="replace") as fh:
                doc = fh.read()
            b64 = base64.b64encode(data).decode("ascii")
            doc = add_listen_card(doc, out, share={"audio": "data:audio/mp4;base64," + b64,
                                                   "transcript": page, "mark": mark})
            if mark not in doc:
                raise ShareError("the Listen card could not be built (the [WARN] above says why)")
            doc = _share_feedback(doc, out)            # the survey link says src=share
            with open(part, "w", encoding="utf-8") as fh:
                fh.write(doc)
            with open(part, encoding="utf-8") as fh:          # what a recipient gets
                doc = fh.read()
            inner = embedded_refs(doc)                        # the transcript: nothing outside
            if inner:
                raise ShareError("the built-in transcript points out of the file: "
                                 + ", ".join(inner)[:300])
            refs = local_refs(doc)
            if [r for r in refs if _POD_REF.search(r)]:
                raise ShareError("it still points at the podcast folder: "
                                 + ", ".join(r for r in refs if _POD_REF.search(r))[:300])
            got = _DATA_AUDIO.search(doc)
            if not got or base64.b64decode(got.group(2), validate=True) != data:
                raise ShareError("the audio read back out of the written file is not the audio "
                                 "put in")
            dec = decode_back(data, tmp)
            if dec and (not _close(dec[0], want or info["duration_s"]) or dec[1] < 300):
                raise ShareError(f"the built-in audio decodes to {dec[0]:.1f} s (peak {dec[1]}); "
                                 f"expected {want or info['duration_s']:.1f} s of sound")
            os.replace(part, path)
    finally:
        if os.path.exists(part):
            os.remove(part)
    size = os.path.getsize(path)
    lines.insert(0, ("OK", f"{SHARE_NAME}: {_mb(size)} ({_mb(_attached(size))} as an email "
                           f"attachment), the report with the audio ({how}, "
                           f"{minutes_label(want)}) and its transcript built in: send this one file"))
    lines.append(("OK", f"built-in audio: MP4 ({info['brand']}) AAC {info['rate']} Hz, "
                        f"{info['duration_s']:.1f} s in its header"
                  + (f"; afconvert/ffmpeg decode it back to {dec[0]:.1f} s" if dec else
                     " (no afconvert or ffmpeg here to decode it back)")
                  + (f" (podcast.json: {want:.1f} s)" if want else "")))
    other = [r for r in refs if not _POD_REF.search(r)]
    if other:
        lines.append(("WARN", f"{len(other)} link(s) in the report point at files beside it and "
                              f"will not open when this file is sent alone: {', '.join(other[:5])}"))
    if play:
        res, why = browser_play(path)
        if res is None or res.get("skip") or res.get("timeout"):
            lines.append(("INFO", "not played in a browser: "
                                  + (why or res.get("skip") or res.get("timeout"))))
        elif res.get("error") or not res.get("duration") or not _close(res["duration"], info["duration_s"]):
            os.remove(path)
            raise ShareError("a headless browser could not play the built-in audio: "
                             + str(res.get("error") or f"duration {res.get('duration')}"))
        elif res.get("players"):
            os.remove(path)
            raise ShareError(f"the built-in transcript has {res['players']} player(s) of its own, "
                             "which cannot play from inside the file")
        else:
            played = (res.get("played") or 0) > 0.5
            lines.append(("OK" if played else "WARN",
                          f"a headless browser loaded it ({res['duration']:.1f} s) and "
                          + (f"played {res['played']:.1f} s" if played
                             else f"did not play it ({res.get('playErr') or 'no progress'})")
                          + f"; the transcript shows {res.get('turns')} of {len(s.turns())} "
                            f"turns and no player; files fetched: {len(res.get('fetched') or [])}"))
            if res.get("fetched"):
                lines.append(("WARN", "it fetched: " + ", ".join(res["fetched"][:5])))
    note = (f"{_mb(size)}: the report with the audio ({how}) and transcript built in; "
            "send this one file")
    if size > SHARE_WARN_MB * 1e6:
        est = size - len(b64) + int(lo * 125 * (want or info["duration_s"]) * 4 / 3)
        lines.append(("WARN", f"{_mb(size)} is {_mb(_attached(size))} as an email attachment "
                              f"(base64, about 4/3 larger), more than the ~{EMAIL_LIMIT_MB:.0f} MB "
                              "many email systems accept: share it through Bioshare instead "
                              "(core_submission.py deliver puts it in the share)"
                      + (f", or re-run with --kbps {lo} (about {_mb(est)}, "
                         f"{_mb(_attached(est))} attached)" if kbps > lo else "")))
        note = (f"{_mb(size)} ({_mb(_attached(size))} attached), more than many email systems "
                "accept: share it through Bioshare. " + note)
    return {"path": path, "bytes": size, "lines": lines, "note": note}


def _share_feedback(doc, out):
    """The shareable report's survey line: the report's own, with src=share (the copy made to be
    forwarded), for a Core run only (submission_report.core_run) and only when the line is there."""
    try:
        import core_submission
        import submission_report
    except ImportError:
        return doc
    if not submission_report.core_run(_session_of(out)):
        return doc
    return core_submission.refresh_feedback(doc, "html", has_podcast(out), src="share")[0]


def ensure_share(outdir):
    """finalize's step: with a podcast, the shareable report current -- left alone when it is,
    built when it is missing or out of date. A copy that cannot be brought up to date is renamed
    Analysis_Report_with_audio.stale.html (scratch_files: never zipped or catalogued) so it is
    never sent; once a current one exists the stale one is deleted. -> None without a podcast,
    else (level, note) for MANIFEST.txt."""
    out = os.path.abspath(outdir)
    current, why = share_state(out)
    if current is None:
        return None
    stale = os.path.join(out, SHARE_STALE)
    if current:
        level, note = "OK", why
    else:
        try:
            level, note = "OK", build_share(out)["note"]
        except Exception as e:                      # recorded in MANIFEST.txt, never swallowed
            level, note = "SKIPPED", f"not built: {e if isinstance(e, ShareError) else repr(e)}"
            path = os.path.join(out, SHARE_NAME)
            if os.path.isfile(path):
                os.replace(path, stale)
                # first, so MANIFEST's 200-character cap on a reason never cuts it off
                note = (f"the older copy no longer matches and was renamed {SHARE_STALE}: never "
                        f"send it -- {note}")
    if level == "OK" and os.path.isfile(stale):
        os.remove(stale)
    return level, note


def cmd_share(a):
    try:
        res = build_share(a.outdir, kbps=a.kbps, unchecked=a.unchecked, play=not a.no_browser)
    except ShareError as e:
        log(f"[share] {e}")
        return e.code
    for level, text in res["lines"]:
        print(f"[{level}] {text}")
    return 0


# ----------------------------------------------------------------------------- the job (HIVE)
def _vtuple(v):
    return tuple(int(x) for x in re.findall(r"\d+", str(v or ""))[:3])


def pipeline_home():
    """The skill's per-user folder: $PROTEOMICS_PIPELINE_HOME (activate.sh), else
    ~/.proteomics-pipeline -- setup.json and the private jobs/ folder live there."""
    return os.environ.get("PROTEOMICS_PIPELINE_HOME") or os.path.join(
        os.path.expanduser("~"), ".proteomics-pipeline")


def private_dir(path):
    """Make `path` a folder only this account can use -> None, or why it cannot be. Core session
    folders on HIVE are group-writable with no sticky bit, so a job script there could be
    swapped for another member's before sbatch reads it (the job then runs as you), and a link
    planted at its log's name would be written through. So the job and its log live here."""
    try:
        os.makedirs(path, mode=0o700, exist_ok=True)
        st = os.lstat(path)
        parent = os.lstat(os.path.dirname(os.path.abspath(path)))
    except OSError as e:
        return f"{path}: {e}"
    uid = os.getuid() if hasattr(os, "getuid") else None
    if not stat.S_ISDIR(st.st_mode):
        return f"{path} is not a folder (a link?)"
    if uid is not None and st.st_uid != uid:
        return f"{path} belongs to another account"
    if uid is not None and (parent.st_uid != uid or parent.st_mode & 0o002 or (
            parent.st_mode & 0o020 and parent.st_gid != os.getgid())):
        return (f"{os.path.dirname(os.path.abspath(path))} can be changed by other accounts, so "
                f"{path} could be put back as theirs")
    if stat.S_IMODE(st.st_mode) & 0o077:
        try:
            os.chmod(path, 0o700)
        except OSError as e:
            return f"{path}: cannot make it private ({e})"
    return None


def job_python(given=None):
    """The python the podcast job runs -> (path, its google-genai version, None), or (None, None,
    why). --python, else the skill's env here (setup.json's `python`, where setup.sh installs
    google-genai), else this python: the first whose google-genai is GENAI_MIN or newer. The
    python that writes the job is often not it -- on HIVE `python3` is the system's."""
    cands = [given]
    if not given:
        try:
            with open(os.path.join(pipeline_home(), "setup.json"), encoding="utf-8") as fh:
                cands = [json.load(fh).get("python")]
        except (OSError, ValueError, AttributeError):
            cands = []
        cands.append(sys.executable)
    tried, seen = [], set()
    for py in cands:
        if not (isinstance(py, str) and py) or py in seen:
            continue
        seen.add(py)
        v = None
        try:
            r = subprocess.run([py, "-c", "import google.genai as g; print(g.__version__)"],
                               capture_output=True, text=True, timeout=120)
            if r.returncode == 0 and r.stdout.strip():
                v = r.stdout.strip().splitlines()[-1]
        except (OSError, subprocess.SubprocessError):
            pass
        if v and _vtuple(v) >= _vtuple(GENAI_MIN):
            return py, v, None
        tried.append(f"{py}: " + (f"google-genai {v}" if v else "no google-genai"))
    return None, None, (f"no python here has google-genai >= {GENAI_MIN}, which Gemini 3 TTS "
                        f"needs ({'; '.join(tried) or 'none found'}). `bash "
                        f"{os.path.join(HERE, 'setup.sh')}` installs it into the pipeline env; "
                        "or name a python with --python")


def write_job(path, py, script, out, consent_args, models, queue, log_path, name):
    """The job script: render --tts gemini (verify follows), link, share -> its text, written
    mode 700 into the private jobs folder (private_dir). The log path is DOUBLE-quoted: sbatch
    reads `--output="a b/c"` as one path in a real run (HIVE, 2026-10-15); an unquoted space
    ends the directive, and single quotes pass --test-only but leave no log in a real run.
    --open-mode=append: a requeued job adds to its log instead of starting it over."""
    if re.search(r'["\n]', log_path):
        raise RenderError(f"the job's log path cannot hold a quote or a line break: {log_path!r}")
    part, acct, qos, requeue = queue
    q = shlex.quote
    L = ["#!/bin/bash -l", f"#SBATCH --job-name=podcast_{name}",
         f'#SBATCH --output="{log_path}"', "#SBATCH --open-mode=append",
         "#SBATCH --cpus-per-task=1", f"#SBATCH --mem={JOB_MEM}", f"#SBATCH --time={JOB_HOURS}:00:00"]
    L += [f"#SBATCH --{k}={v}" for k, v in (("partition", part), ("account", acct), ("qos", qos))
          if v]
    if requeue:
        L.append("#SBATCH --requeue")
    L += ["# Written by make_podcast.py job: the episode in Gemini voices (render; verify, the ASR",
          "# round trip, follows by itself), linked into the report (link) and built into",
          "# Analysis_Report_with_audio.html (share). Minutes of network calls: a job, never the",
          "# login node. Submitting it again resumes -- rendered chunks are cached in podcast/.cache.",
          "set -uo pipefail",
          "# AAC needs ffmpeg, a module on HIVE (this is a login shell, so `module` exists)",
          f"if command -v module >/dev/null 2>&1; then module load {HIVE_FFMPEG_MODULE} >/dev/null "
          "2>&1 || module load ffmpeg >/dev/null 2>&1; fi",
          "command -v ffmpeg >/dev/null 2>&1 || command -v afconvert >/dev/null 2>&1 || echo "
          f"\"[WARN] no ffmpeg here (module load {HIVE_FFMPEG_MODULE} failed): the episode stays "
          "podcast.wav, about 58 MB for 20 min, and share cannot build the report with the "
          "audio\" >&2",
          f"PY={q(py)}", f"MP={q(os.path.abspath(__file__))}",
          " ".join(['"$PY" "$MP" render', q(script), "--tts", "gemini"]
                   + [q(x) for x in consent_args + models]) + " || exit $?",
          f'"$PY" "$MP" link {q(out)} || exit $?',
          f'"$PY" "$MP" share {q(out)}', ""]
    text = "\n".join(L)
    if os.path.islink(path):
        raise RenderError(f"{path} is a link: not writing through it. Remove it and run job again.")
    # a new file moved into place (safe_write.write_file, the one way a job script is written):
    # never through a link planted after the check above, mode set on the file itself
    from safe_write import write_file
    write_file(path, text, 0o700)
    return text


def cmd_job(a):
    """Write the SLURM job that renders a checked script in Gemini voices, then link and share.
    What the job needs is checked here first -- the check, consent, google-genai, a key -- so a
    missing piece is one line now, not a failed job later."""
    if not os.path.isfile(a.script):
        log(f"[job] refusing: {a.script} is not there: put podcast_script.md in the session's "
            "output/podcast/ first.")
        return 2
    s = parse_script(a.script)
    if s.problems:
        for line, msg in s.problems:
            log(f"[FAIL] {'line %d: ' % line if line else ''}{msg}")
        return 2
    out = output_dir_of(a.script)
    bad = output_dir_problem(out) if out else "the script is not in a podcast/ folder"
    if bad:
        log(f"[job] refusing: {bad}. The job links and shares into the session's output/, so "
            "the script must be output/podcast/podcast_script.md.")
        return 2
    _status, _sources, why = read_check(s)
    if why:
        log(f"[job] refusing: {why}. Run `make_podcast.py check {a.script} --source <the report "
            "files>` and fix every flagged line first.")
        return 2
    m_why = read_meaning(s)[1]
    m_override = " ".join((getattr(a, "no_meaning_check", None) or "").split())
    if m_why and not m_override:
        log(f"[job] refusing: {m_why}. A reader who did not write the script checks every "
            f"sentence first: `make_podcast.py meaning {a.script} --brief`, then `--verdict "
            "PASS --reviewer ... --findings ...` (or --no-meaning-check \"<why>\").")
        return 2
    consent, core_rec, refusal = consent_of(a)
    if not consent:
        log("[job] refusing: " + (refusal or "the job sends the transcript and the audio to "
                                  "Google, so it needs consent: --core-run for a Core run, or "
                                  "--cloud-ok \"<who agreed, when>\" from the person whose data "
                                  "it is."))
        return 2
    py, version, why = job_python(a.python)
    if not py:
        log(f"[job] refusing: {why}")
        return 2
    key, src, why = find_key(bool(core_rec))
    if not key:
        log(f"[job] refusing: no Gemini API key for the job: set GEMINI_API_KEY or put it in "
            f"{KEY_FILE} (chmod 600)" + (f"; {why}" if why else ""))
        return 2
    del key                                          # only where it is was needed
    try:
        import diann_parallel                      # needs_requeue: the one rule
        import run_search                          # slurm_queue: never a hard-coded queue
    except ImportError as e:
        log(f"[job] cannot choose a SLURM queue: {e}. The skill's scripts are incomplete here; "
            "put them again (hive_exec.sh --put-skill).")
        return 2
    part, acct, qos = run_search.slurm_queue(peak_cpus=1, preemptible_ok=True)
    requeue = diann_parallel.needs_requeue(part, qos)
    session = _session_of(out)
    jobs = os.path.join(pipeline_home(), "jobs")
    bad = private_dir(jobs)
    if bad:
        log(f"[job] refusing: {bad}. The job script must sit where only you can change it.")
        return 2
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", (core_rec or {}).get("prot")
                  or os.path.basename(session))[:40]
    stem = f"podcast_{name}_{sha256_bytes(os.path.realpath(session).encode('utf-8'))[:8]}"
    path, log_path = os.path.join(jobs, stem + ".sh"), os.path.join(jobs, stem + "_%j.log")
    if core_rec:
        consent_args = ["--core-run"] + (["internal"] if core_rec["kind"] == "internal" else [])
    else:
        consent_args = ["--cloud-ok"] + ([consent] if isinstance(consent, str) else [])
    models = [x for m in (a.model or []) for x in ("--model", m)]
    if m_why and m_override:                       # the render records the override
        models += ["--no-meaning-check", m_override]
    try:
        write_job(path, py, os.path.abspath(a.script), out, consent_args, models,
                  (part, acct, qos, requeue), log_path, name)
    except (RenderError, OSError) as e:
        log(f"[job] refusing: {e}")
        return 2
    print(json.dumps({"job": path, "submit": f"sbatch {shlex.quote(path)}",
                      "log": log_path.replace("%j", "<jobid>"),
                      "python": py, "google_genai": version, "key": src,
                      "consent": (f"Core run ({core_rec['kind']}"
                                  + (f", {core_rec['prot']}" if core_rec.get("prot") else "")
                                  + "): CORE_CLOUD_OK") if core_rec else "--cloud-ok",
                      "queue": {"partition": part, "account": acct, "qos": qos,
                                "requeue": requeue},
                      "steps": ["render --tts gemini (then verify)", "link", "share"]}, indent=2))
    return 0


def cmd_key(a):
    """Where the Gemini key would come from, and whether the Core's shared key can be used.
    Never the key."""
    key, src, why = find_key(a.core_run)
    core, core_why = core_key()
    print(json.dumps({"found": bool(key), "source": src, "core_run": a.core_run,
                      "core_key": ("usable, for a Core run (--core-run)" if core else
                                   core_why or f"none here ({CORE_KEY_FILE})"),
                      "note": why if not key else None}, indent=2))
    return 0 if key else 1


# ----------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    c = sub.add_parser("check", help="check a script against the report (writes check.txt)")
    c.add_argument("script")
    c.add_argument("--source", nargs="+", default=[], action="extend",
                   help="the report files the script is from")
    c.add_argument("--output-dir", help="the session's output folder (default: the folder "
                                        "above the script's podcast/ folder)")
    c.add_argument("--extra-source", action="append", metavar="FILE",
                   help="a source the lab was NOT given (notes, a draft); needs its own --label, "
                        "and anything only in it must be a Claims bullet naming it")
    c.add_argument("--label", action="append", help="what the matching --extra-source is")
    c.add_argument("--forbid-name", nargs="+", default=[], action="extend",
                   help="names that must not appear (the PI, staff, collaborators)")
    def consent_args(p, cloud_help):
        g = p.add_mutually_exclusive_group()
        g.add_argument("--cloud-ok", nargs="?", const=True, default=False, help=cloud_help)
        g.add_argument("--core-run", nargs="?", const="submission", choices=CORE_RUN_KINDS,
                       help="a Core run: the consent is the Core director's standing approval "
                            "(CORE_CLOUD_OK, recorded in podcast.json), and the Core's shared "
                            "key on HIVE may be used. 'submission' (the default) needs the "
                            "session to answer a CoreOmics submission; 'internal' needs its "
                            "session.json to say internal Core work. Both need you on the Core "
                            "staff list")

    r = sub.add_parser("render", help="render a checked script to audio")
    r.add_argument("script")
    r.add_argument("--out", help="output folder (default: the script's folder)")
    r.add_argument("--tts", required=True, choices=sorted(BACKENDS))
    consent_args(r, "consent to send the transcript to the cloud TTS service; optionally a "
                    "note of who agreed and when (recorded in podcast.json)")
    r.add_argument("--model", action="append", help="Gemini TTS model(s) to use, in order "
                                                    "(default: read from the API)")
    r.add_argument("--rate-wait-budget", type=float, default=RATE_WAIT_BUDGET_MIN,
                   metavar="MIN", help="minutes of rate-limit waiting allowed in one render "
                                       "before it stops with the resume message (default "
                                       "%(default)s)")
    r.add_argument("--min-interval", type=float, default=MIN_INTERVAL,
                   help="seconds between requests to one Gemini model (default %(default)s)")
    r.add_argument("--request-timeout", type=float, default=REQUEST_TIMEOUT, metavar="SEC",
                   help="seconds one TTS request's answer is waited for; asked once more, then "
                        "the next candidate model makes every chunk again (with none, the render "
                        "stops with the resume message) (default %(default)s)")
    r.add_argument("--request-timeout-budget", type=float, default=REQUEST_TIMEOUT_BUDGET_MIN,
                   metavar="MIN", help="minutes lost to unanswered requests allowed in one render "
                                       "before it stops with the resume message (default "
                                       "%(default)s)")
    r.add_argument("--no-meaning-check", metavar="WHY",
                   help="render although meaning_check.md has no PASS for this script; the reason "
                        "goes into podcast.json")
    r.add_argument("--no-preflight", action="store_true",
                   help="skip the two-line test each candidate model gets before the first chunk "
                        "this render makes (recorded in podcast.json)")
    r.add_argument("--keep-wav", action="store_true")
    r.add_argument("--unchecked", action="store_true",
                   help="render without a passing check.txt (recorded in podcast.json)")
    r.add_argument("--redo", type=int, nargs="+", metavar="SEGMENT",
                   help="re-make these segments (numbered from 1) although their text is "
                        "unchanged: audio that verify or a listener found dropped or garbled")
    r.add_argument("--no-verify", action="store_true",
                   help="skip the ASR round trip that follows a --tts gemini render")
    r.add_argument("--verify-model", action="append", help="text model(s) for verify")
    r.add_argument("--dry-run", action="store_true",
                   help="print the chunks as they would be spoken; send nothing")
    mn = sub.add_parser("meaning", help="the meaning check: --brief for the reader who did not "
                                        "write the script; --verdict to record what they found")
    mn.add_argument("script")
    mn.add_argument("--brief", action="store_true",
                    help="print the reviewer's brief (and write meaning_brief.md beside the script)")
    mn.add_argument("--source", nargs="+", default=[], action="extend",
                    help="the report files to name in the brief")
    mn.add_argument("--output-dir", help="the session's output folder (default: the folder above "
                                         "the script's podcast/ folder)")
    mn.add_argument("--verdict", choices=("PASS", "FAIL"))
    mn.add_argument("--reviewer", help="who read it: an agent or person that did not write it")
    mn.add_argument("--findings", help="the reviewer's findings, one line per problem (or `None`)")
    v = sub.add_parser("verify", help="transcribe the rendered audio and compare it with the "
                                      "spoken script (sends the audio to Google)")
    v.add_argument("script")
    v.add_argument("--out", help="the folder with podcast.json (default: the script's folder)")
    consent_args(v, "consent to send the audio to Gemini for transcription")
    v.add_argument("--model", action="append", help="text model(s) to transcribe with")
    k = sub.add_parser("link", help="link the podcast into the report, README and AGENTS.md")
    k.add_argument("outdir", help="the folder with Analysis_Report.html and podcast/")
    k.add_argument("--unchecked", action="store_true",
                   help="link although the script's check no longer holds (the report or the "
                        "script changed since)")
    h = sub.add_parser("share", help="write Analysis_Report_with_audio.html: the report with the "
                                     "audio and transcript built in, one file to send")
    h.add_argument("outdir", help="the folder with Analysis_Report.html and podcast/")
    h.add_argument("--kbps", type=int, default=SHARE_KBPS,
                   help="AAC bitrate for the built-in speech, %d-%d (default %%(default)s)"
                        % SHARE_KBPS_RANGE)
    h.add_argument("--unchecked", action="store_true",
                   help="build it although the script's check no longer holds")
    h.add_argument("--no-browser", action="store_true",
                   help="skip the headless-browser playback test")
    j = sub.add_parser("job", help="write the SLURM job (HIVE) that renders the checked script in "
                                   "Gemini voices, verifies, links and shares it")
    j.add_argument("script", help="<session>/output/podcast/podcast_script.md, checked")
    consent_args(j, "consent to send the transcript and audio to Gemini: a note of who agreed "
                    "and when (recorded in podcast.json)")
    j.add_argument("--model", action="append", help="Gemini TTS model(s) for render, in order")
    j.add_argument("--no-meaning-check", metavar="WHY",
                   help="render although no PASS meaning check is recorded for this script; the "
                        "reason goes into podcast.json")
    j.add_argument("--python", help="the python the job runs (default: the pipeline env's, from "
                                    "setup.json, else this one; it needs google-genai >= "
                                    + GENAI_MIN + ")")
    y = sub.add_parser("key", help="say where the Gemini key would come from (never the key)")
    y.add_argument("--core-run", action="store_true",
                   help="for a Core run, which may use the Core's shared key on HIVE")
    a = ap.parse_args(argv)
    if not a.cmd:
        ap.print_help()
        return 2
    try:
        return {"check": cmd_check, "meaning": cmd_meaning, "render": cmd_render,
                "verify": cmd_verify,
                "link": cmd_link, "share": cmd_share, "job": cmd_job, "key": cmd_key}[a.cmd](a)
    except RenderError as e:
        log(f"[{a.cmd}] {e}")
        return 1
    except KeyboardInterrupt:
        log(f"[{a.cmd}] interrupted; finished chunks are cached")
        return 130
    except Exception:                                   # never let a traceback carry the key
        log(f"[{a.cmd}] unexpected error:\n{traceback.format_exc()}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
