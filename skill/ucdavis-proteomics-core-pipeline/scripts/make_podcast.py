#!/usr/bin/env python3
"""
make_podcast.py -- an OPTIONAL audio discussion of a finished report: two synthetic hosts, a
cell biologist and a statistician, talk through the results ("Signal to Noise").

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
  render SCRIPT.md [--out DIR] --tts gemini|say [--cloud-ok [NOTE]] [--model M] [--keep-wav]
         -> DIR/podcast.m4a (podcast.wav when neither afconvert nor ffmpeg is present),
            transcript.html, podcast.json, podcast_script.md. Refuses unless check.txt says PASS
            for this exact script and unchanged sources (--unchecked overrides; podcast.json
            then says so). Every chunk is cached as DIR/.cache/<sha256>.wav, so a re-run resumes
            and editing one line re-synthesizes only its chunk; unused chunks are pruned after
            a successful render. --dry-run prints what would be spoken and exits.
  verify SCRIPT.md --cloud-ok NOTE [--out DIR] [--model M]
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
            check does not hold (--unchecked overrides); reprints an older Analysis_Report.pdf.
            The hooks the report calls never raise: a bad podcast.json is one [WARN].

Privacy: only the final transcript (render: the turns, after pronunciation substitutions) and
the rendered audio (verify, downsampled) ever leave the machine -- never the report -- and only
with --cloud-ok (explicit consent, recorded in podcast.json). --tts say (macOS) is offline. The Gemini key is read from
GEMINI_API_KEY or ~/.config/ucdavis-proteomics/gemini_key; it is never printed or logged and is
scrubbed from every error message (notify_slack.redact: the skill's one list of secret patterns).

Stdlib only at import time: google-genai is imported by the gemini backend alone. No numpy and
no ffmpeg requirement -- audio is assembled with array/wave; afconvert (macOS) or ffmpeg, when
present, encodes the AAC.
"""
import argparse
import array
import base64
import binascii
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
import shutil
import struct
import subprocess
import sys
import tempfile
import time
import traceback
import urllib.parse
import wave

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
from notify_slack import REDACTED, redact as _redact    # noqa: E402  the ONE secret list

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
TERMS_URL = "https://ai.google.dev/gemini-api/terms"
KEY_FILE = os.path.join("~", ".config", "ucdavis-proteomics", "gemini_key")

RATE = 24000                                        # Hz, 16-bit mono throughout
WPM = 150                                           # spoken words per minute, for estimates
WORDS_LO, WORDS_HI = 2200, 3600
SEGMENTS_LO, SEGMENTS_HI = 6, 10
CHUNK_WORDS = 500                                   # a segment longer than this is split
GAP_SEGMENT, GAP_CHUNK, GAP_TURN = 0.7, 0.3, 0.25   # seconds of silence

START, END = "<!-- podcast:start -->", "<!-- podcast:end -->"
T_START, T_END = "<!-- TRANSCRIPT START -->", "<!-- TRANSCRIPT END -->"
DISCLOSURE = ("AI-generated discussion of these results; voices are synthetic; the report is "
              "the record — check anything important against it.")
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
NEXT_STEPS = r"\.(?:html|csv|md|pdf)\b|Analysis[_ ]Report|PropObs|\bvalidat|\btier"

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
    Listen line to the report does not make the check stale."""
    with open(path, "rb") as fh:
        return strip_block(fh.read().decode("utf-8", errors="replace"))


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
def check(s, sources, forbid=()):
    """-> dict with fails / warns / infos (lists of str) and stats. `sources`: [(path, sha,
    text)]."""
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
    if not sources:
        fails.append("no --source given: nothing to check the numbers and symbols against")

    book = NumberBook()
    hay = []
    for path, sha, text in sources:
        book.add_text(text)
        hay.append(normalize_numbers(text).lower())
    hay = "\n".join(hay)
    phrases = _words_norm(hay)
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
            if re.search(r"\b" + re.escape(m.group(0).lower()) + r"\b", s.claims_text.lower()):
                infos.append(f"{where}: '{m.group(0)}' is listed under Claims beyond the report")
                continue
            fails.append(f"{where}: quantity in words '{ph or m.group(0)}' is not in the sources "
                         "-- use the report's own phrase or its digits (the pronunciation step "
                         "handles speech), or list it under Claims beyond the report if it is "
                         "your own gloss")
        for m in SPECIALTY.finditer(t.text):
            fails.append(f"{where}: '{m.group(0)}' -- a host must not claim a real research "
                         "specialty (it lends a synthetic voice false authority); say \"I'm the "
                         "biologist of the pair\" or \"I'm the statistician\"")

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
        warns.append("the last two segments never say what to do with this: which file to open, "
                     "how to tier hits by PropObs, what to validate first")

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
            "infos": infos, "stats": stats, "pron": pron}


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
         f"script: {s.path}",
         f"script_sha256: {s.sha256}",
         f"checked: {now_iso()}"]
    L += [f"source: {p} sha256={sha}" for p, sha, _ in sources]
    L.append("forbidden names: " + (", ".join(forbid) if forbid else "none given"))
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
                  "only re-spells its written form. It does NOT mean the script is right."]
    L += ["", "What check cannot catch -- read the script against the report for these:"]
    L += [f"- {x}" for x in CANNOT_CATCH]
    return "\n".join(L) + "\n"


def cmd_check(a):
    s = parse_script(a.script)
    sources = []
    missing = []
    for p in a.source:
        if not os.path.isfile(p):
            missing.append(p)
            continue
        sha, text = load_source(p)
        sources.append((os.path.abspath(p), sha, text))
    res = check(s, sources, a.forbid_name or [])
    for p in missing:
        res["fails"].insert(0, f"source not found: {p}")
        res["status"] = "FAIL"
    out = os.path.join(os.path.dirname(s.path), "check.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(check_report(s, res, sources, a.forbid_name or []))
    print(json.dumps({"status": res["status"], "fails": len(res["fails"]),
                      "warnings": len(res["warns"]), "check": out, **res["stats"]}, indent=2))
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
    srcs = [{"file": m.group(1), "sha256": m.group(2)}
            for m in re.finditer(r"(?m)^source: (.+) sha256=(\w+)$", text)]
    if not status or status.group(1) != "PASS":
        return (status.group(1) if status else None), srcs, "check.txt does not say PASS"
    if not sha or sha.group(1) != s.sha256:
        return "STALE", srcs, "the script changed after check.txt was written"
    if not srcs:
        return "STALE", srcs, "check.txt names no source"
    for src in srcs:
        try:
            now = source_sha(source_text(src["file"]))
        except OSError:
            return "STALE", srcs, f"source {src['file']} is missing"
        if now != src["sha256"]:
            return "STALE", srcs, (f"source {os.path.basename(src['file'])} changed after "
                                   "check.txt was written")
    return "PASS", srcs, None


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


class GeminiTTS(object):
    """Gemini multi-speaker TTS, one request per chunk. Model names are read from the API
    (models.list), never assumed. Default order: a pro TTS model, then gemini-3.8-flash-tts,
    then gemini-2.5-flash-preview-tts. A model missing (404), without quota, or rejecting the
    request before it has made a chunk hands over to the next. Once a model has made a chunk
    the episode stays on it, because a second model's voices sound different.

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
        self.key = read_key()
        self.types = types
        self.client = genai.Client(api_key=self.key)
        self.user_models = list(a.model or [])
        self.candidates = self.user_models or self.discover()
        self.model = self.candidates[0]
        self.min_interval = float(getattr(a, "min_interval", MIN_INTERVAL))
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

    def prepare(self, chunks, speak, cache_dir):
        """Stay on the model the cache was made with: resuming must not mix voices."""
        if self.user_models:
            return
        have = {m: sum(os.path.isfile(os.path.join(cache_dir, cache_key(self, c, speak, m) +
                                                   ".wav")) for c in chunks)
                for m in self.candidates}
        best = max(self.candidates, key=lambda m: have[m])
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

    def _has_next(self):
        return (self.model in self.candidates
                and self.candidates.index(self.model) + 1 < len(self.candidates))

    def _advance(self):
        if self.produced or not self._has_next():
            return False
        self.model = self.candidates[self.candidates.index(self.model) + 1]
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
        errors, waits, empty, retried_len = 0, 0, 0, False
        while True:
            self._pace()
            try:
                r = self.request(chunk, speak)
            except Exception as e:
                err = error_info(e)
                kind = err["kind"]
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


def error_info(e):
    """-> {kind, msg, delay, quota}. kind: missing | exhausted (daily or zero quota) | rate
    (per-minute, or any 429 that gives a retryDelay) | rejected (400: the request does not fit
    this model) | fatal (key, permission) | transient. The retry delay and the quota come from
    google.rpc.RetryInfo / QuotaFailure in the error details when present, else from the text."""
    msg = scrub(f"{type(e).__name__}: {e}")
    code = next((c for c in (getattr(e, "code", None), getattr(e, "status_code", None))
                 if isinstance(c, int)), None)
    found = {"ids": [], "values": [], "delay": None}

    def walk(x):
        if isinstance(x, dict):
            for k, v in x.items():
                if k == "quotaId" and isinstance(v, str):
                    found["ids"].append(v)
                elif k == "quotaValue":
                    found["values"].append(str(v))
                elif k == "retryDelay" and found["delay"] is None:
                    found["delay"] = _seconds(v)
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
    info = {"msg": msg, "delay": found["delay"], "quota": quota}
    if code == 404 or "NOT_FOUND" in text:
        return dict(info, kind="missing")
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


def read_key():
    """The Gemini key: GEMINI_API_KEY, else ~/.config/ucdavis-proteomics/gemini_key. Never
    printed."""
    key = (os.environ.get("GEMINI_API_KEY") or "").strip()
    path = os.path.expanduser(KEY_FILE)
    if not key and os.path.isfile(path):
        mode = os.stat(path).st_mode & 0o777
        if mode & 0o077:
            log(f"[WARN] {KEY_FILE} is readable by other accounts (mode {oct(mode)[2:]}): "
                f"chmod 600 {KEY_FILE}")
        with open(path, encoding="utf-8") as fh:
            key = fh.read().strip()
    if not key:
        raise RenderError(f"no Gemini API key: set GEMINI_API_KEY or put it in {KEY_FILE} "
                          "(chmod 600)")
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


# ----------------------------------------------------------------------------- render
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
    consent = a.cloud_ok
    refused = isinstance(consent, str) and consent.strip().lower() in NO_CONSENT
    if refused:
        consent = False
    if backend_cls.cloud and not consent:
        if refused:
            log(f"[render] not sending anything: --cloud-ok {a.cloud_ok!r} is a refusal.")
        log(f"[render] not sending anything. --tts {a.tts} sends the transcript -- {words:,} "
            f"words of the turns, after pronunciation substitutions; never the report -- to "
            f"Google's Gemini API, and verify then sends the rendered AUDIO (downsampled) for "
            f"transcription. On a free-tier key Google may use both to improve its products "
            f"and human reviewers may read it ({TERMS_URL}). Once the user has agreed, re-run "
            f"with --cloud-ok (optionally --cloud-ok \"who agreed, when\"), or use --tts say "
            f"(offline, macOS).")
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

    # Each chunk: <key>.wav, and <key>.json with the warnings made when it was synthesized, so
    # a resumed render still reports them.
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
        chk = os.path.join(os.path.dirname(s.path), "check.txt")
        if os.path.isfile(chk) and os.path.abspath(chk) != os.path.join(out, "check.txt"):
            shutil.copyfile(chk, os.path.join(out, "check.txt"))
    man = {
        "show": s.show, "title": s.title,
        "hosts": [{"name": h, "role": r or None, "voice": backend.voices.get(h)} for h, r in s.hosts],
        "tts": {"backend": backend.name, "model": (backend.models_used or [backend.model])[0],
                "models_used": backend.models_used, "sample_rate": RATE,
                "sent_to_cloud": bool(backend.cloud),
                "what_was_sent": ("the transcript turns after pronunciation substitutions, plus "
                                  "one style line; then, for verify (unless --no-verify), the "
                                  "rendered audio downsampled to 16 kHz; never the report"
                                  if backend.cloud else "nothing (offline)")},
        "cloud_tts_consent": (consent if backend.cloud else False),
        "script": "podcast_script.md", "script_sha256": s.sha256,
        "script_author": s.author,
        "check": {"status": status or "not run", "file": "check.txt",
                  "overridden": bool(why and a.unchecked), "reason": why},
        "sources": sources,
        "words": words, "turns": len(s.turns()), "segments": len(s.segments),
        "chunks": len(chunks), "chunks_cached": n_cached, "chunks_synthesized": n_new,
        "duration_s": duration, "created": now_iso(), "ai_generated": True,
        "audio": audio_name, "encoder": tool, "transcript": "transcript.html",
        "warnings": warnings,
        "made_by": "make_podcast.py (UC Davis Proteomics Core pipeline skill)",
    }
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
    if backend.cloud and consent and not getattr(a, "no_verify", False):
        try:
            run_verify(s.path if os.path.dirname(s.path) == out else script_copy, out,
                       getattr(a, "verify_model", None), consent)
        except Exception as e:
            log(f"[WARN] verify skipped: {scrub(e)} (the render is fine; run "
                f"`make_podcast.py verify {script_copy} --cloud-ok ...` to retry)")
    return 0


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


def verify_tokens(text):
    """-> [(word, start, end)]: text as comparable words. Lowercase; numbers without thousands
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

    heard = {n.value for n in numbers_in(transcript)}
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

    def __init__(self, models=None):
        try:
            from google import genai
            from google.genai import types
        except ImportError:
            raise VerifyError("the verify step needs the google-genai package: python3 -m pip "
                              "install google-genai")
        self.types = types
        self.client = genai.Client(api_key=read_key())
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
                    if err["kind"] == "fatal":
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


def run_verify(script_path, out=None, models=None, consent=True):
    """Transcribe the rendered audio and compare it with the spoken script. Writes
    verify.txt, verify_transcript.txt and a `verify` block in podcast.json. -> the result dict.
    Raises VerifyError when it cannot run; never touches the report."""
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
    text, model = Transcriber(models).transcribe(data, mime)
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
    consent = a.cloud_ok
    if isinstance(consent, str) and consent.strip().lower() in NO_CONSENT:
        consent = False
    if not consent:
        log(f"[verify] not sending anything. verify sends the rendered AUDIO (downsampled) to "
            f"Google's Gemini API for transcription. On a free-tier key Google may use it to "
            f"improve its products and human reviewers may hear it ({TERMS_URL}). Once the user "
            f"has agreed, re-run with --cloud-ok \"who agreed, when\".")
        return 2
    try:
        res = run_verify(a.script, a.out, a.model, consent)
    except (VerifyError, RenderError) as e:
        log(f"[verify] could not verify: {scrub(e)}")
        return 1
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


def transcript_html(s, man, out):
    sys.path.insert(0, HERE)
    import report_style as rs
    title = f"{s.show} — {s.title or 'Untitled'}"
    report = os.path.join(os.path.dirname(out), "Analysis_Report.html")
    what = ('<a href="../Analysis_Report.html">the analysis report</a>' if os.path.isfile(report)
            else "the analysis report")
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
                                    else "") + " — check.txt")],
            ["Made", esc(man["created"])],
            ["Script sha256", f"<code>{esc(s.sha256[:16])}…</code>"]]
    pron = (rs.table(["Written", "Spoken"], [[esc(w), esc(sp)] for w, sp in s.pronunciation])
            if s.pronunciation else "<p>Only the built-in substitutions.</p>")
    body = (f"<header><h1>{esc(title)}</h1><p>An AI-generated audio discussion of {what}: "
            f"two synthetic hosts talk through the results.</p></header>\n"
            + rs.callout("info", f"<p>{esc(DISCLOSURE)} The hosts are fictional, their voices are "
                                 "synthetic, and the script was written by an AI from the report; "
                                 "anything said beyond the report is listed below.</p>",
                         title="AI-generated")
            + f'<audio controls preload="none" src="{href(man["audio"])}"></audio>'
            + f"<style>{TRANSCRIPT_CSS}</style>"
            + '<section><h2 id="transcript">Transcript</h2>' + "".join(segs) + "</section>"
            + '<section><h2 id="claims">Claims beyond the report</h2><p>Statements the script '
              "makes that do not come from the report (general knowledge, metaphors, "
              "speculation). Check these before relying on them.</p>" + claims + "</section>"
            + '<section><h2 id="made">How this was made</h2>' + rs.table(["", ""], rows)
            + "<p>Pronunciation substitutions change only what the voices say; the transcript "
              "above keeps the written form.</p>" + pron + "</section>")
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
    "@media print{.pc-card audio{display:none}.pc-card .pc-print{display:block}"
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


def _paths(outdir, base, man):
    pdir = os.path.join(outdir, "podcast")
    return (rel(os.path.join(pdir, man["audio"]), base),
            rel(os.path.join(pdir, man.get("transcript") or "transcript.html"), base))


def _name(man):
    return f"{man.get('show') or SHOW} — {man.get('title') or 'Untitled'}"


@_never_raise("")
def listen_card_html(outdir, base_dir=None):
    """The "Listen" card (with its markers) for a page in `base_dir` (default: `outdir`, where
    Analysis_Report.html sits). '' when `outdir` has no podcast/podcast.json."""
    man = load_manifest(outdir)
    if not man:
        return ""
    audio, tr = _paths(outdir, base_dir or outdir, man)
    return (f"{START}<style>{CARD_CSS}</style>"
            f'<aside class="pc-card" aria-label="Audio discussion of these results">'
            f'<div class="pc-h">{_HEADPHONES}<span class="pc-k">Listen</span>'
            f"<span>{esc(_name(man))}</span>"
            f'<span class="pc-d">{esc(minutes_label(man.get("duration_s")))}</span></div>'
            f'<audio controls preload="none" src="{href(audio)}">Open '
            f'<a href="{href(audio)}">{esc(audio)}</a> to listen.</audio>'
            f"<p>{esc(DISCLOSURE)} <a href=\"{href(tr)}\">Read the transcript</a> (with the claims "
            f"that go beyond the report). The audio is a separate file in the podcast folder "
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
                             f"[{audio}]({href(audio)}). {DISCLOSURE} Transcript: [{tr}]({href(tr)}).",
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


@_never_raise([])
def agents_md_lines(outdir, base_dir):
    """AGENTS.md's section on the podcast. [] without a podcast."""
    man = load_manifest(outdir)
    if not man:
        return []
    pdir = rel(os.path.join(outdir, "podcast"), base_dir)
    return ["## Audio discussion (podcast): a derivative, not a record", "",
            f"`{pdir}/` holds an AI-generated audio discussion of this analysis "
            f"(*{man.get('show') or SHOW}*: two synthetic hosts, a script written by an AI from "
            "the report). It is a derivative of the report and NOT authoritative: never cite it or "
            "take a number from it -- use the report and the tables. Its transcript is "
            f"`{pdir}/transcript.html`; the script, with its \"Claims beyond the report\" ledger "
            f"(statements not in the report), is `{pdir}/podcast_script.md`; the automated "
            f"fidelity check is `{pdir}/check.txt`; how it was made (TTS service, model, voices, "
            f"consent) is `{pdir}/podcast.json`."]


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
def add_listen_card(doc, outdir):
    """`doc` with the Listen card near the top; unchanged when `outdir` has no podcast. The
    hook make_analysis_html.py calls, so regenerating the report keeps the card."""
    card = listen_card_html(outdir)
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
    # The Listen card vouches for the episode: link only while the check still holds.
    script = os.path.join(out, "podcast", man.get("script") or "podcast_script.md")
    stale = (read_check(parse_script(script))[2] if os.path.isfile(script) else
             f"no {os.path.basename(script)} beside podcast.json")
    overridden = bool((man.get("check") or {}).get("overridden"))
    if stale and not (a.unchecked or overridden):
        log(f"[link] refusing: {stale}. Re-run `make_podcast.py check` (and `render`, if the "
            "script changed) so the episode matches the report, or pass --unchecked.")
        return 2
    if stale:
        log(f"[WARN] linking an episode whose check does not hold: {stale}"
            + (" (rendered with --unchecked)" if overridden else " (--unchecked)"))
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
    done.append(refresh_pdf(os.path.join(out, "Analysis_Report.html")))
    for level, path, note in done:
        print(f"[{level}] {rel(path, os.path.dirname(out))}: {note}")
    return 0


def refresh_pdf(html_path):
    """Analysis_Report.pdf is printed from the HTML (html_to_pdf.py, where the skill has it):
    once link has put the card in the HTML, an older PDF lacks it. Reprint it; the card's print
    style shows the audio's file name in place of the player. -> (level, path, note)"""
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
    ok, note = html_to_pdf.convert(html_path, pdf)
    return (("OK", pdf, f"reprinted with the Listen card ({note})") if ok else
            ("INFO", pdf, f"older than the HTML and NOT reprinted: {note}"))


# ----------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    c = sub.add_parser("check", help="check a script against the report (writes check.txt)")
    c.add_argument("script")
    c.add_argument("--source", nargs="+", default=[], action="extend",
                   help="the report files the script is from")
    c.add_argument("--forbid-name", nargs="+", default=[], action="extend",
                   help="names that must not appear (the PI, staff, collaborators)")
    r = sub.add_parser("render", help="render a checked script to audio")
    r.add_argument("script")
    r.add_argument("--out", help="output folder (default: the script's folder)")
    r.add_argument("--tts", required=True, choices=sorted(BACKENDS))
    r.add_argument("--cloud-ok", nargs="?", const=True, default=False,
                   help="consent to send the transcript to the cloud TTS service; optionally a "
                        "note of who agreed and when (recorded in podcast.json)")
    r.add_argument("--model", action="append", help="Gemini TTS model(s) to use, in order "
                                                    "(default: read from the API)")
    r.add_argument("--rate-wait-budget", type=float, default=RATE_WAIT_BUDGET_MIN,
                   metavar="MIN", help="minutes of rate-limit waiting allowed in one render "
                                       "before it stops with the resume message (default "
                                       "%(default)s)")
    r.add_argument("--min-interval", type=float, default=MIN_INTERVAL,
                   help="seconds between requests to one Gemini model (default %(default)s)")
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
    v = sub.add_parser("verify", help="transcribe the rendered audio and compare it with the "
                                      "spoken script (sends the audio to Google)")
    v.add_argument("script")
    v.add_argument("--out", help="the folder with podcast.json (default: the script's folder)")
    v.add_argument("--cloud-ok", nargs="?", const=True, default=False,
                   help="consent to send the audio to Gemini for transcription")
    v.add_argument("--model", action="append", help="text model(s) to transcribe with")
    k = sub.add_parser("link", help="link the podcast into the report, README and AGENTS.md")
    k.add_argument("outdir", help="the folder with Analysis_Report.html and podcast/")
    k.add_argument("--unchecked", action="store_true",
                   help="link although the script's check no longer holds (the report or the "
                        "script changed since)")
    a = ap.parse_args(argv)
    if not a.cmd:
        ap.print_help()
        return 2
    try:
        return {"check": cmd_check, "render": cmd_render, "verify": cmd_verify,
                "link": cmd_link}[a.cmd](a)
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
