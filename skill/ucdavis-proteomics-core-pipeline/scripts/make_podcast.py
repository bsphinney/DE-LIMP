#!/usr/bin/env python3
"""
make_podcast.py -- an OPTIONAL audio discussion of a finished report: two synthetic hosts, a
cell biologist and a statistician, talk through the results ("Signal to Noise").

The agent running the skill WRITES the script itself, from the finished report it has read,
following references/podcast.md. This script never writes a word of the conversation: it
CHECKS the script, RENDERS it to audio and LINKS the audio into the outputs.

  check  SCRIPT.md --source FILE [FILE ...] [--forbid-name NAME ...]
         -> check.txt beside the script; exit 1 on any FAIL. Every number in the transcript
            must be in the sources; every gene/protein-symbol-like token must be in the sources
            or listed under "Claims beyond the report"; the AI disclosure is in the first 3
            turns; the required sections are there; no forbidden name appears; length is
            reported (words, turns, segments, minutes at 150 wpm).
  render SCRIPT.md [--out DIR] --tts gemini|say [--cloud-ok [NOTE]] [--model M] [--keep-wav]
         -> DIR/podcast.m4a (podcast.wav when neither afconvert nor ffmpeg is present),
            transcript.html, podcast.json, podcast_script.md. Refuses unless check.txt says PASS
            for this exact script (--unchecked overrides; podcast.json then says so). Every chunk
            is cached as DIR/.cache/<sha256>.wav, so a re-run resumes and editing one line
            re-synthesizes only its chunk. --dry-run prints what would be spoken and exits.
  link   OUTDIR
         -> a "Listen" card near the top of OUTDIR/Analysis_Report.html, a line near the top of
            the Markdown report, an entry in README.html / README.md and AGENTS.md. Idempotent:
            the card sits between <!-- podcast:start --> and <!-- podcast:end --> and is replaced
            on a re-run. A file that is not there is skipped with [INFO].

Privacy: only the final transcript (the turns, after pronunciation substitutions) ever leaves
the machine -- never the report -- and only with --tts gemini AND --cloud-ok (explicit consent,
recorded in podcast.json). --tts say (macOS) is offline. The Gemini key is read from
GEMINI_API_KEY or ~/.config/podcast/gemini_key; it is never printed or logged and is scrubbed
from every error message.

Stdlib only at import time: google-genai is imported by the gemini backend alone. No numpy and
no ffmpeg requirement -- audio is assembled with array/wave; afconvert (macOS) or ffmpeg, when
present, encodes the AAC.
"""
import argparse
import array
import base64
import datetime
import decimal
import hashlib
import html
import io
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
import urllib.parse
import wave

HERE = os.path.dirname(os.path.abspath(__file__))

SHOW = "Signal to Noise"
HOSTS = (("Maya", "cell biologist"), ("Leo", "statistician"))
GEMINI_VOICES = ("Kore", "Charon")                 # first host, second host
SAY_VOICES = ("Samantha", "Daniel")                # installed on macOS by default (say -v '?')
FALLBACK_TTS_MODEL = "gemini-2.5-flash-preview-tts"
KNOWN_TTS_MODELS = ("gemini-2.5-pro-preview-tts", "gemini-3.8-flash-tts", FALLBACK_TTS_MODEL)
HOST_STYLES = ("curious and energetic", "calm, precise and dryly funny")   # first, second host
MIN_INTERVAL = 20.0                                 # s between requests to one Gemini model
RATE_WAITS = 10                                     # per-minute 429 waits allowed per chunk
TERMS_URL = "https://ai.google.dev/gemini-api/terms"
KEY_FILE = os.path.join("~", ".config", "podcast", "gemini_key")

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

SPELLED_NUMBER = re.compile(
    r"\b(eleven|twelve|thirteen|fourteen|fifteen|sixteen|seventeen|eighteen|nineteen|twenty|"
    r"thirty|forty|fifty|sixty|seventy|eighty|ninety|hundred|thousand|million|billion)\b", re.I)

_SECRETS = []
_sleep = time.sleep                                 # tests replace these two
_now = time.monotonic


# ----------------------------------------------------------------------------- small helpers
def log(msg):
    print(scrub(msg), file=sys.stderr, flush=True)


def scrub(text):
    """Remove every API key from a string: the one read at run time, and anything shaped like a
    Google key or a key= URL parameter."""
    s = str(text)
    for k in _SECRETS:
        if k:
            s = s.replace(k, "[key]")
    s = re.sub(r"AIza[0-9A-Za-z_\-]{20,}", "[key]", s)
    return re.sub(r"(?i)\b(key|api_key|x-goog-api-key)=([^&\s\"'()]+)", r"\1=[key]", s)


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
    """kind: plain | sci | pow10. value: the absolute value. dec: decimals written (the
    mantissa's, for sci)."""

    def __init__(self, text, value, kind, dec):
        self.text, self.value, self.kind, self.dec = text, value, kind, dec


def _decimals(s):
    return len(s.split(".", 1)[1]) if "." in s else 0


def numbers_in(text):
    out = []
    for m in _NUM.finditer(normalize_numbers(text)):
        g = m.groupdict()
        try:
            if g["p"] is not None:
                p = g["p"]
                out.append(Num(p, abs(D(p)), "plain", _decimals(p)))
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
    to fewer decimals (never to an integer), or -- for a p-value -- its order of magnitude."""

    def __init__(self):
        self.exact, self.rounded, self.sci_rounded, self.exponents = set(), {}, {}, set()
        self.mantissas = {}                     # 5.65 -> "5.65e-10": for a helpful message

    def add_text(self, text):
        for n in numbers_in(text):
            self.add(n)

    def add(self, n):
        v = n.value
        self.exact.add(v)
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
                self.rounded.setdefault(d, set()).update(_rounded(v, d))

    def verdict(self, n):
        """-> 'trivial' | 'exact' | 'rounded' | 'magnitude' | None (not in the sources)."""
        v = n.value
        if n.kind == "plain" and n.dec == 0 and (v <= 10 or (v <= 100 and v % 10 == 0)):
            return "trivial"
        if v in self.exact:
            return "exact"
        if n.kind == "plain" and n.dec >= 1 and v in self.rounded.get(n.dec, ()):
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
_SYM = re.compile(r"(?<![A-Za-z0-9])[A-Za-z][A-Za-z0-9]*(?:\.\d+)*")


def is_symbol(tok):
    if tok.lower() in GENERIC_TOKENS:
        return False
    if any(c.isdigit() for c in tok):
        return True                                   # Jph3, Kv2.1, FKBP12.6, C1qa, RYR2
    if len(tok) >= 3 and tok.isupper():
        return True                                   # SERCA, VAPA, BSA
    return bool(re.search(r"[a-z][A-Z]", tok))        # IgG, timsTOF, PropObs


def symbol_found(tok, hay):
    """Case-insensitive, as a whole token on the left; a token ending in a digit must not run
    on into more digits (Kcnb2 is not Kcnb20), one ending in a letter may (RyR -> RyR2)."""
    t = tok.lower()
    tail = r"(?![0-9])" if t[-1].isdigit() else ""
    if re.search(r"(?<![a-z0-9])" + re.escape(t) + tail, hay):
        return True
    return len(t) > 3 and t.endswith("s") and symbol_found(tok[:-1], hay)


# ----------------------------------------------------------------------------- sources
def load_source(path):
    """-> (sha256, text) with the parts that are not prose removed: embedded images (a base64
    blob is full of digit runs), script/style, tags, and any podcast block."""
    with open(path, "rb") as fh:
        raw = fh.read()
    t = strip_block(raw.decode("utf-8", errors="replace"))
    t = re.sub(r"data:[\w/+.\-]+;base64,[A-Za-z0-9+/=]+", " ", t)
    if os.path.splitext(path)[1].lower() in (".html", ".htm"):
        t = re.sub(r"(?is)<(script|style)\b.*?</\1\s*>", " ", t)
        t = re.sub(r"(?s)<!--.*?-->", " ", t)
        t = re.sub(r"<[^>]+>", " ", t)
        t = html.unescape(t)
    return sha256_bytes(raw), t


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
    claims_hay = normalize_numbers(s.claims_text).lower()
    disclosed = NumberBook()
    disclosed.add_text(s.claims_text)
    if len(book.exact) > 20000:
        warns.append(f"the sources hold {len(book.exact):,} distinct numbers: a large numeric "
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
            fails.append(f"{where}: spelled-out number '{m.group(0)}' -- write numbers as digits "
                         "so they can be checked (the pronunciation step handles speech)")

    seen = {}
    for i, t in enumerate(turns, 1):
        for m in _SYM.finditer(t.text):
            tok = m.group(0)
            if is_symbol(tok) and tok not in seen:
                seen[tok] = (i, t)
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

    if turns and not any(DISCLOSURE_RE.search(t.text) for t in turns[:3]):
        fails.append("no AI disclosure in the first 3 turns: one host must say, early and plainly, "
                     "that this is an AI-generated discussion (e.g. 'AI-generated')")

    for name in forbid:
        name = name.strip()
        if not name:
            continue
        rx = re.compile(r"(?<![A-Za-z])" + re.escape(name) + r"(?![A-Za-z])", re.I)
        hits = [f"turn {i} (line {t.line})" for i, t in enumerate(turns, 1) if rx.search(t.text)]
        if s.title and rx.search(s.title):
            hits.insert(0, "the title")
        if rx.search(s.claims_text):
            hits.append("Claims beyond the report")
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
            "infos": infos, "stats": stats}


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
    if res["status"] == "PASS":
        L += ["", "Every number and symbol in the transcript is in the sources or disclosed. This "
                  "checks tokens, not meaning: a sentence built from real numbers can still say "
                  "something the report does not. Read it against the report once more."]
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
    """-> (status, sources, problem) from the check.txt beside the script."""
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
    return "PASS", srcs, None


# ----------------------------------------------------------------------------- audio
def _arr(data):
    a = array.array("h")
    a.frombytes(data[: len(data) // 2 * 2])
    if sys.byteorder == "big":
        a.byteswap()
    return a


def _bytes(a):
    if sys.byteorder == "big":
        a = array.array("h", a)
        a.byteswap()
    return a.tobytes()


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


def decode_audio(data, mime=""):
    """PCM or WAV bytes -> array('h') at RATE. Gemini returns raw 16-bit little-endian PCM
    ("audio/L16;codec=pcm;rate=24000")."""
    if data[:4] == b"RIFF":
        with wave.open(io.BytesIO(data), "rb") as w:
            return _from_wave(w)
    m = re.search(r"rate=(\d+)", mime or "")
    return resample(_arr(data), int(m.group(1)) if m else RATE)


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


def write_wav(path, a):
    tmp = path + ".part"
    with wave.open(tmp, "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(RATE)
        w.writeframes(_bytes(a))
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
        self.last_call = {}
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
        raise RenderError(f"the response held no audio (finish reason: {', '.join(why) or 'none'})")

    @staticmethod
    def interaction_audio(r):
        """interaction.output_audio.data: base64 audio/wav (24 kHz mono 16-bit) by default;
        the RIFF header is read, not spliced, so the rate comes from the file itself."""
        out = getattr(r, "output_audio", None)
        data = getattr(out, "data", None) if out is not None else None
        if not data:
            raise RenderError(f"the interaction held no audio (status: "
                              f"{getattr(r, 'status', None) or 'not given'})")
        if isinstance(data, str):
            data = base64.b64decode(data)
        elif isinstance(data, (bytes, bytearray)) and not bytes(data[:4]) == b"RIFF":
            try:
                data = base64.b64decode(data, validate=True)
            except (ValueError, TypeError):
                pass
        mime = getattr(out, "mime_type", None) or ""
        rate = getattr(out, "sample_rate", None)
        if rate and "rate=" not in str(mime):
            mime = f"{mime};rate={rate}"
        return decode_audio(bytes(data), str(mime))

    def request(self, chunk, speak):
        if self.api(self.model) == "interactions":
            r = self.client.interactions.create(
                model=self.model, input=[{"type": "user_input", "content": self.parts(chunk, speak)}],
                response_format={"type": "audio"},
                generation_config={"speech_config": self.speech_config()})
            return self.interaction_audio(r)
        r = self.client.models.generate_content(model=self.model, contents=self.prompt(chunk, speak),
                                                config=self.config())
        return self.audio(r)

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

    def synth(self, chunk, speak, label):
        errors, waits, retried_len = 0, 0, False
        while True:
            self._pace()
            try:
                pcm = self.request(chunk, speak)
            except Exception as e:
                err = error_info(e)
                kind, msg = err["kind"], err["msg"]
                if kind in ("missing", "exhausted", "rejected") and not self.produced:
                    was = self.model
                    if self._advance():
                        why = {"missing": "not available", "exhausted": "no quota",
                               "rejected": "request rejected"}[kind]
                        log(f"[render] {was}: {why} ({msg[:200]}); trying {self.model}")
                        raise SwitchModel()
                    raise RenderError(
                        f"{label}: no Gemini TTS model could be used (tried "
                        f"{', '.join(self.candidates)}; last: {msg[:200]}). Quotas are per Cloud "
                        "project (see https://aistudio.google.com/rate-limit): re-run later "
                        "(finished chunks are cached), use a paid-tier key, or --tts say.")
                if kind == "exhausted":
                    raise RenderError(
                        f"{label}: {self.model} is out of its daily quota ({err['quota'] or msg[:200]}"
                        f"). The chunks already made are cached: re-run the same command later "
                        f"to resume with {self.model}, or pass --model <another> to render EVERY "
                        "chunk with that model (voices differ between models, so one episode "
                        "never mixes them).")
                if kind == "rate":
                    waits += 1
                    if waits > RATE_WAITS:
                        raise RenderError(
                            f"{label}: still rate-limited on {self.model} after {RATE_WAITS} "
                            f"waits ({err['quota'] or msg[:200]}). The chunks already made are "
                            "cached; re-run the same command to resume.")
                    wait = (err["delay"] if err["delay"] is not None else 60.0) + 5.0
                    log(f"[render] rate-limited, waiting {wait:.0f}s ({label}, {self.model}"
                        + (f", {err['quota']}" if err["quota"] else "") + f"; wait {waits} of "
                        f"{RATE_WAITS})")
                    _sleep(wait)
                    continue
                errors += 1
                if kind in ("fatal", "rejected") or errors >= 4:
                    raise RenderError(f"{label} failed after {errors} error(s) on "
                                      f"{self.model}: {msg[:300]}. The chunks already made are "
                                      "cached; re-run the same command to resume.")
                wait = err["delay"] if err["delay"] is not None else min(60.0, 5.0 * 2 ** (errors - 1))
                log(f"[render] {label}: {kind} error on {self.model} ({msg[:120]}); retrying in "
                    f"{wait:.0f} s")
                _sleep(wait)
                continue
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
    """The Gemini key: GEMINI_API_KEY, else ~/.config/podcast/gemini_key. Never printed."""
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
    if backend_cls.cloud and not a.cloud_ok:
        log(f"[render] not sending anything. --tts {a.tts} sends the transcript -- {words:,} "
            f"words of the turns, after pronunciation substitutions; never the report -- to "
            f"Google's Gemini API. On a free-tier key Google may use it to improve its products "
            f"and human reviewers may read it ({TERMS_URL}). Once the user has agreed, re-run "
            f"with --cloud-ok (optionally --cloud-ok \"who agreed, when\"), or use --tts say "
            f"(offline, macOS).")
        return 2
    out = os.path.abspath(a.out or os.path.dirname(s.path))
    cache = os.path.join(out, ".cache")
    os.makedirs(cache, exist_ok=True)
    backend = backend_cls(s, a)
    backend.prepare(chunks, speak, cache)
    log(f"[render] {len(chunks)} chunk(s), {words:,} words -> about {words / float(WPM):.0f} min "
        f"with {backend.name}")

    pieces, n_cached, n_new = [], 0, 0
    for i, c in enumerate(chunks, 1):
        label = f"chunk {i}/{len(chunks)} (segment {c.seg + 1})"
        while True:
            path = os.path.join(cache, cache_key(backend, c, speak) + ".wav")
            if os.path.isfile(path):
                try:
                    pcm = read_wav(path)
                    n_cached += 1
                    backend.produced += 1
                    if backend.model not in backend.models_used:
                        backend.models_used.append(backend.model)
                    log(f"[render] {label}: cached")
                    break
                except (wave.Error, EOFError, OSError, ValueError):
                    os.remove(path)
            try:
                pcm = backend.synth(c, speak, label)
            except SwitchModel:
                continue
            write_wav(path, pcm)
            n_new += 1
            backend.produced += 1
            log(f"[render] {label}: {len(pcm) / float(RATE):.0f} s from {backend.model}")
            break
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
                                  "one style line; not the report" if backend.cloud else
                                  "nothing (offline)")},
        "cloud_tts_consent": (a.cloud_ok if backend.cloud else False),
        "script": "podcast_script.md", "script_sha256": s.sha256,
        "script_author": s.author,
        "check": {"status": status or "not run", "file": "check.txt",
                  "overridden": bool(why and a.unchecked), "reason": why},
        "sources": sources,
        "words": words, "turns": len(s.turns()), "segments": len(s.segments),
        "chunks": len(chunks), "chunks_cached": n_cached, "chunks_synthesized": n_new,
        "duration_s": duration, "created": now_iso(), "ai_generated": True,
        "audio": audio_name, "encoder": tool, "transcript": "transcript.html",
        "warnings": backend.warnings,
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
                      "warnings": backend.warnings}, indent=2))
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


def load_manifest(outdir):
    try:
        with open(os.path.join(outdir, "podcast", "podcast.json"), encoding="utf-8") as fh:
            man = json.load(fh)
        return man if isinstance(man, dict) and man.get("audio") else None
    except (OSError, ValueError):
        return None


def _paths(outdir, base, man):
    pdir = os.path.join(outdir, "podcast")
    return (rel(os.path.join(pdir, man["audio"]), base),
            rel(os.path.join(pdir, man.get("transcript") or "transcript.html"), base))


def _name(man):
    return f"{man.get('show') or SHOW} — {man.get('title') or 'Untitled'}"


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


def listen_md(outdir, base_dir=None):
    man = load_manifest(outdir)
    if not man:
        return ""
    audio, tr = _paths(outdir, base_dir or outdir, man)
    return "\n".join([START, f"> **Listen:** *{man.get('show') or SHOW}* — "
                             f"{man.get('title') or 'Untitled'} ({minutes_label(man.get('duration_s'))}): "
                             f"[{audio}]({href(audio)}). {DISCLOSURE} Transcript: [{tr}]({href(tr)}).",
                      END])


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


def add_listen_card(doc, outdir):
    """`doc` with the Listen card near the top; unchanged when `outdir` has no podcast. The
    hook make_analysis_html.py calls, so regenerating the report keeps the card."""
    card = listen_card_html(outdir)
    return _upsert(doc, card, _html_slot)[0] if card else doc


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
    man = load_manifest(out)
    if not man:
        log(f"[link] no podcast/podcast.json under {out}: render first, and point link at the "
            "folder that holds Analysis_Report.html and podcast/")
        return 2
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
    r.add_argument("--min-interval", type=float, default=MIN_INTERVAL,
                   help="seconds between requests to one Gemini model (default %(default)s)")
    r.add_argument("--keep-wav", action="store_true")
    r.add_argument("--unchecked", action="store_true",
                   help="render without a passing check.txt (recorded in podcast.json)")
    r.add_argument("--dry-run", action="store_true",
                   help="print the chunks as they would be spoken; send nothing")
    k = sub.add_parser("link", help="link the podcast into the report, README and AGENTS.md")
    k.add_argument("outdir", help="the folder with Analysis_Report.html and podcast/")
    a = ap.parse_args(argv)
    if not a.cmd:
        ap.print_help()
        return 2
    try:
        return {"check": cmd_check, "render": cmd_render, "link": cmd_link}[a.cmd](a)
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
