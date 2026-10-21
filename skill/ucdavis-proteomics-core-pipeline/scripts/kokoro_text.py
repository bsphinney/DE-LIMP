#!/usr/bin/env python3
"""
kokoro_text.py -- the text rules for the local podcast voice (Kokoro-82M on HIVE), applied to a
turn AFTER the script's Pronunciation table, for `make_podcast.py render --tts local` and the
automatic fallback when Gemini cannot be used. They change only what the voice reads; the
transcript keeps the written form.

    from kokoro_text import apply
    text, changes = apply("the subunit k c n b two, analysed on 2026-10-14")
    # -> "the subunit K C N B two, analysed on October 14, 2026"

    python3 kokoro_text.py "PROT_0807 on 2026-10-14"     # what the voice would read

Each rule fixes a misreading found in Core episodes rendered with Kokoro (2026-10), each checked
by two Whisper transcriptions that agreed; ordinary sentences pass through unchanged:

1. spelled-out letters  "c s g a" -> "C S G A"   (a lone lowercase "a" was read as the word "uh");
                        a leading "a" stays the article: "for a r y r effect" -> "for a R Y R effect"
2. ISO dates            "2026-10-14" -> "October 14, 2026"   (was "twenty twenty sixteen fourteen")
3. version numbers      "version 2.13.0" -> "version 2 point 13 point 0"   (was "two thirteen zero")
4. PROT numbers         "PROT_0807" -> "P R O T zero eight zero seven"   (a leading zero was dropped)
5. signs                "+0.61" -> "plus 0.61" (was "plus zero sixty-one"); a minus sign before a
                        number, ASCII or Unicode, -> "minus"
6. number ranges        "91–97%" -> "91 to 97%" (an en or em dash between numbers was dropped or
                        garbled); a hyphen between two plain numbers too ("68-72%" -> "68 to 72%"),
                        never a year-month (2026-10), a chain (1-2-3) or an ID (X1234-C1, Lys8-labelled)

Standard library only.
"""
import re
import sys

MONTHS = ("January February March April May June July August September October November "
          "December").split()
DIGITS = "zero one two three four five six seven eight nine".split()

# 4. PROT_0807, PROT 0807, PROT-0807 (any case); the digits are read one by one, zeros kept
_PROT = re.compile(r"(?<![A-Za-z0-9])PROT[_ -]?(\d+)(?![0-9])", re.I)
# 2. a full ISO date standing alone (not inside a longer name such as 2026-10-15_PROT_0807)
_DATE = re.compile(r"(?<![\w-])(\d{4})-(\d{2})-(\d{2})(?![\w-])")
# 3. three or more dotted parts anywhere (2.13.0); two parts only after "version" or "v" (2.13 is
#    otherwise an ordinary decimal and stays as it is)
_VER3 = re.compile(r"(?<![\w.])[vV]?(\d+(?:\.\d+){2,})(?![\w]|\.\d)")
_VER2 = re.compile(r"\b(?:([Vv]ersion)\s+|[vV](?=\d))(\d+)\.(\d+)(?![\w]|\.\d)")
# 5. a plus sign before a number, not inside an expression (1+1, x+1, C++)
_PLUS = re.compile(r"(?<![\w.)+])\+(?=\d|\.\d)")
# 5b. a minus sign before a number (ASCII or Unicode), not inside a word or an expression
_MINUS = re.compile(r"(?<![\w.)+−-])[−-](?=\d|\.\d)")
# 6. an en or em dash between two numbers
_RANGE = re.compile(r"(?<=[\d%])\s*[–—]\s*(?=[\d+−-])")
# 6b. a hyphen between two plain numbers ("68-72%"), not in a chain (1-2-3) or an ID (X1234-C1)
_HRANGE = re.compile(r"(?<![\w.,%-])(\d[\d,]*(?:\.\d+)?%?)-(\d[\d,]*(?:\.\d+)?)(?![\w,-]|\.\d)")
# 1. two or more single letters in a row, separated by spaces ("c s g a", "delta h n s"); a lone
#    letter is left alone, so the article "a" and the pronoun "I" never change
_LETTERS = re.compile(r"(?<![\w'’.-])[A-Za-z](?: [A-Za-z])+(?![\w'’-])")


def _hrange(m):
    a, b = m.group(1), m.group(2)
    if re.fullmatch(r"(?:19|20)\d\d", a) and re.fullmatch(r"0[1-9]|1[0-2]", b):
        return m.group(0)                        # a year-month (2026-10), not a range
    return f"{a} to {b}"


def _prot(m):
    return "P R O T " + " ".join(DIGITS[int(d)] for d in m.group(1))


def _date(m):
    y, mo, d = int(m.group(1)), int(m.group(2)), int(m.group(3))
    if not (1900 <= y <= 2100 and 1 <= mo <= 12 and 1 <= d <= 31):
        return m.group(0)
    return f"{MONTHS[mo - 1]} {d}, {y}"


def _ver3(m):
    return " point ".join(str(int(p)) if p.isdigit() else p for p in m.group(1).split("."))


def _ver2(m):
    word = m.group(1) or "version"
    return f"{word} {m.group(2)} point {m.group(3)}"


def _letters(m):
    s = m.group(0)
    toks = s.split(" ")
    if toks[0] in ("a", "i"):
        # a leading lowercase a/i is the word ("looking for a r y r effect", "a T cell"): keep
        # it, and capitalise only the letters after it when two or more remain. A symbol that
        # really starts with a ("AKT1" spelled "a k t one") needs fix_table on the Pronunciation
        # table, where there is no article to confuse it with.
        rest = toks[1:]
        if len(rest) >= 2 and any(c.islower() for c in "".join(rest)):
            return toks[0] + " " + " ".join(rest).upper()
        return s
    return s.upper() if any(c.islower() for c in s) else s


RULES = (("prot", _PROT, _prot), ("date", _DATE, _date), ("version", _VER3, _ver3),
         ("version", _VER2, _ver2), ("range", _RANGE, lambda m: " to "),
         ("range", _HRANGE, _hrange), ("plus", _PLUS, lambda m: "plus "),
         ("minus", _MINUS, lambda m: "minus "), ("letters", _LETTERS, _letters))


def apply(text):
    """-> (new_text, [(rule, before, after), ...]) with every rule applied in order."""
    changes = []
    for name, rx, fn in RULES:
        def sub(m, name=name, fn=fn):
            new = fn(m)
            if new != m.group(0):
                changes.append((name, m.group(0), new))
            return new
        text = rx.sub(sub, text)
    return text, changes


def fix_table(rows):
    """The letters rule for a Pronunciation table's spoken column, where there is no article to
    confuse a leading "a" with: [(written, spoken)] -> same shape."""
    return [(w, _LETTERS.sub(lambda m: m.group(0).upper(), sp)) for w, sp in rows]


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print(__doc__.strip())
        sys.exit(0 if len(sys.argv) > 1 else 2)
    new, changed = apply(" ".join(sys.argv[1:]))
    print(new)
    for rule, before, after in changed:
        print(f"  {rule}: {before!r} -> {after!r}", file=sys.stderr)
