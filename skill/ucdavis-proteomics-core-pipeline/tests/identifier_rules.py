"""The ONE definition of a real identifier in the published skill, for the identifier guards.

The repository is public. What may not appear in its text: a real CoreOmics submission number,
a Core staff login or name, and the name of a client -- a PI, a lab, a company, or anyone else
outside the Core. What may appear instead:

* submission numbers: PROT_0807 (the one public example, also as 807 / PROT_807), PROT_08070
  (the malformed test value), PROT_0000 (the placeholder) and the reserved fake ranges
  PROT_9000-PROT_9999 and PROT_90000-PROT_99999, in any spelling (PROT_, PROT-, PROT, prot_ ...);
  "submission <n>" follows the same rule. In docs, say "a submission" or "an earlier submission"
  when a second one is meant;
* people: `<your HIVE user>` in docs, generic logins in tests (bstaff, mstaff, gstaff, staff1,
  tester, ...), "a Core member" in comments;
* clients: "a client lab", "the bat-plate study", "a coral study" -- the study, never the person.

`problems(rel, line)` is copied unchanged from the 2.15 slim-down's guard
(tests/test_no_real_identifiers.py on refactor/skill-215-slim), so both guards read one rule.
`client_names(line)` adds the client names that guard does not know. They are NOT in this
repository in any form: a name split into string parts still spells it out, and a digest of a
name can be reversed with a word list in seconds (2.14.2 C2 security review, B1: 7 of 32 salted
digests fell to the stock macOS dictionary), and most of these names were never public. The list
is a private file, one name per line ('#' starts a comment), read from $SKILL_PRIVATE_NAMES, else
~/.config/ucdavis-proteomics/client_names.txt (mode 600; the Core keeps a group copy on HIVE).
Without it the client-name part is SKIPPED, and says so (client_list); the submission-number and
staff parts still run. A hit names the list's entry number, never the name, so a test log does
not spell it out either. To add a name, add a line to the private file -- never to this one. The
other patterns are built from parts so this file does not match itself.
"""
import hashlib
import os
import re
import unicodedata

TEXT = (".py", ".md", ".R", ".sh", ".json", ".txt", ".yaml", ".yml", ".tsv", ".csv", ".html",
        ".cfg", ".toml", ".js", ".css")

# ---- copied from the slim-down's tests/test_no_real_identifiers.py (keep the two identical) ----
# A submission number written with its PROT prefix (any case, separator or none), and one
# introduced by the word "submission" -- the two forms core_submission.py's identify reads.
PROT = re.compile(r"(?<![A-Za-z])" + "pr" + r"ot[ _#-]?([0-9]{3,5})(?![0-9])", re.I)
SUBMISSION = re.compile(r"\b" + "submis" + r"sion\s*(?:number|no\.?)?\s*#?\s*([0-9]{3,5})(?![0-9])", re.I)
FAKE_NUMBER = re.compile(r"0807|807|08070|0000|9[0-9]{3}|9[0-9]{4}")
# The podcast fixtures' study is synthetic and says so ("no such study"); its number stays.
SYNTHETIC = {(os.path.join("tests", "fixtures", "podcast"), "04" + "21")}

# Staff logins and names, the client company, a client lab and dataset, a person outside the
# Core, and the real CoreOmics id of a client submission.
NAMES = re.compile("|".join((
    "brett" + "sp", "msal" + "emi", "sal" + "emi", r"(?<![a-z])gab" + "rig", "gabri" + "ela",
    "mich" + "elle", "frag" + "matics", "dick" + "son", "sil" + r"va(?=[0-9_])", "adam " + "schaal",
    "0022066" + "cd85f", "dupan" + "loup", r"\bga" + r"rg\b")), re.I)

# (path under the skill, the line exactly, why it must stay)
ALLOWED = (
    # core_admins.txt is the trust anchor: skill_version.py/.sh, notes.py, slack_collab.py and
    # make_podcast.py trust only files owned by a HIVE account listed here, so the real admin
    # login is data the code reads (tests pin it by hash, not by name).
    (os.path.join("scripts", "core_admins.txt"), "brett" + "sp"),
)


def _number_ok(rel, digits):
    if FAKE_NUMBER.fullmatch(digits):
        return True
    return any(rel.startswith(folder + os.sep) and digits == n for folder, n in SYNTHETIC)


def problems(rel, line):
    """What is wrong with one line of the file at `rel` (a path under the skill): [] when none."""
    found = []
    for rx in (PROT, SUBMISSION):
        for m in rx.finditer(line):
            if not _number_ok(rel, m.group(1)):
                found.append(f"submission number {m.group(0)!r}")
    if NAMES.search(line) and (rel, line.strip()) not in ALLOWED:
        found.append(f"a real name or login {NAMES.search(line).group(0)!r}")
    return found
# ---- end of the copy ----


# Client and outside-person names the slim-down's NAMES does not hold (2.14.2, the C2 scrub):
# the PIs, labs, companies and outside people named in the Core's issue reports, session names
# and project notes. They live ONLY in a private file outside the repository (C2 security B1).
# Words that are also ordinary words or the names of cited authors (a grant acknowledgment, a
# published contaminant library) are left out of that list on purpose: they cannot be told apart
# from a client by a pattern, and were checked by hand instead.
PRIVATE_NAMES_ENV = "SKILL_PRIVATE_NAMES"
PRIVATE_NAMES_DEFAULT = os.path.join("~", ".config", "ucdavis-proteomics", "client_names.txt")
SKIPPED = "client-name check SKIPPED"


def private_list_path():
    """Where the private client-name list is read from: $SKILL_PRIVATE_NAMES, else the default."""
    return os.path.expanduser(os.environ.get(PRIVATE_NAMES_ENV) or PRIVATE_NAMES_DEFAULT)


def read_private_list(path):
    """The private list at `path` -> {key: entry number}, a key being the name's words folded and
    joined ("Van Dam" -> "vandam"). Raises OSError when it cannot be read and ValueError for an
    entry with no letters or a list with no names: a list that is there but unusable is an error,
    never a silent skip."""
    keys = {}
    with open(path, encoding="utf-8") as fh:
        for n, raw in enumerate(fh, 1):
            entry = raw.split("#", 1)[0].strip()
            if not entry:
                continue
            words = _words(entry)
            if not words:
                raise ValueError(f"{path}: entry on line {n} holds no letters")
            if len(words) > NAME_WORDS_MAX:
                raise ValueError(f"{path}: entry on line {n} has more than {NAME_WORDS_MAX} words")
            keys.setdefault("".join(words), n)
    if not keys:
        raise ValueError(f"{path}: the private client-name list holds no names")
    return keys


_LISTS = {}


def client_list():
    """-> (keys or None, message). None when the private list is absent: the client-name part is
    then skipped, and the message says so and where to put the list. Cached per path and mtime."""
    path = private_list_path()
    try:
        st = os.stat(path)
    except FileNotFoundError:
        return None, (f"{SKIPPED}: no private list at {path} (set ${PRIVATE_NAMES_ENV}, or put "
                      f"the Core's list there, one name per line); the submission-number and "
                      f"staff checks still ran")
    stamp = (path, st.st_mtime_ns, st.st_size)
    if stamp not in _LISTS:
        _LISTS.clear()
        _LISTS[stamp] = read_private_list(path)
    keys = _LISTS[stamp]
    return keys, f"client-name check: {len(keys)} names from {path}"


def _unaccent(text):
    """Accents dropped, case kept ('Müller' -> 'Muller'), so a word is split as it is written."""
    text = unicodedata.normalize("NFKD", text)
    return "".join(c for c in text if not unicodedata.combining(c))


# One word: a run of capitals (an acronym, or a capital before the next word's capital), a
# capitalised or lower-case run, or digits -- so "SmithBat", "Smith08172026" and "VanDam"
# split into words (a name of several words is also tried joined, "vandam"), and a name inside
# an ordinary word is not a hit.
_WORD = re.compile(r"[A-Z]+(?![a-z])|[A-Z]?[a-z]+|[0-9]+")


def _words(text):
    """The lower-case, unaccented words of `text` that are letters only."""
    return [w for w in (x.lower() for x in _WORD.findall(_unaccent(text))) if w.isalpha()]


# A cited author who shares a surname with a client (the QuantUMS paper, references/de-analysis.md).
_CITED = re.compile("sil" + r"va campitelli", re.I)
# The longest name tried, in words ("Van Dam" is two); read_private_list holds every entry to it.
NAME_WORDS_MAX = 4


def client_names(line, names=None):
    """The private list's entry numbers found on one line: [] when none, or when the list is
    absent (client_list says so). `names`: a list from read_private_list, else client_list()'s."""
    if names is None:
        names, _why = client_list()
    if not names:
        return []
    words = _words(_CITED.sub(" ", line))
    found = []
    for i in range(len(words)):
        for k in range(1, NAME_WORDS_MAX + 1):
            n = names.get("".join(words[i:i + k])) if i + k <= len(words) else None
            if n and n not in found:
                found.append(n)
    return found


# Public text the slim-down's NAMES matches that names no client: (file, phrase, why). Each phrase
# is taken out of that file's lines before problems() reads them, so anything else on the line is
# still checked, and the phrase anywhere else is still a hit. (The slim-down's guard needs the
# same allowance when 2.15 lands: these lines are on main.)
_LICENSOR = "frag" + "matics"
_LICENSOR_WHY = ("the commercial licensor of MSFragger, IonQuant and diaTracer, named with its "
                 "contact in the tools' public licence text: a user who asks whether their work "
                 "needs a commercial licence is sent there (C2 docs N2, Ruling 82)")
PUBLIC_TEXT = tuple((rel, phrase, _LICENSOR_WHY) for rel, phrase in (
    ("SKILL.md", f"via {_LICENSOR}.com"),
    ("SKILL.md", f"info@{_LICENSOR}.com"),
    (os.path.join("references", "fragpipe-diatracer.md"), f"via {_LICENSOR}.com"),
    (os.path.join("references", "search-engines.md"),
     f"via **{_LICENSOR.title()}** ({_LICENSOR}.com, info@{_LICENSOR}.com)"),
    (os.path.join("references", "search-engines.md"), f"via {_LICENSOR.title()}."),
    (os.path.join("references", "search-engines.md"), f"with {_LICENSOR.title()}."),
    (os.path.join("references", "search-engines.md"), f"info@{_LICENSOR}.com"),
))


def _public_text_out(rel, line):
    """`line` with the public phrases allowed in `rel` taken out (PUBLIC_TEXT), longest first."""
    for phrase in sorted((p for r, p, _why in PUBLIC_TEXT if r == rel), key=len, reverse=True):
        line = line.replace(phrase, " ")
    return line


def line_hits(rel, line, names=None):
    """Everything wrong with one line: the slim-down's problems() (after the public text allowed
    in that file, PUBLIC_TEXT) plus the client names, each named by its entry number in the
    private list -- never by the name itself."""
    return problems(rel, _public_text_out(rel, line)) + [
        f"a client name (private list entry {n})" for n in client_names(line, names)]


def baseline_key(rel, text):
    """THE key of a hit in tests/identifier_baseline.json: the sha256 (hex) of "<path>\\0<line>",
    the path relative to the skill with "/" separators and the line stripped (for a file or
    folder name, the name). The baseline holds only these keys and counts, never the text: a
    list of every already-public staff and submission line in one new public file would
    republish them (2.14.2, before the push). A key names a line already in the tree; it is not
    a digest of a name, and client names are never in the baseline."""
    return hashlib.sha256((rel.replace(os.sep, "/") + "\0" + text.strip()).encode("utf-8")).hexdigest()


def walk(skill, skip=()):
    """Yield (rel, line_number, line) for every text line under `skill`, and (rel, 0, name) for
    every file and folder name. `skip` holds paths (relative to `skill`) left out."""
    for dp, dns, fns in os.walk(skill):
        dns[:] = sorted(d for d in dns if d not in ("__pycache__", ".pytest_cache"))
        for name in dns + sorted(fns):
            rel = os.path.relpath(os.path.join(dp, name), skill)
            if rel not in skip:
                yield rel, 0, name
        for fn in sorted(fns):
            path = os.path.join(dp, fn)
            rel = os.path.relpath(path, skill)
            if rel in skip:
                continue
            with open(path, "rb") as fh:
                data = fh.read()
            if not fn.endswith(TEXT) and b"\0" in data[:8192]:
                continue                          # binary; any text file is read, named or not
            for i, line in enumerate(data.decode("utf-8", "replace").splitlines(), 1):
                yield rel, i, line
