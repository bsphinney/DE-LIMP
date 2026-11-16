"""The ONE definition of a real identifier in the published skill, for the identifier guard
(tests/test_no_real_identifiers.py) and the export's check (export_new_repo.py).

The repository is public. What may not appear in its text: a real CoreOmics submission number,
a Core staff login or name, and the name of a client -- a PI, a lab, a company, or anyone else
outside the Core. What may appear instead:

* submission numbers: PROT_0807 (the one public example, also as 807 / PROT_807), PROT_08070
  (the malformed test value), PROT_0000 (the placeholder) and the reserved fake ranges
  PROT_9000-PROT_9999 and PROT_90000-PROT_99999, in any spelling (PROT_, PROT-, PROT, PROT #,
  PROT., PROT–, prot_ ...); "submission <n>" and "submission: <n>" follow the same rule, and so
  does a bare 4-digit number with a leading 0 and at most NEAR_WORDS word between it and
  "submission(s)" or "CoreOmics" ("the 0807 submission", "submission ID 0807"). In docs, say "a submission" or "an earlier submission"
  when a second one is meant;
* people: `<your HIVE user>` in docs, generic logins in tests (bstaff, mstaff, gstaff, staff1,
  tester, ...), "a Core member" in comments;
* clients: "a client lab", "the bat-plate study", "a coral study" -- the study, never the person.

`problems(rel, line)` is the submission-number and fixed-name rule (NAMES: names already public
in the history, built from parts). Two PRIVATE lists add what may never be in this repository in
any form -- a name split into string parts still spells it out, and a digest of a name can be
reversed with a word list in seconds (2.14.2 C2 security review, B1: 7 of 32 salted digests fell
to the stock macOS dictionary):

* the client names (client_names): $SKILL_PRIVATE_NAMES, else
  ~/.config/ucdavis-proteomics/client_names.txt -- one name per line, '#' starts a comment;
* the Core staff (staff_names): $SKILL_PRIVATE_STAFF, else
  ~/.config/ucdavis-proteomics/core_staff.txt -- the format of the Core's staff file that
  scripts/staff.py reads on HIVE (/quobyte/proteomics-grp/.config/core_staff.txt): one HIVE login
  per line, '#' starts a comment; after the login, optionally, the names the person goes by,
  separated by commas ("jdoe  Jane, Jane Doe"), which staff.py ignores.

Both live on the Core director's Mac only (Ruling 93: no GitHub secret, no copy on HIVE). Without
one, its part is SKIPPED and says so (client_list / staff_list); the submission-number and fixed-
name parts still run. With SKILL_REQUIRE_CLIENT_NAMES=1 (REQUIRE_ENV) a missing list FAILS the
guard instead -- the controller sets it on the Mac at every merge into the slim or release branch
and before the release push. A hit names the list's entry number, never the name, so a test log
does not spell it out either. To add a name, add a line to the private file -- never to this one.
The other patterns are built from parts so this file does not match itself.
"""
import os
import re
import unicodedata

TEXT = (".py", ".md", ".R", ".sh", ".json", ".txt", ".yaml", ".yml", ".tsv", ".csv", ".html",
        ".cfg", ".toml", ".js", ".css")

# ---- the submission-number and fixed-name rule (test_no_real_identifiers imports these) ----
# A submission number written with its PROT prefix (any case; a separator or none: _ - # . : a
# space, an en or em dash, up to three of them -- "PROT #0807", "PROT.0807", "PROT–0807"), and
# one introduced by the word "submission" ("submission 807", "submission: 807", "submission ID
# 807") -- the forms core_submission.py's identify reads, and the ones people write.
PROT = re.compile(r"(?<![A-Za-z])" + "pr" + r"ot[ _#.:\u2010-\u2015-]{0,3}([0-9]{3,5})(?![0-9])", re.I)
SUBMISSION = re.compile(r"\b" + "submis" + r"sion\s*(?:number|no\.?|id)?\s*[:#]?\s*#?\s*([0-9]{3,5})(?![0-9])",
                        re.I)
FAKE_NUMBER = re.compile(r"0807|807|08070|0000|9[0-9]{3}|9[0-9]{4}")
# A bare number beside a submission word ("the 0807 submission", "CoreOmics 0807", "submission ID
# 0807", "submissions 0807 and 9756"): 4 digits with a leading 0 -- how the Core writes one --
# with at most NEAR_WORDS words between it and "submission(s)" / "CoreOmics" (other numbers of a
# list in between do not count). Measured on the 2.15 tree: any 4-digit number within 3 tokens
# hit 32 lines (years, dates, paths, a count), none of them a submission; 19xx/20xx left out, 4;
# a leading 0, 1 (a umask: two words away); a leading 0 with at most one word between, none.
BARE_NUMBER = re.compile(r"0[0-9]{3}")
SUBMISSION_WORDS = frozenset(("submis" + "sion", "submis" + "sions", "core" + "omics"))
NEAR_WORDS = 1
_TOKEN = re.compile(r"[A-Za-z]+|[0-9]+")
# The podcast fixtures' study is synthetic and says so ("no such study"); its number stays.
SYNTHETIC = {(os.path.join("tests", "fixtures", "podcast"), "04" + "21")}

# Staff logins and names, the client company, a client lab and dataset, a person outside the
# Core, and the real CoreOmics id of a client submission -- already public in the history; the
# rest of the Core's staff are in the private staff list (staff_names), never here.
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


def near_numbers(line):
    """The bare numbers (BARE_NUMBER) on `line` with at most NEAR_WORDS words between them and a
    submission word; numbers in between (a list of them) are not counted as words."""
    toks = _TOKEN.findall(line)
    found = []

    def reaches(i, step):
        words, j = 0, i + step
        while 0 <= j < len(toks):
            t = toks[j]
            if t.lower() in SUBMISSION_WORDS:
                return True
            if not t.isdigit():
                words += 1
                if words > NEAR_WORDS:
                    return False
            j += step
        return False
    for i, t in enumerate(toks):
        if BARE_NUMBER.fullmatch(t) and (reaches(i, -1) or reaches(i, 1)):
            found.append(t)
    return found


def problems(rel, line):
    """What is wrong with one line of the file at `rel` (a path under the skill): [] when none."""
    found, said = [], set()
    for rx in (PROT, SUBMISSION):
        for m in rx.finditer(line):
            if not _number_ok(rel, m.group(1)):
                found.append(f"submission number {m.group(0)!r}")
                said.add(m.group(1).lstrip("0"))
    for n in near_numbers(line):
        if not _number_ok(rel, n) and n.lstrip("0") not in said:
            found.append(f"submission number {n!r} beside a submission word")
            said.add(n.lstrip("0"))
    if NAMES.search(line) and (rel, line.strip()) not in ALLOWED:
        found.append(f"a real name or login {NAMES.search(line).group(0)!r}")
    return found
# ---- end of the rule ----


# Client and outside-person names the slim-down's NAMES does not hold (2.14.2, the C2 scrub):
# the PIs, labs, companies and outside people named in the Core's issue reports, session names
# and project notes. They live ONLY in a private file outside the repository (C2 security B1).
# Words that are also ordinary words or the names of cited authors (a grant acknowledgment, a
# published contaminant library) are left out of that list on purpose: they cannot be told apart
# from a client by a pattern, and were checked by hand instead.
PRIVATE_NAMES_ENV = "SKILL_PRIVATE_NAMES"
PRIVATE_NAMES_DEFAULT = os.path.join("~", ".config", "ucdavis-proteomics", "client_names.txt")
SKIPPED = "client-name check SKIPPED"
# The Core's staff list (scripts/staff.py's format), private like the client names (Ruling 93).
PRIVATE_STAFF_ENV = "SKILL_PRIVATE_STAFF"
PRIVATE_STAFF_DEFAULT = os.path.join("~", ".config", "ucdavis-proteomics", "core_staff.txt")
STAFF_SKIPPED = "staff-name check SKIPPED"
# "1": a private list that is not there FAILS the guard instead of skipping it -- the controller
# sets it on the Mac at every merge into the slim or release branch and before the release push
# (Ruling 93). CI and HIVE never have the lists, and never set it.
REQUIRE_ENV = "SKILL_REQUIRE_CLIENT_NAMES"


def required():
    """Must every private list be there (REQUIRE_ENV is "1")?"""
    return os.environ.get(REQUIRE_ENV, "").strip() == "1"


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


def private_staff_path():
    """Where the private staff list is read from: $SKILL_PRIVATE_STAFF, else the default."""
    return os.path.expanduser(os.environ.get(PRIVATE_STAFF_ENV) or PRIVATE_STAFF_DEFAULT)


_STAFF_LOGIN = re.compile(r"^[a-z_][a-z0-9_.-]{0,31}$")    # a HIVE username (staff.py's shape)


def read_staff_list(path):
    """The private staff list at `path` -> {"logins": {login: entry number}, "names": {key: entry
    number}}: the login that starts each line, and each comma-separated name after it (a key as
    in read_private_list). Raises OSError when it cannot be read and ValueError for a line whose
    first word is not a login, a name with no letters, or a list with no logins."""
    logins, names = {}, {}
    with open(path, encoding="utf-8") as fh:
        for n, raw in enumerate(fh, 1):
            entry = raw.split("#", 1)[0].strip()
            if not entry:
                continue
            first, rest = (entry.split(None, 1) + [""])[:2]       # as staff.py splits a line
            login = first.lower()
            if not _STAFF_LOGIN.match(login):
                raise ValueError(f"{path}: line {n} does not start with a HIVE login")
            logins.setdefault(login, n)
            for name in (x.strip() for x in rest.split(",")):
                if not name:
                    continue
                words = _words(name)
                if not words:
                    raise ValueError(f"{path}: a name on line {n} holds no letters")
                if len(words) > NAME_WORDS_MAX:
                    raise ValueError(f"{path}: a name on line {n} has more than {NAME_WORDS_MAX} "
                                     "words")
                names.setdefault("".join(words), n)
    if not logins:
        raise ValueError(f"{path}: the private staff list holds no logins")
    return {"logins": logins, "names": names}


_STAFF = {}


def staff_list():
    """-> (lists or None, message), as client_list(): None when the private staff list is absent
    (the staff-name part is then skipped, and the message says so). Cached per path and mtime."""
    path = private_staff_path()
    try:
        st = os.stat(path)
    except FileNotFoundError:
        return None, (f"{STAFF_SKIPPED}: no private staff list at {path} (set "
                      f"${PRIVATE_STAFF_ENV}, or put a copy of the Core's staff file there: one "
                      f"HIVE login per line, then the names the person goes by, comma-separated); "
                      f"the fixed names and the submission-number checks still ran")
    stamp = (path, st.st_mtime_ns, st.st_size)
    if stamp not in _STAFF:
        _STAFF.clear()
        _STAFF[stamp] = read_staff_list(path)
    lists = _STAFF[stamp]
    return lists, (f"staff-name check: {len(lists['logins'])} logins and {len(lists['names'])} "
                   f"names from {path}")


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


def staff_names(rel, line, staff=None):
    """The private staff list's entry numbers found on one line: a login as a whole word (not
    inside a longer login or word: "jdoe" is not in "jdoes"), a name as client_names matches one.
    [] when the list is absent, and for a line ALLOWED (the trust anchor's login)."""
    if staff is None:
        staff, _why = staff_list()
    if not staff or (rel, line.strip()) in ALLOWED:
        return []
    found = []
    low = line.lower()
    for login, n in staff["logins"].items():
        if n not in found and re.search(r"(?<![a-z0-9])" + re.escape(login) + r"(?![a-z0-9])",
                                        low):
            found.append(n)
    for n in client_names(line, staff["names"]) if staff["names"] else []:
        if n not in found:
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
    # 2.15: SKILL.md's "send them to info@..." line (old step 7) lives here in the slim tree
    (os.path.join("references", "fragpipe-diatracer.md"), f"info@{_LICENSOR}.com"),
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


def line_hits(rel, line, names=None, staff=None):
    """Everything wrong with one line: problems() (after the public text allowed in that file,
    PUBLIC_TEXT), the client names and the staff logins and names, each named by its entry
    number in its private list -- never by the name itself. `names` / `staff`: lists as read
    (read_private_list / read_staff_list); None reads the private files; {} = none."""
    return problems(rel, _public_text_out(rel, line)) + [
        f"a client name (private list entry {n})" for n in client_names(line, names)] + [
        f"a Core staff login or name (private staff list entry {n})"
        for n in staff_names(rel, line, staff)]


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


# ---- what the new repository ships outside skill/ (export_new_repo.py; review LOW 15) --------
# The export takes, besides the skill: .claude-plugin/marketplace.json, LICENSE, the CHANGELOG's
# "## [Skill ...]" sections, the CI workflow, the docs pages the skill uses (its two seed pages,
# every page the skill links as a Pages URL or its tests read, and what those link relatively) and
# a README written from a template outside this repository. All are public. The README draft and
# any other file can be added with $SKILL_SCAN_EXTRA (paths, os.pathsep-separated: files or
# folders), as the controller does at a merge; export_new_repo.py runs the same scan on its tree.
SCAN_EXTRA_ENV = "SKILL_SCAN_EXTRA"
SKILL_REL = "skill/ucdavis-proteomics-core-pipeline"
SHIPPED_FILES = (".claude-plugin/marketplace.json", "LICENSE", ".github/workflows/skill-checks.yml")
DOCS_SEED = ("skill-install.html", "how-reports-are-made.html")      # as export_new_repo.py
# the skill's Pages (from 2.15) and DE-LIMP's (reports made before it), as export_new_repo.py
_PAGES = re.compile(r"bsphinney\.github\.io/(?:DE-LIMP|UCD-Proteomics-AI-Skill)/([A-Za-z0-9_./-]+)")
_TEST_DOC = re.compile(r"""["']docs/([A-Za-z0-9_./-]+)["']|["']docs["']\s*,\s*["']([A-Za-z0-9_./-]+)["']""")
_LINK = re.compile(r"""(?:href|src)\s*=\s*["']([^"'#?]+)|\]\(([^)#?\s]+)""")
# A data: URI's payload (an embedded image, base64 or percent-encoded) is not text: its letters
# split into "words" that can spell anything. It is taken out before a line is checked.
DATA_URI = re.compile(r"""data:[^,\s"'()]*,[^\s"'()<>]*""")
# The public Hugging Face Space's URL holds the owner's public username; it stays (Ruling 89b).
OUTSIDE_PUBLIC = ("huggingface.co/spaces/" + "brett" + "sp" + "/",)


def _read_text(path):
    """A file's text, or None for a binary file (a NUL in its first 8 KiB)."""
    with open(path, "rb") as fh:
        data = fh.read()
    return None if b"\0" in data[:8192] else data.decode("utf-8", "replace")


def shipped_outside_skill(repo, skill=None):
    """The public files outside the skill that the new repository ships -> [(rel, path, lines)],
    `lines` the line numbers to check (None: all of them; the CHANGELOG: its skill sections)."""
    skill = skill or os.path.join(repo, SKILL_REL)
    out = [(rel, os.path.join(repo, rel), None) for rel in SHIPPED_FILES
           if os.path.isfile(os.path.join(repo, rel))]
    cl = os.path.join(repo, "CHANGELOG.md")
    if os.path.isfile(cl):
        keep, on = set(), False
        with open(cl, encoding="utf-8", errors="replace") as fh:
            for i, line in enumerate(fh, 1):
                if line.startswith("## "):
                    on = line.startswith("## [Skill")
                if on:
                    keep.add(i)
        out.append(("CHANGELOG.md", cl, keep))
    docs = os.path.join(repo, "docs")
    take = [d for d in DOCS_SEED if os.path.isfile(os.path.join(docs, d))]

    def add(rel):
        rel = os.path.normpath(rel)
        if rel not in take and not rel.startswith("..") and os.path.isfile(os.path.join(docs, rel)):
            take.append(rel)
    for rel, _i, line in walk(skill):
        for m in _PAGES.finditer(line):
            add(m.group(1).rstrip("."))
        if rel.startswith("tests" + os.sep):
            for m in _TEST_DOC.finditer(line):
                add(m.group(1) or m.group(2))
    i = 0
    while i < len(take):
        rel = take[i]
        i += 1
        if rel.endswith((".html", ".md")):
            text = _read_text(os.path.join(docs, rel)) or ""
            for m in _LINK.finditer(text):
                link = m.group(1) or m.group(2)
                if "://" not in link and not link.startswith(("mailto:", "data:", "/")):
                    add(os.path.join(os.path.dirname(rel), link))
    out += [(os.path.join("docs", d), os.path.join(docs, d), None) for d in take]
    for extra in filter(None, os.environ.get(SCAN_EXTRA_ENV, "").split(os.pathsep)):
        extra = os.path.expanduser(extra)
        paths = [extra] if os.path.isfile(extra) else sorted(
            os.path.join(dp, f) for dp, _d, fns in os.walk(extra) for f in fns)
        out += [(os.path.abspath(x), x, None) for x in paths]
    return out


def outside_hits(rel, line, names=None, staff=None):
    """line_hits() for a file outside the skill: its data: URIs and the public Space URL taken
    out first."""
    line = DATA_URI.sub(" ", line)
    for phrase in OUTSIDE_PUBLIC:
        line = line.replace(phrase, " ")
    return line_hits(rel, line, names, staff)
