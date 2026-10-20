#!/usr/bin/env python3
"""
podcast_local.py -- the podcast's local voice and its check, on HIVE: Kokoro-82M (Apache-2.0)
speaks the episode and Whisper (MIT) hears it back, each in the Core's pinned container, so
nothing leaves campus and no consent is needed. make_podcast.py uses it for `render --tts local`,
and by itself when Gemini fails in a way waiting cannot fix (402 prepay credits, the key refused,
no quota on any model, every model stalled): then EVERY chunk is made again with Kokoro, so one
episode never mixes voices (references/podcast.md, "The local voice").

Where the containers run (find_runner):
  * a HIVE compute node (a job, or a step of one): here, `apptainer exec`, sized from the job --
    Kokoro workers of THREADS_PER_WORKER threads and about GB_PER_WORKER GB each;
  * a HIVE login node: each container call on a compute node, in ONE blocking srun with the queue
    written out (on_compute.srun_argv: HPC@UCD's ucd-hpc rules), its files beside the episode's
    cache on group storage;
  * another computer with the skill's HIVE login (record_run.hive_login): the lines (or the
    episode) go to a private folder in your HIVE home, a HIVE job (hive_exec.sh --job) speaks or
    hears them, the results come back, and the folder is removed;
  * nowhere else: Unavailable, saying why (no HIVE login here; this account cannot read the
    Core's containers, which are for members of proteomics-grp; the skill is not on HIVE).
Each image is pinned by its sha256 (KOKORO_SHA256, WHISPER_SHA256), checked in the job before it
runs: another file at that path is refused, never used.

    python3 podcast_local.py where          # which of these applies here, and why (JSON)

Standard library only.
"""
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import uuid

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

# The Core's containers on HIVE (/quobyte/proteomics-grp/tts, README.md in each folder). A new
# image is a NEW versioned file with its own pin here; the old ones are kept.
TTS_ROOT = "/quobyte/proteomics-grp/tts"
KOKORO_IMAGE = TTS_ROOT + "/kokoro/kokoro-1.2.sif"
KOKORO_SHA256 = "9877b44195fbffd81327d85744b68c49807318e19869a459e72cda5d8c16b355"
KOKORO_SYNTH = "1.2"            # /opt/kokoro/synth.py: a dead worker fails the render at once
KOKORO_MODEL = "hexgrad/Kokoro-82M@f3ff3571791e39611d31c381e3a41a3af07b4987"
KOKORO_LABEL = "Kokoro-82M (on HIVE)"
WHISPER_IMAGE = TTS_ROOT + "/whisper/whisper-1.0.sif"
WHISPER_SHA256 = "778806d931e59b6322dd5cee14f0bf14245959a54c57b751d26eb42497aae349"
# the main check, then a different model as the second opinion: a spot counts only where both
# transcriptions show it (Whisper sometimes invents or skips words)
WHISPER_MODELS = ("large-v3-turbo", "small.en")

# Voices chosen by Brett Phinney (Core director), 2026-10-20, from eight stock pairs: Emma, then
# Heart. Both are women's voices, so a second host named Leo is SPOKEN as Lena in a Kokoro
# render; the script, the transcript and podcast.json keep the name Leo (spoken_names).
VOICES = ("bf_emma", "af_heart")                       # first host, second host
VOICES_IN_IMAGE = ("af_heart", "af_bella", "am_michael", "am_fenrir", "am_puck", "bf_emma",
                   "bm_george", "bm_fable")
WOMEN = re.compile(r"^[ab]f_")                         # Kokoro's names: a/b = US/UK, f/m
SPOKEN_NAMES = {"Leo": "Lena"}                         # a man's name, in a woman's voice

# Kokoro on CPUs (measured 2026-10-19, a 24-min episode): 8 workers x 2 threads on 16 CPUs speak
# it in 2.1 min. A fixed 2 threads per worker keeps the audio byte-identical whatever the number
# of workers. Each worker holds its own model, about 2 GB.
THREADS_PER_WORKER = 2
GB_PER_WORKER = 2.5
RESERVE_GB = 3.0                # the main process and the container's start-up
MIN_GB = GB_PER_WORKER + RESERVE_GB
# A container call sent to a compute node (from a login node, or from another computer)
JOB_CPUS, JOB_MEM = 16, "24G"
SPEAK_MINUTES, HEAR_MINUTES = 30, 45
REMOTE_DIR = "~/.proteomics-pipeline/podcast-local"    # on HIVE; private (700), emptied after
EXIT_WORKER_DIED, EXIT_NOT_PINNED = 3, 4               # synth.py 1.2's; the image check's


class Unavailable(Exception):
    """The local voice cannot run from here; the message says why, in plain words."""


class LocalError(Exception):
    """The local voice ran and failed; the message says what happened."""


# ----------------------------------------------------------------------------- voices
def voices_for(hosts, header=None):
    """{host: Kokoro voice}: VOICES by position, then a script's `Local voices:` (or `Kokoro
    voices:`) header, e.g. "Maya=bf_emma, Leo=af_heart". ValueError for a voice not in the image."""
    v = {h: VOICES[min(i, len(VOICES) - 1)] for i, h in enumerate(hosts)}
    names = {h.lower(): h for h in hosts}
    for m in re.finditer(r"([A-Za-z][A-Za-z .'\-]*?)\s*[=:]\s*([^,;]+)", header or ""):
        h = names.get(m.group(1).strip().lower())
        if h:
            v[h] = m.group(2).strip()
    bad = sorted(set(v.values()) - set(VOICES_IN_IMAGE))
    if bad:
        raise ValueError(f"Kokoro voice(s) {', '.join(bad)} are not in the Core's image (it has "
                         f"{', '.join(VOICES_IN_IMAGE)})")
    return v


def spoken_names(voices):
    """{written name: spoken name} for the hosts whose Kokoro voice is a woman's and whose name
    is a man's (SPOKEN_NAMES): Leo in af_heart is said as Lena. Nothing else is renamed."""
    return {h: SPOKEN_NAMES[h] for h, v in voices.items() if h in SPOKEN_NAMES and WOMEN.match(v)}


# ----------------------------------------------------------------------------- sizing
def _gb(size):
    m = re.match(r"^(\d+(?:\.\d+)?)([KMGT])$", str(size))
    if not m:
        raise ValueError(f"not a size with a unit: {size!r}")
    return float(m.group(1)) * {"K": 1.0 / 1024 ** 2, "M": 1.0 / 1024, "G": 1.0, "T": 1024.0}[m.group(2)]


def job_resources(env=None):
    """(CPUs, memory in GB or None) this process may use: its SLURM allocation, else this
    machine (its CPU affinity; its physical memory)."""
    env = os.environ if env is None else env
    cpus = next((int(env[k]) for k in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE")
                 if str(env.get(k) or "").isdigit()), None)
    if cpus is None:
        try:
            cpus = len(os.sched_getaffinity(0))
        except AttributeError:
            cpus = os.cpu_count() or 1
    mem = None
    if str(env.get("SLURM_MEM_PER_NODE") or "").isdigit():
        mem = int(env["SLURM_MEM_PER_NODE"]) / 1024.0
    elif str(env.get("SLURM_MEM_PER_CPU") or "").isdigit():
        mem = int(env["SLURM_MEM_PER_CPU"]) * cpus / 1024.0
    elif not env.get("SLURM_JOB_ID"):
        try:
            mem = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / float(1 << 30)
        except (ValueError, OSError, AttributeError):
            mem = None
    return max(1, cpus), mem


def layout(cpus, mem_gb):
    """-> (workers, threads) for synth.py: workers of THREADS_PER_WORKER threads, as many as the
    CPUs hold and the memory holds at GB_PER_WORKER each after RESERVE_GB; at least one."""
    by_cpu = max(1, cpus // THREADS_PER_WORKER)
    by_mem = by_cpu if mem_gb is None else int((mem_gb - RESERVE_GB) // GB_PER_WORKER)
    workers = max(1, min(by_cpu, by_mem))
    return workers, max(1, min(cpus, workers * THREADS_PER_WORKER))


# ----------------------------------------------------------------------------- where it runs
def _on_hive():
    try:
        from core_submission import on_hive
    except ImportError:
        return os.path.isdir("/quobyte")
    return on_hive()


def _image_problem(path):
    if not os.path.exists(path):
        return (f"{path} is not here (the Core's pinned image: {os.path.dirname(path)}/README.md "
                "says how it is built)")
    if not os.access(path, os.R_OK):
        return (f"this HIVE account cannot read {path}: the Core's voice containers are for "
                "members of proteomics-grp")
    return None


def _apptainer_setup():
    return ("command -v apptainer >/dev/null 2>&1 || "
            "{ [ -r /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh && "
            "module load apptainer >/dev/null 2>&1; }\n"
            "command -v apptainer >/dev/null 2>&1 || "
            "{ echo '[podcast_local] no apptainer here (module load apptainer failed)' >&2; exit 5; }")


def _image_check(path, sha):
    q = shlex.quote(path)
    return (f"got=$(sha256sum {q} | cut -c1-64) || exit 5\n"
            f"[ \"$got\" = {sha} ] || {{ echo \"[podcast_local] {path} is not the pinned image: "
            f"sha256 $got, expected {sha}\" >&2; exit {EXIT_NOT_PINNED}; }}")


class Runner(object):
    """Runs one shell script that calls the containers, where they can run. A call's folder is
    `local` here and `run_dir` (a shell word) where the script runs; send/fetch move files
    between them (nothing to do on HIVE itself)."""
    remote = False
    scratch = None

    def resources(self):
        return self.cpus, self.mem_gb

    def open(self):
        d = tempfile.mkdtemp(prefix="podcast-local-", dir=self._base())
        os.chmod(d, 0o700)
        return d, shlex.quote(d)

    def _base(self):
        return os.environ.get("TMPDIR") or None

    def send(self, local, names):
        return

    def fetch(self, local, names):
        return

    def close(self, local):
        shutil.rmtree(local, ignore_errors=True)

    def describe(self):
        return {"where": self.where, "cpus": self.cpus, "memory_gb": self.mem_gb}


class Here(Runner):
    """A HIVE compute node: the containers run in this job's own allocation."""

    def __init__(self, env=None):
        self.cpus, self.mem_gb = job_resources(env)
        mem = f"{self.mem_gb:.0f} GB" if self.mem_gb else "memory not known"
        self.where = f"this HIVE job ({self.cpus} CPUs, {mem})"

    def run(self, script, minutes):
        if self.mem_gb is not None and self.mem_gb < MIN_GB:
            raise LocalError(
                f"this job has {self.mem_gb:.1f} GB of memory; the local voice needs at least "
                f"{MIN_GB:g} GB. Write the podcast job again (`make_podcast.py job`, which asks "
                f"for {JOB_CPUS} CPUs and {JOB_MEM}) and submit it: finished chunks are cached")
        return _run(["bash", "-c", script], minutes)


class LoginSrun(Runner):
    """A HIVE login node: each container call goes to a compute node in one blocking srun (the
    queue written out). Its folder must be one the node can see -- never this node's own /tmp:
    the episode's cache (`scratch`, set by render and verify), else the private
    podcast-local/ folder in the pipeline home."""

    def __init__(self):
        self.cpus, self.mem_gb = JOB_CPUS, _gb(JOB_MEM)
        self.where = (f"a HIVE compute node (one srun per step, {JOB_CPUS} CPUs and {JOB_MEM}, "
                      "sent from this login node)")

    def _base(self):
        if self.scratch:
            return self.scratch
        home = os.environ.get("PROTEOMICS_PIPELINE_HOME") or os.path.join(
            os.path.expanduser("~"), ".proteomics-pipeline")
        d = os.path.join(home, "podcast-local")
        os.makedirs(d, mode=0o700, exist_ok=True)
        return d

    def run(self, script, minutes):
        import on_compute
        try:
            argv = on_compute.srun_argv(["bash", "-c", script], minutes, JOB_MEM, JOB_CPUS,
                                        name="pp_podcast_local")
        except (ValueError, on_compute.QueueChoice) as e:
            raise LocalError(f"no compute node could be asked for: {e}")
        return _run(argv, minutes + 15)


def _remote(sub=""):
    """REMOTE_DIR/<sub> as a word for a command on HIVE (~ left for HIVE's shell or rsync to
    expand), and as a shell word for the job script ("$HOME/..." inside double quotes)."""
    path = REMOTE_DIR.rstrip("/") + (f"/{sub}" if sub else "")
    word = ("$HOME/" + path[2:]) if path.startswith("~/") else path
    return path, f'"{word}"'


class Remote(Runner):
    """Another computer with the skill's HIVE login: a private folder on HIVE (REMOTE_DIR/<id>,
    mode 700, in the HIVE home), one HIVE job per step (hive_exec.sh --job), files both ways."""
    remote = True

    def __init__(self, user):
        self.user = user
        self.cpus, self.mem_gb = JOB_CPUS, _gb(JOB_MEM)
        self.where = (f"a HIVE job sent from this computer (hive_exec.sh --job, {JOB_CPUS} CPUs "
                      f"and {JOB_MEM})")
        self._ids = {}

    def _hx(self, *args, timeout=600):
        return _run(["bash", os.path.join(HERE, "hive_exec.sh")] + list(args), timeout / 60.0)

    def open(self):
        local = tempfile.mkdtemp(prefix="podcast-local-")
        os.chmod(local, 0o700)
        rid = uuid.uuid4().hex[:12]
        self._ids[local] = rid
        base, here = _remote()[0], _remote(rid)[0]
        try:
            rc, out, err = self._hx(f"mkdir -p {base} && chmod 700 {base} && mkdir -m 700 {here}",
                                    timeout=120)
        except LocalError as e:
            rc, out, err = 1, "", str(e)
        if rc:
            self.close(local)                    # here and, if it was made, there
            raise LocalError(f"could not make a private folder on HIVE: {_tail(err or out)}")
        return local, _remote(rid)[1]

    def send(self, local, names):
        rid = self._ids[local]
        for n in names:
            rc, out, err = self._hx("--put", os.path.join(local, n), _remote(rid)[0] + "/",
                                    timeout=1800)
            if rc:
                raise LocalError(f"could not copy {n} to HIVE: {_tail(err or out)}")

    def run(self, script, minutes):
        return self._hx("--job", "--time", str(minutes), "--mem", JOB_MEM, "--cpus",
                        str(JOB_CPUS), script, timeout=(minutes + 30) * 60)

    def fetch(self, local, names):
        """`names` ending in "/" are folders: each lands as local/<name> (never its contents
        spilled into local/, which a trailing slash means to rsync)."""
        rid = self._ids[local]
        for n in names:
            folder = n.endswith("/")
            n = n.rstrip("/")
            rc, out, err = self._hx("--get", _remote(f"{rid}/{n}")[0],
                                    local + os.sep if folder else os.path.join(local, n),
                                    timeout=1800)
            if rc:
                raise LocalError(f"could not copy {n} back from HIVE: {_tail(err or out)}")

    def close(self, local):
        """Never raises: it runs in `finally` blocks, where its own error would hide the real
        one. A folder it could not remove is named, so a person can."""
        rid = self._ids.pop(local, None)
        shutil.rmtree(local, ignore_errors=True)
        if not rid:
            return
        try:
            rc, out, err = self._hx(f"rm -rf {_remote(rid)[0]}", timeout=120)
        except LocalError as e:
            rc, err = 1, str(e)
        if rc:
            sys.stderr.write(f"[podcast_local] [WARN] could not remove {_remote(rid)[0]} on HIVE "
                             f"({_tail(err, 200)}): remove it by hand (it holds the script's "
                             "lines or the episode)\n")


# Never handed to a container call, an srun or hive_exec.sh: nothing there needs them
PRIVATE_ENV = ("GEMINI_API_KEY", "GOOGLE_API_KEY")


def _run(argv, minutes):
    """-> (exit status, stdout, stderr); a run past `minutes` is killed and is a LocalError. The
    child gets this environment without PRIVATE_ENV."""
    env = {k: v for k, v in os.environ.items() if k not in PRIVATE_ENV}
    try:
        r = subprocess.run(argv, capture_output=True, text=True, timeout=max(60.0, minutes * 60.0),
                           stdin=subprocess.DEVNULL, env=env)
    except subprocess.TimeoutExpired:
        raise LocalError(f"no answer in {minutes:g} min: stopped")
    except OSError as e:
        raise LocalError(f"cannot run {argv[0]}: {e}")
    return r.returncode, r.stdout, r.stderr


def _tail(text, n=600):
    t = " ".join((text or "").split())
    return t[-n:] if t else "(no message)"


def find_runner(probe=True):
    """The Runner for here, or Unavailable with the reason. `probe`: off HIVE, ask HIVE (one
    hive_exec.sh call) whether this account can read the containers and has the skill there."""
    if _on_hive():
        for path in (KOKORO_IMAGE, WHISPER_IMAGE):
            why = _image_problem(path)
            if why:
                raise Unavailable(why)
        import on_compute
        return LoginSrun() if on_compute.on_login_node() else Here()
    try:
        from record_run import hive_login
        user = hive_login()
    except ImportError:
        user = None
    if not user:
        raise Unavailable("this computer has no HIVE login saved (check_access.sh saves it), so "
                          "the Core's local voice on HIVE cannot be reached")
    r = Remote(user)
    if probe:
        test = " && ".join(f"test -r {shlex.quote(p)}" for p in (KOKORO_IMAGE, WHISPER_IMAGE))
        try:
            rc, out, err = r._hx(f"if {test}; then echo PODCAST_LOCAL_IMAGES_OK; fi; "
                                 "test -d ~/proteomics-pipeline/scripts && "
                                 "echo PODCAST_LOCAL_SKILL_THERE; "
                                 "test -f ~/proteomics-pipeline/scripts/on_compute.py && "
                                 "echo PODCAST_LOCAL_SKILL_OK; true", timeout=120)
        except LocalError as e:
            raise Unavailable(f"HIVE did not answer ({e})")
        if rc or "PODCAST_LOCAL_" not in (out or ""):
            raise Unavailable(f"HIVE could not be asked ({_tail(err or out, 200)})")
        if "PODCAST_LOCAL_IMAGES_OK" not in out:
            raise Unavailable(f"the HIVE account {user} cannot read the Core's voice containers "
                              f"({TTS_ROOT}): they are for members of proteomics-grp")
        if "PODCAST_LOCAL_SKILL_OK" not in out:
            raise Unavailable(
                ("the skill on HIVE is older than this one: it has no on_compute.py, which runs "
                 "the HIVE job" if "PODCAST_LOCAL_SKILL_THERE" in out else
                 "the skill is not on HIVE, which runs the HIVE job")
                + ": bash scripts/hive_exec.sh --put-skill")
    return r


def _open(runner):
    """runner.open(), a folder that cannot be made being a LocalError like any other failure."""
    try:
        return runner.open()
    except OSError as e:
        raise LocalError(f"cannot make a work folder: {e}")


def _failed(what, rc, err, workers=None):
    """A container call's failure in plain words."""
    said = _tail(err)
    if rc == EXIT_WORKER_DIED:
        return f"{what} failed: a worker process died ({said})"
    if rc == EXIT_NOT_PINNED:
        return f"{what} was not run: {said}"
    if rc in (137, -9):
        return (f"{what} was killed (exit {rc}), most often for using more memory than the job "
                f"has{f' ({workers} worker(s) at about {GB_PER_WORKER:g} GB each)' if workers else ''}")
    return f"{what} failed (exit {rc}): {said}"


# ----------------------------------------------------------------------------- speak
def speak(runner, voices, lines, log=None):
    """Kokoro speaks each (speaker, text) line, all in ONE container call (start-up once).
    -> (a folder to close with runner.close(), [each line's WAV in it, in order], the call's
    timing record). LocalError when it fails."""
    log = log or (lambda m: None)
    local, run_dir = _open(runner)
    try:
        with open(os.path.join(local, "job.json"), "w", encoding="utf-8") as fh:
            json.dump({"voices": voices, "speed": 1.0,
                       "turns": [{"speaker": sp, "text": tx} for sp, tx in lines]},
                      fh, ensure_ascii=False)
        workers, threads = layout(*runner.resources())
        img = shlex.quote(KOKORO_IMAGE)
        script = "\n".join([
            "set -eo pipefail", f"D={run_dir}", _apptainer_setup(),
            _image_check(KOKORO_IMAGE, KOKORO_SHA256),
            f'apptainer exec --bind "$D" {img} python3 /opt/kokoro/synth.py "$D/job.json" - '
            f'--per-turn "$D/turns" --timing "$D/timing.json" --threads {threads} '
            f"--workers {workers}"])
        runner.send(local, ["job.json"])
        log(f"[render] Kokoro: {len(lines)} line(s) in one container call, on {runner.where}: "
            f"{workers} worker(s) x {THREADS_PER_WORKER} threads")
        t0 = time.time()
        rc, out, err = runner.run(script, SPEAK_MINUTES)
        if rc:
            raise LocalError(_failed("Kokoro", rc, err, workers))
        runner.fetch(local, ["turns/", "timing.json"])
        with open(os.path.join(local, "timing.json"), encoding="utf-8") as fh:
            timing = json.load(fh)
        timing.pop("per_turn", None)
        timing["wall_s"] = round(time.time() - t0, 1)
        paths = [os.path.join(local, "turns", f"turn-{k + 1:04d}.wav") for k in range(len(lines))]
        gone = [p for p in paths if not os.path.isfile(p)]
        if gone:
            raise LocalError(f"Kokoro wrote {len(paths) - len(gone)} of {len(paths)} line(s)")
        log(f"[render] Kokoro made {timing.get('turns_audio_s', 0):.0f} s of speech in "
            f"{timing['wall_s']:.0f} s (start-up {timing.get('startup_in_process_s') or 0:.0f} s)")
        return local, paths, timing
    except BaseException as e:               # Ctrl-C too: the folders hold the script's lines
        runner.close(local)
        if isinstance(e, LocalError) or not isinstance(e, Exception):
            raise
        raise LocalError(f"Kokoro: {type(e).__name__}: {e}")


# ----------------------------------------------------------------------------- hear
def hear(runner, audio, log=None):
    """Whisper transcribes the episode with both WHISPER_MODELS at once, in ONE container call.
    -> {model: its transcribe.py record (text, words with times, ...)}. LocalError when it
    fails."""
    log = log or (lambda m: None)
    local, run_dir = _open(runner)
    try:
        # A copy in the work folder, the only folder bound into the container: never the
        # episode's own folder, which on a login node holds the work folder (a read-only bind
        # of it could hide the writable one) and whose name may hold a , or : that --bind splits
        name = "episode" + (os.path.splitext(audio)[1].lower() or ".m4a")
        shutil.copyfile(audio, os.path.join(local, name))
        runner.send(local, [name])
        cpus = runner.resources()[0]
        each = max(1, cpus // len(WHISPER_MODELS))
        img = shlex.quote(WHISPER_IMAGE)
        script = "\n".join([
            "set -eo pipefail", f"D={run_dir}", f'A="$D"/{name}',
            _apptainer_setup(), _image_check(WHISPER_IMAGE, WHISPER_SHA256),
            f'one() {{ apptainer exec --bind "$D" {img} python3 '
            f'/opt/whisper/transcribe.py "$A" "$D/asr_$1.json" --model "$1" --threads {each} '
            '> "$D/asr_$1.log" 2>&1; }']
            + [f"one {m} & p{k}=$!" for k, m in enumerate(WHISPER_MODELS)]
            + ["r=0"] + [f"wait $p{k} || r=$?" for k in range(len(WHISPER_MODELS))]
            + ['[ "$r" -eq 0 ] || { tail -n 15 "$D"/asr_*.log >&2; exit "$r"; }'])
        log(f"[verify] Whisper ({' and '.join(WHISPER_MODELS)}) hears the episode back, on "
            f"{runner.where}; nothing is sent anywhere")
        t0 = time.time()
        rc, out, err = runner.run(script, HEAR_MINUTES)
        if rc:
            raise LocalError(_failed("Whisper", rc, err))
        names = [f"asr_{m}.json" for m in WHISPER_MODELS]
        runner.fetch(local, names)
        got = {}
        for m, n in zip(WHISPER_MODELS, names):
            with open(os.path.join(local, n), encoding="utf-8") as fh:
                got[m] = json.load(fh)
        log(f"[verify] heard back in {time.time() - t0:.0f} s")
        return got
    except (OSError, ValueError, KeyError) as e:
        raise LocalError(f"Whisper: {type(e).__name__}: {e}")
    finally:
        runner.close(local)


# ----------------------------------------------------------------------------- the agreement rule
# A spot where the audio departs from the text counts only when BOTH transcriptions show it at
# the same place (the 2026-10 trial: two of one model's apparent errors were Whisper's own).
def same_value(tok):
    """4.10 and 4.1 are the same number (the skill's number check agrees): decimals compare by
    value. Whole numbers stay as written, because a leading zero matters (0807 is not 807)."""
    import decimal
    if "." not in tok or not tok[:1].isdigit():
        return tok
    try:
        return format(decimal.Decimal(tok).normalize(), "f")
    except decimal.InvalidOperation:
        return tok


def join_letters(tokens):
    """"d i a n n" and "diann" are the same spelled-out symbol: runs of single letters are
    joined (verify_tokens joins only capital ones), so how a transcriber writes it does not
    matter. [(word, offset)] -> the same, a joined run keeping its first offset."""
    out = []
    for w, at in tokens:
        single = len(w) == 1 and w.isalpha()
        if single and out and out[-1][2]:
            out[-1] = (out[-1][0] + w, out[-1][1], True)
        else:
            out.append((w, at, single))
    return [(w, at) for w, at, _ in out]


_SENTENCE = re.compile(r"(?<=[.!?;:])\s+")


def tokens(mp, text):
    """make_podcast.verify_tokens, sentence by sentence: number words must not join across a full
    stop ("...up to thirty. Three cultures" is 30 and 3, not 33). -> [(word, offset)]"""
    out, pos = [], 0
    for piece in _SENTENCE.split(text or ""):
        start = text.find(piece, pos)
        out += [(w, start + a) for w, a, _ in mp.verify_tokens(piece)]
        pos = start + len(piece)
    return out


def departures(mp, spoken, heard):
    """-> [(kind, i1, i2, script, heard, offset in `heard`)] where the heard text departs from
    the spoken text: a number not heard as that number, 3+ heard words with nothing behind them
    ("invented"), or a gap (make_podcast's: VERIFY_GAP+ spoken words not heard)."""
    import difflib
    h = join_letters([(same_value(w), a) for w, a in tokens(mp, heard)])
    sp = [w for w, _ in join_letters([(same_value(w), a) for w, a in tokens(mp, spoken)])]
    words = [w for w, _ in h]
    found = []
    for tag, i1, i2, j1, j2 in difflib.SequenceMatcher(None, sp, words, autojunk=False).get_opcodes():
        if tag == "equal":
            continue
        s_, h_ = sp[i1:i2], words[j1:j2]
        at = h[j1][1] if j1 < len(h) else len(heard)
        nums = [w for w in s_ if w[0].isdigit()]
        if nums and not all(n in h_ for n in nums):
            found.append(("number", i1, i2, " ".join(s_), " ".join(h_), at))
        elif len(h_) - len(s_) >= 3:
            found.append(("invented", i1, i2, " ".join(s_), " ".join(h_), at))
        elif tag in ("delete", "replace") and i2 - i1 >= getattr(mp, "VERIFY_GAP", 6):
            found.append(("gap", i1, i2, " ".join(s_), " ".join(h_), at))
    return found


def time_at(asr, offset):
    """The start time (s) of the transcript word at character `offset` of asr["text"]."""
    text, pos, best = asr.get("text") or "", 0, None
    for w in asr.get("words") or []:
        t = (w.get("w") or "").strip()
        i = text.find(t, pos) if t else -1
        if i < 0:
            continue
        if i > offset:
            break
        best, pos = w.get("start"), i + len(t)
    return best


def compare_two(mp, spoken, asr1, asr2):
    """-> (agreed, one_only): the spots where BOTH transcripts depart from the spoken text at the
    same place (the audio, most likely), and those only one shows (most likely that
    transcriber). Each: {kind, script, heard, by, at_s}."""
    asrs = (asr1, asr2)
    deps = [departures(mp, spoken, a.get("text") or "") for a in asrs]
    agreed, one_only = [], []
    for k in (0, 1):
        for kind, i1, i2, st, h, at in deps[k]:
            both = any(kk == kind and k1 <= i2 and i1 <= k2 for kk, k1, k2, *_ in deps[1 - k])
            if both and k == 1:
                continue                                   # listed once, from the first
            rec = {"kind": kind, "script": st, "heard": h, "by": "both" if both else k,
                   "at_s": time_at(asrs[k], at)}
            (agreed if both else one_only).append(rec)
    return agreed, one_only


def clock(s):
    return "?" if s is None else f"{int(s // 60)}:{int(s % 60):02d}"


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv[:1] != ["where"]:
        print(__doc__.strip())
        return 0 if argv[:1] in (["-h"], ["--help"]) else 2
    try:
        r = find_runner()
        print(json.dumps(dict(r.describe(), available=True, kokoro=KOKORO_IMAGE,
                              whisper=WHISPER_IMAGE), indent=2))
        return 0
    except Unavailable as e:
        print(json.dumps({"available": False, "why": str(e)}, indent=2))
        return 1


if __name__ == "__main__":
    sys.exit(main())
