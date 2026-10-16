"""
safe_write.py -- the ONE way the skill writes a file another process will run or read: a job
script, a submit.sh, a file list or a cfg a job reads.

Core session folders are group-writable (umask 0002, no sticky bit, every proteomics-grp member
can write into them). `open(path, "w")` then `os.chmod(path, ...)` both FOLLOW a symlink, so a
group member who planted `step2_firstpass.sbatch -> ~victim/.bashrc` had the generator overwrite
the victim's file with a job script and make it executable (security review of 2.12, reproduced
on make_podcast.py). write_file() never writes through a link at `path`:

  1. it creates a NEW file beside `path` -- os.open(O_CREAT | O_EXCL | O_NOFOLLOW): EXCL fails on
     anything already there, a planted link included, and the name is random;
  2. writes the text and sets the mode on that open file (fchmod: the umask plays no part);
  3. os.replace()s it onto `path`, which swaps the directory entry -- a link at `path` is
     replaced, never followed -- and is atomic, so a reader sees the old file or the new one.

Text is written with "\\n" line ends on every platform: these are bash scripts and lists bash
reads. Stdlib only.
"""
import io
import os
import secrets


def write_file(path, text, mode=0o644):
    """Write `text` to `path` as a new file with permissions `mode` (see the module docstring).
    Returns `path`. Nothing is left behind if the write fails."""
    d = os.path.dirname(os.path.abspath(path))
    flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY | getattr(os, "O_NOFOLLOW", 0)
    while True:
        tmp = os.path.join(d, f".{os.path.basename(path)}.{secrets.token_hex(6)}.tmp")
        try:
            fd = os.open(tmp, flags, 0o600)
            break
        except FileExistsError:
            continue
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
            fh.flush()
            if hasattr(os, "fchmod"):
                os.fchmod(fh.fileno(), mode)
        if not hasattr(os, "fchmod"):
            os.chmod(tmp, mode)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise
    return path


class open_text:
    """`with open_text(path, mode) as fh: fh.write(...)` -- the text is collected and handed to
    write_file() when the block ends cleanly; nothing is written if it raises. For a writer that
    streams (json.dump, a loop of fh.write)."""

    def __init__(self, path, mode=0o644):
        self.path, self.mode, self._buf = path, mode, None

    def __enter__(self):
        self._buf = io.StringIO()
        return self._buf

    def __exit__(self, exc_type, exc, tb):
        if exc_type is None:
            write_file(self.path, self._buf.getvalue(), self.mode)
        return False
