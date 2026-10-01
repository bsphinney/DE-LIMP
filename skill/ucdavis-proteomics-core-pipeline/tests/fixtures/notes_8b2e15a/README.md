# Notes written by notes.py at 8b2e15a

The inbox went live on HIVE before merge, holding files written by `notes.py` at commit 8b2e15a.
Every later `notes.py` must read such files as they are: `tests/test_notes.py` (`Compatibility`)
does.

These three files are that code's own output, byte for byte: `Inbox("…", "brettsp").send(...)`
from `git show 8b2e15a:skill/ucdavis-proteomics-core-pipeline/scripts/notes.py`, then
`Inbox("…", "msalemi").ack(..., reply=...)`, with the clock fixed. They follow the live note's
shape (a submission number in the subject, a service folder as the session), with a placeholder
number (`PROT_0000`) and folder (`PI_Example/SET1-28`); the body is illustrative. Never edit them:
regenerate them from 8b2e15a instead.
