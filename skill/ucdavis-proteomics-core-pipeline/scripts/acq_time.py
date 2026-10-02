"""
acq_time.py -- WHEN a run was acquired, read from the run itself, and the one place the skill
decides which time zone that time is in.

Two vendors, two different records of the same moment:

  Bruker .d    analysis.tdf GlobalMetadata AcquisitionDateTime, ISO 8601 WITH the instrument
               PC's UTC offset ("2026-09-30T13:04:03.063-07:00"). Unambiguous.
  Thermo .raw  ThermoRawFileParser's metadata JSON, FileProperties NCIT:C69199 "Content
               Creation Date" (MetadataWriter.cs: FileHeader.CreationDate.ToString()). NO
               offset, and the text is in the parser's culture: "09/30/2026 03:06:18" from the
               Core's TRFP 2.0.0.0 on HIVE (invariant culture), "9/30/2026 3:06:18 AM" from a
               US-English Windows.

What the Thermo value IS, measured on HIVE (srun job 24264231, 2026-10-01): for an Orbitrap
Exploris 480 QC run on a 120-min gradient, TRFP printed 03:06:18, and the .raw was last written
at 05:06:20 PDT -- the end of the run. A Fusion Lumos QC run: 13:39:43, last written 15:39:46
PDT. TRFP printed the same value with TZ=UTC and TZ=Asia/Tokyo, so the parser does not convert
it to the reading machine's zone: it is the instrument PC's own wall clock at the START of
acquisition. The Core's instrument PCs keep Pacific time, so it is read as CORE_TZ, and every
record that carries such a time says so (an assumption, not something the file states).

Times are compared in UTC; they are shown to people as the Core's local time (CORE_TZ), the
clock on the instrument and in STAN's dashboard. zoneinfo is used when the computer has the
time-zone database; Windows Pythons often do not (it needs the tzdata package), and then the US
Pacific rule (DST from the second Sunday of March to the first Sunday of November, 2 a.m.,
since 2007) stands in for CORE_TZ -- and the note says so. Any other zone without zoneinfo is
not guessed: the time is unreadable, with that reason.

Stdlib only; imported by detect_acquisition.py and qc_bracket.py.
"""
import os
import re
import sqlite3
from datetime import datetime, timedelta, timezone, tzinfo

# The Core's instrument PCs (UC Davis): what a Thermo creation date's wall clock is read as.
CORE_TZ = "America/Los_Angeles"
# TRFP metadata: FileProperties, keyed by ACCESSION (TRFP's `name` strings are not stable --
# detect_acquisition.py's CV_* terms are matched the same way).
CV_CREATION = "NCIT:C69199"      # "Content Creation Date"
BRUKER_FIELD = "AcquisitionDateTime"
THERMO_SOURCE = (f"ThermoRawFileParser metadata {CV_CREATION} (Content Creation Date): the "
                 "instrument PC's clock at the start of acquisition, with no time zone")
BRUKER_SOURCE = f"analysis.tdf GlobalMetadata {BRUKER_FIELD} (start of acquisition, with offset)"
# .NET DateTime.ToString() in the two cultures the Core's parsers run in. A day-first culture
# (en-GB: 30/09/2026 03:06:18) prints the same shape as the first and cannot be told apart from
# it; check_against_mtime() catches the impossible readings that produces.
THERMO_FORMATS = (("%m/%d/%Y %H:%M:%S", "invariant culture, MM/dd/yyyy HH:mm:ss"),
                  ("%m/%d/%Y %I:%M:%S %p", "US English, M/d/yyyy h:mm:ss AM/PM"))
# A file cannot be created after it was last written. An acquisition time later than the file's
# modification time by more than this was read wrongly (or the clock was wrong).
AFTER_MTIME_SLACK = timedelta(hours=1)
_FRACTION = re.compile(r"(\.\d+)(?=([+-]\d\d:?\d\d|Z)?$)")


def _zone(name):
    """(tzinfo or None, how). zoneinfo when the computer has the database; for CORE_TZ without
    it, the US Pacific rule below."""
    try:
        from zoneinfo import ZoneInfo
        return ZoneInfo(name), "zoneinfo"
    except Exception:                      # ImportError, ZoneInfoNotFoundError (no tzdata)
        pass
    if name == CORE_TZ:
        return _USPacific(), "built-in US Pacific rule (no time-zone database on this computer)"
    return None, f"no time-zone database on this computer to read {name!r} with"


def _nth_sunday(year, month, n):
    d = datetime(year, month, 1)
    d += timedelta(days=(6 - d.weekday()) % 7)        # first Sunday
    return d + timedelta(weeks=n - 1)


class _USPacific(tzinfo):
    """America/Los_Angeles by the US rule in force since 2007: PDT (UTC-7) from 02:00 on the
    second Sunday of March to 02:00 on the first Sunday of November, PST (UTC-8) otherwise.
    A wall-clock time in November's repeated hour reads as its first (PDT) occurrence, as
    zoneinfo's fold=0 does."""

    @staticmethod
    def _bounds(year):
        """(DST start, DST end) as wall-clock times: 02:00 on each Sunday."""
        return (_nth_sunday(year, 3, 2) + timedelta(hours=2),
                _nth_sunday(year, 11, 1) + timedelta(hours=2))

    def utcoffset(self, dt):
        if dt is None:
            return None
        start, end = self._bounds(dt.year)
        return timedelta(hours=-7 if start <= dt.replace(tzinfo=None) < end else -8)

    def dst(self, dt):
        return None if dt is None else self.utcoffset(dt) + timedelta(hours=8)

    def tzname(self, dt):
        return "PDT" if self.utcoffset(dt) == timedelta(hours=-7) else "PST"

    def fromutc(self, dt):
        start, end = self._bounds(dt.year)
        u = dt.replace(tzinfo=None)
        dst = start + timedelta(hours=8) <= u < end + timedelta(hours=7)   # the same instants, UTC
        return (u + timedelta(hours=-7 if dst else -8)).replace(tzinfo=self)


def parse_iso(text):
    """An ISO 8601 timestamp -> datetime (aware when it carries an offset or Z), or None.
    Python 3.9's fromisoformat takes neither 'Z' nor 7 fractional digits (.NET's "o" format);
    both are normalised first."""
    if not isinstance(text, str) or not text.strip():
        return None
    s = text.strip().replace(" ", "T", 1) if re.match(r"^\d{4}-\d\d-\d\d \d", text.strip()) \
        else text.strip()
    if s.endswith(("Z", "z")):
        s = s[:-1] + "+00:00"
    s = _FRACTION.sub(lambda m: (m.group(1) + "000000")[:7], s)
    s = re.sub(r"([+-]\d\d)(\d\d)$", r"\1:\2", s)
    try:
        return datetime.fromisoformat(s)
    except ValueError:
        return None


def localize(naive, tz_name=CORE_TZ):
    """(aware datetime, how) for a wall-clock time in `tz_name`, or (None, why)."""
    tz, how = _zone(tz_name)
    if tz is None:
        return None, how
    return naive.replace(tzinfo=tz), how


def to_utc(dt):
    return dt.astimezone(timezone.utc)


def iso_utc(dt):
    """'2026-09-30T20:04:03Z' -- the form every record carries."""
    return to_utc(dt).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def local_text(dt, tz_name=CORE_TZ, fmt="%Y-%m-%d %H:%M"):
    """A time as the Core's clock shows it, e.g. '2026-09-30 13:04 PDT'. UTC when the zone
    cannot be read here (said in the text)."""
    tz, _how = _zone(tz_name)
    if tz is None:
        return to_utc(dt).strftime(fmt) + " UTC"
    loc = to_utc(dt).astimezone(tz)
    return f"{loc.strftime(fmt)} {loc.tzname() or tz_name}"


def local_date(dt, tz_name=CORE_TZ):
    tz, _how = _zone(tz_name)
    return (to_utc(dt).astimezone(tz) if tz else to_utc(dt)).date()


def result(when=None, source=None, note=None):
    """The record every reader returns: acquired_at (UTC ISO, or None), its source, and a note
    (the time-zone assumption, or why there is no time)."""
    return {"acquired_at": iso_utc(when) if when is not None else None,
            "acquired_at_source": source if when is not None else None,
            "acquired_at_note": note}


def from_bruker_value(value, tz_name=CORE_TZ):
    """GlobalMetadata AcquisitionDateTime text -> result()."""
    dt = parse_iso(value)
    if dt is None:
        return result(note=f"{BRUKER_FIELD} {value!r} is not an ISO 8601 time")
    if dt.tzinfo is None:
        loc, how = localize(dt, tz_name)
        if loc is None:
            return result(note=f"{BRUKER_FIELD} {value!r} has no UTC offset and {how}")
        return result(loc, BRUKER_SOURCE.replace("with offset", "no offset recorded"),
                      f"{BRUKER_FIELD} carries no UTC offset: read as {tz_name} wall-clock time "
                      f"({how}) -- an assumption")
    return result(dt, BRUKER_SOURCE)


def read_bruker(d, tz_name=CORE_TZ):
    """result() for a Bruker .d, read from analysis.tdf (opened read-only and immutable)."""
    tdf = os.path.join(d, "analysis.tdf")
    if not os.path.exists(tdf):
        return result(note="no analysis.tdf in the .d")
    from bruker_tdf import connect_tdf      # the one way the skill opens a tdf
    try:
        con = connect_tdf(tdf)
        try:
            row = con.execute("SELECT Value FROM GlobalMetadata WHERE Key = ?",
                              (BRUKER_FIELD,)).fetchone()
        finally:
            con.close()
    except sqlite3.Error as e:
        return result(note=f"analysis.tdf could not be read ({e})")
    if not row or not row[0]:
        return result(note=f"analysis.tdf GlobalMetadata has no {BRUKER_FIELD}")
    return from_bruker_value(str(row[0]), tz_name)


def thermo_creation_text(meta):
    """The raw Content Creation Date text from TRFP metadata, or None."""
    for term in (meta or {}).get("FileProperties") or []:
        if isinstance(term, dict) and term.get("accession") == CV_CREATION:
            v = term.get("value")
            return v.strip() if isinstance(v, str) and v.strip() else None
    return None


def from_thermo_value(value, tz_name=CORE_TZ):
    """TRFP's Content Creation Date text -> result(): read as `tz_name` wall-clock time."""
    if not value:
        return result(note=f"ThermoRawFileParser metadata has no {CV_CREATION} (Content Creation "
                           f"Date)")
    dt = parse_iso(value)
    form = "ISO 8601"
    if dt is None:
        for fmt, label in THERMO_FORMATS:
            try:
                dt, form = datetime.strptime(value, fmt), label
                break
            except ValueError:
                continue
    if dt is None:
        return result(note=f"Content Creation Date {value!r} is in no format this reads "
                           f"({'; '.join(lab for _f, lab in THERMO_FORMATS)}; ISO 8601)")
    if dt.tzinfo is not None:
        return result(dt, THERMO_SOURCE.replace("with no time zone", "with an offset"))
    loc, how = localize(dt, tz_name)
    if loc is None:
        return result(note=f"Content Creation Date {value!r} has no time zone and {how}")
    return result(loc, THERMO_SOURCE,
                  f"read {value!r} ({form}) as {tz_name} wall-clock time ({how}): the .raw "
                  f"records no time zone, and the Core's instrument PCs keep Pacific time -- an "
                  f"assumption for a .raw acquired anywhere else")


def check_against_mtime(rec, path):
    """A warning (str) when the acquisition time is later than the file was last written --
    impossible for a correct reading (a day-first date read month-first, a wrong PC clock).
    None when it is plausible or either time is missing."""
    when = parse_iso(rec.get("acquired_at"))
    try:
        mtime = datetime.fromtimestamp(os.path.getmtime(path), tz=timezone.utc)
    except OSError:
        return None
    if when is not None and when - mtime > AFTER_MTIME_SLACK:
        return (f"the acquisition time read ({iso_utc(when)}) is after the file was last "
                f"written ({iso_utc(mtime)}): the date was probably read wrongly (a day-first "
                f"date read month-first?) or the instrument PC's clock was wrong")
    return None
