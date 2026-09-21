# An analysis.tdf is read so SQLite never replays a stale -wal and never writes a
# -wal/-shm anywhere. On HIVE a stale mid-acquisition analysis.tdf-wal beside a
# finished tdf, replayed by a read-write open, truncated 342 runs (see
# skill/ucdavis-proteomics-core-pipeline/scripts/bruker_tdf.py). The readers here
# already work on a temp COPY of the tdf, so the original is safe; the immutable open
# stops a WAL-mode copy leaving <tmp>.tdf-shm / -wal orphans in tempdir() on every read,
# and keeps the open safe if the copy is ever dropped.

source(file.path(project_root, "R", "helpers_instrument.R"))

N_FRAMES <- 1000L
N_STALE <- 10L

# A finished run with a stale mid-acquisition -wal beside it -- the HIVE state.
make_stale_wal_tdf <- function(dir_name = "run.d") {
  d <- file.path(tempfile("tdf_fixture_"), dir_name)
  dir.create(d, recursive = TRUE)
  tdf <- file.path(d, "analysis.tdf")
  con <- DBI::dbConnect(RSQLite::SQLite(), tdf)
  DBI::dbGetQuery(con, "PRAGMA journal_mode=WAL")
  DBI::dbGetQuery(con, "PRAGMA wal_autocheckpoint=0")
  keys <- c(InstrumentName = "timsTOF HT", InstrumentSerialNumber = "1",
            AcquisitionSoftware = "timsControl", AcquisitionSoftwareVersion = "5",
            SampleName = "s", AcquisitionDateTime = "2025-08-09T00:00:00",
            OperatorName = "op", MzAcqRangeLower = "100", MzAcqRangeUpper = "1700",
            OneOverK0AcqRangeLower = "0.7", OneOverK0AcqRangeUpper = "1.3",
            MethodName = "dia-PASEF.m")
  DBI::dbExecute(con, "CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT)")
  DBI::dbExecute(con, "INSERT INTO GlobalMetadata VALUES (?, ?)",
                 params = list(names(keys), unname(keys)))
  DBI::dbExecute(con, paste("CREATE TABLE Frames (Id INTEGER PRIMARY KEY, Time REAL,",
                            "MsMsType INTEGER, TimsId INTEGER, NumScans INTEGER,",
                            "SummedIntensities INTEGER)"))
  ins <- function(ids) {
    n <- length(ids)
    DBI::dbExecute(con, "INSERT INTO Frames VALUES (?, ?, ?, ?, ?, ?)",
                   params = list(ids, ids * 0.1, rep(0L, n), (ids - 1L) * 64L,
                                 rep(5L, n), rep(1000L, n)))
  }
  ins(seq_len(N_STALE))
  wal <- paste0(tdf, "-wal")
  stale <- readBin(wal, "raw", file.size(wal))
  ins(seq.int(N_STALE + 1L, N_FRAMES))
  DBI::dbGetQuery(con, "PRAGMA wal_checkpoint(TRUNCATE)")
  DBI::dbDisconnect(con)
  writeBin(stale, wal)
  tdf
}

fingerprint <- function(dir) {
  f <- list.files(dir, full.names = TRUE)
  setNames(unname(tools::md5sum(f)), basename(f))
}

count_frames <- function(con) DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM Frames")$n

test_that("fixture: a plain SQLITE_RO open replays the stale -wal", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  tdf <- make_stale_wal_tdf()
  con <- DBI::dbConnect(RSQLite::SQLite(), tdf, flags = RSQLite::SQLITE_RO)
  on.exit(DBI::dbDisconnect(con))
  # proves the fixture has the hazard, so the tests below cannot pass vacuously
  expect_equal(count_frames(con), N_STALE)
})

test_that("tdf_dbconnect_ro reads the finished file and writes nothing beside it", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  tdf <- make_stale_wal_tdf()
  before <- fingerprint(dirname(tdf))
  con <- tdf_dbconnect_ro(tdf)
  expect_equal(count_frames(con), N_FRAMES)
  DBI::dbDisconnect(con)
  expect_identical(fingerprint(dirname(tdf)), before)
})

test_that("tdf_dbconnect_ro refuses writes", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  con <- tdf_dbconnect_ro(make_stale_wal_tdf())
  on.exit(DBI::dbDisconnect(con))
  expect_error(DBI::dbExecute(con, "CREATE TABLE scribble (x)"))
})

test_that("the immutable URI survives characters that end or escape a URI path", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  tdf <- make_stale_wal_tdf("odd #1? 50%.d")
  expect_match(sqlite_immutable_uri(tdf), "^file:///.*\\?mode=ro&immutable=1$")
  con <- tdf_dbconnect_ro(tdf)
  on.exit(DBI::dbDisconnect(con))
  expect_equal(count_frames(con), N_FRAMES)
})

test_that("the tdf readers see the finished file and leave nothing behind", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  tdf <- make_stale_wal_tdf()
  orig <- fingerprint(dirname(tdf))
  tmp_before <- list.files(tempdir())

  meta <- parse_timstof_from_tdf(tdf)
  expect_null(meta$parse_error)
  expect_equal(meta$n_frames, N_FRAMES)
  expect_equal(meta$instrument_model, "timsTOF HT")

  tic <- extract_tic_timstof(tdf)
  expect_equal(nrow(tic), N_FRAMES)

  expect_identical(fingerprint(dirname(tdf)), orig)
  leftovers <- setdiff(list.files(tempdir()), tmp_before)
  expect_equal(leftovers[grepl("\\.tdf", leftovers)], character(0))
})
