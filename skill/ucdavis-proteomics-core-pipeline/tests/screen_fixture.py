#!/usr/bin/env python3
"""A finished limpa (DPC-Quant) analysis on a synthetic DIA-NN report, for the protein screen's tests.

make_session(tmp) writes <tmp>/session/output/search/{report.parquet,report.stats.tsv} and
<tmp>/session/output/tables/ (run_de.R --method dpc --contrasts B-A,C-A) and returns the paths.

3 groups x 3 runs named S-01.a ... S-09.a. P00001-P00005 are up 1.5 log2 in B; P00051-P00055 are in C
only (not measured in A or B, so B-A does not test them). One precursor of P00006 has Q.Value 0.02
in S-04.a. `P00090;Cont_P99999` is a contaminant group.

Run names: run_de.R matches the conditions' File.Name to the report's Run exactly (no extension is
stripped), so the conditions carry the same odd names (S-01.a), not S-01.a.raw.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, HERE)

from test_run_de_contaminants import r_has, rscript  # noqa: E402,F401

RUN_DE = os.path.join(SCRIPTS, "run_de.R")
R_NEEDS = ("limpa", "limma", "arrow", "nanoparquet", "jsonlite", "statmod")

SYNTH_R = r'''
set.seed(11)
runs <- sprintf("S-%02d.a", 1:9); grp <- rep(c("A", "B", "C"), each = 3)
rows <- list()
add <- function(pg, prec, base, only = NULL, effB = 0, ids = pg) {
  for (r in seq_along(runs)) {
    if (!is.null(only) && !(grp[r] %in% only)) next
    lv <- base + rnorm(1, 0, 0.2) + (grp[r] == "B") * effB
    int <- 2^(lv + rnorm(1, 0, 0.3))
    if (base < 17 && runif(1) < plogis(-(lv - 13) * 2)) next
    rows[[length(rows) + 1]] <<- data.frame(Run = runs[r], Precursor.Id = prec,
      Protein.Group = pg, Protein.Ids = ids, Protein.Names = paste0(pg, "_X"), Genes = sub("P0*", "G", pg),
      Proteotypic = 1L, Precursor.Normalised = int, Precursor.Quantity = int * 1.1,
      Normalisation.Factor = 0.9, Quantity.Quality = 0.95, Empirical.Quality = 0.8,
      Ms1.Area = int / 3, RT = 10 + (match(prec, unique(prec)) %% 7),
      Q.Value = 0.001, Lib.Q.Value = 0.001, Lib.PG.Q.Value = 0.001, PG.Q.Value = 0.001,
      Global.Q.Value = 0.001, Global.PG.Q.Value = 0.001, stringsAsFactors = FALSE)
  }
}
for (i in 1:55) {
  pg <- sprintf("P%05d", i)
  base <- if (i == 6) 19 else if (i > 50) runif(1, 17, 20) else runif(1, 14, 20)
  only <- if (i > 50) "C" else NULL
  for (k in seq_len(NPREC)) add(pg, sprintf("PEP%dK%d2", i, k), base + rnorm(1, 0, 0.5), only = only,
                     effB = if (i <= 5) 1.5 else 0)
}
for (k in seq_len(min(3, NPREC))) add("P00090;Cont_P99999", sprintf("CONT%d2", k), 19, ids = "P00090;Cont_P99999")
d <- do.call(rbind, rows)
hit <- d$Run == "S-04.a" & d$Precursor.Id == "PEP6K12"
stopifnot(sum(hit) == 1)
d$Q.Value[hit] <- 0.02
pgm <- aggregate(Precursor.Normalised ~ Protein.Group + Run, d, sum)
names(pgm)[3] <- "PG.MaxLFQ"
d <- merge(d, pgm, by = c("Protein.Group", "Run"))
d$PG.MaxLFQ.Quality <- 0.9
dir.create(file.path(OUT, "output", "search"), recursive = TRUE, showWarnings = FALSE)
arrow::write_parquet(d, file.path(OUT, "output", "search", "report.parquet"))
# DIA-NN writes File.Name as the raw file's full path (the last run here is a Windows path)
fn <- paste0("/quobyte/lab/raw/", runs, ".raw"); fn[length(fn)] <- paste0("D:", intToUtf8(92), "raw", intToUtf8(92), runs[length(runs)], ".d")
write.table(data.frame(File.Name = fn, Precursors.Identified = as.integer(table(d$Run)[runs]),
                       Normalisation.Instability = round(runif(length(runs), 0.01, 0.2), 3)),
            file.path(OUT, "output", "search", "report.stats.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
write.csv(data.frame(File.Name = runs, Group = grp), file.path(OUT, "conditions.csv"), row.names = FALSE)
'''


def make_session(tmp, method="dpc", precursors=4, contrasts=None, de_args=(), block_col=None):
    """Build the session in <tmp>/session and run the DE into output/tables. Raises on failure.

    `contrasts` (default B-A,C-A) are run_de.R's --contrasts, any limma expression; `de_args` are more run_de.R
    arguments; `block_col` adds a metadata column of that name (three blocks, one sample of each group in each:
    a paired design) for --block.
    `precursors` per protein (contaminant group: at most 3); with 1, limpa quantifies row by row.
    """
    root = os.path.join(tmp, "session")
    os.makedirs(root, exist_ok=True)
    rscript(f'OUT <- "{root}"\nNPREC <- {int(precursors)}L\n' + SYNTH_R)
    output = os.path.join(root, "output")
    tables = os.path.join(output, "tables")
    contrasts = list(contrasts or ["B-A", "C-A"])
    if block_col:
        cf = os.path.join(root, "conditions.csv")
        with open(cf) as fh:
            lines = fh.read().splitlines()
        with open(cf, "w") as fh:
            fh.write(lines[0] + "," + block_col + "\n")
            fh.write("\n".join(l + ",%d" % (i % 3 + 1) for i, l in enumerate(lines[1:])) + "\n")
    p = subprocess.run(["Rscript", RUN_DE, "--input", os.path.join(output, "search", "report.parquet"),
                        "--metadata", os.path.join(root, "conditions.csv"), "--method", method,
                        "--contrasts", ",".join(contrasts), "--outdir", tables, *de_args],
                       capture_output=True, text=True, cwd=tmp)
    if p.returncode != 0:
        raise AssertionError("run_de.R failed:\n" + p.stderr[-3000:])
    return dict(session=root, output=output, tables=tables, contrasts=contrasts)
