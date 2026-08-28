# Environment & acquisition reference

How the skill adapts to where it's running, plus FASTA resolution.

## Platform classes (`detect_env.sh`)

| class | how detected | engine acquisition | execution |
|---|---|---|---|
| `hpc` | `sbatch` on PATH, or `/quobyte/proteomics-grp` visible | reuse existing Apptainer `.sif` (DIA-NN); download Sage binary | **submit via sbatch** (never login-node) |
| `mac` | `uname` = darwin | DIA-NN via **Docker** (no native mac build); Sage native | inline |
| `linux` | everything else | native binaries | inline |

`uc_davis_hive` is true when `/quobyte/proteomics-grp` exists → enables FASTA reuse
and the existing `.sif`.

Container runtime preference: hpc→apptainer, mac→docker, linux→native.

## DIA-NN by environment

- **No native macOS build exists.** On mac you must run DIA-NN through Docker. Set
  `DIANN_DOCKER_IMAGE` to a built image, or build one from the Academia Linux zip's
  bundled Dockerfile. `acquire_tools.sh` writes a note when this is unresolved.
- **HIVE (Proteomics Core):** DIA-NN is kept under `/quobyte/proteomics-grp/dia-nn/`.
  Recent versions are **native builds**, e.g. DIA-NN **2.6.0** at
  `build_260/diann-2.6.0/diann-linux`; older 2.3.0 is a `.sif`. `acquire_tools.sh`
  resolves the **pinned** version by looking for `build_*/diann-<version>/diann-linux`
  first, then a version-matched `.sif`, and **never silently substitutes a different
  version** (reproducibility). The facility's `run_diann_*.sbatch` in that folder is
  the reference invocation. AlphaDIA is also on HIVE at
  `/quobyte/proteomics-grp/apptainers/alphadia.sif` (auto-reused).
- **Linux native:** needs glibc ≥ Linux Mint 21.2 and .NET 8. If missing, prefer
  Docker/Apptainer.
- DIA-NN reads `.raw`/`.d` natively from 2.1+.

## Version pinning (reproducibility)

`acquire_tools.sh` honors `PIN_ENGINE`/`PIN_VERSION` from the workflow bundle and
caches under `~/.proteomics-pipeline/tools/<engine>/<version>/`. Different pinned
versions coexist. The written `tools.json` records `pinned` and `versions` so the
report can state exactly what ran. **Always pass the bundle's engine+version** — a
result from "latest" is not reproducible.

## FASTA resolution (`fetch_fasta.py`)

### `resolve` — organism → proteome (always run this; never guess a UPID)
`fetch_fasta.py resolve --organism "mouse"` returns ranked `candidates` +
`selected` + `needs_menu` + `notes`. Accepts a common name, a scientific name, an
NCBI taxid, or a `UP…` accession (`--taxid` also works). Show the user `selected`
and confirm. When `needs_menu` is true, present the menu instead of auto-picking.

`scripts/resolve_organism.py` is a thin alias for this same resolver, kept for
older call sites. **One implementation** — don't add a third.

How it picks, and why each rule exists:
- **Curated table first.** 18 organisms a core facility actually sees map
  name/alias → *taxid*, and the proteome is then resolved live from that taxid.
  Free-text search answers some queries badly: "Escherichia coli K-12" returns
  five Non-Reference MG1655 assemblies and never surfaces the real reference
  `UP000000625` at all. Only the taxid is curated — a pinned `UP…` accession goes
  stale (the dog reference moved from `UP000002254` to `UP000805418`), so the
  table's accession is used only as an offline fallback and a staleness check,
  which is reported in `notes`.
- **Exact proteome-type match.** `"reference" in "Non Reference proteome"` is
  True, so a substring test ranks strain assemblies as references.
- **Strain qualifier stripped when name-matching.** `scientificName` is
  `Saccharomyces cerevisiae (strain ATCC 204508 / S288c)`, so the bare species
  name matched nothing and ranking fell through to protein count — which put
  *S. pastorianus* first.
- **"Excluded" proteomes dropped** unless nothing else matches.
- **Auto-picks only when unambiguous**: one reference proteome the user actually
  named. "mouse" also matches *Myotis myotis* and mouse-ear cress; "baker's yeast"
  matches two S. cerevisiae strains → menu.

Common IDs (still confirm with `resolve`): human `UP000005640`, mouse
`UP000000589`, rat `UP000002494`, yeast `UP000002311`, E. coli `UP000000625`.

### `fetch` — proteome → search FASTA
Priority, cheapest/most-trusted first:
1. `--path` override → used verbatim (pre-staged proteome).
2. **HIVE** (`--hive`): reuse `/quobyte/proteomics-grp/MRS/`. Matches only files
   whose name starts with the proteome ID and skips `*_plus_*contam*` /
   `*decoy*` / `*predicted*` variants — appending contaminants to a database that
   already contains them would duplicate them.
3. **UniProt.**

### Database type (`--content`, default `one_per_gene`)
| value | what you get | human size |
|---|---|---|
| `one_per_gene` | canonical, one protein per gene — **the default** | 20,652 |
| `reviewed` | Swiss-Prot only | ~20,400 |
| `full` | + unreviewed TrEMBL | **147,506** |
| `*_isoforms` | + splice isoforms | larger still |

**`one_per_gene` only comes from the reference-proteome FTP tree**
(`{Kingdom}/{UPID}/{UPID}_{TAXID}.fasta.gz`), because UniProt's REST
`&onePerGene=true` is **silently ignored** — verified 2026-07-29: the yeast stream
returns byte-identical output (6,067 entries) with and without it. A plain REST
`(proteome:X)` stream is therefore always the *full* set. This is why DE-LIMP
(`R/helpers_search.R`) uses FTP, and why the Core's staged HIVE database is
`UP000005640_9606.fasta` = 20,663 sequences, not 147k.

If no FTP file exists (non-reference proteome), the script warns loudly, records
the warning in its output, and falls back to the REST full set — it never swaps
databases silently. If REST also fails it exits with the reason.

### Contaminants (`--contaminants`, default `universal`)
Sets: `universal` (default; what the Core stages on HIVE), `cell_culture`,
`mouse_tissue`, `rat_tissue`, `neuron_culture`, `stem_cell_culture`, `none`.
Source order: `--contaminants-path` → HIVE (matching the *requested* set) →
a DE-LIMP checkout's `contaminants/` → the Hao lab GitHub repo (public;
Frankenfield et al. 2022, JPR 21(9):2104-2113, doi:10.1021/acs.jproteome.2c00145).

Headers are `Cont_`-tagged, so DIA-NN's `--cont-quant-exclude Cont_` keeps them
out of quantification and normalisation. The DIA-NN **Linux binary ships no
contaminant FASTA of its own** (the GUI's "Contaminants" checkbox is a
Windows-side asset), so they must be appended here. `fetch_fasta.py` reports the
tag as `diann_cont_quant_exclude`; pass the sidecar to `estimate_params.py
--fasta-meta` and the flag lands in the cfg automatically.

**A failure to fetch contaminants is fatal, not a warning** — the GPM cRAP URL
this script used previously now 404s, and the old warn-and-continue behaviour
meant searches silently ran with no contaminants at all, which also made the
contaminant-dominance QC check meaningless. Override deliberately with
`--contaminants none` or `--allow-missing-contaminants`.

### Output
Refuses to proceed on 0 sequences, and warns if a full-set download comes back
>5% short of UniProt's declared count (truncated stream). Writes
`<out>.meta.json` alongside the FASTA — sha256, source URL, organism, taxid,
content type, UniProt release, per-part sequence counts, contaminant set +
citation, and any warnings. Pass it to `provenance.py --fasta-info` so
`reproduce.sh` rebuilds the database that was *actually searched*.

## SLURM submission (hpc)

`run_search.py --sbatch job.sh` emits a login-node-safe script:
`--partition=high --qos=genome-center-grp-high-qos`, 64G, 12h. If that queue is
full, the DE-LIMP fallback is `publicgrp/low` (the only fallback — there is no
`genome-center-grp` LOW partition). Submit with `sbatch job.sh`, poll the
`<job>_<id>.log`, then run `run_search.py --adapt-only` for Sage/FragPipe to build
`report.parquet`.
