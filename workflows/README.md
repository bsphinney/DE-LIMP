# Workflow registry — RETIRED (2026-08-14)

This directory used to hold per-organism "validated workflow" bundles that the
`ucdavis-proteomics-core-pipeline` skill fetched over the network at run time.
It no longer does. **Search parameters now ship with the skill.**

Everything that used to live here is now:

| Was | Is now |
|---|---|
| `workflows/<id>/workflow.yaml` | `scripts/resolve_defaults.py` — one defaults table keyed on data type |
| Instrument → mass accuracy | `scripts/estimate_params.py` — already the only real source of these numbers |
| `DIA_SpecLib_Quant*.workflow`, `default.radiantConfig` | `presets/` templates + `scripts/make_presets.py`, which writes the config per run |
| `fetch_workflows.py match` / `pull` | `resolve_defaults.py` (same `workflow.manifest.json` output contract) |

## Why

**Species never affected a search parameter.** Only the FASTA depends on the
organism, and that always came from the user's own organism answer — never from
the bundle. The organism dimension bought nothing and forced a near-duplicate
bundle per species, so a mouse timsTOF run and a human timsTOF run needed two
entries that differed only in a FASTA ID nobody read.

**A remote, mutable registry broke reproducibility.** One skill version could
produce different parameters on different days depending on what `main` said. It
did: two bundles pinned two different DIA-NN versions (2.3.0 and 2.6.0), and a
run silently got whichever one its instrument happened to match. Parameters that
ship with the skill mean the skill version fully determines the search.

**It added a network dependency that can fail for unrelated reasons.** Every run
hit the GitHub API. Unauthenticated calls rate-limit at 60/hour per IP, so on a
shared network a first run could fail with nothing wrong with the data.

**The validation gate fired on everything.** Every bundle was unvalidated or
placeholder-validated, so the "is this validated?" check tripped on essentially
every run and turned a first interaction into a menu of ⚠ warnings.

## Site-specific parameters

Overriding a default is now an explicit flag rather than a repo edit:

```bash
python3 scripts/resolve_defaults.py --acquisition DIA --instrument "timsTOF HT" \
    --ms1-ppm 12 --ms2-ppm 12 --dest ./wf
```

Overrides are tagged as such in `workflow.manifest.json`, so a run record always
distinguishes a site SOP value from the shipped default.

## Old skill installs

`index.json` is kept as a tombstone with an empty `workflows` list so an
un-updated skill degrades to "no workflow matched → estimate parameters" instead
of crashing on a 404. Update the skill and it stops looking here at all.
