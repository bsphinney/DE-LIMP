# The no-DE route: an analysis with no differential expression

SKILL.md points here from step 3 (the design questions), step 4 and the top of step 8.

## When it applies
The session has no study design to test, or the question is not a comparison of groups:
- a detection question ("is protein X there?"), an identification list, a species check;
- a QC run, or one sample per condition;
- a **re-analysis** of a published or earlier dataset, a **benchmark**, or a comparison
  against a published result that asks its own question (not a volcano of two groups).

Say so in step 3, confirm it with the defaults in step 4, and record it
(`log_decision.py`). Do not file the steps after step 8 as "after the DE / later": they are
the finish of THIS analysis. Three sessions that did ended with a hand-made report, a hand-made
reproducibility zip and no AGENTS.md, README, bundle or finalize, and one delivered raw files
and search results with its own copy script, past the QC check.

## The route, in order
1. **The search (step 7). One rule: a no-DE DDA search with Sage is identifications only
   (`estimate_params.py --ids-only`, step 6b) unless the user -- for a Core submission, the
   client -- wants amounts.** Sage then writes `identifications.tsv` and no `report.parquet`. On
   a timsTOF `.d` it reads the `.d` itself, which identified more than the mzML conversion a
   quantifying search needs (896 peptides against 753 on two timsTOF HT DDA-PASEF runs,
   `references/search-engines.md`).
   - **Quantify instead** (no `--ids-only`) when the question needs amounts: abundance compared
     between runs or samples, an intensity table in the delivery, or a DE that may follow
     later. An identifications-only search cannot feed a DE; it would be searched again.
   - A detection question ("is protein X there, and on what evidence?"), an identification
     list or a species check is identifications only.
   - When the request does not say, ask once. Record the choice and why with
     `log_decision.py`.
   - A DIA search (DIA-NN) always quantifies (`--ids-only` is refused with DIA-NN). That is
     fine: there is just no step 8.
2. **Skip steps 8–8d and 8f** (the DE, set tests, DE figures, audit, expert review, peak check).
   **Core runs still do step 8e** (`qc_bracket.py --session <S>`): `deliver` holds a no-DE
   delivery on the QC check exactly as it holds an analysis delivery.
3. **The report (step 9):** `make_analysis_html.py` works without DE tables. With no DE record
   and no DE table, `analysis_prompt.py` writes the **no-DE brief**, titled "identification and
   QC (no differential expression)". Its body is written for that session, not only its title.
   Pass `--session <S>`, as step 9 says: it finds the session's search (`session.search_folder`;
   else `--report`'s folder).
   - **It names the search's own files to read.** For an identifications-only search:
     `identifications.tsv` and its record, `sage_identifications.json` or
     `fragpipe_identifications.json` (the rule, the protein-level statement, the row unit). For a
     quantifying one: the search report, `report.parquet` and `report.stats.tsv`. For both:
     `search_provenance.json`. It says "identifications only: no quantities" when the search ran
     none, and carries the first-pass caveat of a DIA-NN `--dda` fallback.
   - **It makes no DE demands:** no significance rule, contrasts, normalisation check, volcano,
     DE tables or DE Methods. Its sections are Overview, Identification depth per run, QC
     Assessment, Notable proteins (presence only, each with its evidence per run), How This
     Analysis Works, Data Quality Notes, Methods and Next steps (plus Instrument & Acquisition
     when `--instrument` is given).
   - **With two or more conditions it carries the caveat** that differences between them were
     not tested. With one run per condition (n = 1) they cannot be separated from run-to-run
     variation; in DDA, a protein missing from one run is not shown to be absent. The README and
     AGENTS.md say it in the same words (`session_docs.no_de_condition_caveat`).

   Write about what was identified and how well. Never write about fold changes, or about a
   protein "unique to" a condition.
4. **Methods (step 9d):** `make_methods.py` without `--de-dir`, or leave it to finalize.
5. **The podcast (step 9e)**, for a CoreOmics submission (by default).
6. **The reproducibility bundle (step 10):** `provenance.py` as in step 10, or leave it to
   finalize, which builds it from the session's records when it is missing or out of date.
   **The file catalogue (step 11):** in the session, `make_report.py --out
   output/OUTPUT_FILES.md --session . --repro output/reproducibility --extra
   input/conditions.csv output/figures output/AUDIT.md output/Analysis_Report.html
   output/methods.md` (without `--de-dir`). `--session` catalogues the session's search folder:
   with no DE, the search output IS the result, and without `--session` (or `--search-out`) no
   search folder is listed. finalize only adds to an `OUTPUT_FILES.md` that exists, and step 12
   points the user at it. A benchmark of two searches of the same files runs **11a**
   (`compare_searches.py`): that comparison is the answer to its question.
7. **Finalize (step 12):** `session.py finalize --dir <S> --zip`. It is the last step of
   **every** analysis, with or without a DE: it writes AGENTS.md, README.md/README.html,
   MANIFEST.txt (each DE-only part `[SKIPPED]` with the reason), the Methods and the bundle when
   they are missing, and the finalize stamp. Then `session.py status --session <S>`, and show
   its result before you say the analysis is finished.
8. **Deliver (step 12b), Core submissions only: `deliver --mode no-de`**
   (`core_submission.py deliver --summary <...> --session "$S" --mode no-de`, dry run, then
   `--apply`). It needs the report, the Methods, the bundle and the finalize stamp, and refuses
   a session that has a DE record or DE files. Internal sessions stop at step 12.
   - Raw files and the search results **before** the report: `deliver --mode interim`
     (`core_submission.py deliver ... --mode interim`). It needs no report or finalize, and its
     `owed` record keeps `session.py status` open until a later `--mode no-de` (or analysis)
     delivery. Never copy search folders to a share by hand.
   - **Two or more searches** (a benchmark, 11a): with no DE to choose by, `--mode no-de`
     carries **no** search folder and says why, and the Methods' Database-search paragraph is
     tagged blanks (`make_methods.py --session`: none chosen). Deliver **`--mode interim`
     first** (every search folder, each under its own name), **then `--mode no-de`** (the
     report, Methods and bundle), which closes `owed`. To deliver one search only, rename each
     other folder to `output/superseded_<name>/` before finalize (SKILL.md 12b, "Several
     searches").
   - "I only require raw data": `--mode raw-only`.

## Work done outside the skill
- Analysis code you write outside the session (a species-composition script, a notebook) goes
  in with `session.py finalize --code-dir <dir>` (repeatable; `--extra-dir <dir>` for other
  files). It is copied to `output/custom/`, checksummed, and listed in the README and AGENTS.md.
- A project's own reproducibility package (for example a CLAUDE.md "collaborator package") is
  IN ADDITION to the skill's bundle, never instead of it.
- A README.md or AGENTS.md written by hand is never overwritten: finalize writes
  `README.skill.md` / `AGENTS.skill.md` beside it.
