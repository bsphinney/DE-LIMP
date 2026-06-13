# DE-LIMP PG Farm setup — runbook for a Claude session

> **Goal.** Stand up non-interactive (service-account) access to DE-LIMP's
> central Postgres on UC Davis Library's **PG Farm**, the same way STAN's was
> wired. After this, DE-LIMP code + the `scripts/import_coreomics_*.py`
> importers can read/write `uc-davis-genome-center-proteomics-core/delimp`
> with no human login.
>
> **Status when this doc was written (2026-06-13):** database not yet created,
> service account not yet provisioned, `scripts/migrate_pg_v1.sql` staged but
> never run. PG Farm admin contact: **Justin Merz** (UC Davis Library).
>
> The pattern mirrors STAN exactly — STAN's equivalent doc is
> `~/Documents/STAN/docs/PG_FARM.md` + `docs/PG_FARM_ACCESS.md`. Read those if
> anything here is unclear; the auth model is identical.

---

## Do you need a SEPARATE service account? (the actual question)

**Short answer: technically no, but you probably want one.**

PG Farm service accounts are **organization-scoped**, and the token they mint
is a JWT that identifies the *account*, not a database. Postgres then enforces
access per-database via GRANTs. Because `delimp` lives in the **same org** as
`stan` (`uc-davis-genome-center-proteomics-core`), the existing
`genome-proteomics-service-account` *could* simply be granted on the `delimp`
database and reused — same secret, same token, one less thing to manage.

**But the recommendation is a separate account** (e.g.
`genome-proteomics-delimp-service-account`), because:

- **Blast radius / least privilege.** `delimp` holds the *internal* layer
  (customer/facility data per `docs/HISTORY_DB_DESIGN_V2.md`); STAN holds only
  aggregate QC. A leaked DE-LIMP token must not also open the STAN DB, and
  vice-versa.
- **The repo already assumes it.** `scripts/import_coreomics_api.py` reads its
  password from a **separate** token file, `.pgfarm_delimp_token` — distinct
  from STAN's `.pgfarm_token`. Keeping a distinct account behind that distinct
  file is the clean version of what the code already expects.
- **Independent rotation.** Rotating one secret doesn't disrupt the other
  system's dispatch/cron.

The two paths below are labelled **[SEPARATE]** (recommended) and
**[SHARED]** (the shortcut). Pick one before you start. Everything else is the
same.

---

## 0. Prerequisites

- PG Farm UI access (Brett's CAS login) to the
  `uc-davis-genome-center-proteomics-core` organization.
- Justin Merz available to create the database and/or the service account if
  the UI doesn't expose it to Brett directly.
- The shared Quobyte FS mounted (Mac: `/Volumes/proteomics-grp/...`,
  Hive: `/quobyte/proteomics-grp/...` — same bytes, two mount points).
- `psycopg2` available in whatever Python you run (Mac: `/opt/anaconda3/bin/python3`;
  Hive: `/quobyte/proteomics-grp/brett/stan_venv/bin/python`). For R, the
  `RPostgres` package.

---

## 1. Create the database on PG Farm

In the PG Farm UI (or ask Justin):

1. Organization: `uc-davis-genome-center-proteomics-core`.
2. Create database **`delimp`** → full name
   `uc-davis-genome-center-proteomics-core/delimp`.
3. Confirm Brett's CAS user is the **owner** (you need owner to run the
   migration + grants).

---

## 2. Provision the service account

### [SEPARATE] — recommended

1. In the PG Farm UI for the org, create a service account, e.g.
   **`genome-proteomics-delimp-service-account`** (if the UI can't, Justin
   creates it).
2. Add that service account as a **user on the `delimp` database**.
3. Click **"rotate"** to download `service-account.json`
   (`{username, secret}`, ~512-char secret). This is the only time you see the
   secret — admins can't read it afterward.

### [SHARED] — shortcut (reuse STAN's account)

1. Add the existing **`genome-proteomics-service-account`** as a user on the
   `delimp` database (UI or Justin).
2. No new secret needed — you already have it at `.pgfarm_secret.json`. Skip
   to step 4 (grants) using that account name.

---

## 3. Store the secret + wire token refresh

Store the secret chmod 600 on the shared FS (both mount points see the same
file):

```
# [SEPARATE]
/quobyte/proteomics-grp/brett/.pgfarm_delimp_secret.json   # Hive
/Volumes/proteomics-grp/brett/.pgfarm_delimp_secret.json   # Mac (same file)
```

```bash
chmod 600 /Volumes/proteomics-grp/brett/.pgfarm_delimp_secret.json
```

Mint the 7-day token. **Reuse STAN's refresher** — it's generic, just point it
at the DE-LIMP files (the script lives in the STAN repo):

```bash
python ~/Documents/STAN/scripts/pgfarm_refresh_token.py \
    --secret-file /Volumes/proteomics-grp/brett/.pgfarm_delimp_secret.json \
    --token-file  /Volumes/proteomics-grp/brett/.pgfarm_delimp_token \
    --max-age-days 0          # 0 = mint now
```

That writes the token (chmod 600) to `.pgfarm_delimp_token` — exactly the path
`scripts/import_coreomics_api.py` already reads.

**[SHARED] variant:** skip the new secret; point `--secret-file` at the
existing `.pgfarm_secret.json` but still write a DE-LIMP-named token file if you
want the import scripts to keep using `.pgfarm_delimp_token`. (Or just symlink
`.pgfarm_delimp_token -> .pgfarm_token`.)

### Auto-refresh cron

Add a tick to the same Hive cron that refreshes STAN's token (it already runs
`pgfarm_refresh_token.py`). One extra invocation with the DE-LIMP paths and
`--max-age-days 5`. **Do not put cron on a Hive login node** (STAN got burned
by exactly that — see STAN's `docs/PG_FARM.md` history section). Use the
established cron host / launchd agent.

---

## 4. Run the migration + grants (as the OWNER, `brettsp`)

The schema is already written: `scripts/migrate_pg_v1.sql`. Run it **once**
against the fresh DB. First grant yourself schema privileges if needed (the SQL
header notes this).

```bash
PGPASSWORD=$(cat /Volumes/proteomics-grp/brett/.pgfarm_token) \
psql "host=pgfarm.library.ucdavis.edu port=5432 \
      dbname=uc-davis-genome-center-proteomics-core/delimp \
      user=genome-proteomics-service-account sslmode=require" \
  -c "SELECT 1;"   # sanity: can the OWNER connect? (use brettsp creds if SA isn't owner)
```

> Connect as the **owner** to run DDL. If Brett's CAS user is owner, use those
> creds for the migration; the service account only needs DML grants.

```bash
# as owner
psql "<owner connection string for .../delimp>" -f ~/Documents/claude/scripts/migrate_pg_v1.sql
```

Then grant the service account DML on everything the app touches, and make
future tables auto-granted (mirrors what STAN did once):

```sql
-- run as owner, against .../delimp
GRANT USAGE ON SCHEMA public TO "genome-proteomics-delimp-service-account";
GRANT SELECT, INSERT, UPDATE, DELETE
  ON ALL TABLES IN SCHEMA public
  TO "genome-proteomics-delimp-service-account";
GRANT USAGE, SELECT ON ALL SEQUENCES IN SCHEMA public
  TO "genome-proteomics-delimp-service-account";
ALTER DEFAULT PRIVILEGES IN SCHEMA public
  GRANT SELECT, INSERT, UPDATE, DELETE ON TABLES
  TO "genome-proteomics-delimp-service-account";
ALTER DEFAULT PRIVILEGES IN SCHEMA public
  GRANT USAGE, SELECT ON SEQUENCES
  TO "genome-proteomics-delimp-service-account";
```

(For **[SHARED]**, swap the account name to
`genome-proteomics-service-account`.)

---

## 5. Connect — copy/paste

Connection params (identical to STAN except the database name and the user if
you went [SEPARATE]):

```
host      pgfarm.library.ucdavis.edu
port      5432
database  uc-davis-genome-center-proteomics-core/delimp
user      genome-proteomics-delimp-service-account   # or ...-service-account if [SHARED]
sslmode   require            # NOT verify-full — cert path breaks on Mac
password  = contents of .pgfarm_delimp_token
```

Driver: **psycopg2** (Python) / **RPostgres** (R). NOT psycopg v3.

### R (DBI + RPostgres) — for the Shiny app

```r
library(DBI)
pw <- trimws(readLines("/Volumes/proteomics-grp/brett/.pgfarm_delimp_token", warn = FALSE)[1])
con <- dbConnect(
  RPostgres::Postgres(),
  host     = "pgfarm.library.ucdavis.edu",
  port     = 5432,
  dbname   = "uc-davis-genome-center-proteomics-core/delimp",
  user     = "genome-proteomics-delimp-service-account",
  password = pw,
  sslmode  = "require"
)
dbGetQuery(con, "SELECT version FROM delimp_schema_version;")
dbDisconnect(con)
```

### Python (psycopg2) — for the importers

```python
import psycopg2
pwd = open('/Volumes/proteomics-grp/brett/.pgfarm_delimp_token').read().strip()
with psycopg2.connect(
    host='pgfarm.library.ucdavis.edu', port=5432,
    database='uc-davis-genome-center-proteomics-core/delimp',
    sslmode='require',
    user='genome-proteomics-delimp-service-account', password=pwd,
) as c:
    cur = c.cursor()
    cur.execute("SELECT version FROM delimp_schema_version")
    print(cur.fetchall())
```

The existing importers already read `PGPASSWORD` from
`.pgfarm_delimp_token` — so just:

```bash
PGPASSWORD=$(cat /Volumes/proteomics-grp/brett/.pgfarm_delimp_token) \
python scripts/import_coreomics_api.py --database uc-davis-genome-center-proteomics-core/delimp ...
```

> Note: `import_coreomics_api.py` uses **two** tokens — `COREOMICS_TOKEN`
> (the Coreomics REST API) and `PGPASSWORD`/`.pgfarm_delimp_token` (the DB).
> They're unrelated; this doc only covers the PG Farm DB token.

---

## 6. Verify

```sql
-- service account can connect and see the schema
SELECT current_user, current_database();
SELECT version FROM delimp_schema_version;          -- should be '1.0.0'
SELECT count(*) FROM information_schema.tables
  WHERE table_schema = 'public';                    -- schema present
-- write check (then roll back)
BEGIN; INSERT INTO delimp_schema_version (version, notes)
  VALUES ('test', 'access check'); ROLLBACK;
```

If the write check succeeds inside the transaction, grants are correct.

---

## 7. Update STAN's reference (optional but tidy)

If you went [SEPARATE], add a one-liner to STAN's
`~/Documents/STAN/docs/PG_FARM.md` noting that `delimp` uses its own
`genome-proteomics-delimp-service-account` + `.pgfarm_delimp_token`, so a
future Claude doesn't assume one account spans both DBs.

---

## Gotchas (same family as STAN's)

- **`sslmode=verify-full` is broken on Mac** — use `require`.
- **psycopg2 not in system/cron Python** — use the venv interpreter explicitly.
- **No compute / cron on Hive login nodes** — SLURM or the designated cron host.
- **Don't open a connection per row** in loops — reuse one connection.
- **Never commit the secret or token** to git.
- **Token age ≠ DB grant.** A fresh token still fails if the account was never
  GRANTed on `delimp` (step 4) — that's the most likely first-run error.
```
