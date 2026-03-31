# Auto-Queue Switching Logic

## Overview

DE-LIMP automatically moves SLURM jobs between partitions to optimize resource usage:
- **genome-center-grp/high** → priority queue, per-user CPU limit (64 CPUs)
- **publicgrp/low** → preemptible, large capacity (1000+ idle CPUs), jobs can be killed and requeued

## Decision Flow

```
Job submitted on genome-center-grp/high
  ↓
Monitoring observer polls every 15 seconds
  ↓
After N minutes pending (default: 5 min):
  → Query sinfo -p low for idle CPUs
  → If publicgrp has idle CPUs → move array steps (2, 4) to low
  → Assembly steps (1, 3, 5) stay on high (not preemptible)
  ↓
If job pending reason is InvalidQOS or QOSMaxCpuPerUserLimit:
  → Move ALL pending steps (force_move_all)
  ↓
When genome-center-grp has capacity again:
  → Move queued publicgrp jobs back to high
```

## What Gets Moved

| Step | Type | Moved to low? | Reason |
|------|------|---------------|--------|
| Step 1 | Library prediction | Only if nothing started yet | Single job, fast |
| Step 2 | Per-file quant (array) | YES | Embarrassingly parallel, safe to preempt |
| Step 3 | Empirical library assembly | NO | Single job, uses all quant files, can't restart mid-way |
| Step 4 | Final quant (array) | YES | Same as Step 2 |
| Step 5 | Cross-run report | NO | Single job, needs all Step 4 output |

## SLURM Attributes Set on Move

When `slurm_move_job()` moves a job:
1. **Account** → `publicgrp` (or `genome-center-grp` on move-back)
2. **Partition** → `low` (or `high`)
3. **QOS** → `{account}-{partition}-qos` (e.g., `publicgrp-low-qos`)
4. **Requeue=1** → only when moving to `low` (preemptible — auto-restart on preemption)

## Known Issues & Fixes (March 2026)

### Preemption on Low Partition
- Jobs on low can be killed (PREEMPTED) when higher-priority users need nodes
- `Requeue=1` makes SLURM auto-restart preempted tasks
- `check_slurm_status()` maps PREEMPTED → "queued" (not "failed")

### Partial Move Tracking
- When only steps 2 & 4 are moved, `slurm_account` stays "genome-center-grp"
- `partially_on_public` flag tracks these split-partition jobs
- Move-back logic checks both `slurm_account == "publicgrp"` AND `partially_on_public == TRUE`

### Retry Dependency Chain (Critical Fix)
When step 2 tasks fail on low and are retried:
1. Retry creates a NEW SLURM job ID for the failed tasks
2. Step 3 had `--dependency=afterany:ORIGINAL_STEP2_ID`
3. Original step 2 was already "complete" (some tasks succeeded, some failed)
4. Step 3 starts immediately → quant_verify fails (retry .quant files don't exist yet)

**Fix**: After submitting retry, use `scontrol update` to change Step 3's dependency to `afterany:RETRY_STEP2_ID`. Same for Step 4 retries → Step 5.

### pending_reason Detection
- Must fetch pending reason for the first QUEUED step, not the first non-completed step
- If current step is RUNNING, it has no pending reason
- `QOSMaxCpuPerUserLimit` triggers force_move_all (not just `InvalidQOS`)

## Key Functions

| Function | File | Purpose |
|----------|------|---------|
| `slurm_move_job()` | helpers_search.R:3123 | scontrol update Account/Partition/QOS/Requeue |
| `check_slurm_status()` | helpers_search.R:1961 | Maps SLURM states to app states |
| `select_best_partition()` | helpers_search.R:3158 | Decides initial partition at submission |
| `check_cluster_resources()` | helpers_search.R:2926 | Queries sacctmgr/squeue/sinfo for capacity |
| Auto-switch observer | server_search.R:5735 | Monitors pending jobs and moves them |
| Move-back observer | server_search.R:5862 | Returns jobs to high when capacity available |

## SLURM State Mapping

| SLURM State | App State | Notes |
|-------------|-----------|-------|
| PENDING | queued | |
| RUNNING, COMPLETING | running | |
| COMPLETED | completed | |
| FAILED, TIMEOUT, OUT_OF_MEMORY | failed | |
| CANCELLED | cancelled | |
| PREEMPTED, REQUEUED | queued | Job will auto-restart if Requeue=1 |
| NODE_FAIL, BOOT_FAIL | failed | Infrastructure failure |
| SUSPENDED, STOPPED | running | Can be resumed |

## Debugging

Check the console log for these messages:
```
[Auto-queue] pub_idle=1206 CPUs, wait_min=5
[DE-LIMP] Auto-switched 'Bovine_Liver_' [step2, step4] to publicgrp/low after 5 min pending
[DE-LIMP] Updated Step 3 (12345) dependency to wait for retry Step 2 (12346)
[DE-LIMP] Partial retry: step 2 tasks [25,27] with 128 GB (was 96 GB)
```

Queue switch events are also appended to `search_info.md` in the output directory.
