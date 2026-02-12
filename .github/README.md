# GitHub Actions Workflows

## Sync to Hugging Face (`sync-to-hf.yml`)

**Purpose**: Automatically sync code changes from GitHub to Hugging Face Spaces while properly managing platform-specific README files.

### How It Works

1. **Trigger**: Runs on every push to `main` branch (except changes to README.md or .github)
2. **Process**:
   - Checks out GitHub repository
   - Clones Hugging Face repository
   - Syncs all files **except** README.md and .github
   - Copies `README_HF.md` → `README.md` (with YAML frontmatter for HF)
   - Commits and pushes to Hugging Face

### Setup Required

**⚠️ IMPORTANT**: Add your Hugging Face token as a GitHub secret:

1. Go to https://github.com/bsphinney/DE-LIMP/settings/secrets/actions
2. Click "New repository secret"
3. Name: `HF_TOKEN`
4. Value: Your Hugging Face write token from https://huggingface.co/settings/tokens
5. Click "Add secret"

### What This Solves

**Problem**: The two repositories need different README.md files:
- GitHub: Documentation-focused (no YAML)
- Hugging Face: Configuration-focused (requires YAML frontmatter)

**Previous Solution**: Manual sync with recovery script (failed 9 times)

**New Solution**: Automated sync that:
- ✅ Never overwrites HF README with GitHub README
- ✅ Runs automatically on every push
- ✅ Eliminates manual "git push hf main" step
- ✅ Prevents README conflicts entirely

### Workflow After Setup

**New simplified workflow:**
```bash
# 1. Make changes locally
# 2. Commit to GitHub
git add DE-LIMP.R app.R Dockerfile  # or whatever files
git commit -m "Your changes"
git push origin main

# 3. GitHub Action automatically syncs to HF!
# (No manual push to hf remote needed)
```

### Monitoring

- Check workflow status: https://github.com/bsphinney/DE-LIMP/actions
- View HF build logs: https://huggingface.co/spaces/brettsp/de-limp-proteomics/logs

### Notes

- The workflow ignores changes to README.md and .github to prevent loops
- Uses `rsync` to copy files, preserving structure and timestamps
- Hugging Face builds still take 5-10 minutes (or 30-45 for Dockerfile changes)
- The HF remote in your local repo is no longer needed (can remove it)
