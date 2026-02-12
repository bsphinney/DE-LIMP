# DE-LIMP Project Context for Claude

## Working Preferences (for Claude)
**IMPORTANT**: Proactively suggest updates to this file when you notice:
- New important patterns or gotchas that emerge during work
- Solutions to non-trivial problems that should be documented
- Architectural decisions or implementation details worth preserving
- Common mistakes or confusing aspects that need clarification

Don't wait to be asked - suggest adding these to CLAUDE.md to help future sessions!

## Project Overview
DE-LIMP is a Shiny proteomics data analysis pipeline using the LIMPA R package for differential expression analysis of DIA-NN data.

**GitHub Repository**: https://github.com/bsphinney/DE-LIMP

## System Requirements
**CRITICAL**: limpa package requires specific R/Bioconductor versions:
- **R version**: 4.5 or newer (required)
- **Bioconductor**: 3.22+ (automatically set with R 4.5+)
- **Installation**: Download R from https://cloud.r-project.org/

**Note**: The app will automatically detect if R version is too old and provide upgrade instructions.

## Key Files
- **DE-LIMP.R** - Main Shiny app (1923 lines) - For GitHub releases and local users
- **app.R** - Copy of DE-LIMP.R for Hugging Face Spaces (HF requires this naming)
- **README.md** - Full documentation for GitHub
- **README_HF.md** - Hugging Face version with YAML frontmatter (source for HF's README.md)
- **Main URL**: http://localhost:3838 when running locally

## Deployment & Release Management

### Dual Remote Setup
This project is deployed to **two platforms** with separate git remotes:
- **`origin`** → GitHub (https://github.com/bsphinney/DE-LIMP)
- **`hf`** → Hugging Face Spaces (https://huggingface.co/spaces/brettsp/de-limp-proteomics)

**IMPORTANT**: Changes to GitHub do NOT automatically sync to Hugging Face. You must push to both remotes manually.

**CRITICAL WARNING**: The two remotes have **different README.md files**:
- GitHub: Full documentation (no YAML)
- Hugging Face: YAML frontmatter required for Spaces configuration

**⚠️ DO NOT push README.md to both remotes!** If you push to both, you'll overwrite one of them and break HF configuration.

### 🚨 CRITICAL: Avoiding README Conflicts

**THE PROBLEM**: Every time you push to HF, if README.md is in your recent commits OR working directory, it will overwrite HF's YAML README and break deployment.

**THE SOLUTION**: When making HF-specific changes (Dockerfile, etc.), follow this EXACT workflow:

```bash
# 1. Make your Dockerfile changes
# Edit Dockerfile...

# 2. Add ONLY the specific file (NOT README.md!)
git add Dockerfile  # Be specific! Don't use git add .

# 3. Commit
git commit -m "Fix Dockerfile"

# 4. Push to HF only
git push hf main

# 5. IF you get "Missing configuration in README" error:
#    Run the recovery script below
```

**Recovery Script** (copy-paste this when README breaks):
```bash
# Quick fix for broken HF README
cp README.md README_GITHUB_BACKUP.md
cp README_HF.md README.md
git add README.md
git commit -m "Fix HF: Restore README with YAML frontmatter"
git push hf main
cp README_GITHUB_BACKUP.md README.md
git add README.md
git commit -m "Restore GitHub README"
git push origin main
rm README_GITHUB_BACKUP.md
```

**WHY THIS HAPPENS**:
- Git pushes ALL commits in your branch history, not just the files you changed
- If any recent commit includes README.md, it gets pushed to HF
- HF needs YAML frontmatter in README.md, GitHub doesn't
- When GitHub's README replaces HF's README → deployment breaks

**PREVENTION RULES**:
1. ✅ **DO**: Use `git add <specific-file>` before pushing to both remotes
2. ❌ **DON'T**: Use `git add .` or `git add -A` before pushing to both remotes
3. ✅ **DO**: Keep README.md in separate commits for each remote
4. ❌ **DON'T**: Include README.md in commits that go to both remotes

### When to Update Each Platform

#### For Code Changes (app.R, CLAUDE.md, etc.):
```bash
# Make your changes to app.R, CLAUDE.md, or other shared files
# ⚠️ DO NOT include README.md in this commit!
git add app.R CLAUDE.md
git commit -m "Description of changes"

# Push to BOTH remotes
git push origin main && git push hf main
```

**⚠️ WARNING**: Never use `git add .` or `git add -A` before pushing to both remotes! This will include README.md and break HF configuration.

#### For GitHub-Only Changes (DE-LIMP.R, full README.md):
```bash
# Changes to DE-LIMP.R or documentation
git add DE-LIMP.R README.md
git commit -m "Description"

# Push to GitHub only
git push origin main
```

#### For Hugging Face-Only Changes (Dockerfile, README YAML):
```bash
# Make changes to Dockerfile or README_HF.md
git add Dockerfile
git commit -m "Description"

# Push to HF only
git push hf main
```

### Key Differences Between Platforms

| File | GitHub | Hugging Face | Notes |
|------|--------|--------------|-------|
| Main app | `DE-LIMP.R` | `app.R` | HF requires `app.R` naming |
| README | Full documentation (no YAML) | YAML frontmatter required | HF needs config metadata |
| Deployment | Manual download | Docker container | HF auto-builds from Dockerfile |

### Release Checklist (for Major Versions)

When creating a new release (e.g., v2.0):

1. **Update Documentation**:
   - [ ] Update `CLAUDE.md` with new features/fixes
   - [ ] Update `USER_GUIDE.md` with new workflows
   - [ ] Update `README.md` with feature highlights
   - [ ] Update `README_HF.md` if HF deployment changes

2. **Sync Code**:
   - [ ] Ensure `app.R` and `DE-LIMP.R` are identical (or copy: `cp app.R DE-LIMP.R`)
   - [ ] Test app locally: `shiny::runApp('DE-LIMP.R', port=3838, launch.browser=TRUE)`

3. **Commit & Tag**:
   ```bash
   git add .
   git commit -m "Release vX.X: Description"
   git tag -a vX.X -m "Version X.X release"
   git push origin main
   git push origin vX.X
   ```

4. **Update Hugging Face**:
   ```bash
   # Push code changes to HF
   git push hf main

   # Update HF README with YAML if needed
   cp README_HF.md README.md
   git add README.md
   git commit -m "Update HF README with YAML frontmatter"
   git push hf main

   # Restore GitHub README
   git checkout HEAD~1 README.md  # or use README_GITHUB.md backup
   git add README.md
   git commit -m "Restore GitHub README"
   git push origin main
   ```

5. **Create GitHub Release**:
   - [ ] Go to https://github.com/bsphinney/DE-LIMP/releases/new
   - [ ] Select tag vX.X
   - [ ] Write comprehensive release notes (features, fixes, requirements, links)
   - [ ] Attach `DE-LIMP.R` file for download
   - [ ] Link to Hugging Face Space in release notes
   - [ ] Publish release

6. **Verify Deployments**:
   - [ ] GitHub release page shows correctly
   - [ ] Hugging Face Space builds successfully (check Logs tab)
   - [ ] Test both local download and HF web version

### If You Accidentally Break HF Configuration

If you get "Missing configuration in README" error on Hugging Face, it means the GitHub README overwrote the HF README. Fix it:

```bash
# 1. Backup current README
cp README.md README_GITHUB_BACKUP.md

# 2. Replace with HF version (has YAML)
cp README_HF.md README.md

# 3. Commit and push to HF ONLY
git add README.md
git commit -m "Fix HF: Restore README with YAML frontmatter"
git push hf main

# 4. Restore GitHub version and push to GitHub ONLY
cp README_GITHUB_BACKUP.md README.md
git add README.md
git commit -m "Restore GitHub README"
git push origin main
rm README_GITHUB_BACKUP.md
```

### Quick Commands Reference

```bash
# Check which remote you're on
git remote -v

# Push to both remotes (⚠️ ONLY for files other than README.md)
git add app.R CLAUDE.md  # Be specific!
git commit -m "Description"
git push origin main && git push hf main

# Check Hugging Face build status (in browser)
# Visit: https://huggingface.co/spaces/brettsp/de-limp-proteomics/logs

# Test local version
shiny::runApp('DE-LIMP.R', port=3838, launch.browser=TRUE)
```

## Recent Changes & Important Fixes

### 2026-02-11: v2.0 Release Preparation & HF Deployment Fixes

**CRITICAL LESSON LEARNED**: README.md management with dual remotes

**The Recurring Problem:**
HF deployment broke **multiple times** during v2.0 release because GitHub README (without YAML) kept overwriting HF README (with YAML frontmatter). This happened even when we:
1. Fixed Dockerfile and pushed to HF → README broke
2. Updated CLAUDE.md and pushed to both remotes → README broke again
3. Fixed README, then pushed any commit to HF → README broke again

**Root Cause - Git Commit History:**
Git doesn't just push the files you changed - it pushes your **entire branch history**. When we:
- Made commit A: Changed README.md (GitHub version)
- Made commit B: Changed CLAUDE.md
- Pushed commit B to HF → Git also pushed commit A (with README changes)

**Why This Is Tricky:**
- We documented the workflow correctly in CLAUDE.md
- We followed the "use `git add <specific-file>`" rule
- But the problem persisted because **previous commits** in the branch history included README changes
- Every push to HF brought those README commits along

**The Real Solution:**
1. ✅ **DO**: Keep README changes in SEPARATE commits for each remote
2. ✅ **DO**: After fixing HF README, never push those commits back to origin
3. ✅ **DO**: Accept that the two branches will diverge on README.md
4. ✅ **DO**: Use the recovery script IMMEDIATELY when "Missing configuration" error appears
5. ❌ **DON'T**: Push to both remotes if ANY recent commit includes README.md changes
6. ❌ **DON'T**: Try to keep both remotes in perfect sync - they MUST diverge on README.md

**Pattern That Works:**
```bash
# Making app changes (safe to push to both):
git add app.R CLAUDE.md  # Specific files only
git commit -m "Update app"
git push origin main
git push hf main  # OK if no README in recent commits

# After this, if README broke on HF:
# Run recovery script IMMEDIATELY
# Don't continue working until README is fixed
```

**Lesson:** The two remotes MUST have different README.md files permanently. This is not a bug - it's a feature of having platform-specific configurations. The branches will diverge on README.md and that's OK.

**Recovery Count:** Fixed README on HF 4+ times during v2.0 release session

1. **Updated USER_GUIDE.md for v2.0**
   - Updated prerequisites to R 4.5+ requirement
   - Added Load Example Data section (Option A vs Option B)
   - Documented streamlined "Assign Groups & Run Pipeline" workflow
   - Added customizable covariate support section
   - Added QC enhancements (group averages, fullscreen view)
   - Added Reproducibility & Code Export section
   - Added Education & Resources section
   - Added Save & Load Analysis Sessions section
   - Updated volcano plot section (raw p-values vs FDR)
   - Added multi-protein violin plot details
   - Added "Accessing DE-LIMP" section with Hugging Face link
   - Enhanced troubleshooting with 10 updated entries

2. **Created Comprehensive v2.0 Release**
   - Tagged v2.0 in git
   - Created GitHub release with detailed notes (12+ features)
   - Included system requirements and quick start
   - Links to Hugging Face Space, GitHub Pages, documentation
   - Attached DE-LIMP.R for local download

3. **Fixed Hugging Face Deployment Issues**
   - Fixed README.md configuration error (missing YAML frontmatter)
   - Separated README.md for GitHub vs Hugging Face
   - HF now has README with proper YAML metadata
   - GitHub keeps full documentation README
   - Both remotes working correctly

4. **File Management for Dual Deployment**
   - Created `DE-LIMP.R` for GitHub releases (copied from app.R)
   - `app.R` used exclusively for Hugging Face
   - Both files identical in content, different in purpose/location
   - Added to git on GitHub only (not pushed to HF)

### 2026-02-11: Added Save/Load Session Feature
1. **Save/Load Analysis Sessions as RDS** (sidebar UI + server handlers)
   - Feature: Save and Load buttons in sidebar under new "3. Session" section
   - **Save**: Downloads an .rds file containing all analysis state (data, results, settings)
   - **Load**: Opens modal with file picker and confirmation warning
   - Saved state includes: raw_data, metadata, fit, y_protein, dpc_fit, design, qc_stats, gsea_results, repro_log, covariate names, UI settings (contrast, logfc cutoff, q cutoff)
   - Load validates the session file has required fields before restoring
   - Restores UI state (contrast selector choices/selection, slider values)
   - Appends load event to reproducibility log with timestamp of original save
   - Uses `%||%` (null coalescing) for backwards compatibility with older session files
   - File naming: `DE-LIMP_session_YYYYMMDD_HHMMSS.rds`

### 2026-02-10: Added Reproducibility Log Download Feature
1. **Download Button for Reproducibility Log** (lines 460, 1085-1113)
   - Feature: Added download button in "Reproducibility > Code Log" tab
   - Downloads as timestamped .R file: `DE-LIMP_reproducibility_log_YYYYMMDD_HHMMSS.R`
   - Includes full analysis log + session info for complete reproducibility
   - Button styled with success class (green) and download icon for visibility

2. **Fixed CRAN Mirror Popup Issue** (line 7)
   - Problem: VS Code and non-interactive terminals couldn't display CRAN mirror selection popup
   - Solution: Set default CRAN mirror (`https://cloud.r-project.org`) before any package installations
   - Prevents "Please select a CRAN mirror" hanging during startup

3. **Fixed Missing Package Installation** (lines 92-95)
   - Problem: `readr`, `dplyr`, `ggplot2`, `rhandsontable`, `arrow`, `shiny`, `bslib` not in auto-install list
   - Error: "there is no package called 'readr'" after installation completed
   - Solution: Added all library() packages to required_pkgs list for auto-installation
   - Ensures complete package installation before loading any libraries

4. **Fixed Invalid Icon Warning** (line 415)
   - Problem: "chart-scatter" is not a valid Font Awesome icon name
   - Warning: "The `name` provided ('chart-scatter') does not correspond to a known icon"
   - Solution: Changed to "chart-line" (valid Font Awesome icon)

5. **Enhanced Methodology Documentation** (lines 467, 1116-1200)
   - Feature: Comprehensive, well-formatted methodology text in "Reproducibility > Methodology" tab
   - Changed from `textOutput` to `verbatimTextOutput` for proper formatting with line breaks
   - Covers: Data input, DPC-CN normalization, maxLFQ quantification, limma statistical framework
   - Includes: Detailed explanation of empirical Bayes, FDR correction, and references
   - Provides proper citations for limpa, limma, and DIA-NN
   - Publication-ready text with clear section breaks and bullet points

6. **Added Fullscreen View for QC Trend Plot** (lines 404, 937-977)
   - Feature: "🔍 View Fullscreen" button in QC Trends tab
   - Opens plot in a large modal dialog (700px height, extra-large width)
   - Better for viewing on small monitors or detailed inspection
   - Modal has interactive controls (zoom, pan, export via plotly toolbar)
   - Click outside modal or "Close" button to dismiss

7. **Added Group Average Lines to QC Trend Plot** (lines 937-977)
   - Feature: Dashed horizontal lines showing average value for each group
   - Lines span the full range of each group on the x-axis
   - Color-matched to their respective groups for easy identification
   - Helps quickly identify group-level trends and compare averages
   - Works in both regular view and fullscreen modal

8. **Streamlined Workflow: Run Pipeline from Assign Groups Modal** (lines 359-362, 649-780)
   - Improvement: Removed standalone "Run Pipeline" button from sidebar
   - "Assign Groups" button now labeled "Assign Groups & Run Pipeline"
   - Modal footer has "▶ Run Pipeline" button instead of "Save & Close"
   - Workflow: Upload data → Open modal → Assign groups → Click "Run Pipeline" (all in one flow)
   - Validates groups (need 2+ groups) before running pipeline
   - Automatically closes modal and navigates to QC Plots tab when complete
   - More intuitive: users assign groups and immediately run analysis without extra steps

9. **Added Customizable Multiple Covariate Support** (lines 563-920)
   - Feature: Three covariate columns: Batch, Covariate1, Covariate2
   - **Customizable names**: Text inputs let you rename Covariate1/2 to anything (e.g., "Sex", "Diet", "Age")
   - Column headers update in real-time based on your custom names
   - Individual checkboxes: Select which covariates to include in the model
   - Design matrix: Dynamically builds formula using your custom names
   - Formula examples:
     - No covariates: `~ 0 + groups`
     - One covariate: `~ 0 + groups + batch`
     - Custom names: `~ 0 + groups + batch + sex + diet` (if you named them "Sex" and "Diet")
   - Validation: Only includes covariates with >1 unique non-empty values
   - Smart logging: Reproducibility log uses your custom covariate names
   - Use cases: Batch effects, Sex, Age, Diet, Instrument, Time_Point, etc.
   - Fully flexible: Name and use covariates for any categorical variable

10. **Improved Volcano Plot to Show Raw P-Values** (lines 1047-1051, 1302, 1634-1643)
   - Change: Volcano plot y-axis now uses non-adjusted P.Value instead of adj.P.Val
   - Rationale: Best practice in proteomics - visualize raw p-values, determine significance by FDR
   - **Significance coloring unchanged**: Red points still indicate FDR-corrected significance (adj.P.Val < 0.05)
   - Y-axis label: Explicitly labeled as "-log10(P-Value)" for clarity
   - **Results table**: Now shows both P.Value and adj.P.Val columns
   - **Grid view**: Also includes both P.Value and adj.P.Val for consistency
   - Horizontal line at -log10(0.05) represents unadjusted p-value threshold
   - Users can see both raw and adjusted p-values to understand statistical stringency

### 2026-02-10: Fixed Package Installation for First-Time Users
1. **Fixed Installation Conflicts** (lines 7-82)
   - Problem: Auto-installation tried to update already-loaded packages → ggplot2 unload errors
   - Solution: Check for missing packages BEFORE loading any libraries
   - Uses `update = FALSE` to prevent conflicts with loaded packages

2. **Fixed limpa Installation & R Version Check** (lines 16-62)
   - Problem: limpa requires R 4.5+, but error messages were unclear
   - Solution:
     - Check R version and provide clear upgrade instructions
     - Try Bioconductor release, then devel as fallback
     - Show exact version requirements (R 4.5+, Bioc 3.22+)
   - Key insight: limpa has been in Bioconductor since 3.21 (R 4.5 requirement)

3. **Improved Error Messages** (lines 46-58, 82-94)
   - Clear step-by-step instructions for upgrading R
   - Shows current vs required versions
   - Links to download page and Bioconductor docs

### 2026-02-09: Major Updates
1. **Fixed Startup Issue** (line 70)
   - Problem: dplyr `summarise()` syntax error in `get_diann_stats_r()`
   - Solution: Added explicit braces and `.groups = 'drop'`

2. **Enhanced Reproducibility Logging** (lines 399-426, 832-853)
   - Cumulative logging (appends instead of overwrites)
   - Timestamps for every action
   - Logs: data upload, pipeline runs, contrast changes, GSEA, exports
   - Helper function: `add_to_log(action_name, code_lines)`

3. **Fixed Violin Plot Button** (lines 943-954)
   - Problem: Button didn't recognize table row selections
   - Solution: Added observer to sync `input$de_table_rows_selected` → `values$plot_selected_proteins`

4. **Fixed Multi-Select in Results Table** (lines 763-807)
   - Problem: Reactive loop - table called `volcano_data()` which re-rendered on each selection
   - Solution: Made table build its own data independently, breaking the reactive loop
   - Result: Ctrl+Click and Shift+Click now work properly for multi-select

5. **Multi-Protein Violin Plots** (lines 956-1004)
   - Feature: Violin plot button now shows ALL selected proteins, not just first one
   - Layout: 2-column grid with facet_wrap
   - Dynamic height: Adjusts based on number of proteins (200px per row)
   - Individual scales: Each protein gets its own Y-axis (`scales = "free_y"`)

## Architecture Notes

### Key Reactive Values
- `values$plot_selected_proteins` - Selected proteins from table/volcano (used by all viz)
- `values$fit` - limma fit object from DE analysis
- `values$y_protein` - Protein-level quantification matrix
- `values$repro_log` - Cumulative R code for reproducibility

### LIMPA Pipeline Flow
1. `readDIANN()` - Load DIA-NN parquet file
2. `dpcCN()` - Data Point Correspondence normalization
3. `dpcQuant()` - Peptide → protein quantification
4. `dpcDE()` - Differential expression model
5. `contrasts.fit()` + `eBayes()` - Pairwise comparisons

### Selection System
- Volcano plot click/box-select → updates `values$plot_selected_proteins`
- Table row selection → updates `values$plot_selected_proteins`
- Both sync bidirectionally for highlighting

## Development Workflow

### Running the App
**In VS Code R Terminal (recommended):**
```r
shiny::runApp('/Users/brettphinney/Documents/claude/DE-LIMP.r', port=3838, launch.browser=TRUE)
```

**From command line (background mode):**
```bash
Rscript -e "shiny::runApp('/Users/brettphinney/Documents/claude/DE-LIMP.r', port=3838)" &
```

**Important Notes:**
- **DO NOT use** `source()` to launch the app - it doesn't work properly in VS Code
- **Always use** `shiny::runApp()` instead
- Shiny apps don't hot-reload - must restart after every code change
- Stop the app with `Ctrl+C` in the terminal, or `pkill -f "DE-LIMP.r"` from command line

### Useful Commands
- **Check if app is running**: `lsof -i :3838`
- **Stop the app**: `pkill -f "DE-LIMP.r"`
- **Restart after changes**: Stop the app, then run the `shiny::runApp()` command again

## Key Patterns & Gotchas

### Package Installation & Loading
- **CRITICAL**: Check and install packages BEFORE any `library()` calls
- Attempting to install/update packages after they're loaded causes unload errors
- Pattern from this project (lines 16-81):
  ```r
  # 1. Check R/Bioc versions
  # 2. Install missing packages
  # 3. THEN load libraries with library() calls
  ```
- Use `update = FALSE` in `BiocManager::install()` to prevent updating already-loaded packages
- Use `upgrade = "never"` in `remotes::install_*()` for same reason

### R Shiny Reactivity
- **DT table row indices** refer to CURRENTLY DISPLAYED data (filtered or not)
- **Avoid circular reactivity**: Don't filter data based on a reactive value that the filter updates
- **CRITICAL Reactive Loop Pattern**:
  - If `renderDT` depends on a reactive that uses `values$selection`
  - AND selecting rows updates `values$selection`
  - → Creates infinite loop: selection → reactive update → table re-render → selection reset
  - **Solution**: Make `renderDT` build data independently without reactive dependencies on selection state
- Use `isolate()` to break reactive dependencies when needed
- DT tables with `selection = "multiple"` support Ctrl+Click and Shift+Click

### dplyr/tidyverse Best Practices
- In `summarise()`, use explicit `{}` braces for multi-line if statements
- Always set `.groups = 'drop'` to avoid grouping warnings
- Can't use bracket subsetting `[condition]` directly on column names inside summarise - need proper context
- Example from this project (line 70):
```r
summarise(
  Proteins = if(has_pg_q) {
    n_distinct(Protein.Group[PG.Q.Value <= 0.01])
  } else {
    n_distinct(Protein.Group)
  },
  .groups = 'drop'
)
```

### Reproducibility Logging Pattern
- Log actions with timestamps, not just final code
- Use cumulative logging (append) rather than overwrite
- Include context: which button was clicked, which parameters changed
- Add human-readable action names before code blocks
- Implementation: `add_to_log(action_name, code_lines)` helper function

## Common Issues & Solutions

### Installation Issues
1. **"Please select a CRAN mirror" popup not showing (VS Code)**:
   - Root cause: VS Code's R terminal doesn't display interactive popups
   - Fixed in current version: CRAN mirror set automatically (line 7)
   - Manual fix: Add `options(repos = c(CRAN = "https://cloud.r-project.org"))` at top of script

1a. **"there is no package called 'readr'" (or other package) after installation**:
   - Root cause: Package not included in auto-installation list (required_pkgs)
   - Fixed in current version: All library() packages now auto-installed (lines 92-95)
   - Manual fix: Add missing package to required_pkgs vector

2. **"limpa package not found" or "R 4.5 required"**:
   - Root cause: limpa requires R 4.5+ (in Bioconductor 3.22+)
   - Solution: Upgrade R from https://cloud.r-project.org/
   - Then run: `BiocManager::install('limpa')`
   - The app will auto-detect version mismatch and show instructions

3. **"Package ggplot2 cannot be unloaded" during installation**:
   - Root cause: Installation happening after packages already loaded
   - Fixed in current version: packages checked/installed BEFORE library() calls (lines 16-81)
   - Workaround: Restart R session completely before running app

4. **"limpa is not available for Bioconductor version X.XX"**:
   - limpa requires Bioc 3.22+ (which requires R 4.5+)
   - Check your version: `BiocManager::version()`
   - Upgrade R to 4.5+ first, then Bioconductor will auto-update

### Runtime Issues
5. **App doesn't start when using `source()` in VS Code**:
   - Problem: Using `source("/path/to/DE-LIMP.r")` doesn't launch the app properly
   - Solution: Use `shiny::runApp('/Users/brettphinney/Documents/claude/DE-LIMP.r', port=3838, launch.browser=TRUE)`
   - This is the correct way to launch Shiny apps in VS Code

6. **App won't start after code changes**:
   - Shiny apps don't hot-reload
   - Must stop (`Ctrl+C`) and rerun the `shiny::runApp()` command

7. **Violin plot shows "select protein first"**:
   - Make sure table/volcano selection observer is working (line ~1115)

8. **Can't select multiple proteins**:
   - Ensure table is NOT filtered by selections
   - Check for reactive loops - table should NOT call `volcano_data()` directly

9. **Selections disappear after clicking**:
   - Reactive loop! Table re-renders on selection changes
   - Solution: Make table build data independently without reactive dependencies on selection state

### Git/Version Control Issues
10. **GitHub shows old/wrong version despite local file being correct**:
   - Problem: Local file (1715 lines) correct, but GitHub shows old version (459 lines)
   - Root cause: Modified file was NOT staged/committed to git (`git status` shows "Changes not staged")
   - Symptoms: `git diff HEAD origin/main` shows no difference, but GitHub is wrong
   - How to check: `git show HEAD:DE-LIMP.R | wc -l` to see what's actually in the git commit
   - Solution: `git add DE-LIMP.R && git commit -m "message" && git push origin main`
   - **Key lesson**: Always check `git status` and commit working files before assuming they're in the repo
   - Fixed in commit b2fbf11 (2026-02-10)

## AI Chat Feature
- Uses Google Gemini API
- Uploads top 800 proteins via File API
- Bi-directional: User selects proteins → AI analyzes, AI finds proteins → highlights in plots
- Selection format: `[[SELECT: P12345; P67890]]`

## Next Steps / TODO
- [x] Add download button for reproducibility log ✅ (2026-02-10)
- [x] Add option to save/load analysis sessions (RDS) ✅ (2026-02-11)
- [x] Documentation: Re-write README.md and USER_GUIDE.md to incorporate new features and workflows ✅ (2026-02-11, v2.0 release)
- [ ] Consider adding plot theme customization
- [ ] Volcano plot: Add annotation/legend indicating FDR threshold used for coloring
- [ ] GSEA: Add KEGG and Reactome enrichment to the pathway analysis functionality
- [ ] GSEA: Clarify which DE results (contrast) the GSEA analysis is being performed for
- [ ] Consider auto-pushing to both remotes with git hooks or GitHub Actions
