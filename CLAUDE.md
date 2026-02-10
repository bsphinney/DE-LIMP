# DE-LIMP Project Context for Claude

## Project Overview
DE-LIMP is a Shiny proteomics data analysis pipeline using the LIMPA R package for differential expression analysis of DIA-NN data.

## Key Files
- **DE-LIMP.r** - Main Shiny app (948 lines)
- **Main URL**: http://localhost:3838 when running

## Recent Changes & Important Fixes

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
- **Always restart app after code changes**: `pkill -f "DE-LIMP.r" && Rscript -e "shiny::runApp('DE-LIMP.r', port=3838)" &`
- **Check status**: `lsof -i :3838`
- **View logs**: `tail -f /tmp/shiny-app.log`
- **Note**: Shiny apps don't hot-reload - must restart after every code change

## Key Patterns & Gotchas

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
1. **App won't start**: Check line 70 syntax, ensure all packages installed
2. **Violin plot shows "select protein first"**: Make sure table/volcano selection observer is working (line 944)
3. **Can't select multiple proteins**:
   - Ensure table is NOT filtered by selections (line 764)
   - Check for reactive loops - table should NOT call `volcano_data()` (line 767)
4. **Selections disappear after clicking**: Reactive loop! Table re-renders on selection changes
   - Solution: Make table build data independently without reactive dependencies on selection state

## AI Chat Feature
- Uses Google Gemini API
- Uploads top 800 proteins via File API
- Bi-directional: User selects proteins → AI analyzes, AI finds proteins → highlights in plots
- Selection format: `[[SELECT: P12345; P67890]]`

## Next Steps / TODO
- [ ] Add download button for reproducibility log
- [ ] Consider adding plot theme customization
- [ ] Add option to save/load analysis sessions (RDS)
