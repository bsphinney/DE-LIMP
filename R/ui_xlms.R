# ==============================================================================
#  UI HELPERS — XL-MS mode stub panels
#  Called from R/ui.R for XL-MS conditionalPanel content
# ==============================================================================

# Stub: XL-MS setup panel for the Search tab sidebar
xlms_setup_ui <- function() {
  div(
    class = "alert alert-info", role = "alert",
    style = "margin-top: 10px;",
    icon("diagram-project"),
    tags$strong(" Cross-Linking Mass Spectrometry"),
    tags$br(),
    tags$span(
      style = "font-size: 13px; color: #555;",
      "MeroX + xiSearch pipeline for protein-protein interaction mapping. ",
      "Coming soon."
    )
  )
}
