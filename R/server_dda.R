# ==============================================================================
#  SERVER MODULE — DDA Search (Sage + Casanovo pipeline)
#  Called from app.R as: server_dda(input, output, session, values, add_to_log)
#  Status: Stub — no pipeline logic yet
# ==============================================================================

server_dda <- function(input, output, session, values, add_to_log) {

  # --- Mode observer: sync input to reactive values ---
  observeEvent(input$acquisition_mode, {
    values$acquisition_mode <- input$acquisition_mode
  })

  # --- Contextual label for the mode switcher ---
  output$mode_context_label <- renderUI({
    mode <- input$acquisition_mode %||% "dia"
    switch(mode,
      "dia" = tags$span(
        style = "font-size: 12px; color: #0d6efd; font-style: italic;",
        icon("circle-check", style = "color: #198754;"),
        " DIA-NN + limpa pipeline"
      ),
      "dda" = tags$span(
        style = "font-size: 12px; color: #6c757d; font-style: italic;",
        icon("circle-info", style = "color: #0d6efd;"),
        " Sage + Casanovo pipeline"
      ),
      "xlms" = tags$span(
        style = "font-size: 12px; color: #6c757d; font-style: italic;",
        icon("diagram-project", style = "color: #6f42c1;"),
        " MeroX + xiSearch + network"
      )
    )
  })

}
