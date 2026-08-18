make_fig4 <- function() {
  make_temporal_effect_figure(
    data_name = "quantum_yield",
    resp_var = "qy",
    panel_suffix = "Fv/Fm",
    output_file = "fig4_quantum_yield_effects.pdf",
    height_mm = 125
  )
}
