make_fig3 <- function() {
  make_temporal_effect_figure(
    data_name = "growth",
    resp_var = "diameter_inc_t0",
    panel_suffix = "diameter incr.",
    output_file = "fig3_diameter_increment_effects.pdf",
    height_mm = 125
  )
}
