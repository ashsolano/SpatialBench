# Purpose:  Shared plotting helper functions used across manuscript figure
#           scripts. Includes wrappers for common panel types (spatial maps,
#           violin plots, dot plots) so figure scripts stay concise.
# Inputs:   none (functions only — source into figure scripts)
# Outputs:  none (defines functions in the calling environment)

# Planned helpers (to be added):
#
#   plot_spatial_cells()    — ggplot wrapper for spatial scatter of cell coords,
#                             coloured by a metadata variable
#
#   plot_umap_panel()       — standard UMAP panel with theme_sb applied
#
#   plot_violin_panel()     — violin + jitter for per-cell metrics by method
#
#   plot_dot_panel()        — dot plot for cell-type marker expression
#
#   save_figure()           — thin wrapper around ggsave using dims constants
#                             from theme.R

message("plot_helpers.R — in development")
