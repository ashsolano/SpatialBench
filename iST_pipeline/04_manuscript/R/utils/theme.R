# Purpose:  Shared ggplot2 theme, colour palettes, and figure dimension
#           constants for all SpatialBench manuscript figures.
# Inputs:   none
# Outputs:  pal_muted, pal_fill, scale_*_platform(), theme_sb(),
#           theme_bg_panels, dims (in calling environment)

library(ggplot2)
library(patchwork)

# Platform palette
pal_muted <- c(
  VisiumHD = "#6DAA2C",   # olive-green
  Xenium   = "#4477AA",   # indigo
  MERSCOPE = "#AA3377",   # plum
  snRNAseq = "#D4A017"
)

# Extend with grey for any "Background" fill
pal_fill <- c(pal_muted, Background = "grey85")

scale_colour_platform <- function(...) scale_colour_manual(values = pal_muted, ...)
scale_fill_platform   <- function(...) scale_fill_manual(values = pal_fill, ...)

# Core publication theme (base_size = 7 pt matches Nature/Cell figure guidelines)
theme_sb <- function(base_family = "Arial", base_size = 7) {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid   = element_blank(),
      axis.title   = element_text(size   = base_size,
                                  colour = "black",
                                  family = base_family),
      axis.text    = element_text(size   = base_size - 1,
                                  colour = "black",
                                  family = base_family),
      legend.title = element_text(size   = base_size - 1,
                                  colour = "black",
                                  family = base_family),
      legend.text  = element_text(size   = base_size - 1,
                                  colour = "black",
                                  family = base_family),
      plot.margin  = margin(1, 1, 1, 1, "mm")
    )
}

# Additions for the background / Moran's I panels (add after theme_sb());
# shared by fig2_background.R and extended/fig2ext_background_persample.R
theme_bg_panels <- theme(
  axis.text.x      = element_text(angle = 45, hjust = 1),
  legend.position  = "none",
  axis.line        = element_line(colour = "black"),
  strip.background = element_rect(colour = "black", fill = "white"),
  strip.text       = element_text(face = "italic"),
  panel.background = element_blank(),
  panel.border     = element_blank(),
  plot.background  = element_blank(),
  panel.spacing    = unit(2, "mm")
)

# A4 paper size and margins (millimetres)
paper_mm   <- list(width = 210, height = 297)
margins_mm <- list(horiz = 20, vert = 20)

# Figure dimensions (millimetres)
dims <- list(
  full_w  = paper_mm$width  - 2 * margins_mm$horiz,        # ~170 mm — double-column
  half_w  = (paper_mm$width - 2 * margins_mm$horiz) / 2,   # ~85 mm  — single-column
  height  = paper_mm$height - 2 * margins_mm$vert,          # ~257 mm — full page
  line    = 0.5,
  point   = 1,
  tagface = "bold",
  tagsize = 7
)

# Saving convention:
#   ggsave(
#     filename = "figures/figX.pdf",
#     plot     = p,
#     width    = dims$full_w,
#     height   = dims$height / 2,
#     units    = "mm",
#     device   = cairo_pdf
#   )
