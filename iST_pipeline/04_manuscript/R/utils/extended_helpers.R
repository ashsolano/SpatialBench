# Purpose:  Shared helpers for 04_manuscript/R/extended/ scripts: animal order,
#           VisiumHD animal map, missing-animal check, theme_ext, facet_animals(),
#           layout measuring and fixed-height saving.
#           Source after theme.R (uses dims and patchwork).
# Inputs:   none (functions only)
# Outputs:  none (defines objects in the calling environment)

# ---------------------------------------------------------------------------
# Animal IDs and column order
# ---------------------------------------------------------------------------
# Animal ID = sample name without its batch suffix (e.g. ko167_batch10 -> ko167)
sample_to_animal <- function(sample) sub("_batch[0-9]+$", "", sample)

# Order by condition (reference column order), then numeric ID
order_animals <- function(animals, condition_order = c("ko", "wt", "ctrl")) {
  condition <- sub("[0-9]+$", "", animals)
  number    <- as.integer(sub("^[a-z]+", "", animals))
  animals[order(match(condition, condition_order), number)]
}

# Ordered animal levels from the config sample lists; stops if the two
# platforms do not cover the same animals
animal_levels_from_config <- function(cfg) {
  config_animals <- list(
    MERSCOPE = sample_to_animal(cfg$spatial_analysis$merscope_samples),
    Xenium   = sample_to_animal(cfg$spatial_analysis$xenium_default_samples)
  )
  if (!setequal(config_animals$MERSCOPE, config_animals$Xenium)) {
    stop("Animal IDs differ between platforms:\n  MERSCOPE: ",
         paste(sort(config_animals$MERSCOPE), collapse = ", "), "\n  Xenium:   ",
         paste(sort(config_animals$Xenium), collapse = ", "))
  }
  order_animals(unique(config_animals$Xenium))
}

# VisiumHD sample names carry no condition (e.g. batch33_709), so map them to
# animal IDs via config gene_comparison.animals (WT709 -> visiumhd: batch33_709).
# Stops if a config VisiumHD sample has no animal.
visiumhd_animal_map <- function(cfg) {
  animals <- cfg$gene_comparison$animals
  map     <- setNames(tolower(names(animals)), purrr::map_chr(animals, "visiumhd"))
  unmapped <- setdiff(names(cfg$visiumhd$samples), names(map))
  if (length(unmapped) > 0) {
    stop("VisiumHD samples missing from config gene_comparison.animals: ",
         paste(unmapped, collapse = ", "))
  }
  map
}

# Samples deliberately excluded upstream (no background signal), per platform
excluded_samples_from_config <- function(cfg) {
  list(
    MERSCOPE = unlist(cfg$qc_backgrounds$exclude_samples$merscope),
    Xenium   = unlist(cfg$qc_backgrounds$exclude_samples$xenium)
  )
}

# Every animal absent from a platform must be an upstream exclusion
# (config qc_backgrounds.exclude_samples); anything else is a bug
check_missing_animals <- function(df, sample_col, what, animal_levels,
                                  excluded_samples, platform_levels) {
  purrr::walk(platform_levels, function(plat) {
    present          <- unique(sample_to_animal(df[[sample_col]][df$platform == plat]))
    missing          <- setdiff(animal_levels, present)
    expected_missing <- sample_to_animal(excluded_samples[[plat]])
    if (!setequal(missing, expected_missing)) {
      stop(what, " / ", plat, ": animals missing from data (", paste(missing, collapse = ", "),
           ") do not match config exclude_samples (", paste(expected_missing, collapse = ", "), ")")
    }
    if (length(missing) > 0) {
      message(what, " / ", plat, ": not shown (config exclude_samples): ",
              paste(missing, collapse = ", "))
    }
  })
}

# ---------------------------------------------------------------------------
# Theme and faceting
# ---------------------------------------------------------------------------
# Extended figures only (main figures keep axis lines): boxed panels and an
# italic row title. Borders are unclipped (strip.clip here, clip = "off" in each
# coord) so the shared strip/panel edge draws as one line; linewidth is halved
# so unclipped edges match the clipped default thickness.
theme_ext <- theme(
  panel.border     = element_rect(colour = "black", fill = NA, linewidth = rel(0.5)),
  strip.background = element_rect(colour = "black", fill = "white", linewidth = rel(0.5)),
  strip.clip       = "off",
  axis.line        = element_blank(),
  plot.title       = element_text(face = "italic", size = 7, hjust = 0,
                                  margin = margin(0, 0, 1, 0, "mm"))
)

# One facet per animal present on this platform, packed across the row; with
# fixed scales y tick labels appear on the leftmost panel only
facet_animals <- function(...) {
  facet_wrap(~ animal, nrow = 1, labeller = as_labeller(toupper), ...)
}

# ---------------------------------------------------------------------------
# Saving
# ---------------------------------------------------------------------------
# Page height (mm) of a patchwork with absolute panel heights. Measured on a
# cairo device opened first: grobs measure text on the current device, and the
# default pdf device lacks Arial
layout_height_mm <- function(p_fixed, width_mm = dims$full_w) {
  cairo_pdf(tempfile(fileext = ".pdf"), width = width_mm / 25.4, height = 10)
  on.exit(dev.off())
  gt <- patchwork::patchworkGrob(p_fixed)
  grid::convertHeight(sum(gt$heights), "mm", valueOnly = TRUE)
}

# Page width (mm) of a patchwork with absolute panel widths (as layout_height_mm)
layout_width_mm <- function(p_fixed, width_mm = dims$full_w) {
  cairo_pdf(tempfile(fileext = ".pdf"), width = width_mm / 25.4, height = 10)
  on.exit(dev.off())
  gt <- patchwork::patchworkGrob(p_fixed)
  grid::convertWidth(sum(gt$widths), "mm", valueOnly = TRUE)
}

# Panel height (mm) that makes a stacked patchwork fill page_h_mm: non-panel
# height is measured once, the rest split equally across rows
panel_h_for_page <- function(p, page_h_mm, n_rows, width_mm = dims$full_w) {
  probe_h  <- 10
  overhead <- layout_height_mm(
    p + plot_layout(heights = unit(rep(probe_h, n_rows), "mm")), width_mm
  ) - n_rows * probe_h
  panel_h <- (page_h_mm - overhead) / n_rows
  if (panel_h <= 0) stop("page_h_mm = ", page_h_mm, " leaves no room for panels (overhead ",
                         round(overhead, 1), " mm)")
  panel_h
}

# Save a stacked patchwork with panels panel_h_mm tall; page height follows
# from the layout
save_fixed_panels <- function(p, file, panel_h_mm, n_rows, width_mm = dims$full_w) {
  p_fixed <- p + plot_layout(heights = unit(rep(panel_h_mm, n_rows), "mm"))
  page_h  <- layout_height_mm(p_fixed, width_mm)
  ggsave(
    file, p_fixed,
    width  = width_mm,
    height = page_h,
    units  = "mm",
    device = cairo_pdf,
    bg     = "white"
  )
  message("Saved: ", basename(file), sprintf(" (%.1f x %.1f mm, panels %.1f mm tall)",
                                             width_mm, page_h, panel_h_mm))
}
