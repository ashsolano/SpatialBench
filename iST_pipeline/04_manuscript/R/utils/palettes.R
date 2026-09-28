# Purpose:  Colour palettes for manuscript figures beyond the platform-level
#           palette in theme.R. Includes cell-type, condition, and GC zone
#           colour mappings shared across all figure scripts, and light
#           platform variants for density-gradient low ends.
# Inputs:   none
# Outputs:  pal_cell_type, pal_condition, pal_gc_zone, pal_segmentation,
#           pal_muted_light, pal_density (in calling environment)

# --- Density contour palette ---
# Gradient endpoints from the reference SFig1 (ColorBrewer Greens / RdPu /
# Blues steps 2 and 9); interpolate in Lab (scales::seq_gradient_pal())
pal_density <- list(
  VisiumHD = c(low = "#e5f5e0", high = "#00441b"),   # near-white -> deep green
  MERSCOPE = c(low = "#fde0dd", high = "#7a0177"),   # near-white -> plum
  Xenium   = c(low = "#deebf7", high = "#08306b")    # near-white -> navy
)

# --- Light platform palette ---
# Light variants of pal_muted (theme.R) for the Fig 1c density gradient
pal_muted_light <- c(
  VisiumHD = "#e5f5d6",   # near-white olive-green
  MERSCOPE = "#f5d9ec",   # near-white plum
  Xenium   = "#d9e8f5"    # near-white indigo
)

pal_cell_type <- c(
  "Erythrocytes"   = "#00A087FF",
  "Naive B cells"  = "#E9967A",
  "GC B cells"     = "#525ecc99",
  "Plasma B cells" = "#480607",
  "T cells"        = "#4DBBD5FF",
  "NK cells"       = "#F0E685FF",
  "ILC"            = "#802268FF",
  "Monocytes"      = "#7E6148FF",
  "Macrophages"    = "#CCEBC5",
  "DC"             = "#386CB0",
  "Granulocytes"   = "#cc52c099",
  "Stem cells"     = "#F0027F"
)

pal_gc_zone <- c(
  "Dark zone"  = "#FF8C00",
  "Light zone" = "#68228B"
)

pal_condition <- c(
  "wt"   = "#647dc9",
  "ko"   = "#5C5C5C",
  "ctrl" = "#b1b1b1"
)

pal_segmentation <- c(
  "Default"   = "#666666",
  "Cellpose"  = "#0072B2",
  "Proseg"    = "#009E73"
)

