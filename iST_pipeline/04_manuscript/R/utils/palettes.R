# Purpose:  Colour palettes for manuscript figures beyond the platform-level
#           palette in theme.R. Includes cell-type, condition, and GC zone
#           colour mappings shared across all figure scripts.
# Inputs:   none
# Outputs:  pal_cell_type, pal_condition, pal_gc_zone, pal_segmentation (in calling environment)

# --- Cell type palette ---
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

# --- GC zone palette ---
pal_gc_zone <- c(
  "Dark zone"  = "#FF8C00",
  "Light zone" = "#68228B"
)

# --- Condition palette ---
pal_condition <- c(
  "wt"   = "#647dc9",
  "ko"   = "#5C5C5C",
  "ctrl" = "#b1b1b1"
)

# --- Segmentation palette ---
pal_segmentation <- c(
  "Default"   = "#666666",
  "Cellpose"  = "#0072B2",
  "Proseg"    = "#009E73"
)

