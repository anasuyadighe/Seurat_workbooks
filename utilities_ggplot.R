# Customization of ggplot themes and geoms. 

# Dependencies:
#   - Custom fonts:
#     - Source Sans Pro (download from https://fonts.google.com/specimen/Source+Sans+Pro)
#   - Packages (don't load here, but require them to be loaded into the global env 
#               wherever this script is being sourced):
#     - ggplot2
#     - cowplot
#     - ggrepel

# define variables ============================================================
.use_font <- "Source Sans Pro"

# Update geoms ================================================================
## text -----------------------------------------------------------------------
ggplot2::update_geom_defaults("text", list(family = .use_font))
ggplot2::update_geom_defaults(geom = "text_repel", list(family = .use_font)) # requires ggrepel

ggplot2::update_geom_defaults("label", list(family = .use_font))

message("Updated default font for geoms `text` and `text_repel` to: ", .use_font)

# colour palettes =============================================================
pal_okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

# Themes ======================================================================
## half open ------------------------------------------------------------------
### with minimal grid
theme_half_open_grid <- function(font_size = 12,
                                 line_size = 0.25,
                                 font_family = .use_font,
                                 ...) {
  cowplot::theme_half_open(
    font_family = font_family,
    font_size   = font_size,
    line_size   = line_size,
    ...
  ) %+replace%
    theme(
      strip.text.x     = element_text(hjust = 0),
      strip.background = element_blank(),
      panel.grid.major = element_line(size = 0.25, color = "grey85")
    )
}

### no grid
theme_half_open_nogrid <- function(font_size = 12,
                                   line_size = 0.25,
                                   font_family = .use_font,
                                   ...) {
  cowplot::theme_half_open(
    font_family = font_family,
    font_size   = font_size,
    line_size   = line_size,
    ...
  ) %+replace%
    theme(
      strip.text.x = element_text(hjust = 0),
      strip.background = element_blank(),
      panel.grid = element_blank()
    )
}

## bw variants ----------------------------------------------------------------
theme_bw_2 <- function(base_family = .use_font,
                       ...) {
  theme_bw(base_family = base_family, ...) %+replace%
    theme(
      strip.background = element_blank(),
      panel.grid.major = element_line(size = 0.25, 
                                      color = "grey85"),
      panel.grid.minor = element_blank()
    )
}

## for heatmaps ---------------------------------------------------------------
theme_tileplot <- function(base_family = .use_font,
                           ...) {
  cowplot::theme_minimal_grid(font_family = base_family, ...) %+replace%
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90, 
                                     hjust = 1, 
                                     vjust = 0.5),
          axis.text.y = element_text(angle = 0, 
                                     hjust = 1, 
                                     vjust = 0.5),
          axis.title = element_blank(),
          axis.ticks = element_blank())
}

# =============================================================================
# message("Loaded the following custom themes: \n\t", 
#         paste0(objects(pattern = "theme_"), 
#                collapse = "\n\t"))

# =============================================================================
