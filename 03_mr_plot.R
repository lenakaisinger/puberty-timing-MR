# =============================================================================
# Title: Likely Causal Effects of Circulating Hormones on Puberty Timing:
#        a Mendelian Randomisation Study
# Description: Code to create Figure 2 of the manuscript. 
# Author: Lena Kaisinger
# =============================================================================

# Dependencies
library(ggplot2)
library(dplyr)
library(readxl)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
# Update this to your input file path
# Expected column names are: Puberty Timing, beta, se, pvalue, name, Group

df1 <- read_xlsx(
  "/your/file/path",
  col_names = TRUE, 
  sheet = "Sheet1"
)

# ---------------------------------------------------------------------------
# Data prep
# ---------------------------------------------------------------------------

# Rename spaced column
df1 <- df1 %>% rename(Puberty_Timing = `Puberty Timing`)

# Fix order of "Group" variable and preserve original spreadsheet row order
df1 <- df1 %>%
  mutate(
    Group = factor(Group, levels = unique(Group)),
    orig_order = row_number()
  ) %>%
  group_by(Group, name) %>%
  mutate(
    row_within = if_else(Puberty_Timing == "AAM", 1L, 2L),
    name_order = min(orig_order)
  ) %>%
  ungroup()

# Compute confidence intervals and specify significance threshold
df1 <- df1 %>%
  mutate(lci = beta - 1.96 * se,
         uci = beta + 1.96 * se,
         is_header = FALSE,
         significant = pvalue < (0.05 / 10),  # Bonferroni threshold
         shape_key = paste0(Puberty_Timing, "_", ifelse(significant, "sig", "ns"))
         )

# Create an empty "dummy" row for each group 
header_rows <- df1 %>%
  distinct(Group) %>%
  mutate(
    name           = NA_character_,
    Puberty_Timing = NA_character_,
    beta           = NA_real_,
    se             = NA_real_,
    pvalue         = NA_real_,
    lci            = NA_real_,
    uci            = NA_real_,
    row_within     = 0L,       
    name_order     = 0L,       
    orig_order     = 0L,
    is_header      = TRUE,
    significant    = NA, 
    shape_key      = NA
  )

# Combine the headers with the data
# Maintaining the order of: Group → name_order → row_within
df_all <- bind_rows(df1, header_rows) %>%
  arrange(Group, name_order, row_within) %>%
  mutate(y_pos = rev(seq_len(n()))) # reverse so first group plots at the top

# Separate back out
df_data   <- df_all %>% filter(!is_header)
df_header <- df_all %>% filter(is_header)

# Midpoint of each exposure pair for name labels
name_labels <- df_data %>%
  group_by(Group, name) %>%
  summarise(y_mid = mean(y_pos), .groups = "drop")

# ---------------------------------------------------------------------------
# Calculate axis limits
# ---------------------------------------------------------------------------

# y limits
group_ranges <- df_all %>%
  group_by(Group) %>%
  summarise(
    y_min = min(y_pos) - 0.5,
    y_max = max(y_pos) + 0.5,
    .groups = "drop"
  )

y_total_min <- min(group_ranges$y_min)
y_total_max <- max(group_ranges$y_max)

# x limits
x_data_min  <- min(df_data$lci, na.rm = TRUE)
x_data_max  <- max(df_data$uci, na.rm = TRUE)
x_range     <- x_data_max - x_data_min
x_axis_min  <- -0.4
x_axis_max  <- max(x_data_max + x_range * 0.05, 0.1)
label_space <- x_range * 0.60
x_plot_min  <- x_axis_min - label_space
label_x     <- x_axis_min - x_range * 0.04

# ---------------------------------------------------------------------------
# Specify the color, shapes and theme used for the plot
# ---------------------------------------------------------------------------

# Colourblind-friendly colours 
puberty_colours <- c("AAM" = "#E07B39", "AVB" = "#0072B2")

# Shapes:
# significant — filled circle (AAM=16) / filled triangle (AVB=17)
# non-significant — open circle  (AAM=1)  / open triangle   (AVB=2)
shape_values <- c(
  "AAM_sig" = 16, "AAM_ns" = 1,
  "AVB_sig" = 17, "AVB_ns" = 2
)

# Theme adhering to JCEM's art guidelines
journal_theme <- theme_classic(base_family = "Arial", base_size = 8) +
  theme(
    axis.line        = element_blank(),
    axis.ticks.x     = element_line(linewidth = 0.5, colour = "black"),
    axis.ticks.y     = element_blank(),
    axis.text.x      = element_text(size = 8, family = "Arial", colour = "black"),
    axis.text.y      = element_blank(),
    axis.title.x     = element_text(size = 8, family = "Arial", colour = "black"),
    axis.title.y     = element_blank(),
    legend.text      = element_text(size = 8, family = "Arial"),
    legend.title     = element_blank(),
    legend.key.size  = unit(3, "mm"),
    legend.background = element_rect(fill = "white", linewidth = 0.5,
                                     linetype = "solid", colour = "black"),
    legend.position  = "right",
    panel.border     = element_blank(),
    plot.margin      = margin(3, 3, 3, 3, "mm"),
    plot.title       = element_blank()
  )

# ---------------------------------------------------------------------------
# Build plot
# ---------------------------------------------------------------------------

p1 <- ggplot(df_data, aes(x = beta, y = y_pos,
                          colour = Puberty_Timing,
                          shape  = shape_key)) +
  
  # Alternating group shading 
  geom_rect(
    data        = group_ranges %>% filter(row_number() %% 2 == 0),
    aes(ymin = y_min, ymax = y_max, xmin = x_axis_min, xmax = x_axis_max),
    inherit.aes = FALSE,
    fill = "#F5F5F5", alpha = 0.7
  ) +
  
  # Manual panel border
  annotate("rect",
           xmin = x_axis_min, xmax = x_axis_max,
           ymin = y_total_min, ymax = y_total_max,
           fill = NA, colour = "black", linewidth = 0.5) +
  
  # Zero reference line
  annotate("segment", x = 0, xend = 0,
           y = y_total_min, yend = y_total_max,
           colour = "black", linetype = "dashed",
           linewidth = 0.5, alpha = 0.7) +
  
  # CI bars
  geom_linerange(aes(xmin = lci, xmax = uci),
                 linewidth = 0.5, orientation = "y") +
  
  # Point estimates
  geom_point(size = 1.8, stroke = 0.4) +
  
  # Separator between groups
  geom_segment(
    data = group_ranges %>% filter(row_number() < n()),
    aes(x = x_axis_min, xend = x_axis_max,
        y = y_min, yend = y_min),
    inherit.aes = FALSE,
    colour = "black", linewidth = 0.4
  ) +
  
  # Exposure labels
  geom_text(
    data = name_labels,
    aes(x = label_x, y = y_mid, label = name),
    inherit.aes = FALSE,
    hjust  = 1, 
    size = 8 / .pt
  ) +
  
  # Group header labels
  geom_text(
    data = df_header,
    aes(x = label_x, y = y_pos, label = as.character(Group)),
    inherit.aes = FALSE,
    hjust    = 1, vjust = 0.5,
    size     = 8 / .pt,
    fontface = "bold"
  ) +
  
  scale_x_continuous(
    limits = c(x_plot_min, x_axis_max),
    breaks = scales::pretty_breaks()(c(x_axis_min, x_axis_max)),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(y_total_min, y_total_max),
    expand = expansion(add = c(0, 0))
  ) +
  
  # Give colour and shape identical name/breaks/labels so ggplot merges into one legend box
  scale_colour_manual(
    name   = "Puberty timing",
    breaks = c("AAM", "AVB"),
    values = puberty_colours,
    labels = c("AAM \u2014 age at menarche", "AVB \u2014 age at voice breaking")
  ) +
  scale_shape_manual(
    name   = "Puberty timing",
    breaks = c("AAM_sig", "AVB_sig"),
    values = shape_values,
    labels = c("AAM \u2014 age at menarche", "AVB \u2014 age at voice breaking")
  ) +
  
  labs(x = "Puberty timing (years)") +
  journal_theme

# ---------------------------------------------------------------------------
# Save the plot
# ---------------------------------------------------------------------------

# Specify the height of the plot
n_rows <- nrow(df_all)
fig_h  <- min(9.25, n_rows * 0.15 + 0.5)

# Save the plot
ggsave(
  filename = "/your/path/plot.pdf",
  plot     = p1,
  width    = 5.8,
  height   = fig_h,
  units    = "in",
  device   = cairo_pdf
)
