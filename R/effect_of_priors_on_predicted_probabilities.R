rm(list = ls())
library(ggplot2)
library(dplyr)
library(patchwork)

# Parameters

N       <- 2000
n.draws <- 1000

set.seed(4555632)
x1 <- rnorm(N, 0, 1)
x2 <- rnorm(N, 0, 1)
x3 <- rnorm(N, 0, 1)
X  <- cbind(1, x1, x2, x3)

prior.sds  <- c(0.5, 1, 2, 10, 100, 1000)

# Compute induced predicted probabilities for each prior SD

collect <- bind_rows(lapply(prior.sds, function(sd) {
  Beta  <- matrix(rnorm(n.draws * 4, mean = 0, sd = sd),
                  nrow = n.draws, ncol = 4)
  P     <- plogis(X %*% t(Beta))   # N x n.draws
  data.frame(p = as.vector(P), sd = sd)
}))


# Shared theme

base_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "#FAFAF8", colour = NA),
    panel.background = element_rect(fill = "#FAFAF8", colour = NA),
    panel.border     = element_rect(colour = "#CCCCCC", fill = NA, linewidth = 0.6),
    panel.grid.major = element_line(colour = "#E8E8E4", linewidth = 0.35),
    panel.grid.minor = element_line(colour = "#F0F0EC", linewidth = 0.2),
    axis.title       = element_text(face = "bold", colour = "#333333", size = 9.5),
    axis.text        = element_text(colour = "#555555", size = 8.5),
    plot.title       = element_text(face = "bold", colour = "#1A1A1A", size = 10,
                                    hjust = 0.5, margin = margin(b = 4)),
    legend.position  = "none",
    plot.margin      = margin(6, 10, 6, 8)
  )

# Colour per panel
sd.cols <- c(
  "0.5"  = "#2A9D8F",
  "1"    = "#0077BB",
  "2"    = "#F4A500",
  "10"   = "#E63946",
  "100"  = "#7B2D8B",
  "1000" = "#111111"
)


# One plot per prior SD

LINE_COL <- "#0077BB"

make.panel <- function(sd) {
  d   <- collect %>% filter(.data$sd == !!sd)
  # Compute density manually so we can close the polygon to the x-axis
  dens <- density(d$p, from = 0, to = 1, n = 512, adjust = 1.2)
  dens.df <- data.frame(x = dens$x, y = dens$y)
  
  ggplot(dens.df, aes(x = x, y = y)) +
    geom_area(fill = "transparent", colour = LINE_COL, linewidth = 0.8) +
    scale_x_continuous(
      name   = "Predicted probability",
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25)
    ) +
    scale_y_continuous(name = "Density", expand = expansion(mult = c(0, 0.05))) +
    labs(title = bquote(sigma == .(sd))) +
    base_theme
}

panels <- lapply(prior.sds, make.panel)


# Combine — 2 rows of 3

combined <- (panels[[1]] | panels[[2]] | panels[[3]]) /
  (panels[[4]] | panels[[5]] | panels[[6]]) +
  plot_annotation(
    title   = "Induced prior on predicted probability",
    theme   = theme(
      plot.title      = element_text(face = "bold", size = 13, hjust = 0.5,
                                     colour = "#1A1A1A", margin = margin(b = 3)),
      plot.subtitle   = element_text(size = 8.5, hjust = 0.5,
                                     colour = "#666666", margin = margin(b = 8)),
      plot.caption    = element_text(colour = "#888888", size = 8,
                                     hjust = 0.5, margin = margin(t = 8)),
      plot.background = element_rect(fill = "#FAFAF8", colour = NA)
    )
  )

print(combined)
