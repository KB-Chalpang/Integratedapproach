################################################################################

library(pROC)
library(ggplot2)
library(dplyr)
library(patchwork)

# Shared parameters

N               <- 2000
tol.prev        <- 0.0001
tol.cstat       <- 0.0001
max.iter.inner  <- 1000
max.iter.outer  <- 1000
n.rep.cstat     <- 5

target.cstat  <- .999999
target.prevs  <- c(0.10, 0.50, 0.8, 0.90)

# Fixed covariates

set.seed(4555632)
x1 <- rnorm(N, 0, 1)
x2 <- rnorm(N, 0, 1)
x3 <- rnorm(N, 0, 1)
X  <- cbind(1, x1, x2, x3)
B.fixed <- c(log(2), log(1.2), log(2.8))

# Helpers

det.prev  <- function(beta, X) mean(plogis(X %*% beta))
get.cstat <- function(Y, p) as.numeric(pROC::auc(Y, p, quiet = TRUE))

find.intercept <- function(sigma, X, B.fixed, target.prev, tol, max.iter) {
  B.scaled <- sigma * B.fixed
  lo <- -20; hi <- 20; mid <- 0
  for (i in seq_len(max.iter)) {
    mid  <- (lo + hi) / 2
    prev <- det.prev(c(mid, B.scaled), X)
    if (prev < target.prev) lo <- mid else hi <- mid
    if (abs(prev - target.prev) <= tol) break
  }
  list(intercept = mid, B.scaled = B.scaled,
       prev = det.prev(c(mid, B.scaled), X))
}


# Algorithm 3 — AUC = 1 special case: push sigma very large

algo3 <- function(target.prev, target.cstat, seed = 42) {
  set.seed(seed)
  
  if (target.cstat == 1.0) {
    inner <- find.intercept(1000, X, B.fixed, target.prev, tol.prev, max.iter.inner)
    beta  <- c(inner$intercept, inner$B.scaled)
    return(list(p = as.numeric(plogis(X %*% beta)), beta = beta))
  }
  
  lo <- 0; hi <- 20
  int.mid <- 0; B.scaled <- B.fixed
  
  for (outer in seq_len(max.iter.outer)) {
    sigma.mid <- (lo + hi) / 2
    inner     <- find.intercept(sigma.mid, X, B.fixed,
                                target.prev, tol.prev, max.iter.inner)
    int.mid  <- inner$intercept
    B.scaled <- inner$B.scaled
    p        <- plogis(X %*% c(int.mid, B.scaled))
    
    cstats <- numeric(n.rep.cstat)
    for (r in seq_len(n.rep.cstat)) {
      set.seed(outer * 100 + r + seed)
      Y         <- rbinom(N, 1, p)
      cstats[r] <- get.cstat(Y, p)
    }
    cstat.emp <- mean(cstats)
    
    if (cstat.emp < target.cstat) lo <- sigma.mid else hi <- sigma.mid
    if (abs(cstat.emp  - target.cstat) <= tol.cstat &&
        abs(inner$prev - target.prev)  <= tol.prev) break
  }
  
  beta <- c(int.mid, B.scaled)
  list(p = as.numeric(plogis(X %*% beta)), beta = beta)
}

# Run scenarios

cat(sprintf("Running scenarios: AUC = %.2f, prevalence = %s\n",
            target.cstat, paste(target.prevs, collapse = ", ")))

plot.df <- bind_rows(lapply(target.prevs, function(tp) {
  cat(sprintf("  prev = %.2f\n", tp))
  out <- algo3(tp, target.cstat, seed = 42)
  data.frame(p = out$p, target.prev = tp)
}))

# Shared theme

LINE_COL <- "#0077BB"

base_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = "#FAFAF8", colour = NA),
    panel.background = element_rect(fill = "#FAFAF8", colour = NA),
    panel.border     = element_rect(colour = "#CCCCCC", fill = NA, linewidth = 0.6),
    panel.grid.major = element_line(colour = "#E8E8E4", linewidth = 0.35),
    panel.grid.minor = element_line(colour = "#F0F0EC", linewidth = 0.2),
    axis.title       = element_text(face = "bold", colour = "#333333", size = 9.5),
    axis.text        = element_text(colour = "#555555", size = 8.5),
    plot.title       = element_text(face = "bold", colour = "#1A1A1A", size = 11,
                                    hjust = 0.5, margin = margin(b = 4)),
    legend.position  = "none",
    plot.margin      = margin(6, 10, 6, 8)
  )


# One panel per prevalence level

make.panel <- function(tp) {
  d    <- plot.df %>% filter(near(target.prev, tp))
  dens <- density(d$p, from = 0, to = 1, n = 512, adjust = 1.5)
  dens.df <- data.frame(x = dens$x, y = dens$y)
  
  ggplot(dens.df, aes(x = x, y = y)) +
    geom_area(fill = "transparent", colour = LINE_COL, linewidth = 0.8) +
    scale_x_continuous(name   = "Predicted probability",
                       limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.25)) +
    scale_y_continuous(name   = "Density",
                       expand = expansion(mult = c(0, 0.05))) +
    labs(title = bquote("Prevalence" == .(tp))) +
    base_theme
}

panels <- lapply(target.prevs, make.panel)

# Combine — 2 rows x 2 columns

combined <- (panels[[1]] | panels[[2]]) /
  (panels[[3]] | panels[[4]]) +
  plot_annotation(
    title    = "Algorithm 3 — Predicted probabilities under perfect seperation",
    theme    = theme(
      plot.title      = element_text(face = "bold", size = 13, hjust = 0.5,
                                     colour = "#1A1A1A", margin = margin(b = 3)),
      plot.subtitle   = element_text(size = 9, hjust = 0.5,
                                     colour = "#666666", margin = margin(b = 8)),
      plot.background = element_rect(fill = "#FAFAF8", colour = NA)
    )
  )

print(combined)
