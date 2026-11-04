################################################################################
library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)


# Shared parameters

N               <- 2000
tol.prev        <- 0.0001
tol.cstat       <- 0.0001
max.iter.inner  <- 1000
max.iter.outer  <- 1000
n.rep.cstat     <- 5

target.prevs  <- c(0.10, 0.50, 0.80)
target.cstats <- c(0.50, 0.70, 0.85)

# Fixed covariates

set.seed(4555632)
x1 <- rnorm(N, 0, 1)
x2 <- rnorm(N, 0, 1)
x3 <- rnorm(N, 0, 1)
X  <- cbind(1, x1, x2, x3)
B.fixed <- c(log(2), log(1.2), log(2.8))


det.prev <- function(beta, X) mean(plogis(X %*% beta))
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

# Algorithm 1 — bisect intercept only (sigma fixed at 1)

algo1 <- function(target.prev, target.cstat, seed = 42) {
  set.seed(seed)
  sigma    <- 1
  B.scaled <- sigma * B.fixed
  lo <- -20; hi <- 20; mid <- 0
  
  for (i in seq_len(max.iter.inner)) {
    mid  <- (lo + hi) / 2
    beta <- c(mid, B.scaled)
    prev <- det.prev(beta, X)
    if (prev < target.prev) lo <- mid else hi <- mid
    if (abs(prev - target.prev) <= tol.prev) break
  }
  
  beta <- c(mid, B.scaled)
  list(p = as.numeric(plogis(X %*% beta)), beta = beta)
}

# Algorithm 2 — bisect scalar only (intercept fixed at 0)

algo2 <- function(target.prev, target.cstat, seed = 42) {
  set.seed(seed)
  int.fixed <- 0
  lo <- 0; hi <- 20; sigma <- 1
  
  if (target.cstat == 0.5) {
    B.scaled <- 0 * B.fixed
    beta     <- c(int.fixed, B.scaled)
    return(list(p = as.numeric(plogis(X %*% beta)), beta = beta))
  }
  
  for (outer in seq_len(max.iter.outer)) {
    sigma    <- (lo + hi) / 2
    B.scaled <- sigma * B.fixed
    p        <- plogis(X %*% c(int.fixed, B.scaled))
    
    cstats <- numeric(n.rep.cstat)
    for (r in seq_len(n.rep.cstat)) {
      set.seed(outer * 100 + r + seed)
      Y         <- rbinom(N, 1, p)
      cstats[r] <- get.cstat(Y, p)
    }
    cstat.emp <- mean(cstats)
    
    if (cstat.emp < target.cstat) lo <- sigma else hi <- sigma
    if (abs(cstat.emp - target.cstat) <= tol.cstat) break
  }
  
  beta <- c(int.fixed, sigma * B.fixed)
  list(p = as.numeric(plogis(X %*% beta)), beta = beta)
}

# Algorithm 3 — nested bisection

algo3 <- function(target.prev, target.cstat, seed = 42) {
  set.seed(seed)
  
  if (target.cstat == 0.5) {
    inner <- find.intercept(0, X, B.fixed, target.prev, tol.prev, max.iter.inner)
    beta  <- c(inner$intercept, rep(0, length(B.fixed)))
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

# Run all scenarios for all algorithms

run.scenarios <- function(algo.fn, algo.label) {
  cat(sprintf("\nRunning %s...\n", algo.label))
  bind_rows(lapply(target.prevs, function(tp) {
    bind_rows(lapply(target.cstats, function(tc) {
      cat(sprintf("  prev=%.2f, cstat=%.2f\n", tp, tc))
      out <- algo.fn(tp, tc, seed = 42)
      data.frame(p = out$p, target.prev = tp, target.cstat = tc,
                 algo = algo.label)
    }))
  }))
}

df1 <- run.scenarios(algo1, "Algorithm 1")
df2 <- run.scenarios(algo2, "Algorithm 2")
df3 <- run.scenarios(algo3, "Algorithm 3")

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


# Panel builder

make.panel <- function(df, tp, tc) {
  d    <- df %>% filter(target.prev == tp, target.cstat == tc)
  dens <- density(d$p, from = 0, to = 1, n = 512, adjust = 1.5)
  dens.df <- data.frame(x = dens$x, y = dens$y)
  
  ggplot(dens.df, aes(x = x, y = y)) +
    geom_area(fill = "transparent", colour = LINE_COL, linewidth = 0.8) +
    scale_x_continuous(name   = "Predicted probability",
                       limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.25)) +
    scale_y_continuous(name   = "Density",
                       expand = expansion(mult = c(0, 0.05))) +
    labs(title = sprintf("Prev = %.2f  |  AUC = %.2f", tp, tc)) +
    base_theme
}


# Build 3x3 grids per algorithm and combine into one figure

make.grid <- function(df, title) {
  panels <- lapply(target.prevs, function(tp) {
    lapply(target.cstats, function(tc) make.panel(df, tp, tc))
  })
  
  (panels[[1]][[1]] | panels[[1]][[2]] | panels[[1]][[3]]) /
    (panels[[2]][[1]] | panels[[2]][[2]] | panels[[2]][[3]]) /
    (panels[[3]][[1]] | panels[[3]][[2]] | panels[[3]][[3]]) +
    plot_annotation(
      title = title,
      theme = theme(
        plot.title      = element_text(face = "bold", size = 13, hjust = 0.5,
                                       colour = "#1A1A1A", margin = margin(b = 6)),
        plot.background = element_rect(fill = "#FAFAF8", colour = NA)
      )
    )
}

grid1 <- make.grid(df1, "Algorithm 1 — Distributions of predicted probabilities")
grid2 <- make.grid(df2, "Algorithm 2 — Distributions of predicted probabilities")
grid3 <- make.grid(df3, "Algorithm 3 — Distributions of predicted probabilities")

# Print each separately so they appear as three distinct figures
print(grid1)
print(grid2)
print(grid3)
