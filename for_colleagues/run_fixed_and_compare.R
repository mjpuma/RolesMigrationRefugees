# =============================================================================
# run_fixed_and_compare.R
#
# Step 1: re-runs the ORIGINAL role method and checks that it reproduces the
#         saved roles in legacy/results/ (adjusted Rand should be 1.0).
# Step 2: runs FIXED_role_assessment_section.R on the migrant (UNDESA) and
#         refugee (UNHCR) networks, with binary and with valued ties, and
#         compares each result with the original roles.
# Nothing in legacy/results/ is overwritten.
#
# Run from the repository root:   Rscript for_colleagues/run_fixed_and_compare.R
# Needs the raw data in data/raw/ (see README) and these packages:
#   install.packages(c("igraph", "blockmodeling", "cluster"))
# Outputs: output/fixed_section/ (role files, k diagnostics, comparison.csv)
# =============================================================================

suppressPackageStartupMessages({
  library(igraph)                                  # network objects
  library(blockmodeling)                           # REGE.for / REGE.nm.for
  library(cluster)                                 # silhouette()
})

out_dir <- "output/fixed_section"                  # all results of this script go here
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Edge lists, built exactly as in the original scripts --------------------
States <- read.csv("data_ref/country_list.csv")    # Abel & Cohen country list (name <-> iso3c)

# Migrants (FINAL_undesa_role_assessment.R, lines 49-59)
edges.mig <- read.csv("data/raw/bilat_mig.csv")[, c(1:3, 9)]   # year0, orig, dest, da_pb_closed
edges.mig$state.origin.name      <- States$name[match(edges.mig$orig, States$iso3c)]
edges.mig$state.destination.name <- States$name[match(edges.mig$dest, States$iso3c)]
edges.mig$year5 <- edges.mig$year0 + 5             # label each 5-year period by its end year
edges.mig <- data.frame(from = paste0(edges.mig$state.origin.name, ".", edges.mig$year5),
                        to   = paste0(edges.mig$state.destination.name, ".", edges.mig$year5),
                        flow = edges.mig$da_pb_closed)

# Refugees (FINAL_unhcr_role_assessment.R, lines 52-61)
edges.ref <- read.csv("data/raw/UNHCR_Flow_Data_REF.csv")[, -c(2, 4)]  # names + refflows.YYYY columns
colnames(edges.ref)[1:2] <- c("from", "to")
edges.ref.long <- reshape(edges.ref, idvar = c("from", "to"),
                          varying = colnames(edges.ref)[3:ncol(edges.ref)],
                          direction = "long")      # one row per dyad-year; column "time" = year
edges.ref.long <- data.frame(from = paste(edges.ref.long$from, edges.ref.long$time, sep = "."),
                             to   = paste(edges.ref.long$to,   edges.ref.long$time, sep = "."),
                             flow = edges.ref.long$refflows)

# ---- Adjusted Rand index (agreement between two partitions; 1 = identical) ---
ari <- function(a, b) {
  tab <- table(a, b); c2 <- function(n) n * (n - 1) / 2
  s_ij <- sum(c2(tab)); s_a <- sum(c2(rowSums(tab))); s_b <- sum(c2(colSums(tab)))
  e <- s_a * s_b / c2(sum(tab))
  (s_ij - e) / (0.5 * (s_a + s_b) - e)
}

# ---- Run the FIXED section with given settings, in a clean environment --------
fixed_src <- readLines("for_colleagues/FIXED_role_assessment_section.R")
set_line  <- function(src, name, value)            # replace one "NAME <- ..." settings line
  sub(paste0("^", name, "\\s*<-.*"), paste(name, "<-", value), src)

run_fixed <- function(dataset, edges, use_weights) {
  tag <- if (use_weights) "_FIXED_valued" else "_FIXED_binary"
  src <- fixed_src
  src <- set_line(src, "DATASET",     deparse(dataset))
  src <- set_line(src, "USE_WEIGHTS", use_weights)
  src <- set_line(src, "RESULTS_DIR", deparse(out_dir))
  src <- set_line(src, "OUT_TAG",     deparse(tag))
  env <- new.env()                                 # keeps the two runs from interfering
  env$edges.mig <- edges                           # the section reads its edge list from 'edges.mig'
  pdf(file.path(out_dir, paste0(dataset, tag, "_plots.pdf")))   # dendrogram + diagnostics plots
  on.exit(dev.off())
  eval(parse(text = src), envir = env)
  env$c                                            # role assignments: id, role
}

# ---- Step 1: reproduce the ORIGINAL method exactly ----------------------------
# Same steps as FINAL_*_role_assessment.R lines 61-134: 0/1 ties, REGE.nm.for
# with default settings, then what hcut(x = a9, method = "average") actually
# did: Euclidean distances between rows of 1 - REGE, Ward.D2 linkage, k = 3.
# If this matches the saved roles (adjusted Rand = 1), the data and setup are
# the same as the original run, and any difference in step 2 comes from the fixes.
run_original <- function(edges) {
  a  <- na.omit(edges[edges[, 3] >= 1, ])          # ties with at least 1 person (original cutoff)
  g  <- graph_from_data_frame(a)                   # directed graph of "Country.Year" nodes
  A  <- as.matrix(as_adjacency_matrix(g, sparse = FALSE))  # unweighted adjacency, as in the original
  A[A > 0] <- 1                                    # 0/1 ties
  a9 <- 1 - REGE.nm.for(A, E = 1)$E                # original REGE call, default iterations
  dimnames(a9) <- dimnames(A)
  b  <- cutree(hclust(dist(a9), method = "ward.D2"), k = 3)  # what hcut() actually computed
  data.frame(id = names(b), role = b)
}

# ---- Step 2: run the FIXED section and compare with the original roles --------
results <- list()
for (dataset in c("UNDESA", "UNHCR")) {
  edges    <- if (dataset == "UNDESA") edges.mig else edges.ref.long
  original <- read.csv(file.path("legacy/results", paste0(dataset, "_full_roles_k3.csv")))[, c("id", "role")]

  cat("\n================", dataset, "original method (reproduction check) ================\n")
  rep <- run_original(edges)
  m   <- merge(original, rep, by = "id", suffixes = c("_original", "_rerun"))
  print(table(original = m$role_original, rerun = m$role_rerun))
  results[[length(results) + 1]] <- data.frame(
    dataset = dataset, ties = "original method (reproduction)",
    n_matched = nrow(m), adjusted_rand_vs_original = round(ari(m$role_original, m$role_rerun), 3),
    role_sizes = paste(table(rep$role), collapse = "/"))

  for (w in c(FALSE, TRUE)) {
    cat("\n================", dataset, if (w) "valued REGE (REGE.for)" else "binary REGE (REGE.nm.for)", "================\n")
    new <- run_fixed(dataset, edges, w)
    m   <- merge(original, new, by = "id", suffixes = c("_original", "_fixed"))
    cat("Cross-tab, original roles (rows) vs fixed roles (columns):\n")
    print(table(original = m$role_original, fixed = m$role_fixed))
    results[[length(results) + 1]] <- data.frame(
      dataset = dataset, ties = if (w) "valued" else "binary",
      n_matched = nrow(m), adjusted_rand_vs_original = round(ari(m$role_original, m$role_fixed), 3),
      role_sizes = paste(table(new$role), collapse = "/"))
  }
}
comparison <- do.call(rbind, results)
print(comparison)
write.csv(comparison, file.path(out_dir, "comparison.csv"), row.names = FALSE)
