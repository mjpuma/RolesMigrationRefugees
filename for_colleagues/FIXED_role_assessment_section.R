##### ==========================================================================
##### FIXED ROLE-ASSESSMENT SECTION
##### Drop-in replacement for lines 61-140 of FINAL_undesa_role_assessment.R
##### (from "##### Migrant networks - strict cutoff - unweighted #####"
#####  down to   clusters <- read.csv("results/UNDESA_full_roles_k3.csv")[,2:3])
#####
##### The same block, with "UNDESA" -> "UNHCR" and edges.mig -> edges.ref.long,
##### replaces lines 61-131 of FINAL_unhcr_role_assessment.R.
#####
##### What was wrong in the original block:
#####  (1) hcut(x = a9, method = link): hcut() has no 'method' argument, so it
#####      was silently ignored and Ward.D2 was used, not "average".
#####  (2) a9 is a plain matrix, so hcut()/fviz_nbclust() treated it as DATA and
#####      clustered Euclidean distances between its rows, not the REGE distances.
#####  (3) The k diagnostics (elbow/silhouette) were run on a different clustering
#####      than the final cut, so they did not justify k = 3.
#####  (4) The network was binary (>= 1 migrant = tie). Almost every country both
#####      sends and receives, so REGE on 0/1 ties carries almost no information
#####      and the resulting roles track degree (number of partners), not flows.
#####
##### Needs: igraph, blockmodeling (already loaded above) + cluster (base R install).
##### Everything downstream (the "Network Metrics" section) still works, because
##### this block produces the same objects: c, clusters, and the results CSV.
##### ==========================================================================

library(cluster)                                   # silhouette(); ships with every R installation

# ---- Settings (change here, nowhere else) -----------------------------------
DATASET     <- "UNDESA"                            # label used in plot titles and file names
USE_WEIGHTS <- TRUE                                # TRUE = REGE on flow values; FALSE = original 0/1 ties
MIN_FLOW    <- 1                                   # a tie exists when flow >= MIN_FLOW people (as before)
REGE_ITER   <- 3                                   # REGE iterations; state this in the Methods (results depend on it)
LINK        <- "average"                           # linkage method the manuscript reports
K_RANGE     <- 2:25                                # k values scanned (say "2-25" in the paper, not "2 to n-1")
K_FINAL     <- 3                                   # final k; set ONLY after inspecting the diagnostics below
RESULTS_DIR <- "results"                           # where outputs are written (as in the original scripts)
OUT_TAG     <- ""                                  # added to file names; "" keeps the original names

# ---- 1. Build the stacked country-period network -----------------------------
a <- na.omit(edges.mig[edges.mig[, 3] >= MIN_FLOW, ])  # keep dyads with at least MIN_FLOW migrants; drop NAs
colnames(a)[3] <- "weight"                         # third column = flow size; igraph reads it as edge weight
a <- a[sub("\\.\\d{4}$", "", a[[1]]) !=           # drop self-loops (origin = destination; columns 1-2
       sub("\\.\\d{4}$", "", a[[2]]), ]            #   are "Country.Year"), which REGE should never see
g <- graph_from_data_frame(a, directed = TRUE)     # directed graph; node names are "Country.Year"
g <- igraph::delete_vertices(g,                    # remove isolated nodes (no ties at all);
       V(g)[igraph::degree(g) == 0])               #   same step as the original a4

if (USE_WEIGHTS) {                                 # valued network:
  M <- as.matrix(as_adjacency_matrix(g, attr = "weight",
                                     sparse = FALSE))  # n x n matrix of flow sizes (0 = no tie)
  M <- log1p(M)                                    # log(1 + flow): flows span 1 to millions; without the log
                                                   #   a handful of huge corridors dominate every comparison
} else {                                           # binary network (original behaviour):
  M <- as.matrix(as_adjacency_matrix(g, sparse = FALSE))  # 1 = tie, 0 = no tie
  M[M > 0] <- 1                                    # guard against multi-edges counting as 2, 3, ...
}
diag(M) <- 0                                       # no self-ties (a country-year is not its own partner)

# ---- 2. Regular equivalence --------------------------------------------------
# REGE.for is White's classic REGGE and is designed for VALUED ties: two ties
# match in proportion to how similar their sizes are. REGE.nm.for (the original
# choice) is a variant that rescales ties by each node's row/column totals; with
# 0/1 input that rescaling is essentially a degree comparison, which is why the
# original migrant roles are predicted ~98% by degree alone.
rege_fun  <- if (USE_WEIGHTS) REGE.for else REGE.nm.for  # valued REGGE, or the original function for comparison
iter_name <- intersect(c("iter", "it"),            # name of the iterations argument differs between
                       names(formals(rege_fun)))[1] #   blockmodeling versions; find the one this version uses
rege_args <- list(M = M, E = 1)                    # network matrix and starting similarity (all 1s)
if (!is.na(iter_name)) rege_args[[iter_name]] <- REGE_ITER  # fixed number of iterations, if settable
rege_out  <- do.call(rege_fun, rege_args)          # run REGE
S <- rege_out$E                                    # REGE similarity matrix, 0 (different) .. 1 (equivalent)
dimnames(S) <- dimnames(M)                         # keep "Country.Year" labels on rows and columns
a9 <- 1 - S                                        # dissimilarity = 1 - similarity (same as original a9)

# Sanity check: if nearly all REGE scores are 1, there is no role structure to find
share_equiv <- mean(S[upper.tri(S)] > 0.999)       # share of pairs scored as (almost) perfectly equivalent
cat("Share of pairs with REGE ~ 1:", round(share_equiv, 3), "\n")  # report it
if (share_equiv > 0.9) warning("REGE is almost all 1s: roles below will be driven by noise/degree, not structure")

# ---- 3. Hierarchical clustering ON THE REGE DISTANCES -----------------------
d  <- as.dist(a9)                                  # tell R these ARE distances (the key fix: not "data" rows)
hc <- hclust(d, method = LINK)                     # average-linkage clustering, as the Methods state
plot(hc, labels = FALSE,                           # dendrogram (labels off: 1,000+ nodes are unreadable)
     main = paste(DATASET, "- REGE dendrogram"), sub = "", xlab = "")

# ---- 4. Choose k using diagnostics computed on THE SAME clustering ----------
Dm <- as.matrix(d)                                 # full distance matrix, needed for the WCSS formula
diag_k <- do.call(rbind, lapply(K_RANGE, function(k) {  # one row of diagnostics per candidate k
  cl  <- cutree(hc, k = k)                         # cut the SAME tree that gives the final roles
  sil <- mean(silhouette(cl, d)[, "sil_width"])    # average silhouette width (higher = better separated)
  wcss <- sum(sapply(unique(cl), function(j) {     # within-cluster sum of squares from distances:
    idx <- which(cl == j)                          #   members of cluster j
    sum(Dm[idx, idx]^2) / (2 * length(idx))        #   sum of squared pairwise distances / (2 n_j)
  }))
  data.frame(k = k, avg_silhouette = sil, wcss = wcss,
             smallest_cluster = min(table(cl)))    # flags "clusters" of 1-2 nodes (outliers, not roles)
}))
print(diag_k)                                      # inspect: silhouette peak, WCSS elbow, cluster sizes
write.csv(diag_k, file.path(RESULTS_DIR, paste0(DATASET, OUT_TAG, "_k_diagnostics.csv")),
          row.names = FALSE)                       # keep the table for the Supplement

par(mfrow = c(2, 1))                               # two stacked panels
plot(diag_k$k, diag_k$avg_silhouette, type = "b",  # silhouette vs k
     xlab = "k", ylab = "Average silhouette", main = paste(DATASET, "- silhouette (higher = better)"))
plot(diag_k$k, diag_k$wcss, type = "b",            # WCSS vs k (look for the elbow)
     xlab = "k", ylab = "WCSS", main = paste(DATASET, "- elbow"))
par(mfrow = c(1, 1))                               # reset plotting layout
# Rule of thumb for the text: average silhouette < 0.25 means no substantial
# cluster structure; report the value for the k you choose.

# ---- 5. Final role assignment -------------------------------------------------
b <- cutree(hc, k = K_FINAL)                       # roles from the SAME average-linkage tree
print(table(b))                                    # cluster sizes; singletons are outliers, not roles
c <- data.frame("id" = names(b), "role" = b)       # same structure as the original output

# Are the roles just degree? (the concern with the original migrant roles)
deg_tab <- data.frame(role = b,                    # role of each country-year
                      degree = igraph::degree(g, v = names(b), mode = "all"))  # its number of ties
print(aggregate(degree ~ role, deg_tab,            # min / median / max degree per role:
                function(x) c(min = min(x), median = median(x), max = max(x))))
# If the degree ranges barely overlap across roles, the roles are degree classes.

roles_file <- file.path(RESULTS_DIR, paste0(DATASET, OUT_TAG, "_full_roles_k", K_FINAL, ".csv"))  # same name pattern as before
write.csv(c, roles_file)                           # save role assignments
clusters <- read.csv(roles_file)[, 2:3]            # used by the Network Metrics section
