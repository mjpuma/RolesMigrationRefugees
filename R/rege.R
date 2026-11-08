# =============================================================================
# rege.R — regular equivalence (REGE) for directed networks, written out in
# plain R so every step can be read and checked.
#
# Idea (Borgatti, Everett & Freeman 1996): two countries play the same role if
# they are "equally related to equivalent others". Country i sends migrants to
# some set of countries; country j sends to another set. i and j are similar if
# every destination of i can be matched to an equally-similar destination of j
# (and vice versa), and the same holds for their origins.
#
# Iteration (E starts as all 1s):
#
#   out-term(i,j) = sum_{k in out(i)} max_{m in out(j)} E[k,m]
#                 + sum_{m in out(j)} max_{k in out(i)} E[k,m]
#   in-term(i,j)  = same with in-neighbours
#   E_new[i,j]    = (out-term + in-term) / (|out(i)|+|out(j)|+|in(i)|+|in(j)|)
#
# E_new = 1 when i and j can match every tie perfectly; it falls towards 0 as
# ties go unmatched. A pure sender compared to a pure receiver scores 0.
#
# With several tie thresholds (e.g. >=1, >=1,000, >=100,000 people), each
# threshold is treated as its own binary relation and the terms are summed;
# ties must then be matched at every magnitude level ("multi-threshold REGE").
#
# Speed trick: for a fixed j, max_{m in out(j)} E[k,m] is a column vector over
# all k (matrix V below). Then sum_{k in out(i)} V[k,j] for all i at once is
# simply the matrix product A %*% V.
#
# Worked example of the speed trick. Let A be the 0/1 adjacency matrix
# (A[i,k] = 1 if i sends to k). Then
#   (A %*% V)[i, j] = sum_k A[i,k] V[k,j] = sum_{k in out(i)} max_{m in out(j)} E[k,m]
# which is exactly the first half of out-term(i,j). The second half is the
# same quantity with i and j swapped, i.e. the transpose t(A %*% V). The
# in-terms use t(A) (columns = senders) in place of A.
#
# Nodes are COUNTRY-PERIODS (e.g. "MEX_1995"), not countries: all periods are
# stacked into one network (see roles.R), so a country can be compared with
# itself or with any other country in any period, and roles are comparable
# over time.
# =============================================================================

# V[k, j] = max over neighbours m of j of E[k, m]   (0 if j has no neighbours)
#   E           current n x n similarity matrix
#   neighbours  list of length n; neighbours[[j]] = indices of j's out- (or in-)
#               neighbours under one threshold
# A node with no neighbours gets a column of 0s: none of its partner's ties can
# be matched, so they contribute nothing to the numerator.
best_match_matrix <- function(E, neighbours) {
  n <- nrow(E)
  V <- matrix(0, n, n)
  for (j in seq_len(n)) {
    nb <- neighbours[[j]]
    # one neighbour: the "max" is simply that column of E
    if (length(nb) == 1L) V[, j] <- E[, nb]
    # several neighbours: row-wise maximum over those columns (fast C code)
    else if (length(nb) > 1L) V[, j] <- matrixStats::rowMaxs(E[, nb, drop = FALSE])
  }
  V
}

#' @param A_list list of binary adjacency matrices (one per threshold), same node order
#' @param iterations number of REGE updates (3 in the original analysis; the
#'        returned `change` shows how much the last update still moved E)
#' @return list(similarity = n x n matrix in [0,1], change = max change per iteration)
rege <- function(A_list, iterations = 3, verbose = TRUE) {
  # a single matrix is accepted and treated as one threshold
  if (is.matrix(A_list)) A_list <- list(A_list)
  n <- nrow(A_list[[1]])
  # For each threshold ("relation"): make the matrix strictly 0/1, remove
  # self-ties, and pre-compute every node's out- and in-neighbour lists.
  rel <- lapply(A_list, function(A) {
    A <- (A > 0) * 1; diag(A) <- 0
    list(A = A,
         out = lapply(seq_len(n), function(i) which(A[i, ] > 0)),   # destinations of i
         inn = lapply(seq_len(n), function(i) which(A[, i] > 0)))   # origins of i
  })
  # denominator: total number of ties the pair has to match, summed over
  # thresholds: |out(i)| + |out(j)| + |in(i)| + |in(j)|. outer(x, x, "+")
  # gives the n x n matrix of x[i] + x[j].
  den <- Reduce(`+`, lapply(rel, function(r) {
    d_out <- lengths(r$out); d_in <- lengths(r$inn)
    outer(d_out, d_out, "+") + outer(d_in, d_in, "+")
  }))

  # start: every pair fully equivalent
  E <- matrix(1, n, n, dimnames = dimnames(A_list[[1]]))
  change <- numeric(0)
  for (it in seq_len(iterations)) {
    num <- matrix(0, n, n)
    for (r in rel) {
      # P_out[i,j] = sum over i's destinations k of the best match among j's destinations
      P_out <- r$A    %*% best_match_matrix(E, r$out)   # i's destinations vs j's
      # P_in[i,j]  = the same for origins (t(A)[i,k] = 1 if k sends to i)
      P_in  <- t(r$A) %*% best_match_matrix(E, r$inn)   # i's origins vs j's
      # add both directions of matching (i->j and j->i) for out and in ties
      num <- num + P_out + t(P_out) + P_in + t(P_in)
    }
    # pairs with no ties at all (den = 0) are defined as fully equivalent;
    # isolated nodes do not occur in practice because nodes come from flows
    E_new <- ifelse(den > 0, num / den, 1)
    dimnames(E_new) <- dimnames(E)
    # largest change of any coefficient in this iteration (convergence check)
    change[it] <- max(abs(E_new - E))
    E <- E_new
    if (verbose) msg(sprintf("    REGE iteration %d: max change %.4f", it, change[it]))
  }
  list(similarity = E, change = change)
}

# -----------------------------------------------------------------------------
# Valued REGE (White's REGGE) — REGE computed on the flow values themselves.
# Ties are matched in proportion to how similar their sizes are, so the
# coefficients reflect both the pattern and the magnitude of flows.
# The C++ core is in src/rege_valued.cpp (needs Rcpp + a C++ compiler; on a Mac
# install Xcode command-line tools once: `xcode-select --install`).
# Used only for the robustness comparison (analysis/migrant_volume_robustness.R);
# the main roles use the multi-threshold binary REGE above.
# -----------------------------------------------------------------------------
rege_valued <- function(W, iterations = 3, transform = c("log1p", "none"),
                        no_self_match = FALSE) {
  transform <- match.arg(transform)
  # compile the C++ function the first time it is needed in an R session
  if (!exists("rege_valued_cpp", mode = "function")) {
    src <- file.path(dirname(sys.frame(1)$ofile %||% "R/rege.R"), "..", "src", "rege_valued.cpp")
    if (!file.exists(src)) src <- "src/rege_valued.cpp"
    Rcpp::sourceCpp(src, env = globalenv())
  }
  # log1p compresses the very skewed flow sizes so that a 100,000 vs 200,000
  # difference counts about as much as 100 vs 200
  X <- if (transform == "log1p") log1p(W) else W
  diag(X) <- 0
  res <- rege_valued_cpp(X, iterations, no_self_match)
  dimnames(res$similarity) <- dimnames(W)
  res
}

# Weighted stacked matrix (country-period x country-period) of flow values.
# Rows are senders, columns receivers; a flow in period p only links nodes of
# period p (e.g. MEX_1995 -> USA_1995), so the matrix is block-diagonal by period.
build_stacked_valued <- function(flows, threshold = 1) {
  e <- flows |> filter(flow >= threshold)
  nodes <- sort(unique(c(node_id(e$origin, e$period), node_id(e$destination, e$period))))
  W <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  # cbind(row index, column index) fills all edges in one vectorised assignment
  W[cbind(match(node_id(e$origin, e$period), nodes),
          match(node_id(e$destination, e$period), nodes))] <- e$flow
  W
}
