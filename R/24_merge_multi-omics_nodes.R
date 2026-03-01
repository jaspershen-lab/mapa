# setwd(r4projects::get_project_wd())
# load("demo_data/demo_multi-omics/multiplex_het.rda")
# diffusion_profile <- compute_diffusion_profiles(multiplex_het)
# save(diffusion_profile, file = "demo_data/demo_multi-omics/diffusion_profile.rda")
# Compute diffusion profiles for all nodes in a MultiplexHet network

compute_diffusion_profiles <- function(
    MultiplexHet_Object,
    TransitionMatrix = NULL,
    r = 0.7,
    eta = 0.5,
    lambda = 0.5,
    delta1 = 0.5,
    delta2 = 0.5,
    verbose = TRUE
) {

  # Compute transition matrix
  if (is.null(TransitionMatrix)) {
    if (verbose) message("[1/3] Computing transition matrix ...")
    TransitionMatrix <- compute.transition.matrix(
      MultiplexHet_Object,
      lambda = lambda,
      delta1 = delta1,
      delta2 = delta2
    )
  } else {
    if (verbose) message("[1/3] Using provided transition matrix.")
  }

  n_nan <- sum(is.nan(TransitionMatrix@x))
  n_inf <- sum(is.infinite(TransitionMatrix@x))
  if (n_nan > 0 || n_inf > 0) {
    if (verbose) message(sprintf("  Cleaning TransitionMatrix: %d NaN and %d Inf replaced with 0.", n_nan, n_inf))
    TransitionMatrix@x[is.nan(TransitionMatrix@x)] <- 0
    TransitionMatrix@x[is.infinite(TransitionMatrix@x)] <- 0
  }

  seeds_to_run  <- list(mol = MultiplexHet_Object$Multiplex1$Pool_of_Nodes,
                        pathway = MultiplexHet_Object$Multiplex2$Pool_of_Nodes)

  total_seeds <- sum(sapply(seeds_to_run, length))

  # Run RWR for each seed (node) -> extract global score vector

  run_single_seed <- function(seed, which_multiplex) {

    if (which_multiplex == "mol") {
      m1_seeds <- seed
      m2_seeds <- character(0)
    } else {
      m1_seeds <- character(0)
      m2_seeds <- seed
    }

    result <- tryCatch(
      Random.Walk.Restart.MultiplexHet.default(
        x = TransitionMatrix,
        MultiplexHet_Object = MultiplexHet_Object,
        Multiplex1_Seeds = m1_seeds,
        Multiplex2_Seeds = m2_seeds,
        r = r,
        eta = eta,
        DispResults = "Alphabetic"
      ),
      error = function(e) {
        warning(sprintf("RWR failed for seed '%s': %s", seed, e$message))
        return(NULL)
      }
    )
    if (is.null(result)) return(NULL)

    # global results
    global_df <- result$RWRMH_GlobalResults
    score_vec <- setNames(global_df$Score, global_df$NodeNames)
    return(score_vec)
  }

  profile_list <- list()
  counter <- 0

  for (net_type in names(seeds_to_run)) {
    seeds_vec <- seeds_to_run[[net_type]]
    for (seed in seeds_vec) {
      counter <- counter + 1
      if (verbose && (counter %% 50 == 0 || counter == total_seeds)) {
        message(sprintf("  ... %d / %d seeds done", counter, total_seeds))
      }
      profile_list[[seed]] <- run_single_seed(seed, net_type)
    }
  }

  failed <- names(which(sapply(profile_list, is.null)))
  if (length(failed) > 0) {
    warning(sprintf("%d seeds failed and were removed: %s",
                    length(failed), paste(failed, collapse = ", ")))
    profile_list <- profile_list[!names(profile_list) %in% failed]
  }

  # build profile matrix
  if (verbose) message("[3/3] Assembling profile matrix ...")

  all_col_names <- unique(unlist(lapply(profile_list, names)))

  profile_matrix <- do.call(rbind, lapply(profile_list, function(v) {
    out <- numeric(length(all_col_names))
    names(out) <- all_col_names
    out[names(v)] <- v
    out
  }))

  if (verbose) message(sprintf(
    "Done. Profile matrix: %d rows (seeds) x %d cols.",
    nrow(profile_matrix), ncol(profile_matrix)
  ))

  # return(list(
  #   profile_matrix   = profile_matrix,
  #   TransitionMatrix = TransitionMatrix
  # ))

  profile_matrix
}

.calculate_sim <- function(m) {

  expected_gb <- (nrow(m)^2 * 8) / 1e9
  if (expected_gb > 2) {
    warning(sprintf("Output matrix will be ~%.1f GB, consider subsetting.", expected_gb))
  }

  norms <- sqrt(rowSums(m^2))
  norms[norms == 0] <- NA

  m_norm <- m / norms
  cosine_sim <- m_norm %*% t(m_norm)
  return(cosine_sim)
}


