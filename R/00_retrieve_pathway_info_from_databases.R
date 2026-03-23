# # setwd(r4projects::get_project_wd())
#
# # This script retrieves pathway/term information from GO, KEGG, and Reactome
# # and outputs a single tab-delimited .txt file where each row is one pathway/term.
# #
# # Output columns:
# #   source_db, id, id_suffix, name, description, category, record_status
#
# library(KEGGREST)
# library(httr2)
# library(jsonlite)
# library(tidyverse)
#
# # 1. GO terms  (from go-basic.obo) ====
# parse_go_obo <- function(obo_file = "demo_data/go-basic.obo") {
#   lines <- readLines(obo_file)
#
#   # Split into stanzas by [Term] header
#   term_starts <- which(lines == "[Term]")
#   if (length(term_starts) == 0) stop("No [Term] entries found in OBO file.")
#
#   parse_stanza <- function(stanza_lines) {
#     get_field <- function(tag) {
#       hits <- grep(paste0("^", tag, ": "), stanza_lines, value = TRUE)
#       if (length(hits) == 0) return(NA_character_)
#       sub(paste0("^", tag, ": "), "", hits[1])
#     }
#
#     id_raw <- get_field("id")          # e.g. "GO:0000001"
#     name <- get_field("name")
#     namespace <- get_field("namespace")   # biological_process / molecular_function / cellular_component
#     def_raw   <- get_field("def")         # "description text" [refs]
#     is_obsolete   <- get_field("is_obsolete")  # "true" or NA
#
#     # Extract description text from def field (between the first pair of quotes)
#     description <- if (!is.na(def_raw)) {
#       m <- regmatches(def_raw, regexpr('"[^"]*"', def_raw))
#       if (length(m) > 0) gsub('^"|"$', "", m) else NA_character_
#     } else NA_character_
#
#     # record_status
#     record_status <- if (!is.na(is_obsolete) && is_obsolete == "true") "TRUE" else NA
#
#     tibble(
#       source_db  = "GO",
#       id = id_raw,
#       id_suffix = id_raw,
#       name = name,
#       description = description,
#       category = namespace,
#       record_status = record_status
#     )
#   }
#
#   # Determine stanza boundaries
#   end_indices <- c(term_starts[-1] - 1, length(lines))
#
#   results <- vector("list", length(term_starts))
#   for (i in seq_along(term_starts)) {
#     stanza <- lines[term_starts[i]:end_indices[i]]
#     results[[i]] <- parse_stanza(stanza)
#   }
#
#   bind_rows(results)
# }
#
# go_data <- parse_go_obo("demo_data/go-basic.obo")
#
# # 2. KEGG pathways  (via KEGGREST, reference pathways only) ====
#
# fetch_kegg_pathways <- function(batch_size = 10) {
#   message("Fetching KEGG pathway IDs ...")
#   pathway_list <- keggList("pathway")
#   pathway_ids  <- names(pathway_list)
#
#   message(sprintf("  -> %d reference KEGG pathways found.", length(pathway_ids)))
#
#   # Batch-fetch (KEGGREST allows up to 10 per call)
#   batches <- split(pathway_ids, ceiling(seq_along(pathway_ids) / batch_size))
#
#   records <- vector("list", length(batches))
#   for (i in seq_along(batches)) {
#     if (i %% 10 == 1) message(sprintf("  Fetching batch %d / %d ...", i, length(batches)))
#     batch <- tryCatch(keggGet(batches[[i]]), error = function(e) NULL)
#     if (is.null(batch)) next
#     records[[i]] <- lapply(batch, function(entry) {
#       entry_id <- if (!is.null(entry$ENTRY)) unname(entry$ENTRY)[1] else NA_character_
#       name_val <- entry$NAME
#       desc_val <- entry$DESCRIPTION
#
#       tibble(
#         source_db = "KEGG",
#         id = entry_id,
#         id_suffix = sub("map", "", entry_id),
#         name = if (!is.null(name_val)) paste(name_val, collapse = "; ") else NA_character_,
#         description = if (!is.null(desc_val)) paste(desc_val, collapse = " ") else "",
#         category = NA,
#         record_status = NA
#       )
#     })
#   }
#
#   bind_rows(unlist(records, recursive = FALSE))
# }
#
# kegg_data <- fetch_kegg_pathways()
# message(sprintf("  -> %d KEGG pathways retrieved.", nrow(kegg_data)))
#
# # 3. Reactome pathways  (from ReactomePathways.txt + REST API) ====
# #
# # Strategy:
# #   a) Query all human (Homo sapiens) pathways first; store description by id_suffix.
# #   b) For every other species, query each pathway:
# #        - isInferred == FALSE -> use retrieved description directly.
# #        - isInferred == TRUE  -> fall back to the human description for the same id_suffix;
# #                                 all other fields still come from the retrieved record.
# #   c) Sleep `sleep_sec` seconds between every API call to respect rate limits.
#
# fetch_reactome_pathways <- function(
#     reactome_file = "demo_data/ReactomePathways.txt",
#     sleep_sec     = 0.5
# ) {
#   # Read pathway list (no header: pathway_id, name, species)
#   rp <- read.table(reactome_file, sep = "\t", header = FALSE,
#                    col.names = c("pathway_id", "name", "species"),
#                    stringsAsFactors = FALSE, quote = "")
#
#   # id_suffix = numeric part of the pathway id (e.g. "73843" from "R-HSA-73843")
#   rp$id_suffix <- sub("^R-[A-Z]+-", "", rp$pathway_id)
#
#   # Helper: query Reactome content service for one pathway id
#   get_pathway_events <- function(pathway_id) {
#     url <- paste0("https://reactome.org/ContentService/data/query/", pathway_id)
#     tryCatch({
#       response <- request(url) |>
#         req_headers(accept = "*/*") |>
#         req_perform()
#       resp_body_json(response, simplifyVector = TRUE)
#     }, error = function(e) NULL)
#   }
#
#   # Helper: extract description text from a query result
#   extract_description <- function(result) {
#     if (is.null(result) || is.null(result$summation)) return(NA_character_)
#     sums <- result$summation
#     if (is.data.frame(sums) && "text" %in% names(sums)) {
#       paste(sums$text, collapse = " ")
#     } else if (is.character(sums)) {
#       paste(sums, collapse = " ")
#     } else NA_character_
#   }
#
#   # ---- Step a: fetch all human pathways ----
#   rp_human <- rp[rp$species == "Homo sapiens", ]
#   message(sprintf("Fetching %d human Reactome pathways ...", nrow(rp_human)))
#
#   human_desc <- list()   # named by id_suffix
#   for (i in seq_len(nrow(rp_human))) {
#     if (i %% 100 == 1) message(sprintf("  Human pathway %d / %d ...", i, nrow(rp_human)))
#     pid   <- rp_human$pathway_id[i]
#     suffix <- rp_human$id_suffix[i]
#     result <- get_pathway_events(pid)
#     human_desc[[suffix]] <- extract_description(result)
#     Sys.sleep(sleep_sec)
#   }
#
#   # Build records for human pathways
#   human_records <- vector("list", nrow(rp_human))
#   for (i in seq_len(nrow(rp_human))) {
#     suffix <- rp_human$id_suffix[i]
#     human_records[[i]] <- tibble(
#       source_db = "Reactome",
#       id = rp_human$pathway_id[i],
#       id_suffix = suffix,
#       name = rp_human$name[i],
#       description = human_desc[[suffix]],
#       category = NA,
#       record_status = NA
#     )
#   }
#
#   # ---- Step b: fetch all non-human pathways ----
#   rp_other <- rp[rp$species != "Homo sapiens", ]
#   message(sprintf("Fetching %d non-human Reactome pathways ...", nrow(rp_other)))
#
#   other_records <- vector("list", nrow(rp_other))
#   for (i in seq_len(nrow(rp_other))) {
#     if (i %% 100 == 1) message(sprintf("  Non-human pathway %d / %d ...", i, nrow(rp_other)))
#     pid    <- rp_other$pathway_id[i]
#     suffix <- rp_other$id_suffix[i]
#     result <- get_pathway_events(pid)
#
#     is_inferred <- if (!is.null(result) && !is.null(result$isInferred)) {
#       isTRUE(result$isInferred)
#     } else FALSE
#
#     description <- if (!is_inferred) {
#       # Use the description retrieved directly for this pathway
#       extract_description(result)
#     } else {
#       # Fall back to the human pathway description with the same id_suffix
#       if (!is.null(human_desc[[suffix]])) human_desc[[suffix]] else NA_character_
#     }
#
#     other_records[[i]] <- tibble(
#       source_db = "Reactome",
#       id = pid,
#       id_suffix = suffix,
#       name = rp_other$name[i],
#       description = description,
#       category = NA,
#       record_status = NA
#     )
#     Sys.sleep(sleep_sec)
#   }
#
#   bind_rows(c(human_records, other_records))
# }
#
# reactome_data <- fetch_reactome_pathways("demo_data/ReactomePathways.txt")
# message(sprintf("  -> %d Reactome pathways retrieved.", nrow(reactome_data)))
#
# save(reactome_data, file = "demo_data/reactome_data.rda")
#
# sum(is.null(reactome_data$description))
# sum(is.na(reactome_data$description))
# which(is.na(reactome_data$description))
# # 4. Combine and write output ====
# all_pathways <- bind_rows(go_data, kegg_data, reactome_data)
#
# output_file <- "demo_data/pathway_database_info.txt"
# write.table(
#   all_pathways,
#   file      = output_file,
#   sep       = "\t",
#   row.names = FALSE,
#   col.names = TRUE,
#   quote     = FALSE,
#   na        = "NA"
# )
#
# message(sprintf("Done. Output written to: %s (%d rows)", output_file, nrow(all_pathways)))
