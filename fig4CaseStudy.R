library(dplyr)
library(readxl)
library(org.Hs.eg.db)
library(reactome.db)
library(AnnotationDbi)
library(KEGGREST)
library(mapa)
load("C:\\Users\\cgxjd\\Documents\\WeChat Files\\wxid_b6mvr687a1e722\\FileStorage\\File\\2025-06\\ora_enriched_functional_module.rda")
control_data <- read_excel("D:/NTU/mapa_manuscript/2_data/control_data.xlsx", sheet = 1)
module_counts <- table(control_data$expected_module)
modules_to_keep <- names(module_counts[module_counts > 1])
control_data_filtered <- control_data[control_data$expected_module %in% modules_to_keep, ]

map_pathway_to_genes <- function(pid) {
  ## 1) GO 通路
  if (grepl("^GO", pid, perl = TRUE)) {
    go_genes <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys     = pid,
      columns  = "ENSEMBL", # 将 SYMBOL 改为 ENSEMBL
      keytype  = "GO"
    )
    return(paste(unique(go_genes$ENSEMBL), collapse = "/"))
  }

  ## 2) Reactome 通路（ID 形如 R-HSA-69306）
  if (grepl("^R-", pid, perl = TRUE)) {
    eid <- unlist(as.list(reactomePATHID2EXTID[pid]))
    gene_df <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys     = eid,
      columns  = "ENSEMBL", # 将 SYMBOL 改为 ENSEMBL
      keytype  = "ENTREZID"
    )
    return(paste(unique(gene_df$ENSEMBL), collapse = "/"))
  }

  ## 3) KEGG 通路（ID 形如 hsa00190）
  if (grepl("^hsa", pid, perl = TRUE)) {
    eid_raw <- keggLink("hsa", paste0("path:", pid))
    eid     <- sub("hsa:", "", unname(eid_raw))
    gene_df <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys     = eid,
      columns  = "ENSEMBL", # 将 SYMBOL 改为 ENSEMBL
      keytype  = "ENTREZID"
    )
    return(paste(unique(gene_df$ENSEMBL), collapse = "/"))
  }

  ## 4) 其它情况：返回 NA
  return(NA_character_)
}

# ── 2. 批量应用到数据框 ────────────────────────────────────────
control_data_filtered <- control_data_filtered %>%
  mutate(
    gene_symbols = vapply(id, map_pathway_to_genes, FUN.VALUE = character(1))
  )


randomize_modules <- function(df, module_col_name = "expected_module") {
  # 检查指定的列是否存在
  if (!module_col_name %in% names(df)) {
    stop(paste("数据框中找不到列:", module_col_name))
  }

  # 1. 统计每个模块的原始数量
  module_counts <- df %>%
    count(!!sym(module_col_name), name = "count")

  # 2. 创建新的带有“Random”前缀的模块名称
  new_modules <- paste0("Random ", module_counts[[module_col_name]])

  # 3. 根据原始数量重新构建一个包含所有新模块的向量
  # 这个向量将用于随机抽样
  all_new_modules_for_sampling <- rep(new_modules, times = module_counts$count)

  # 4. 随机分配这些新模块
  # sample函数默认是不放回抽样，这里我们要保证总数不变
  df[[module_col_name]] <- sample(all_new_modules_for_sampling)

  return(df)
}

control_data_randomized = randomize_modules(control_data_filtered,)
combined_data <- rbind(control_data_filtered, control_data_randomized)

formated_result <- aggregate(cbind(id, name) ~ expected_module,
                     data = combined_data,
                     FUN = function(x) paste(x, collapse = ";"))
formated_result$geneID <- aggregate(gene_symbols ~ expected_module,
                                    data = combined_data,
                                    FUN = function(x) paste(x, collapse = "/"))$gene_symbols

module_content_counts <- combined_data %>%
  group_by(expected_module) %>%
  summarise(module_content_number = n())


formated_result <- left_join(formated_result, module_content_counts, by = "expected_module")
colnames(formated_result)[1:3] = c("module","pathway_id","Description")

#加上一列$module_content_number
openai_key = "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"
enriched_functional_module@merged_module[["functional_module_result"]] = formated_result
processed_module = llm_interpret_module(enriched_functional_module, api_key=openai_key,embedding_output_dir = "D:\\NTU\\mapa_doc\\embedding")
