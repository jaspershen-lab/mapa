variable_info <-
  read.csv("data.csv")

###if the id type is ensembl, then map to uniprot and entrezid
colnames(variable_info) <- c("ensembl")
library(clusterProfiler)
library(org.Hs.eg.db)
other_id <-
  clusterProfiler::bitr(
    variable_info$ensembl,
    fromType = "ENSEMBL",
    toType = c("UNIPROT", "ENTREZID"),
    OrgDb = org.Hs.eg.db
  )
### remove duplicated rows
other_id <-
  dplyr::distinct(other_id, ENSEMBL, .keep_all = TRUE)
variable_info <-
  dplyr::left_join(variable_info,
                   other_id, by = c("ensembl" = "ENSEMBL"))
colnames(variable_info) <-
  c("ensembl", "uniprot", "entrezid")


##if the id type is uniprot, then map to ensembl and entrezid

colnames(variable_info) <- c("uniprot")
library(clusterProfiler)
library(org.Hs.eg.db)
other_id <-
  clusterProfiler::bitr(
    variable_info$uniprot,
    fromType = "UNIPROT",
    toType = c("ENSEMBL", "ENTREZID"),
    OrgDb = org.Hs.eg.db
  )
# remove duplicated rows
other_id <-
  dplyr::distinct(other_id, UNIPROT, .keep_all = TRUE)
variable_info <-
  dplyr::left_join(variable_info,
                   other_id, by = c("uniprot" = "UNIPROT"))
colnames(variable_info) <-
  c("uniprot", "ensembl", "entrezid")

##if the id type is entrezid, then map to ensembl and uniprot
colnames(variable_info) <- c("entrezid")
library(clusterProfiler)
library(org.Hs.eg.db)
other_id <-
  clusterProfiler::bitr(
    variable_info$entrezid,
    fromType = "ENTREZID",
    toType = c("ENSEMBL", "UNIPROT"),
    OrgDb = org.Hs.eg.db
  )
# remove duplicated rows
other_id <-
  dplyr::distinct(other_id, ENTREZID, .keep_all = TRUE)
variable_info <-
  dplyr::left_join(variable_info,
                   other_id,
                   by = c("entrezid" = "ENTREZID"))
colnames(variable_info) <-
  c("entrezid", "ensembl", "entrezid")
