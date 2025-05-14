# Library
library(tidyverse)
library(MotifDb)
library(motifbreakR)
library(universalmotif)
library(TFBSTools)
library(memes)
`%ni%` = Negate(`%in%`)

# Hocomoco motifbreakR (old)
data('hocomoco')
hocomoco.mm <- convert_motifs(hocomoco, "TFBSTools-PFMatrix")
hocomoco.mm.name <- unlist(lapply(hocomoco.mm, function(x) {x@ID}))
hocomoco.mm <- do.call(PFMatrixList, hocomoco.mm)
names(hocomoco.mm) <- hocomoco.mm.name

# Hocomoco v12 (new)
hocomoco_meta <- read.table("homomoco_meta.txt", header = T) %>%
  as_tibble()
hocomoco_meme <- read_meme("../Downloads/H12CORE_meme_format.meme")
hocomoco_mdb <- convert_motifs(hocomoco_meme, "MotifDb-MotifList")
hocomoco_pfml <- convert_motifs(hocomoco_mdb, "TFBSTools-PFMatrix")
hocomoco_pfml_name <- unlist(lapply(hocomoco_pfml, function(x) {x@ID}))
hocomoco_pfml <- do.call(PFMatrixList, hocomoco_pfml)
names(hocomoco_pfml) <- hocomoco_pfml_name

# Cis-BP
df <- data.table::fread("../Downloads/Homo_sapiens_2023_10_25_4_59_pm/TF_Information_all_motifs_plus.txt")
colnames(df)[6] <- "ENSGID"
cisbp_meta <- read.table("cisbp_meta.txt", header = T) %>%
  as_tibble()
cisbp_meta$TF_status <- df$TF_Status
humanTF_meta <- data.table::fread("../Downloads/Human_TF_MotifList_v_1.01.csv") %>%
  dplyr::mutate(Other_ID = `Motif ID`,
                ENSGID = `Ensembl ID`) %>%
  dplyr::filter(`Best Motif(s)? (Figure 2A)` == T)
cisbp_meta %>% 
  filter(TF_status == "D")
files <- list.files("../Downloads/Homo_sapiens_2023_10_25_4_59_pm/pwms_all_motifs/")
cisbp <- lapply(1:length(files), function(x) {
  print(files[x])
  in_motif <- read_cisbp(paste0("../Downloads/Homo_sapiens_2023_10_25_4_59_pm/pwms_all_motifs/", files[x]))
  in_motif@name <- gsub(".txt", "", files[x])
  return(in_motif)
  })
cisbp_mdb <- convert_motifs(cisbp, "MotifDb-MotifList")
cisbp_pfml <- convert_motifs(cisbp_mdb, "TFBSTools-PFMatrix")
cisbp_pfml_name <- unlist(lapply(cisbp_pfml, function(x) {x@ID}))
cisbp_pfml <- do.call(PFMatrixList, cisbp_pfml)
names(cisbp_pfml) <- cisbp_pfml_name

# Vierstra
vierstra_df <- readRDS("Vierstra_Archetype_Motifs_v2.1.rds")
for (x in 1:length(vierstra_df)) {
  vierstra_df[[x]]@strand <- "+"
}
vierstra_pfml <- convert_motifs(vierstra_df, "TFBSTools-PFMatrix")
vierstra_pfml.name <- unlist(lapply(vierstra_pfml, function(x) {x@ID}))
vierstra_pfml <- do.call(PFMatrixList, vierstra_pfml)
names(vierstra_pfml) <- vierstra_pfml.name
vierstra_mdb <- convert_motifs(vierstra_pfml, "MotifDb-MotifList")

# Jaspar 2022
jaspar_meme <- read_meme("../Downloads/JASPAR2022_CORE_non-redundant_pfms_meme.txt")
jaspar_mdb <- convert_motifs(jaspar_meme, "MotifDb-MotifList")
jaspar_pfml <- convert_motifs(jaspar_mdb, "TFBSTools-PFMatrix")
jaspar_pfml_name <- unlist(lapply(jaspar_pfml, function(x) {x@ID}))
jaspar_pfml <- do.call(PFMatrixList, jaspar_pfml)
names(jaspar_pfml) <- jaspar_pfml_name

# Add names
for (x in 1:length(cisbp_pfml)) {
  cisbp_pfml[x]@listData[[1]]@name <- cisbp_pfml[x]@listData[[1]]@ID
  cisbp_pfml[x]@listData[[1]]@strand <- "+-"
}
for (x in 1:length(hocomoco_pfml)) {
  hocomoco_pfml[x]@listData[[1]]@name <- hocomoco_pfml[x]@listData[[1]]@ID
  hocomoco_pfml[x]@listData[[1]]@strand <- "+-"
}
for (x in 1:length(jaspar_pfml)) {
  jaspar_pfml[x]@listData[[1]]@name <- jaspar_pfml[x]@listData[[1]]@ID
  jaspar_pfml[x]@listData[[1]]@strand <- "+-"
}
for (x in 1:length(vierstra_pfml)) {
  vierstra_pfml[x]@listData[[1]]@name <- vierstra_pfml[x]@listData[[1]]@ID
  vierstra_pfml[x]@listData[[1]]@strand <- "+-"
}

# Get meta data (TF labels)
hocomoco_tfs <- bind_cols("motif" = names(hocomoco_pfml), "TF" = names(hocomoco_pfml) %>% as_tibble() %>% inner_join(hocomoco_meta, by = c("value" = "motif")) %>% .$TF, motif_db = "hocomoco")
jaspar_tfs <- bind_cols("motif" = lapply(1:length(jaspar_meme), function(x) {jaspar_meme[[x]]@name}) %>% unlist(), "TF" = lapply(1:length(jaspar_meme), function(x) {jaspar_meme[[x]]@altname}) %>% unlist(), motif_db = "jaspar")
vierstra_tfs <- bind_cols("motif" = names(vierstra_pfml), "TF" = names(vierstra_pfml) %>% as_tibble() %>% separate(value, c("ID", "TF", "description"), ":") %>% .$TF, motif_db = "vierstra")
cisbp_tfs <- cisbp_meta %>% distinct("motif" = Motif_ID, TF, motif_db = "cisbp")
motifs_meta <- bind_rows(hocomoco_tfs, jaspar_tfs, vierstra_tfs, cisbp_tfs) %>%
  dplyr::mutate(TF = toupper(TF)) %>%
  distinct()
motifs_meta <- motifs_meta %>%
  dplyr::filter(motif %in% c(names(cisbp_pfml), names(hocomoco_pfml), names(jaspar_pfml), names(vierstra_pfml)))
motifs_mdb <- c(hocomoco_mdb, jaspar_mdb, cisbp_mdb, vierstra_mdb)
motifs_meta <- motifs_meta %>% 
  inner_join(motifs_mdb %>% 
               to_df() %>% 
               mutate(width = nchar(consensus)) %>% 
               select("motif" = altname, width, icscore))

save(jaspar_mdb, jaspar_pfml, 
     hocomoco_mdb, hocomoco_pfml, 
     cisbp_mdb, cisbp_pfml, 
     vierstra_mdb, vierstra_pfml,
     motifs_meta, 
     file = "motifs_processed_20230225.RData")

# Convert to meme 
load("motifs_processed_20230225.RData")
write_meme(cisbp_pfml, "cisbp.meme", strand = "+ -", overwrite = T)
write_meme(hocomoco_pfml, "hocomoco.meme", strand = "+ -", overwrite = T)
write_meme(jaspar_pfml, "jaspar.meme", strand = "+ -", overwrite = T)
system("sed -i 's/1e+05/10000/g' jaspar.meme")
write_meme(vierstra_pfml, "vierstra.meme", strand = "+ -", overwrite = T)
motifs_pfml <- c(hocomoco_pfml, jaspar_pfml, cisbp_pfml, vierstra_pfml)
write_meme(motifs_pfml, "motifs.meme", strand = "+ -", overwrite = T)
system("sed -i 's/1e+05/100000/g' motifs.meme")


# Merge motifs
universalmotif::compare_motifs(hocomoco_pfml[1:10], method = "PCC", tryRC = T)
