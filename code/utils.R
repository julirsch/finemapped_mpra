# Opposite of %in% operator
"%ni%" <- Negate("%in%")

# T/F if variant is in MHC region
in_MHC <- function(chromosome, position) {
  return(chromosome == "chr6" & 25e6 <= position & position <= 36e6)
}

# Read in MPRA data
read_mpra <- function(fullfile, cell_type, library) {
  read_delim(fullfile, delim  = "\t", col_names = T) %>%
    dplyr::filter(grepl("wC", ID)) %>%
    dplyr::mutate(variant = paste0("chr", SNP)) %>%
    rename_with(recode, Skew_logFDR = "Skew_logPadj") %>%
    dplyr::select(ID, variant,
                  A_Ctrl_Mean, A_Exp_Mean, A_log2FC, A_log2FC_SE, A_logPadj_BF,
                  B_Ctrl_Mean, B_Exp_Mean, B_log2FC, B_log2FC_SE, B_logPadj_BF,
                  Log2Skew, Skew_SE, Skew_logPadj) %>%
    group_by(ID) %>%
    dplyr::mutate(log2FC = ifelse(max(abs(A_log2FC) > abs(B_log2FC), na.rm = T), A_log2FC, B_log2FC),
                  log2FC_SE = ifelse(max(abs(A_log2FC) > abs(B_log2FC), na.rm = T), A_log2FC_SE, B_log2FC_SE),
                  logPadj_BF = max(A_logPadj_BF, B_logPadj_BF, na.rm = T)) %>%
    dplyr::mutate(cell_type = cell_type,
                  library = library) %>%
    ungroup()
}

# Read in MPRA data
read_mpra2 <- function(cell_types, libraries, min_DNA = 100) {
  pmap_dfr(
    crossing(cell_type = cell_types, library = libraries),
    function(cell_type, library) {
      tryCatch({
        read_delim(list.files(path = paste0("../results/"), pattern = paste0(library, "_", cell_type, "_emVAR_", ".*.out"), full.names = T), delim  = "\t", col_names = T) %>%
          dplyr::filter(grepl("wC", ID)) %>%
          dplyr::mutate(variant = paste0("chr", SNP)) %>%
          dplyr::select(ID, variant,
                        A.Ctrl.Mean, A.Exp.Mean, A.log2FC, A.log2FC_SE, A.logPadj_BF,
                        B.Ctrl.Mean, B.Exp.Mean, B.log2FC, B.log2FC_SE, B.logPadj_BF,
                        LogSkew, Skew.logP, Skew.logFDR) %>%
          group_by(ID) %>%
          dplyr::filter(A.Ctrl.Mean > min_DNA & B.Ctrl.Mean > min_DNA) %>%
          dplyr::mutate(log2FC = ifelse(max(abs(A.log2FC) > abs(B.log2FC), na.rm = T), A.log2FC, B.log2FC),
                        log2FC_SE = ifelse(max(abs(A.log2FC) > abs(B.log2FC), na.rm = T), A.log2FC_SE, B.log2FC_SE),
                        logPadj_BF = max(A.logPadj_BF, B.logPadj_BF, na.rm = T)) %>%
          dplyr::mutate(cell_type = cell_type,
                        library = library) %>%
          ungroup()
      }, error = function(e) {
        mpra_input[0,]
      })
    }
  )
}

# Read in annotation control variant info
read_ctrl_annot <- function() {
  read_delim(paste0("../mpra_design/GTEx/GTEx_annotcontrol_variants.df.txt.gz"), delim = " ", col_names = T) %>%
    dplyr::select(variant, variant_hg38, chromosome, position, indexVar, pip, type) %>% 
    dplyr::mutate(variant = gsub("_", ":", variant),
                  indexVar = gsub("_", ":", indexVar))
}

# Read in location control variant info
read_ctrl_loc <- function() {
  read_delim(paste0("../mpra_design/GTEx/GTEx_locationcontrol_variants.df.txt.gz"), delim = " ", col_names = T) %>%
    dplyr::mutate("chromosome" = seqnames, "position" = start) %>%
    dplyr::select(variant, variant_hg38, chromosome, position, indexVar, pip, type) %>%
    dplyr::mutate(variant = gsub("_", ":", variant),
                  indexVar = gsub("_", ":", indexVar))
}

# Change DHS size
GRsize <- function(gr, shift = 250) {
  midpt <- (start(gr) + floor(width(gr)/2))
  start(gr) <- as.integer(midpt - shift)
  end(gr) <- as.integer(midpt + shift)
  return(gr)
}

# Define PIP bin breaks
pip_bin_breaks <- c(0, 0.01, 0.1, 0.5, 0.9, 1.0)


