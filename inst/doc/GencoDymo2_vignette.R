## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.width = 8,
  fig.height = 6
)

## ----install, eval=FALSE------------------------------------------------------
# # Install pak if not already installed
# if (!require("pak")) install.packages("pak")
# # Install from GitHub
# pak::pkg_install("github::monahton/GencoDymo2")
# 
# # Load the package
# library(GencoDymo2)

## ----get_release, eval=FALSE--------------------------------------------------
# # Fetch the most recent human and mouse GENCODE release identifiers
# human_release <- get_latest_release("human", verbose = T)
# mouse_release <- get_latest_release("mouse", verbose = T)

## ----get_release_ex, echo=FALSE-----------------------------------------------
# Get latest human and mouse release
cat("Latest human GENCODE release: release_47")
cat("Latest human GENCODE release: release_M36") 

## ----get-files, eval=FALSE----------------------------------------------------
# 
# # Download latest human long noncoding RNAs GTF
# lnc_47_gtf <- get_gtf(
#   species = "human",
#   release_version = human_release,
#   annotation_type = "long_noncoding_RNAs.gtf.gz",
#   dest_folder = tempdir()
# )
# 
# # Download previous human release (release_46) for comparison
# lnc_46_gtf <- get_gtf(
#   species = "human",
#   release_version = "release_46",
#   annotation_type = "long_noncoding_RNAs.gtf.gz",
#   dest_folder = tempdir()
# )
# 
# # Download latest mouse primary assembly annotations (GFF3)
# mouse_36_gff3 <- get_gff3(
#   species = "mouse",
#   release_version = mouse_release,
#   annotation_type = "primary_assembly.annotation.gff3.gz",
#   dest_folder = tempdir()
# )

## ----annotation-types, echo=FALSE---------------------------------------------
cat("Valid Annotation Types:\n")
valid_annotation_types <- c(
    "annotation",
    "basic.annotation",
    "chr_patch_hapl_scaff.annotation",
    "chr_patch_hapl_scaff.basic.annotation",
    "long_noncoding_RNAs",
    "primary_assembly.annotation",
    "primary_assembly.basic.annotation",
    "tRNAs",
    "polyAs")
valid_annotation_types


## ----load-data, eval=FALSE----------------------------------------------------
# # Loading using the stored paths from previous steps
# lnc_47_df <- load_file(lnc_47_gtf)
# head(lnc_47_df)
# 
# # Alternatively, specify the file path directly
# lnc_46_df <- load_file(file.path(tempdir(), "gencode.v46.long_noncoding_RNAs.gtf.gz"))
# head(lnc_46_df)
# 
# # Load mouse GFF3
# mouse_pri_36 <- load_file(file.path(tempdir(),"gencode.vM36.primary_assembly.annotation.gff3.gz"))
# head(mouse_pri_36)

## ----compare-releases, eval=FALSE---------------------------------------------
# # Compare gene counts between release 47 and 46
# gene_comparison <- compare_release(lnc_47_df, lnc_46_df, type = "gene")
# 
# # Compare exon counts
# exon_comparison <- compare_release(lnc_47_df, lnc_46_df, type = "exon")
# 
# # Compare a specific gene biotype (e.g., TEC) using a custom baseline
# comparison <- compare_release(
#   lnc_47_df,
#   lnc_46_df,
#   type = "gene",
#   gene_biotype = "TEC",
#   baseline = "count1"
# )

## ----introns, eval=FALSE------------------------------------------------------
# # Human lncRNA introns for release 47
# introns_lnc_47 <- extract_introns(lnc_47_df, verbose = T)
# 
# # Mouse introns (filtering to primary chromosomes first)
# mouse_pri_36 <- mouse_pri_36[grepl("^chr", mouse_pri_36$seqnames), ]
# mouse_introns_pri_36 <- extract_introns(mouse_pri_36, verbose = T)
# 

## ----splice-sites, eval=FALSE-------------------------------------------------
# # Human
# library(BSgenome.Hsapiens.UCSC.hg38)
# lnc_47_ss <- assign_splice_sites(
#   introns_lnc_47,
#   genome = BSgenome.Hsapiens.UCSC.hg38,
#   verbose = T
# )
# 
# # Mouse
# library(BSgenome.Mmusculus.UCSC.mm39)
# mouse_pri_36_ss <- assign_splice_sites(
#   mouse_introns_pri_36,
#   genome = BSgenome.Mmusculus.UCSC.mm39,
#   verbose = T
# )

## ----cryptic, eval=FALSE------------------------------------------------------
# # Identify cryptic (non-canonical) splice sites
# cryptic_ss <- find_cryptic_splice_sites(
#   lnc_47_ss,
#   genome = BSgenome.Hsapiens.UCSC.hg38,
#   canonical_donor = "GT",
#   canonical_acceptor = "AG",
#   verbose = TRUE
# )

## ----motifs, eval=FALSE-------------------------------------------------------
# # Donor motifs (5'ss)
# motifs_donor <- extract_ss_motif(
#   input = lnc_47_ss,
#   genome = BSgenome.Hsapiens.UCSC.hg38,
#   type = "5ss",
#   verbose = T,
#   save_fasta = T,
#   output_file = file.path(tempdir(), "lnc_47_5ss_motifs.fa")
# )
# 
# # Acceptor motifs (3'ss)
# motifs_acc <- extract_ss_motif(
#   input = lnc_47_ss,
#   genome = BSgenome.Hsapiens.UCSC.hg38,
#   type = "3ss",
#   verbose = T,
#   save_fasta = T,
#   output_file = file.path(tempdir(), "lnc_47_3ss_motifs.fa")
# )

## ----unspliced, eval=FALSE----------------------------------------------------
# ## identify single exon genes and transcripts
# single_exon_genes <- extract_single_exon(lnc_47_df, level = "gene")
# single_exon_trans <- extract_single_exon(lnc_47_df, level = "transcript")

## ----exon_class, eval=FALSE---------------------------------------------------
# # Assign the ordinal position of exons
# lnc_47_class_exons <- classify_exons(lnc_47_df, verbose = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# # Length of spliced transcript
# lnc_47_spliced_length <- spliced_trans_length(lnc_47_df)
# head(lnc_47_spliced_length)

## ----stat, eval=FALSE---------------------------------------------------------
# # Exon length statistics
# lnc_47_exon_stats <- stat_summary(lnc_47_class_exons, type = "exon")
# 
# # Intron length statistics
# lnc_47_intron_stats <- stat_summary(introns_lnc_47, type = "intron")

## ----gc-content, eval=FALSE---------------------------------------------------
# # Human
# lnc_47_gc <- calculate_gc_content(
#   lnc_47_df,
#   genome = BSgenome.Hsapiens.UCSC.hg38,
#   verbose = TRUE
# )
# # Mouse
# mouse_pri_36_gc <- calculate_gc_content(
#   mouse_pri_36,
#   genome = BSgenome.Mmusculus.UCSC.mm39,
#   verbose = TRUE
# )

## ----cds, eval=FALSE----------------------------------------------------------
# # Convert to GRanges and extract
# library(GenomicRanges)
# mouse_pri_36_granges <- GRanges(mouse_pri_36)
# mouse_cds_seqs <- extract_cds_sequences(
#   mouse_pri_36_granges,
#   BSgenome.Mmusculus.UCSC.mm39,
#   save_fasta = TRUE,
#   output_file = file.path(tempdir(), "mouse_pri_36_CDS.fa.gz")
#   verbose = TRUE
# )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
   devtools::session_info()

