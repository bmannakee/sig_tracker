#' run sigtracker
#' Stub here
#' vcf is a mutect vcf 4.0 file
#' reference needs to be a FAFile object
#' vr <- smartcallr::compute_odds('./data/non_pdx/vcfs/AZ1013T1.vcf',reference=BSgenome.Hsapiens.1000genomes.hs37d5)

run_tracksig <- function(vcf,reference, bin_size = 100, sample_name='TUMOR', file_type = "vcf"){
  if (!file.exists(vcf)){
    stop('The VCF path is not valid.')
  }

  # load the vcf as vranges
  # Rules applied are filters for SNV and sample name
  vars <- load_variants(vcf,reference,sample_name, file_type)
  with_vars <- draw_vafs(vars)
  with_vars <- with_vars %>% dplyr::arrange(desc(rbeta))
  with_times <- order_time(with_vars, bin_size)
  with_proportions <- get_proportions(with_times)
  with_one_hot <- get_one_hot(with_times)
  return(with_one_hot)
}

