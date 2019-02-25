#' run sigtracker
#' Stub here
#' vcf is a mutect vcf 4.0 file
#' reference needs to be a FAFile object
#' vr <- smartcallr::compute_odds('./data/non_pdx/vcfs/AZ1013T1.vcf',reference=BSgenome.Hsapiens.1000genomes.hs37d5)

run_tracksig <- function(vcf,reference,sample_name='TUMOR', file_type = "vcf"){
  if (!file.exists(vcf)){
    stop('The VCF path is not valid.')
  }

  # load the vcf as vranges
  # Rules applied are filters for SNV and sample name
  vars <- load_variants(vcf,reference,sample_name, file_type)
  vars
}
