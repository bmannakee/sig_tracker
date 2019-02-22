#' Load a variant call file.
#' Allowed types are either maflite, vcf (version >=4), call_stats (MuTect1 format),
#' or snp (Varscan2 format) files
#'
#' The function is generally internal, but is exported to help the user in debugging input
#'
#' @param variant_file /path/to/variant/file
#' @param file_type either maf, vcf, call_stats, or snp
#' @return Data frame with 5 columns: chr, start, end, ref_allele, alt_allele
#'
#' @examples
#' # An example of converting oncotator output to maflite and calling load_variants
#' \dontrun{
#'
#'
#' }
load_variants <- function(vcf,reference,sample_name,file_type="vcf"){

  if (! file.exists(vcf)){
    stop('Variant file does not exist - Please check the path')
  }
  if (file_type == "vcf"){
    .load_and_filter_vcf(vcf,reference,sample_name)
  }
  else {
    .get_variant_locs_cs(vcf,reference,sample_name)
  }
}

.load_and_filter_vcf <- function(vcf,reference,sample_name){
  message("Loading VCF")
  # Load vcf as ranges, filter.
  # Need SNVs that pass plus SNVs that only fail tlod
  hard_filters <- S4Vectors::FilterRules(list(snv=VariantAnnotation::isSNV,
                                              sample=function(x) {VariantAnnotation::sampleNames(x)==sample_name}))
  vr <- VariantAnnotation::readVcfAsVRanges(vcf)
  vr <- S4Vectors::subsetByFilter(vr,filter=hard_filters)
  # Need to get rid of the GL sequences, leave just the chromosomes
  GenomeInfoDb::seqlevels(vr) <- GenomicAlignments::seqlevelsInUse(vr)
  GenomeInfoDb::genome(vr) <- GenomeInfoDb::genome(reference)[1:length(GenomeInfoDb::genome(vr))] # These need to have exactly the same name
  # Now get contexts from SomaticSignatures
  message("Loading contexts")
  vr <- SomaticSignatures::mutationContext(vr,reference)
  message("contexts loaded")
  mcols(vr)$TLOD <- as.numeric(mcols(vr)$TLOD)
  # TODO: softfiltermatrix values are NA in many cases, this happens because they put a "." in where all filters pass in GATK4
  VariantAnnotation::softFilterMatrix(vr)[is.na(VariantAnnotation::softFilterMatrix(vr))] <- TRUE # This fixes it
  # GATK 4 changes the filter t_lod_fstar to t_lod. need different functions for mutect2 in gatk3-4
  if (sum(stringr::str_detect(colnames(VariantAnnotation::softFilterMatrix(vr)),'t_lod_fstar'))>0){
    message("Detected a GATK 3 VCF")
    mcols(vr)$tlod_only <- rowSums(VariantAnnotation::softFilterMatrix(vr))==ncol(VariantAnnotation::softFilterMatrix(vr))-1 & !VariantAnnotation::softFilterMatrix(vr)[,"t_lod_fstar"]
  }else{
    message("Detected a GATK 4 VCF")
    mcols(vr)$tlod_only <- rowSums(VariantAnnotation::softFilterMatrix(vr))==ncol(VariantAnnotation::softFilterMatrix(vr))-1 & !VariantAnnotation::softFilterMatrix(vr)[,"t_lod"]
  }
  mcols(vr)$pass_all <- rowSums(VariantAnnotation::softFilterMatrix(vr))==ncol(VariantAnnotation::softFilterMatrix(vr))
  # GATK 4 has lists of probabilities in mcols. For now just get rid of them
  mcols(vr)$SA_MAP_AF <- ""
  mcols(vr)$SA_POST_PROB <- ""
  vars <- tibble::as_tibble(vr)
  vars <- vars %>% dplyr::mutate(cref=stringr::str_sub(alteration,1L,1L),
                                 calt=stringr::str_sub(alteration,2L,2L),
                                 mutect_odds=10**(TLOD - 6.0))
  vars
}