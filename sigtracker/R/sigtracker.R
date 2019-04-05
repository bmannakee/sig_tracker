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
  print("got_vars")
  with_vars <- draw_vafs(vars)
  with_vars <- with_vars %>% dplyr::arrange(desc(rbeta))
  with_times <- order_time(with_vars, bin_size)
  with_proportions <- get_proportions(with_times)
  with_one_hot <- get_one_hot(with_times)
  sig_props <- with_times %>% dplyr::group_by(time) %>% tidyr::nest() %>% dplyr::mutate(pi = purrr::map(data, ~ em_alg(.x)))
  # returns a list with the original data, the time points, and top five estimates with their signature names
  pis <- do.call(rbind,sig_props$pi)
  final <-  pis %>% as_tibble() %>% dplyr::mutate(time=sig_props$time) %>%
    unnest() %>%
    pivot_wider(key = time,
                names_from = sigs,
                values_from = pi,
                values_fill = list(pi=0))
  # We want signatures presented as ordered from sig 1 to sig 30.
  to_sort <- colnames(final %>% dplyr::select(contains("_"))) %>% stringr::str_split("_",simplify=T)
  new_order <- paste0("signature_",sort(as.integer(to_sort[,2])))
  final %>% dplyr::select(time,new_order)
}

