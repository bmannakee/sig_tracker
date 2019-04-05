


# function to draw beta VAFs
draw_vafs <- function(df){
  df$rbeta <- rep(0,nrow(df))
  for(i in 1:nrow(df)){
    beta_i <-rbeta(1, df$altCount[i],df$refCount[i])
    df$rbeta[i]<-beta_i
  }
  return(df)
}

order_time <- function(df, size = 100){
  max_time <- ceiling(nrow(df)/size)
  df$time <- NA
  for(i in 1:nrow(df)){
    df$time[i]<-ceiling(i/size)
  }
  df
}

get_proportions <- function(df){
  out_fr <- df %>% dplyr::group_by(time,trinucleotide_context) %>%
    dplyr::summarize(count=n()) %>%
    dplyr::mutate(proportion=count/sum(count))
  out_fr
}

get_one_hot <- function(df){
 #pull out the first three rows
  df_counts<-unique(df[,1:3])
  #make the 98 columns and fill with 0
  df_counts[,4:(length(unique(df$trinucleotide_context))+3)]<-0
  #name each column
  colnames(df_counts)[4:ncol(df_counts)]<-unique(df$trinucleotide_context)
  for (i in 1:nrow(df)){
    trinuc<- df$trinucleotide_context[i]
    df_counts[i,trinuc]<-df_counts[i,trinuc]+1
  }

  one_hot_encoding<- merge(df,df_counts)
  return(one_hot_encoding)
}

em_alg <- function(our_fr){
  sig_fr <- get_signature_fr('~/Desktop/projects/sig_tracker/sigtracker/cosmic_signatures.txt')
  eps <- .00001 # converges at abs(sum(pi_current - pi_new)) <= eps
  pi_current <- rep(1/30,30) # vector of mixing coefficients
  names(pi_current) <- paste0("signature_",1:30)
  sigs_split <- sig_fr %>% split(.$trinucleotide_context) %>%
    map(~ select(.,contains("signature_")))
  this_fr <- our_fr %>% dplyr::mutate(this_sig = map(trinucleotide_context, ~ as.double(sigs_split[[.x]])))
  this_fr <- this_fr %>% dplyr::mutate(pi_cur = list(pi_current))
  this_fr <- this_fr %>% dplyr::mutate(zi = purrr::map2(pi_cur,this_sig, ~ unlist(.x)*unlist(.y))) # zi carries an expectation step
  pi_fr <- as_tibble(Reduce(rbind,this_fr$zi))
  pi_new <- pi_fr %>% summarise_all(mean)
  while(abs(sum(pi_current - pi_new)) > eps){
    pi_current <- pi_new
    this_fr <- this_fr %>% dplyr::mutate(pi_cur = list(pi_current))
    this_fr <- this_fr %>% dplyr::mutate(zi = purrr::map2(pi_cur,this_sig, ~ unlist(.x)*unlist(.y))) # zi carries an expectation step
    pi_fr <- as_tibble(Reduce(rbind,this_fr$zi))
    pi_fr <- pi_fr %>% summarise_all(mean) # maximization
    pi_new <- pi_fr[1,] %>% as.double
  }
  # Get the top five weighted signatures and redo
  top_5_index <- order(-pi_new)[1:5]
  pi_current <- rep(1/5,5)

  names(pi_current) <- paste0("signature_",top_5_index)

  sigs_split <- sig_fr %>% split(.$trinucleotide_context) %>%
    map(~ select(.,names(pi_current)))

  this_fr <- our_fr %>% dplyr::mutate(this_sig = map(trinucleotide_context, ~ as.double(sigs_split[[.x]])))

  this_fr <- this_fr %>% dplyr::mutate(pi_cur = list(pi_current))

  this_fr <- this_fr %>% dplyr::mutate(zi = purrr::map2(pi_cur,this_sig, ~ unlist(.x)*unlist(.y))) # zi carries an expectation step

  pi_fr <- as_tibble(Reduce(rbind,this_fr$zi))

  pi_new <- pi_fr %>% summarise_all(mean)

  while(abs(sum(pi_current - pi_new)) > eps){
    pi_current <- pi_new
    this_fr <- this_fr %>% dplyr::mutate(pi_cur = list(pi_current))
    this_fr <- this_fr %>% dplyr::mutate(zi = purrr::map2(pi_cur,this_sig, ~ unlist(.x)*unlist(.y))) # zi carries an expectation step
    pi_fr <- as_tibble(Reduce(rbind,this_fr$zi))
    pi_fr <- pi_fr %>% summarise_all(mean) # maximization
    pi_new <- pi_fr[1,] %>% as.double
  }

  list(pi=pi_new/sum(pi_new), sigs=paste0("signature_",top_5_index)) # return the weights for the restricted sample of signatures
  #pi_new
}

get_signature_fr <- function(sig_file){
  fr <- suppressMessages(read_tsv(sig_file) %>% dplyr::select(-contains("X")))
  names(fr) %<>% stringr::str_replace_all("\\s","_") %>% tolower
  # get the contexts in the correct form
  fr <- fr %>% dplyr::mutate(alteration = str_replace(substitution_type,">",""),
                             context = paste0(str_sub(trinucleotide,1L,1L),".",str_sub(trinucleotide,3L,3L)),
                             trinucleotide_context = paste0(alteration,"_",context))
  fr <- fr %>% dplyr::select(trinucleotide_context,contains("signature_"))
  fr
}
