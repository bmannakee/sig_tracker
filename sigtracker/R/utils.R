


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
  out_fr <- df %>% group_by(time,trinucleotide_context) %>%
    summarize(count=n()) %>%
    mutate(proportion=count/sum(count))
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
