library(readr)
library(ggplot2)
library(tibble)
library(dplyr)
library(Biostrings)
library(stringr)

get_bp_info_as_tibble <- function(data){
  bp_info <- matrix(, nrow = 0, ncol = 4)
  for( i in seq(1,nrow(data))){
    base = data$base[i]
    string <- strsplit(data$reads[i], "")[[1]]
    As <- 0
    Cs <- 0
    Ts <- 0
    Gs <- 0
    for( char in string){
      if(char %in% c('A','a') | data$base[i] %in% c('A','a')) {
        As <- As + 1
      }else if(char %in% c('C','c') | (data$base[i] %in% c('C','c')  & char %in% c('.',',','') )){
        Cs <- Cs + 1
      }else if(char %in% c('T','t') | (data$base[i] %in% c('T','t')  & char %in% c('.',',','') )){
        Ts <- Ts + 1
      }else if(char %in% c('G','g') | (data$base[i] %in% c('G','g')  & char %in% c('.',',','') )  ){
        Gs <- Gs + 1
      }
    }
    bp_info <- rbind(bp_info, c(As,Cs,Gs,Ts))
  }
  bp_info <- as.tibble(bp_info)
  colnames(bp_info) <- c('A','C','G','T')
  return(bp_info)
}

get_info <- function(data){
  bp_info <- get_bp_info_as_tibble(data)
  var <- apply(bp_info,1,function(row){
    variance <- 1 - (max(as.numeric(row)) / sum(as.numeric(row)))
    return(variance)
  })
  bp_info$variance <- var
  bp_info$sample <- data$sample
  bp_info$base <- data$base
  bp_info$count <- seq(1,nrow(bp_info))
  return(bp_info)
}

files <- c('PileUps-April-2019/8663_A_CGGCTATG-CAGGACGT_L002_.pileup',
           'PileUps-April-2019/8801_TCCGGAGA-GGCTCTGA_L002_.pileup',
           'PileUps-April-2019/8663_O_CTGAAGCT-GTACTGAC_L002_.pileup',
           'PileUps-April-2019/8802_TCCGGAGA-AGGCGAAG_L002_.pileup',
           'PileUps-April-2019/8663_T_CTGAAGCT-TATAGCCT_L002_.pileup',
           'PileUps-April-2019/8803_GAGATTCC-CCTATCCT_L002_.pileup',
           'PileUps-April-2019/8664_2_TCCGCGAA-TATAGCCT_L002_.pileup',
           'PileUps-April-2019/8664_3_TCCGCGAA-ATAGAGGC_L002_.pileup',
           'PileUps-April-2019/8806_ATTCAGAA-AGGCGAAG_L002_.pileup',
           'PileUps-April-2019/8664_4_TCCGCGAA-CCTATCCT_L002_.pileup',
           'PileUps-April-2019/8664_5_TCTCGCGC-GGCTCTGA_L002_.pileup',
           'PileUps-April-2019/8810_TCCGGAGA-TAATCTTA_L002_.pileup',
           'PileUps-April-2019/8664_6_TCTCGCGC-AGGCGAAG_L002_.pileup',
           'PileUps-April-2019/8865_U_ATTCAGAA-CAGGACGT_L002_.pileup',
           'PileUps-April-2019/8665_N_ATTCAGAA-TAATCTTA_L002_.pileup',
           'PileUps-April-2019/8935_B_TAATGCGC-CCTATCCT_L001_.pileup',
           'PileUps-April-2019/8668_A_CGGCTATG-GTACTGAC_L002_.pileup',
           'PileUps-April-2019/8935_P_ATTACTCG-ATAGAGGC_L002_.pileup',
           'PileUps-April-2019/8668_B_CGCTCATT-CAGGACGT_L002_.pileup',
           'PileUps-April-2019/melanophthalma_REF_.pileup',
           'PileUps-April-2019/8668_G_CTGAAGCT-ATAGAGGC_L002_.pileup',
           'PileUps-April-2019/novomexicana_8684d_.pileup',
           'PileUps-April-2019/8678_ATTACTCG-TATAGCCT_L002_.pileup',
           'PileUps-April-2019/occulta_8807i_.pileup',
           'PileUps-April-2019/8684dLnovomexicana_S1_L001_001_.pileup',
           'PileUps-April-2019/occulta_8807iv_.pileup',
           'PileUps-April-2019/8795_TAATGCGC-GGCTCTGA_L002_.pileup',
           'PileUps-April-2019/parilis_8805_.pileup',
           'PileUps-April-2019/8796_TAATGCGC-AGGCGAAG_L002_.pileup',
           'PileUps-April-2019/polymorpha_8665u_.pileup',
           'PileUps-April-2019/8797_CGGCTATG-TAATCTTA_L002_.pileup',
           'PileUps-April-2019/polymorpha_8807iii_.pileup',
           'PileUps-April-2019/8800_ATTACTCG-CCTATCCT_L002_.pileup')

names <- c('port_8663_A',
           'mela_8801_T',
           'poly_8663_O',
           'mela_8802_T',
           'poly_8663_T',
           'pari_8803_G',
           'shus_8664_2',
           'shus_8664_3',
           'pari_8806_A',
           'shus_8664_4',
           'shus_8664_5',
           'mela_8810_T',
           'shus_8664_6',
           'hayd_8865_U',
           'pari_8665_N',
           'poly_8935_B',
           'port_8668_A',
           'hayd_8935_P',
           'mela_8668_B',
           'mela_REF',
           'poly_8668_G',
           'novo_8684',
           'hayd_8678_A',
           'occulta_8807i',
           'novo_8684dL_S1',
           'occulta_8807iv',
           'port_8795_T',
           'parilis_8805',
           'port_8796_T',
           'poly_8665u',
           'port_8797_C',
           'poly_8807iii',
           'mela_8800_A')

get_sections <- function(files,names,start,span){
  print(files[[1]])
  data <- read_tsv(files[[1]],col_names = FALSE)
  data <- data[start:(start+span),]
  colnames(data) <- c('sample','position','base','count','reads','quality')
  data <- get_info(data)
  data$sample <- names[[1]]
  for(i in seq(2,length(files))){
    temp <- read_tsv(files[[i]],col_names = FALSE)
    if(all(dim(temp) == c(0,0))){
      print(files[[i]])
      print('FAIL!!!')
    }else{
      colnames(temp) <- c('sample','position','base','count','reads','quality')
      temp <- get_info(temp[start:(start+span),])
      temp$sample <- names[[i]]
      data <- bind_rows(data,temp) 
    }
  }
  return(data)
}

get_seq <- function(data,name){
  paste0(data[data$sample==name,]$base,collapse = '')
}

get_consensus <- function(x,ref){
  #assuming that the lenght of the longest string will not be more than double that of the reference sequence
  As <- numeric((nchar(ref)*2))
  Cs <- numeric((nchar(ref)*2))
  Gs <- numeric((nchar(ref)*2))
  Ts <- numeric((nchar(ref)*2))
  spaces <- numeric((nchar(ref)*2))
  counts <- c()
  bp_infos = 'Empty'
  for( samp in unique(its1$sample)){
    temp <- x[x$sample == samp,]
    seq <- get_seq(temp,samp)
    aln <- subject(pairwiseAlignment(ref,seq,gapExtension=2,gapOpening=4))
    count  = 1
    spaces_count <- 0
    for( letter in strsplit(toString(aln),'')[[1]]){
      if(letter == 'a'){
        As[count] <- As[count] + 1
      }else if (letter == 'c'){
        Cs[count] <- Cs[count] + 1
      }else if (letter == 'g'){
        Gs[count] <- Gs[count] + 1
      }else if (letter == 't'){
        Ts[count] <- Ts[count] + 1
      }else if (letter == '-'){
        spaces[count] <- spaces[count] + 1
        #get the index of the current count in the temp tibble
        num <- match((count-spaces_count),as.numeric(temp$count))
        #add an indel to the bp information
        temp = add_row(temp,A=0,C=0,G=0,T=0,variance=0,sample=temp$sample[1],base='-',count=-1,.after=(num-1))
        spaces_count <- spaces_count  + 1
      }else{
        print('Dammit')
      }
      count <- count + 1
    }
    counts <- c(counts,count)
    temp$count <- seq(1,nrow(temp))
    #add the current temp to the bp_infos, this is so that we have bo infos that contain indels
    if(bp_infos == 'Empty'){
      bp_infos <- temp
    }else{
      bp_infos <- bind_rows(bp_infos,temp)
    }
  }
  bases_counts <- tibble(As,Cs,Ts,Gs,spaces)
  bases <- c('a','c','t','g','-')
  str <- apply(bases_counts,1,function(x){
    bases[which.max(x)]
  })
  #cut it off so that the string is only as long as the longest observed sequence
  str <- str[1:max(counts)]
  #this is the consensus sequence
  consensus <- paste(str,collapse = '')
  return(list(consensus,bp_infos))
}

get_simple_consensus <- function(x,ref){
  #assuming that the lenght of the longest string will not be more than double that of the reference sequence
  As <- numeric((nchar(ref)*2))
  Cs <- numeric((nchar(ref)*2))
  Gs <- numeric((nchar(ref)*2))
  Ts <- numeric((nchar(ref)*2))
  spaces <- numeric((nchar(ref)*2))
  counts <- c()
  bp_infos = 'Empty'
  for( samp in unique(x$sample)){
    temp <- x[x$sample == samp,]
    seq <- get_seq(temp,samp)
    count  = 1
    spaces_count <- 0
    for( letter in strsplit(seq,'')[[1]]){
      if(letter == 'a'){
        As[count] <- As[count] + 1
      }else if (letter == 'c'){
        Cs[count] <- Cs[count] + 1
      }else if (letter == 'g'){
        Gs[count] <- Gs[count] + 1
      }else if (letter == 't'){
        Ts[count] <- Ts[count] + 1
      }else if (letter == '-'){
        spaces[count] <- spaces[count] + 1
        #get the index of the current count in the temp tibble
        num <- match((count-spaces_count),as.numeric(temp$count))
        #add an indel to the bp information
        temp = add_row(temp,A=0,C=0,G=0,T=0,variance=0,sample=temp$sample[1],base='-',count=-1,.after=(num-1))
        spaces_count <- spaces_count  + 1
      }else{
        print('Dammit')
      }
      count <- count + 1
    }
    counts <- c(counts,count)
    temp$count <- seq(1,nrow(temp))
    #add the current temp to the bp_infos, this is so that we have bo infos that contain indels
    if(bp_infos == 'Empty'){
      bp_infos <- temp
    }else{
      bp_infos <- bind_rows(bp_infos,temp)
    }
  }
  bases_counts <- tibble(As,Cs,Ts,Gs,spaces)
  bases <- c('a','c','t','g','-')
  str <- apply(bases_counts,1,function(x){
    bases[which.max(x)]
  })
  #cut it off so that the string is only as long as the longest observed sequence
  str <- str[1:max(counts)]
  #this is the consensus sequence
  consensus <- paste(str,collapse = '')
  return(list(consensus,bp_infos))
}



#add a column to the bp_info of where the base matches the consensus sequence or not
compare_seq_to_con <- function(data,con){
  con_vec <- strsplit(toString(con),'')[[1]]
  data$consensus_match <- apply(data,1,function(row){
    base <- row[7]
    pos <- as.numeric(str_trim(row[8]))
    #print(pos)
    #print(length(con_vec))
    if(pos > length(con_vec)){
      return(base)
    }else if(con_vec[[pos]] != base){
      return(base)
    }else{
      return('match')
    }
  })
  return(data)
}

its1_trim <- 600
its1_span <- 195
s58_trim <- its1_trim + its1_span
s58_span <- 161
its2_trim <- s58_trim + s58_span
its2_span <- 157

its1 <- get_sections(files,names,its1_trim,its1_span)
s58 <- get_sections(files,names,s58_trim,s58_span)
its2 <- get_sections(files,names,its2_trim,its2_span)


#plot the deviations from the consensus sequence 
ref_its1_8678 <- get_seq(its1,'hayd_8678_A')
its1_info <- get_consensus(its1,ref_its1_8678)
its1_final_info <- compare_seq_to_con(its1_info[[2]],ref_its1_8678)

its1_con <- get_simple_consensus(its1,ref_its1_8678)
its1_info <- get_consensus(its1_con[[2]],its1_con[[1]])
its1_final_info <- compare_seq_to_con(its1_info[[2]],its1_con[[1]])

#plot the variation

samps <- unique(its1_final_info$sample)

small <- its1_final_info[its1_final_info$sample %in% c("hayd_8678_A",samps[13]),]
p <- ggplot(data = small,aes(x=count,y=sample,color=base),fill=consensus_match) + 
  geom_point(shape=15) + 
  theme_minimal() + 
  scale_color_manual(values=c("#000000", "#ff0000", "#0000ff", "#009900","#e5e5e5","#ffff00"))
ggsave('ITS1-sample.png',width = 12, height = 5, units = c("in"))


#sample, number, base, variance, A, T, C, G