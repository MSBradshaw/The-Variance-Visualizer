library(readr)
library(ggplot2)
library(tibble)
library(dplyr)
library(Biostrings)
library(stringr)
library(plyr)

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

get_majority_base <- function(data){
  bases <- c('a','c','g','t')
  apply(data,1,function(row){
    bases[which.max(row[1:4])]
  })
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
           'PileUps-April-2019/8678_ATTACTCG-TATAGCCT_L002_.pileup',
           'PileUps-April-2019/occulta_8807i_.pileup',
           'PileUps-April-2019/occulta_8807iv_.pileup',
           'PileUps-April-2019/8795_TAATGCGC-GGCTCTGA_L002_.pileup',
           'PileUps-April-2019/parilis_8805_.pileup',
           'PileUps-April-2019/8796_TAATGCGC-AGGCGAAG_L002_.pileup',
           'PileUps-April-2019/polymorpha_8665u_.pileup',
           'PileUps-April-2019/8797_CGGCTATG-TAATCTTA_L002_.pileup',
           'PileUps-April-2019/polymorpha_8807iii_.pileup',
           'PileUps-April-2019/8800_ATTACTCG-CCTATCCT_L002_.pileup')

names <- c('Rhpor A',
           'Rhme A',
           'Rhpol A',
           'Rhme B',
           'Rhpol B',
           'Rhpa A',
           'Rhsu A',
           'Rhsu B',
           'Rhpa B',
           'Rhsu C',
           'Rhsu D',
           'Rhme C',
           'Rhsu E',
           'hayd A',
           'Rhpa C',
           'Rhpol C',
           'Rhpor B',
           'hayd B',
           'Rhme D',
           'Rhme REF',
           'Rhpol D',
           'hayd C',
           'Rhoc A',
           'Rhoc B',
           'Rhpor C',
           'Rhpa D',
           'Rhpor D',
           'Rhpol E',
           'Rhpor E',
           'Rhpol F',
           'Rhme F')

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
  data$base <- get_majority_base(data)
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
    aln <- seq#subject(pairwiseAlignment(ref,seq,gapExtension=4,gapOpening=10))
    count  = 1
    spaces_count <- 0
    for( letter in strsplit(toString(aln),'')[[1]]){
      letter <- tolower(letter)
      if(letter == 'a'){
        As[count] <- As[count] + 1
      }else if (letter == 'c'){
        Cs[count] <- Cs[count] + 1
      }else if (letter == 'g'){
        Gs[count] <- Gs[count] + 1
      }else if (letter == 't'){
        Ts[count] <- Ts[count] + 1
      }else if (letter == '-' || letter == 'n'){
        spaces[count] <- spaces[count] + 1
        #get the index of the current count in the temp tibble
        num <- match((count-spaces_count),as.numeric(temp$count))
        #add an indel to the bp information
        temp = add_row(temp,A=0,C=0,G=0,T=0,variance=0,sample=temp$sample[1],base='-',count=-1,.after=(num-1))
        spaces_count <- spaces_count  + 1
      }else{
        print(letter)
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
      letter <- tolower(letter)
      if(letter == 'a'){
        As[count] <- As[count] + 1
      }else if (letter == 'c'){
        Cs[count] <- Cs[count] + 1
      }else if (letter == 'g'){
        Gs[count] <- Gs[count] + 1
      }else if (letter == 't'){
        Ts[count] <- Ts[count] + 1
      }else if (letter == '-' || letter == 'n'){
        spaces[count] <- spaces[count] + 1
        #get the index of the current count in the temp tibble
        num <- match((count-spaces_count),as.numeric(temp$count))
        #add an indel to the bp information
        temp = add_row(temp,A=0,C=0,G=0,T=0,variance=0,sample=temp$sample[1],base='-',count=-1,.after=(num-1))
        spaces_count <- spaces_count  + 1
      }else{
        print(letter)
        print('Dammit2')
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

add_variance_points <- function(data){
  data <- bind_rows(data,data)
  data$sample[(nrow(data)/2):nrow(data)] <- lapply(data[(nrow(data)/2):nrow(data),]$sample,function(x){
    return(paste(x,'*',sep=''))
  })
  data$sample <- factor(as.character(unlist(data$sample)))
  #data$sample <- revalue(data$sample, c("Rhme A*"=" ",   "Rhpor A*"="  "))
  data$hybrid <- data$consensus_match
  data$variance[is.na(data$variance)] <- 0
  data$hybrid[(nrow(data)/2):nrow(data)] <- apply(data[(nrow(data)/2):nrow(data),],1,function(x){
    #5 is the variance
    var <- x[5]
    if(var == 0 || var  < .001){
      return('0')
    }else if(var > .1){
      return('+10%')
    }else if(var > .01){
      return('1-9%')
    }else if(var > 0){
      return('< 1%')
    }else{
      print('Error')
      return('0')
    }
  })
  return(data)
}

#add a column to the bp_info of where the base matches the consensus sequence or not
compare_seq_to_con <- function(data,con){
  data[is.na(data$base),7] <- '-'
  con_vec <- strsplit(toString(con),'')[[1]]
  i = 1
  data$consensus_match <- apply(data,1,function(row){
    base <- row[7]
    pos <- as.numeric(str_trim(row[8]))
    #print(i)
#    print(pos)
#    print(length(con_vec))
    if(i == 10532 || i == 10343){
      one =1
      #stop here
    }
    if(pos > length(con_vec)){
      i <<- i +1
      return(base)
    }else if(con_vec[[pos]] != base){
      i <<- i +1
      return(base)
    }else{
      i <<- i +1
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

its1 <- get_sections(files,names,its1_trim,513)

#plot the deviations from the consensus sequence 
ref_its1_8678 <- get_seq(its1,'Rhpor A')
its1_info <- get_consensus(its1,ref_its1_8678)
its1_final_info <- compare_seq_to_con(its1_info[[2]],ref_its1_8678)

its1_con <- get_simple_consensus(its1,ref_its1_8678)
its1_info <- get_consensus(its1_con[[2]],its1_con[[1]])
its1_final_info <- compare_seq_to_con(its1_info[[2]],its1_con[[1]])

#plot the variation

samps <- unique(its1_final_info$sample)

#change the order of the consensus_match column's factors so they appear correctly in the legend
its1_final_info$consensus_match <- factor(its1_final_info$consensus_match,levels=c('match','a','c','g','t','0','+10%','1-9%','< 1%'))
its1_final_info <- add_variance_points(its1_final_info)

p <- ggplot(data = its1_final_info,aes(x=count,y=sample,color=hybrid),fill=hybrid) + 
  geom_point(shape=15) + 
  theme_minimal() + 
  scale_color_manual(values=c( "#e5e5e5","#ff0000", "#4d4dff", "#0000ff","#b3b300",'#000000','#ff80ff','#ff00ff','#800080')) +
  xlab('Nucleotide Position') + 
  ylab('Sample') + 
  theme(legend.title=element_blank()) + 
  geom_vline(xintercept=its1_span, linetype="solid", color = "red") +
  geom_vline(xintercept=(its1_span+s58_span), linetype="solid", color = "red") 
ggsave('ITS-Hybrid.png',width = 12, height = 5, units = c("in"))


#sample, number, base, variance, A, T, C, G