library(readr)
library(tibble)
library(dplyr)
library(ggplot2)
#return a matrix of the bp
get_bp_info_as_tibble <- function(data){
  bp_info <- matrix(, nrow = 0, ncol = 4)
  for( i in seq(1,nrow(data))){
    base = data$bp[i]
    string <- strsplit(data$info[i], "")[[1]]
    As <- 0
    Cs <- 0
    Ts <- 0
    Gs <- 0
    for( char in string){
      if(char %in% c('A','a') | data$bp[i] %in% c('A','a')) {
        As <- As + 1
      }else if(char %in% c('C','c') | (data$bp[i] %in% c('C','c')  & char %in% c('.',',','') )){
        Cs <- Cs + 1
      }else if(char %in% c('T','t') | (data$bp[i] %in% c('T','t')  & char %in% c('.',',','') )){
        Ts <- Ts + 1
      }else if(char %in% c('G','g') | (data$bp[i] %in% c('G','g')  & char %in% c('.',',','') )  ){
        Gs <- Gs + 1
      }
    }
    bp_info <- rbind(bp_info, c(As,Cs,Gs,Ts))
  }
  bp_info <- as.tibble(bp_info)
  colnames(bp_info) <- c('A','C','G','T')
  return(bp_info)
}

files <- c('pileups-by-species/pileups-by-species/R.hayd-35.pileup',
           'pileups-by-species/pileups-by-species/R.hayd-78.pileup',
           'pileups-by-species/pileups-by-species/R.mel-00.pileup',
           'pileups-by-species/pileups-by-species/R.mel-01.pileup',
           'pileups-by-species/pileups-by-species/R.mel-02.pileup',
           'pileups-by-species/pileups-by-species/R.mel-10.pileup',
           'pileups-by-species/pileups-by-species/R.mel-68_B.pileup',
           'pileups-by-species/pileups-by-species/R.pari-03.pileup',
           'pileups-by-species/pileups-by-species/R.pari-05.pileup',
           'pileups-by-species/pileups-by-species/R.pari-65.pileup',
           'pileups-by-species/pileups-by-species/R.poly-3O.pileup',
           'pileups-by-species/pileups-by-species/R.poly-3T.pileup',
           'pileups-by-species/pileups-by-species/R.poly-5B.pileup',
           'pileups-by-species/pileups-by-species/R.poly-8G.pileup',
           'pileups-by-species/pileups-by-species/R.port-3A.pileup',
           'pileups-by-species/pileups-by-species/R.port-8A.pileup',
           'pileups-by-species/pileups-by-species/R.port-95.pileup',
           'pileups-by-species/pileups-by-species/R.port-96.pileup',
           'pileups-by-species/pileups-by-species/R.port-96TA.pileup',
           'pileups-by-species/pileups-by-species/R.port-97.pileup',
           'pileups-by-species/pileups-by-species/R.shus-2.pileup',
           'pileups-by-species/pileups-by-species/R.shus-3.pileup',
           'pileups-by-species/pileups-by-species/R.shus-4.pileup',
           'pileups-by-species/pileups-by-species/R.shus-5.pileup',
           'pileups-by-species/pileups-by-species/R.shus-6.pileup')
names <- c('Rhha a','Rhha b',
           'Rhme a','Rhme b','Rhme c','Rhme d','Rhme e',
           'Rhpa a','Rhpa b','Rhpa c',
           'Rhpol a','Rhpol b','Rhpol c','Rhpol d',
           'Rhpor a','Rhpor b','Rhpor c','Rhpor d','Rhpor e','Rhpor f',
           'Rhsu a','Rhsu b','Rhsu c','Rhsu d','Rhsu e')

its1_trim <- 600
its1_span <- 195
s58_trim <- its1_trim + its1_span
s58_span <- 161
its2_trim <- s58_trim + s58_span
its2_span <- 157

data <- read_tsv(files[1],col_names = FALSE)
if(length(colnames(data)) != 6){
  colnames(data) <- c('name','position','bp','depth','info')
}else{
  colnames(data) <- c('name','position','bp','depth','info','quality') 
}



#front_trim <- 600
#span <- 100
#instead of trimming, use the second field as a length parameter
#front trim should be 600bp to account for the buffers attached or 600 + lenght of preluding sequences. 
#Span is the lenght of sequence of interest, ITS1, 5.8S of ITS2
front_trim <- its1_trim
span <- its1_span
data_trimmed <- data[front_trim:(front_trim + span),]

data_bp_info <- get_bp_info_as_tibble(data_trimmed)

variances <- c()
for(i in seq(1,nrow(data_bp_info))){
  variance <- 1 - (max(as.numeric(data_bp_info[i,])) / sum(as.numeric(data_bp_info[i,]))) 
  variances <- c(variances,variance)
}
data_bp_info$variance <- variances
data_bp_info$count <- seq(1,nrow(data_bp_info))

p <- ggplot(data=data_bp_info,aes(x=count,y=variance)) + geom_line()
p

As <- apply(data_bp_info,1,function(row){
  row <- unlist(row)
  return(c('A',as.character(row[1]),as.numeric(unlist(row[5]))))
})
As <- as.tibble(t(as.tibble(As)))
Cs <- apply(data_bp_info,1,function(row){
  row <- unlist(row)
  return(c('C',as.character(row[2]),as.numeric(unlist(row[5]))))
})
Cs <- as.tibble(t(as.tibble(Cs)))
Gs <- apply(data_bp_info,1,function(row){
  row <- unlist(row)
  return(c('G',as.character(row[3]),as.numeric(unlist(row[5]))))
})
Gs <- as.tibble(t(as.tibble(Gs)))
Ts <- apply(data_bp_info,1,function(row){
  row <- unlist(row)
  return(c('T',as.character(row[4]),as.numeric(unlist(row[5]))))
})
Ts <- as.tibble(t(as.tibble(Ts)))

tiddy_bp_info <- bind_rows(As,Ts,Cs,Gs)
colnames(tiddy_bp_info) <- c('nucleotide','count','order')

tiddy_bp_info$count <- as.numeric(tiddy_bp_info$count)
tiddy_bp_info$order <- as.numeric(tiddy_bp_info$order)

p <-  ggplot(data=tiddy_bp_info,aes(x=order,y=count,color=nucleotide)) + geom_line()
p
ggsave('delete-this.png',width = 30, height = 10, units = c("in"))



