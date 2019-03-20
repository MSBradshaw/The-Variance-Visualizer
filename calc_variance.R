library(readr)
library(tibble)
library(dplyr)

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

get_variance <- function(filep_path,front_trim,span){
  data <- read_tsv(filep_path,col_names = FALSE)
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
  data_trimmed <- data[front_trim:(front_trim + span),]
  
  data_bp_info <- get_bp_info_as_tibble(data_trimmed)
  variances <- c()
  for(i in seq(1,nrow(data_bp_info))){
    variance <- 1 - (max(as.numeric(data_bp_info[i,])) / sum(as.numeric(data_bp_info[i,]))) 
    variances <- c(variances,variance)
  }
  return(variances)
}

get_variance_table <- function(path,sample_name,sample_region,front_trim,span){
  vars <- get_variance(path,front_trim,span)
  #sample_name <- 'R.hayd-35'
  #sample_region <- 'ITS1'
  var_table <- tibble(sample=character(),region=character(),bp=character(),percent_variance=numeric())
  for( i in seq(1,length(vars))){
    if( vars[i] > 0 ){
      row <- c(sample_name,sample_region,i,vars[i])
      var_table[(nrow(var_table)+1),] <- row
    }
  }
  return(var_table)
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

regions <- replicate(length(files), "ITS1")

its1_trim <- 600
its1_span <- 195
s58_trim <- its1_trim + its1_span
s58_span <- 161
its2_trim <- s58_trim + s58_span
its2_span <- 157

#variance information for the ITS1 region
output_its1 <- get_variance_table(files[1],names[1],'ITS1',its1_trim,its1_span)
for( i in seq(2,length(files))){
  print(i)
  temp <- get_variance_table(files[i],names[i],'ITS1',its1_trim,its1_span)
  if(nrow(temp) != 0){
    output_its1 <- bind_rows(temp,output_its1) 
  }
}
print(output_its1)
write_csv(output_its1,'percent_variance--IST1.csv')

#variance information for 5.8S region
output_58s <- get_variance_table(files[1],names[1],'ITS1',s58_trim,s58_span)
for( i in seq(2,length(files))){
  print(i)
  temp <- get_variance_table(files[i],names[i],'5.8s',s58_trim,s58_span)
  if(nrow(temp) != 0){
    output_58s <- bind_rows(temp,output_58s) 
  }
}
print(output_58s)
write_csv(output_58s,'percent_variance--5.8s.csv')

#variance information for the ITS2 region
output_its2 <- get_variance_table(files[1],names[1],'ITS1',its2_trim,its2_span)
for( i in seq(2,length(files))){
  print(i)
  temp <- get_variance_table(files[i],names[i],'ITS2',its2_trim,its2_span)
  if(nrow(temp) != 0){
    output_its2 <- bind_rows(temp,output_its2) 
  }
}
print(output_its2)
write_csv(output_its2,'percent_variance--ITS2.csv')

total_output <- bind_rows(output_its1,output_58s,output_its2)
print(total_output)
write_csv(total_output,'percent_variance--all.csv')

filtered_total <- total_output[total_output$percent_variance > .01,]
write_csv(filtered_total,'percent_variance--all-filtered.csv')
