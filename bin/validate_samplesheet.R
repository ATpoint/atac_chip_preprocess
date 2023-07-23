args = commandArgs(trailingOnly=TRUE)
ss_file <-args[1]
ss <- read.delim(ss_file, sep=",", header=TRUE)
ss[is.na(ss)] <- ""

l_error <- list()

length_ss <- length(colnames(ss))
sum_ss    <- sum(colnames(ss) %in% c("sample", "r1", "r2"))
if(length_ss != 3 | sum_ss != 3) 
  l_error$headers <- "samplesheet header must be 'sample,r1,r2' -- 3 columns -- csv!"

unique_sample <- unique(ss$sample)
for(u in unique_sample){
  
  sx <- ss[ss$sample==u,]
  
  # R1 and R2 must not have duplicates
  sx_r1 <- table(sx$r1)
  sx_r1w <- which(sx_r1>1)
  if(length(sx_r1w)>0) 
    l_error[paste0("r1_", u)] <- paste("Duplicated r1:", paste(names(sx_r1)[sx_r1w], collapse = ","))
  
  if(unique(sx$r2)[1]!=""){
    sx_r2 <- table(sx$r2)
    sx_r2w <- which(sx_r2>1)
    if(length(sx_r2w)>0) 
      l_error[paste0("r2_", u)] <- paste("Duplicated r2:", paste(names(sx_r2)[sx_r2w], collapse = ","))
  }
  
  # All tech reps per sample must be all paired or all single
  if(!length(sx$r1) == length(sx$r2))
    l_error[paste0("notsame_", u)] <- paste("Sample:", u, "has a mix of SE and PE files")
  
  rm(list=c(ls(pattern="^sx")))
  
}

if(length(l_error)>0){
  for(i in 1:length(l_error)){
    if(i==1) message("[Validation errors for samplesheet]")
    message(l_error[[i]])
  }
} else {
  write.table(ss, "samplesheet_validated.csv", sep=",", col.names=TRUE, 
              row.names=FALSE, quote=FALSE)
}

