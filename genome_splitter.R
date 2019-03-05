#/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop(paste("Supply an input vcf (1) and an output file name/location (2)",
             sep=""),
       call.=FALSE)
}

# Script to pick the 100bp region with the most SNPs for
# each of the 6206 0.5Mb chunks contained in the genome

hg19 <- read.delim(args[1], header = F)

options(scipen = 999) # Remove scientific notation

binmaker <- function(input, output_name){
  output_name <- data.frame(matrix(ncol = 3, nrow = 0))
  for(i in 1:dim(input)[1]){
    chr <- hg19[i,][1]
    start <- as.integer(hg19[i,][2])
    end <- as.integer(hg19[i,][3])
    bins <- as.integer(floor(end/500000))
    tmp_value <- 0
    for (i in ((1:bins)-1)){
      tmp_value_old <- tmp_value
      tmp_value <- tmp_value + 500000
      names(output_name) <- c("Chr", "Start", "End")
      output_name[nrow(output_name) + 1,] <- list(as.character(chr[1,]), tmp_value_old+1, tmp_value)
    }
    output_name[nrow(output_name) + 1,] <- list(as.character(chr[1,]), tmp_value_old+500001, end)
  }
  return(output_name)
}

df <- binmaker(hg19, output_name)
df <- apply(df, 2, as.character) # Flattens the data frame so you can write to file

write.table(x = df, file = args[2], 
            sep = "\t", quote = F, row.names = F, col.names = F)
