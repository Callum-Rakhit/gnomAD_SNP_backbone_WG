#/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop(paste("Supply an input vcf (1) and an output file name/location (2)",
             sep=""),
       call.=FALSE)
}

input_vcf <- read.delim(args[1], header = F)

SNP_maximiser <- function(example_vcf, output_name){
  output_name <- data.frame(matrix(ncol = 3, nrow = 0))
  chr <- example_vcf[1,][1]
  start <- as.integer(example_vcf[1,][2])
  end <- as.integer(tail(example_vcf[2], 1))
  bins <- as.integer(floor((end-start)/100))
  tmp_value <- start
  for (i in (1:(bins-1))){
    tmp_value_old <- tmp_value
    tmp_value <- tmp_value + 100
    names(output_name) <- c("Chr", "Start", "End")
    output_name[nrow(output_name) + 1,] <- list(
      as.character(chr[1,]), tmp_value_old + 1, tmp_value)
  }
  output_name[nrow(output_name) + 1,] <- list(
    as.character(chr[1,]), tmp_value_old + 101, end)
  return(output_name)
}

output <- SNP_maximiser(input_vcf, output)

write.table(x = output, file = args[2], sep = "\t", quote = F, 
            row.names = F, col.names = F)
