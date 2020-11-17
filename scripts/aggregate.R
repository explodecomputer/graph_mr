# assume that every chunk file has a list of data frames
# the list in each file is called 'l'

library(dplyr)

args <- commandArgs(T)

chunk_filelist <- args[-length(args)]
output <- args[length(args)]

result <- list()
for(i in chunk_filelist)
{
	load(i)
	result[[i]] <- dplyr::bind_rows(l)
}

result <- dplyr::bind_rows(result)

save(result, file=output)
