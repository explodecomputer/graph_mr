# assume that every chunk file has a list of data frames
# the list in each file is called 'l'

library(dplyr)

args <- commandArgs(T)

chunk_filelist <- args[-length(args)]
output <- args[length(args)]

#args <- commandArgs(T)
#chunk_filelist <- args[1]
#output <- args[2]

print(chunk_filelist)
class(chunk_filelist)
print(output)

result <- list()
for(i in chunk_filelist)
{
	load(i)
	result[[i]] <- dplyr::bind_rows(l)
}

result <- dplyr::bind_rows(result)

save(result, file=output)




#chunk_filelist=c("/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out1.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out2.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out3.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out4.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out5.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out6.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out7.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out8.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out9.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out10.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out11.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out12.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out13.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out14.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out15.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out16.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out17.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out18.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out19.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out20.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out21.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out22.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out23.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out24.rdata","/mnt/storage/home/kf19639/repo/graph_mr/data/sim_data/out25.rdata")
