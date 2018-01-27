library(EPOM)
library(rtracklayer)
library(parallel)

library(readr)
bw_list <- read_csv("bw_list.txt", col_names = FALSE)
bed_list <- read_csv("bed_list.txt", col_names = FALSE)
cell_type <- read_csv("group_levels.txt", col_names = FALSE)

bw_list <- as.character(unlist(bw_list))
bed_list <- as.character(unlist(bed_list))
cell_type <- as.character(unlist(cell_type))

histone_mark = "h3k4me1"
state = 10:18

bed_header = FALSE
bed_sep = "\t"
save = TRUE

score <- epom(bw_list, bed_list, cell_type = cell_type, histone_mark = histone_mark, state = state, state_name = "enhancer", m = 12, cores = 8)
print(score)
