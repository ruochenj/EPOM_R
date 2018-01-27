#' Select the corresponding regions of a specific chromatin state.
#' @description Select the corresponding regions of a specific chromatin state.
#'
#' @param gen_mat A matrix generated from matrix_bp_trans.
#' @param bed A character specifying the path of the bed file that contains the bed file specify the corresponding regions for each state.
#' @param state A character specifying the chromatin state of interest.
#' @return A vector indicating the rows inside the input gen_mat matrix that correspond to state-region information of input bed file.
#' @export
index_num4state <- function(gen_mat, bed, state){
  gen_mat_2 <- cbind(rep(0,dim(gen_mat)[1]), gen_mat)
  seqname <- levels(bed[,1])
  if(sum(!state %in% levels(bed[,4]))>0){
    warning("Some of the states input is not included in the bed file")
  }
  bed_select <- bed[bed[,4] %in% state,]
  bed_select[,1] = as.character(bed_select[,1])
  bed_select[,2] = as.numeric(as.character(bed_select[,2]))
  bed_select[,3] = as.numeric(as.character(bed_select[,3]))
  for(seq_i in seqname){
    if(seq_i != "track"){
      mat_trans <- gen_mat_2[gen_mat_2[,2] == seq_i, ]
      print(seq_i)
      bed_chi = bed_select[bed_select[,1] == seq_i,]
      begin_select = bed_chi[,2]/200
      end_select = ceiling(bed_chi[,3]/200)
      select_rows = c()
      k = length(begin_select)
      #this selects all the rows that correspond to certain state in E00i
      for(i in 1:k){
        select_rows = c(select_rows, begin_select[i]:end_select[i])
      }
      mat_trans[select_rows,1] = state[1]
      gen_mat_2[gen_mat_2[,2] == seq_i, ]  <- mat_trans
    }
  }
  rows_selected <- which(gen_mat_2[,1] %in% state)
  return(rows_selected)
}


