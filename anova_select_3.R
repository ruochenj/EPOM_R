#' Use ANOVA to select candidate associated regions that are significant
#' @description Step 2 in the EPOM method: use ANOVA to test whether a histone mark has the same group mean signals across different cell and tissue types.
#'
#' @param cell_type A vector indicating the cell type of each sample.
#' @param row_select A vector indicating the row indicies corresponding to a certain state with respect to the original 200 bp matrix.
#' @param matrix_gen A matrix which containing all the samples where the signal are correponding to the same bp interval. For example, if we have 127 available samples, the dimension of the matrix will be x*(3+127). x is the total length of 200 bp regions. First 3 columns are chrom name (eg. chr1), begin (eg. 201), width (eg. 200) and the other 127 are signals (eg. 0.04) of each sample.
#' @param alpha_1 A number indicating the significance level used in the ANOVA. The default is 1e-10, the larger this number is, the more regions will be selected.
#' @param in_dir (Optional) Path for the .RData file to be read in.
#' @param save A Boolean indicating whether the output should be saved.
#' @param out_dir (Optional) Path for the output of the function.
#' @param cores Indicate the number of cores you want to use in your server. The default value is 16.
#'
#' @return This function returns a list l (lower case "L") including (1) a vector called select_reg which indicates
#' the indices for the regions in the original matrix. (2) a matrix called mat_anova which is the matrix
#' after anova selection.
#' @return  An RData file called "mat_anova4state.RData" (eg. "mat_anova425.RData" indicates the RData file after processing anova for state 25) includes the object l and l$select_reg corresponds to (1) and l$mat_anova corresponds to (2).
#' @export
anova_select3 <- function(cell_type, row_select, matrix_gen, state = "enhancer", alpha_1 = 1e-10, in_dir = NULL, save = TRUE, out_dir = NULL, cores = 16){
  if(is.null(in_dir)){
    in_dir = getwd()
  }
  if(is.null(out_dir)){
    out_dir = getwd()
  }
  flag = eval(parse(text =
                      paste("!file.exists('", in_dir, "/mat_anova4", state, ".RData')", sep = "")
  ))
  if(flag){
    n = dim(matrix_gen)[1]
    mat_state <- matrix_gen[row_select,]
    mat_2 <- as.matrix(mat_state[,-c(1,2,3)])
    start = mat_state[,c(1,2)]
    rm(mat_state)
    gc()
    j = 0

    cell_type <- as.factor(cell_type)
    levels_ori <- levels(cell_type)
    levels(cell_type) <- 1:length(levels(cell_type))
    levels_trans <- levels(cell_type)
    eval(parse(text =
                 paste("write.table(cbind(levels_ori, levels_trans), file = '", out_dir, "/cell_type_correspondence.txt', row.names = FALSE)", sep = "")
    ))

    chosen_group <- c()
    not_chosen <- c()
    for(i in levels(cell_type)){
      i = as.numeric(i)
      if(sum(cell_type %in% i) != 1){
        chosen_group <- c(chosen_group, i)
      }
      else{
        not_chosen <- c(not_chosen, which(cell_type %in% i))
      }
    }
    if(length(not_chosen) != 0 ){
      cell_type <- as.factor(as.character(cell_type[-not_chosen]))
    }else{
      cell_type <- as.factor(as.character(cell_type))
    }


    p_val = mclapply(1:dim(mat_2)[1], function(x){
      return(anova(lm(mat_2[x,]~as.factor(cell_type)))$'Pr(>F)'[1])
    }, mc.cores = cores)
    p_val = unlist(p_val)

    index_select <- which(p_val < alpha_1/n)
    select_reg <- row_select[index_select]
    mat_anova <- mat_2[index_select,]
    start <- start[index_select,]
    mat_anova <- cbind(start, mat_anova)
    l = list(select_reg = select_reg, mat_anova = mat_anova)
    eval(parse(text =
                 paste("save(l, file = '", out_dir, "/mat_anova4", state, ".RData')", sep = "")
    ))
    return(l)
  }else{
    eval(parse(text =
                 paste("load('", in_dir, "/mat_anova4", state, ".RData')", sep = "")
    ))
    return(l)
  }
}
