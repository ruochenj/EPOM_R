#' Use t test to select associated regions of a certian cell type.
#' @description Step 3 in the EPOM method: use t test to identify associated regions of a certain cell type.
#'
#' @param select_reg A numeric vector indicating the index of rows corresponding to the orignal transformed matrix after selecing by ANOVA.
#' @param mat_anova A matrix containing the signals only on regions selected in ANOVA (rows for regions and columns for samples).
#' @param cell_type A character vector indicating the cell type each sample belongs to.
#' @param state_name A character string indicating the name for the state. For example, "enhancer1".
#' @param m An integer indicating the significant values required for a region to be selected as the associated region for a cell type.
#' @param alpha_2 A number indicating the significance level used in the t tests. The default value is 0.01.
#' @param in_dir (Optional) Path of the .RData file to be read in.
#' @param save A Boolean indicating whether the output should be saved.
#' @param out_dir (Optional) Path of the output of the function.
#' @param cores Indicate the number of cores you want to use in your server. The default value is given by function detectCores()-1.
#'
#' @return A list l (lower case "L") that contains selected enhancer region corresponding to the cell types which have more than one corresponding genes.
#' @return An RData file called "enhancer_signal.RData" will be saved corresponding to the result of the t test, each row indicates how t-test results are significant for each cell type.
#' @return An RData file called "mat_t.RData" will be saved which includes the list l, where l[[i]] corresponds to the selected region of cell type i. Some of l[[i]] might be NULL if the cell type has only one gene corresponding to it.
#' @export

t_select <- function(select_reg, mat_anova, cell_type, state_name = "enhancer", alpha_2 = 0.01, m = 13, in_dir = NULL, save = TRUE, out_dir = NULL, cores = detectCores()-1){
  if(is.null(in_dir)){
    in_dir = getwd()
  }
  if(is.null(out_dir)){
    out_dir = getwd()
  }
  start = mat_anova[,c(1,2)]
  mat_t = mat_anova[,-c(1,2)]
  rm(mat_anova)

  cell_type <- as.factor(cell_type)
  levels_ori <- levels(cell_type)
  levels(cell_type) <- 1:length(levels(cell_type))
  levels_trans <- levels(cell_type)
  eval(parse(text =
               paste("write.table(cbind(levels_ori, levels_trans), file = '", out_dir, "/cell_type_correspondence.txt', row.names = FALSE)", sep = "")
  ))

  group_list <- list()
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

  for(i in levels(cell_type)){
    i = as.numeric(i)
    group_list[[i]] = which(cell_type %in% i)
  }

  loop_length = dim(mat_t)[1]

  enhancer_region <- mclapply(1:loop_length, function(k){
    if(k %% 200000 == 0){
      print(k/loop_length)
    }
    detect_matrix <- matrix(1, length(chosen_group), length(chosen_group))
    for(i in 1:(length(chosen_group)-1)){
      for(j in (i+1):length(chosen_group)){
        #k is one iteration for each row in anova selected matrix
        item1 = as.numeric(mat_t[k,group_list[[chosen_group[i]]]])
        item2 = as.numeric(mat_t[k,group_list[[chosen_group[j]]]])
        if(sd(item1) != 0 & sd(item2) != 0){
          detect_matrix[i,j] = t.test(item1, item2, mu = 0, alternative = "greater", var.equal = TRUE, paired = FALSE)$p.value
          detect_matrix[j,i] = 1-detect_matrix[i,j]
        }
      }
    }
    result <- apply(detect_matrix < alpha_2, 1, sum)
    if(sum(result >= m) >0){
      record_vec = rep(0, length(result))
      record_vec[which(result > m)] = 1
      enhancer_region <- c(k, record_vec)
    }
    else{
      enhancer_region = NULL
    }
    return(enhancer_region)
  }, mc.cores = cores)

  enhancer_region = matrix(unlist(enhancer_region), ncol = length(chosen_group)+1, byrow = TRUE)
  print("begin selection")
  eval(parse(text =
               paste("save(enhancer_region, file = '", out_dir, "/enhancer_signal_dopar.RData')", sep = "")
  ))
  print("enhancer_signal_saved")
  l = list()
  for(i in 1:length(chosen_group)){
    row_select = enhancer_region[which(enhancer_region[,i+1] == 1),1]
    l[[chosen_group[i]]] = cbind(select_reg[row_select], start[row_select,])
  }
  eval(parse(text =
               paste("save(l, file = '", out_dir, "/mat_t_dopar.RData')", sep = "")
  ))
  return(l)
}
