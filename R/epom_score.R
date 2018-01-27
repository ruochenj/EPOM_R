#' Calculate epom_score after the vector chosen by t_test
#' @description Step 4, Calculate epom_score after the vector chosen by t_test.
#'
#' @param l_t A list of length of chosen groups which have more than one element. Each element of the list is a matrix which have three columns: column one indicates the selected indicies, column two indicates the chr, colun hte
#' @param cell_type A vector indicates the cell type of the ith gene.
#' @return A matrix output the epom score.
#' @export
epom_score <- function(l_t, cell_type, save = TRUE, out_dir = NULL){
  if(is.null(out_dir)){
    out_dir = getwd()
  }
  cell_type <- as.factor(cell_type)
  levels_ori <- levels(cell_type)
  levels(cell_type) <- 1:length(levels(cell_type))
  levels_trans <- levels(cell_type)
  if(!file.exists("cell_type_correspondence.txt")){
    print("file does not exist")
    write.table(cbind(levels_ori, levels_trans), file = "cell_type_correspondence.txt", row.names = FALSE)
  }else{
    print("file found!")
  }

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

  l_t2 = lapply(l_t, function(x){
    if(!is.null(x)){
      return(x[,-c(2,3)])
    }
    else{
      return(NULL)
    }
  })
  n = length(Reduce(union, l_t2))
  k = length(chosen_group)
  p = matrix(1, k, k)
  for(i in 1:(k-1)){
    for(j in i:k){
      A_obj <- l_t[[chosen_group[i]]][,1]
      B_obj <- l_t[[chosen_group[j]]][,1]
      A_abs = length(A_obj)
      B_abs = length(B_obj)
      print(A_abs)
      print(B_abs)
      intersect_AB = length(intersect(A_obj, B_obj))
      print("intersect got")
      print(intersect_AB)
      min_AB = min(A_abs, B_abs)
      print("minimum got")
      print(min_AB)
      p_ij = lapply(intersect_AB:min_AB, function(x){
        exp(lchoose(n, x)+lchoose(n-x, A_abs-x)+lchoose(n-A_abs, B_abs-x)-lchoose(n, A_abs)-lchoose(n, B_abs))
      })
      p[i,j] = sum(unlist(p_ij))
      print(p[i,j])
      p[j,i] = p[i,j]
    }
  }
  p <- data.frame(p)
  epom_score = -log10(p*16*16/2)
  epom_score[epom_score > 300] = 300
  epom_score[epom_score < 0] = 0
  epom_score <- data.frame(epom_score)
  colnames(epom_score) <- chosen_group
  rownames(epom_score) <- chosen_group
  eval(parse(text =
               paste("write.table(p, file = '", out_dir, "/score table', row.names = TRUE, col.names = TRUE)", sep = "")
  ))
  write.table(p, file = "score table", row.names = TRUE, col.names = TRUE)
  return(epom_score)
}


