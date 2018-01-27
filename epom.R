#' Calculate the epom scores between every pair of samples.
#' @description This function outputs the epom scores for given parameters. It’s package for the paper: Epigenome overlap measure (EPOM) for comparing tissue/cell types based on chromatin states by Wei Vivian Li, Zahra S. Razaee and Jingyi Jessica Li (2016).
#'
#' @param bw_list A vector containing the full path and file name for all the big wig files.The length and order of bw_list, bed_list and bed_sep should be matched.
#' @param bed_list A vector containing the full path and file name for all the bed files. The length and order of bw_list, bed_list and bed_sep should be matched.
#' @param bed_head This is used in reading in the bed file. The program uses read.table to read the bed files. This parameter corresponds to parameter header in read.table. You can use ?read.table for further reference.
#' @param bed_sep This is used in reading in the bed file. The program uses read.table to read the bed files. This parameter corresponds to parameter sep in read.table. You can use ?read.table for further reference.
#' @param cell_type A vector containing all the cell type information for each cell. The length and order of bw_list, bed_list and bed_sep should be matched.
#' @param histone_mark A string indicating the current histone_mark you are processing
#' @param state A vector indicating the chromatin states of interest.
#' @param state_name A string indicating the name for the state. For example, "enhancer1".
#' @param in_dir (Optional) Path of the .RData file to be read in.
#' @param save A Boolean indicating whether the output should be saved.
#' @param out_dir (Optional) Path of the output of the function.
#' @param cores An integer indicating the number of cores to use in parallel computation. The default value is 16.
#' @param m An integer indicating the significant values required for a region to be selected as the associated region for a cell type.
#' @param bed_head This is used in reading in the bed file. The program uses read.table to read the bed files. This parameter corresponds to parameter header in read.table. You can use ?read.table for further reference.
#' @param alpha1 A number indicating the significance level used in ANOVA. The default is 1e-10, the larger this number is, the more regions will be selected.
#' @param alpha2 A number indicating the significance level used in the t tests. The default is 0.01.
#'
#' @return A matrix output the epom score.
#' @references Li, Wei Vivian, Zahra S. Razaee, and Jingyi Jessica Li. "Epigenome overlap measure (EPOM) for comparing tissue/cell types based on chromatin states." BMC genomics 17.1 (2016): S10.

epom <- function(bw_list, bed_list, bd_header=FALSE, bed_sep="\t", cell_type, histone_mark, state, state_name, in_dir = NULL, save = TRUE, out_dir = NULL, alpha_1 = 1e-10, cores = 16, alpha_2 = 0.01, m = 13){
  if(is.null(in_dir)){
    in_dir = getwd()
  }
  if(is.null(out_dir)){
    out_dir = getwd()
  }
  #check input
  len_wig = length(bw_list)
  len_bed = length(bed_list)
  len_cel = length(cell_type)
  if(len_wig != len_bed){
    stop("number of bigwig files and bed files not matched")
  }
  if(len_wig != len_cel){
    stop("number of bigwig files and number of cell types not matched")
  }
  tb <- table(cell_type)
  if(sum(tb %in% 1) > 0){
    stop("cell type with only one cell not excluded")
  }
  eval(parse(text =
               paste("flag = !file.exists('", in_dir, "/matrix_gen.RData') & !file.exists('", in_dir, "/row_select.RData')", sep = "")
  ))
  if(flag){
    i = 1
    eval(parse(text=
                 paste("bw",i,"<- import.bw('", bw_list[i],"')", sep = "")
    ))
    eval(parse(text=
                 paste("bed_",i,"<- as.data.frame(read.table('", bed_list[i],"', header =",bed_header, ", sep = '", bed_sep, "', stringsAsFactors=TRUE, quote =''))", sep = "")
    ))

    matrix_gen <- matrix_bp_trans(bw1, 200)
    index4state <- list()
    index4state[[1]] <- index_num4state(matrix_gen, bed_1, state)

    for(i in 2:len_wig){
      print(paste(i, "th file reading in", sep = ""))
      eval(parse(text=
                   paste("bw",i,"<- import.bw('", bw_list[i],"')", sep = "")
      ))
      eval(parse(text=
                   paste("bed_",i,"<- as.data.frame(read.table('", bed_list[i],"', header =",bed_header, ", sep = '", bed_sep, "', stringsAsFactors=TRUE, quote =''))", sep = "")
      ))

      eval(parse(text =
                   paste("matrix_gen_", i,"= matrix_bp_trans(bw",i,", 200)", sep = "")
      ))
      eval(parse(text =
                   paste("matrix_gen <- cbind(matrix_gen, matrix_gen_",i,"[,4])", sep = "")
      ))
      eval(parse(text =
                   paste("index4state[[",i,"]] <- index_num4state(matrix_gen_",i,", bed_",i,", state)", sep = "")
      ))

      #clean the storage space
      eval(parse(text =
                   paste("rm(bw", i,")", sep = "")
      ))
      eval(parse(text =
                   paste("rm(matrix_gen_",i,")", sep = "")
      ))
      gc()
    }

    row_select <- Reduce(union, index4state)
  }
  else{
    print("begin loading")
    #load("matrix_gen.RData")
    eval(parse(text =
                 paste("load('", in_dir, "/matrix_gen.RData')", sep = "")
    ))
    #load("row_select.RData")
    eval(parse(text =
                 paste("load('", in_dir, "/row_select.RData')", sep = "")
    ))
    print("loading complete")
  }

  if(save){
    eval(parse(text =
                 paste("save(matrix_gen, file = '", out_dir, "/matrix_gen.RData')", sep = "")
    ))
    #save(matrix_gen, file = "matrix_gen.RData")
    eval(parse(text =
                 paste("save(row_select, file = '", out_dir, "/row_select.RData')", sep = "")
    ))
  }


  l = anova_select3(cell_type = cell_type, row_select = row_select, matrix_gen = matrix_gen, state = state_name, alpha_1 = alpha_1, in_dir = in_dir, save = save, out_dir = out_dir, cores = cores)
  print("anova_finished")

  rm(matrix_gen)
  rm(row_select)
  gc()
  print("row_selected_from_state")
  print("begin_anova")

  l = t_select3(select_reg = l$select_reg, mat_anova = l$mat_anova, cell_type = cell_type, state = state_name, alpha_2 = alpha_2, m = m, in_dir = in_dir, save = save, out_dir = out_dir, cores = cores)
  print("t_test_finished")

  score = epom_score(l, cell_type, save = save, out_dir = out_dir)
  eval(parse(text =
               paste("write.csv(score, file = '", out_dir, "/", histone_mark, "_score.csv')", sep = "")
  ))
  return(score)
}
