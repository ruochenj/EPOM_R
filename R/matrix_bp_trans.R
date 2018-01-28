#' Transform BigWig file into a matrix
#' @description Used in step 1 to transform BigWig files into a matrix by averaging the signals within every window with a fixed size (interval).
#'
#' @param bw A BigWig file to be transformed.
#' @param interval An integer specifying the bandwidth used to average the bigwig signals. We set it to 200 in order to match with the state information provied by bed files.
#' @return A matrix with four columns. The four columns includes: chromatin state name (eg. chr1), begin of each region (eg. 201), width of each region (eg. 200) and signal of each region (there is only one column of signal in the output) after transformation.
#' @export
matrix_bp_trans <- function(bw, interval){
  seq_name <- bw@seqnames
  start = bw@ranges@start
  width = bw@ranges@width
  score = bw$score

  #seq = "chr1"
  sig <- c()
  begin <- c()
  new_width <- c()
  seqname <- c()
  for(seq in levels(seq_name)){
    seq_select <- as.character(seq_name) %in% seq
    sig_trans <- score[seq_select]
    width_trans <- width[seq_select]
    one_bp <-  rep(sig_trans, width_trans)
    #account for the last component
    ln <- length(one_bp) %% interval
    len_200 <- length(one_bp)-ln
    last_element <- sum(tail(one_bp, ln)/ln)
    one_bp <- one_bp[1:len_200]
    m1 <- matrix(one_bp, length(one_bp)/interval, interval, byrow = TRUE)
    pre_element <- rowSums(m1)/interval
    sig_i <- c(pre_element, last_element)
    sig <- c(sig, sig_i)
    new_width_i <- rep(interval, length(pre_element))
    new_width_i <- c(new_width_i, ln)
    new_width <- c(new_width, new_width_i)
    new_begin <- seq(1, len_200+1, interval)
    begin <- c(begin, new_begin)
    seqname <- c(seqname ,rep(seq, length(sig_i)))
  }
  return(data.frame(cbind(seqname, begin, new_width, sig)))
}
