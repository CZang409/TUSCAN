# infer copy number profile

####################################################################
## runCNV

#' Inferring copy number profile (CNP).
#'
#' This function performs data filtering and whole genome copy number inferring.
#'
#' @param object TUSCAN object.
#' @param gene_order_file Gene location file, containing chromosome, start, and stop position for each gene.
#'
#' example:
#'          chr  start   stop
#' DDX11L1 chr1  11869  14412
#' WASH7P  chr1  14363  29806
#' FAM138A chr1  34554  36081
#' OR4F5   chr1  69091  70008
#' OR4F29  chr1 367640 368634
#' OR4F16  chr1 621059 622053
#' ...
#'
#' @param annotations_file File containing annotations of all spots, with two columns representing the spot barcode and the spot clustering name (support '.txt' and '.csv' file).
#' If the user do not provide an annotation file, the \code{pre_cluster} result will be used as the annotations.
#'
#' example:
#'   V1   V2
#' 10x13	1
#' 10x14	2
#' 10x15	1
#' 10x16	1
#' ...
#'
#' @param ref_group_names A vector containing the cluster name of the reference (normal) spots used for inferring CNV.
#' @param min_ave_counts_per_gene Minimum average gene expression level across all spots (default: 0.01).
#' @param min_spots_per_gene Minimum number of spots requiring expression measurements to include the corresponding gene (default: 3).
#' @param chr_exclude List of chromosomes in the reference genome annotations that should be excluded from analysis (default: c('chrX', 'chrY', 'chrM').
#' @param max_centered_threshold Maximum value a value can have after centering. Also sets a lower bound of -1 * this value (default: 3).
#' @param min_counts_per_spot Minimum number of counts per cell.
#' @param min_counts_per_spot Minimum number of counts per cell.
#' @param platform Spatial sequencing platform. Choices: "Visium" or "ST".
#' @param normalize_factor If the normalize factor is null (default setting), normalize_factor=median(UMI per spot)/UMI per spot.
#' @param window_length Length of the window for the moving average (smoothing). Should be an odd integer (default: 101).
#' @param denoise Whether perform de-noise step (default: TRUE).
#'
#' @return Return TUSCAN object containing inferred CNV data
#' @export
#'
#' @examples

runCNV <- function(object,
                   gene_order_file=NULL,
                   annotations_file=NULL,
                   ref_group_names,
                   min_ave_counts_per_gene=0.01,
                   min_spots_per_gene=3,
                   min_counts_per_spot=1,
                   chr_exclude=c('chrX', 'chrY', 'chrM'),
                   max_centered_threshold=3,
                   normalization=TRUE,
                   platform=c("Visium","ST"),
                   normalize_factor=NULL,
                   window_length=101,
                   denoise=TRUE

){

  platform <- match.arg(platform)

  print("create pre-CNV analysis object ...")

  #### input expression data
  raw.data <- as.matrix(object@counts.data)

  #### get gene order info
  ## we provide hg38 gencode file
  ## user can use their own gencode file

  if(is.null(gene_order_file)){

    data(hg38_annotation)
    gene_order <- hg38_annotation

  }else if(is.matrix(gene_order_file)||is.data.frame(gene_order_file)){
    gene_order <- gene_order_file
  }else{
    extension <- tools::file_ext(gene_order_file)
    if (extension == "txt") {
      gene_order <- read.table(gene_order_file, header=FALSE, row.names=1, sep="\t")
    } else if (extension == "csv") {
      gene_order <- read.csv(gene_order_file, header=FALSE, row.names=1)
    } else {
      stop("Unsupported file format. Only txt and csv files are supported.")
    }
  }

  names(gene_order) <- c('chr', 'start', 'stop')

  if (!is.null(chr_exclude) && any(gene_order$chr %in% chr_exclude)) {
    gene_order <- gene_order[-which(gene_order$chr %in% chr_exclude),]
  }

  #### read annotations file
  if(is.null(annotations_file)){
    input_classifications <- object@annotation[,2,drop = FALSE]
  }else{
    extension <- tools::file_ext(annotations_file)
    if (extension == "txt") {
      input_classifications <- read.table(annotations_file, header=FALSE, row.names=1, sep="\t", stringsAsFactors=FALSE, colClasses = c('character', 'character'))
    } else if (extension == "csv") {
      input_classifications <- read.table(annotations_file, header=FALSE, row.names=1, sep=",", stringsAsFactors=FALSE, colClasses = c('character', 'character'))
    } else {
      stop("Unsupported file format. Only txt and csv files are supported.")
    }
  }

  #### filtering data
  print("filtering data ...")
  print(paste(nrow(raw.data), "genes and", ncol(raw.data), "spots in raw count matrix", sep=" "))

  ## remove genes that aren't sufficiently expressed
  average_gene <- rowMeans(raw.data)

  # Find averages above a certain threshold
  if(sum(average_gene < min_ave_counts_per_gene)>0){
    raw.data2 <- raw.data[-which(average_gene < min_ave_counts_per_gene),,drop=FALSE]
  }else{
    raw.data2 <- raw.data
  }

  n_count_remove <- nrow(raw.data)-nrow(raw.data2)
  print(paste("filtered out", n_count_remove, "genes with less than", min_ave_counts_per_gene, "average counts; remaining", nrow(raw.data2), "genes", sep=" "))

  ## remove genes pass the min_spots_per_gene threshold
  spots_per_gene <- apply(raw.data2, 1, function(x)(sum(x>0)))

  if(sum(spots_per_gene < min_spots_per_gene)>0){
    raw.data3 <- raw.data2[-which(spots_per_gene < min_spots_per_gene),,drop=FALSE]
  }else{
    raw.data3 <- raw.data2
  }

  n_count_remove <- nrow(raw.data2)-nrow(raw.data3)
  print(paste("filtered out", n_count_remove, "genes that expressed in less than", min_spots_per_gene, "spots; remaining", nrow(raw.data3), "genes", sep=" "))

  #############################################################################################################################

  raw.data.new <- raw.data3
  ## extract the genes indicated in the gene ordering file:
  same_genes <- intersect(row.names(raw.data.new),row.names(gene_order))
  n_gene <- length(same_genes)
  print(paste("find",n_gene,"mathed genes in the count matrix and gene order file", sep=" "))

  # gene expression matrix has same gene names with the gene order file
  raw.data.combine <- raw.data.new[same_genes, , drop=FALSE]
  gene_order.combine <- gene_order[same_genes, , drop=FALSE]

  ## remove spots that have no counts
  cs = colSums(raw.data.combine)
  min_counts_per_spot=1
  if(sum(cs < min_counts_per_spot)>0){
    raw.data.combine2 <- raw.data.combine[,-which(cs < min_counts_per_spot),drop=FALSE]
  }else{
    raw.data.combine2 <- raw.data.combine
  }

  n_count_remove <- ncol(raw.data.combine)-ncol(raw.data.combine2)
  print(paste("filtered out", n_count_remove, "spots that expressed less than", min_counts_per_spot, "counts; remaining", ncol(raw.data.combine2), "spots", sep=" "))

  input_classifications <- input_classifications[rownames(input_classifications) %in% colnames(raw.data.combine2), , drop=FALSE]
  raw.data.combine2 <- raw.data.combine2[,colnames(raw.data.combine2) %in% rownames(input_classifications)]

  # Sort genomic position file and expression file to genomic position file
  # Order genes by genomic region
  # This code generates a vector order_names containing the row names of the gene order data frame.
  # It applies the order() function to sort the rows of gene order based on the values in the chr, start, and stop columns.
  # The with() function is used to specify the columns directly within the order() function without explicitly referencing the data frame each time. The resulting sorted row names are stored in order_names.
  # Assuming that gene_order.combine already exists and contains the columns chr, start, and stop.
  # Extract the numeric part; if no numbers are present, set it to the maximum value.
  gene_order.combine2 <- gene_order.combine
  chr_num <- gsub("chr", "", gene_order.combine2$chr)
  chr_num[chr_num=='X'] <- 23
  chr_num[chr_num=='Y'] <- 24
  chr_num[chr_num=='M'] <- 25
  gene_order.combine2$chr_num <- as.numeric(chr_num)

  # sort
  sorted_indices <- with(gene_order.combine2, order(chr_num, start, stop))
  gene_order.combine <- gene_order.combine[sorted_indices,,drop=FALSE]
  chr_levels <- unique(gene_order.combine$chr)
  gene_order.combine$chr <- factor(gene_order.combine$chr,levels = chr_levels)
  raw.data.combine2 <- raw.data.combine2[match(rownames(gene_order.combine), rownames(raw.data.combine2)),,drop=FALSE]

  ## reorder spot classifications according to expression matrix column names
  input_classifications <- input_classifications[order(match(row.names(input_classifications), colnames(raw.data.combine2))),,drop=FALSE]

  ## get indices for reference spots
  input_ref_group_names = ref_group_names
  ref_group_spot_indices = list()
  for (name_group in input_ref_group_names) {
    spot_indices = which(input_classifications[,1] == name_group)

    if (length(spot_indices) == 0 ) {
      stop(sprintf("Error, not identifying spots with classification %s", name_group))
    }
    ref_group_spot_indices[[ toString(name_group) ]] <- spot_indices
  }

  ## rest of the spots are the 'observed' set.
  all_group_names <- unique(input_classifications[,1])
  obs_group_names <- sort(setdiff(all_group_names, input_ref_group_names))

  ## define groupings according to the observation annotation names
  obs_group_spot_indices = list()
  for (name_group in obs_group_names) {
    spot_indices = which(input_classifications[,1] == name_group)
    obs_group_spot_indices[[ toString(name_group) ]] <- spot_indices
  }


  img <- object@image
  img <- img[-which(cs < min_counts_per_spot),,drop=FALSE]
  object@image <- img

  object@cnv.data <- raw.data.combine2
  object@gene_order <- gene_order.combine
  object@reference_grouped_spot_indices <- ref_group_spot_indices
  object@observation_grouped_spot_indices <- obs_group_spot_indices

  ##########################################
  ## start
  print("start cnv analysis process ...")

  ###########################################

  ## Normalization by sequencing depth
  i=1

  if(normalization==TRUE){

    print(paste("step ",i,": normalize gene expression matrix by sequencing depth ...",sep=""))

    expr_data <- object@cnv.data

    # if the normalize factor is null (default setting), normalize_factor=median(UMI per spot)/UMI per spot
    # normalized_factor can be customized by user
    if (is.null(normalize_factor)) {

      cs <- colSums(expr_data)
      normalize_factor = median(cs)/cs

    }

    normalized_data <- sweep(expr_data, STATS=normalize_factor, MARGIN=2, FUN="*")
    object@cnv.data <- normalized_data

    i=i+1
  }



  ###########################################
  ## Log transformation
  # x <- log(x+1)

  print(paste("step ",i,": perform log transformation ...",sep=""))
  x <- object@cnv.data
  object@cnv.data <- log2(x + 1)
  i=i+1

  ###########################################
  ## Subtract average reference gene expression
  ## Since we're in log space, this now becomes log(fold_change)

  print(paste("step ",i,": subtract average normal spots gene expression ...",sep=""))
  object <- subtract_normal(object)
  i=i+1

  ###########################################
  ## Apply maximum centered expression thresholds to data
  ## Outlier processing

  # For relative expression of genes centered by normal spots (log-fold-change values),
  # any value of abs(log(x+1)) above 'max_centered_threshold' (default=3) will be replaced by 3, -3.
  # If the value is greater than 3, it is replaced with 3; if it is less than -3, it is replaced with -3, so that all values are within the range of plus or minus 3, and there are no abnormal values.
  print(paste("step ",i,": Apply maximum centered expression thresholds to data ...",sep=""))

  threshold = max_centered_threshold
  object@cnv.data[object@cnv.data > threshold] <- threshold
  object@cnv.data[object@cnv.data < (-1 * threshold)] <- -1 * threshold
  i=i+1

  ###########################################
  ## Smooth the data along chromosome with gene windows
  print(paste("step ",i,": smooth the data along chromosome with gene windows ...",sep=""))

  object <- smooth_by_chromosome(object, window_length=window_length, smooth_ends=TRUE)

  i=i+1

  ###########################################
  ## Center spots/observations after smoothing.
  print(paste("step ",i,": center spots/observations after smoothing ...",sep=""))
  object <- center_spot_expr_across_chromosome(object, method="median")

  i=i+1

  ###########################################
  ## Subtract average reference (adjustment)
  print(paste("step ",i,": subtract average reference ...",sep=""))

  object <- subtract_ref_expr_from_obs(object)
  i=i+1

  ###########################################
  ## Invert log transform  (convert from log(FC) to FC)
  print(paste("step ",i,": invert log transform ...",sep=""))

  x <- object@cnv.data
  object@cnv.data <- 2^x

  ###########################################
  ## denoise (optional)

  if(denoise){
    object <- clear_noise_via_ref_mean_sd(object)
  }
  return(object)
}


########################################################################
### internal functions

#' Subtracting the mean of the reference expr distributions from the observed spots.
#'
#' @param object TUSCAN object.
#'
#' @return Return TUSCAN object containing the reference subtracted values.
#'
#' @keywords internal
#' @noRd
#'

subtract_normal <- function(object) {

  ref_groups = object@reference_grouped_spot_indices
  cnv.data <- object@cnv.data

  ref_grp_gene_means = do.call(rbind, apply(cnv.data, 1, function(x) {
    grp_means = vapply(ref_groups, function(ref_group) {
      mean(x[ref_group])
    }, double(1))
    names(grp_means) = names(ref_groups)
    as.data.frame(t(data.frame(grp_means)))
  }))

  rownames(ref_grp_gene_means) = rownames(cnv.data)

  my.rownames = rownames(object@cnv.data)
  my.colnames = colnames(object@cnv.data)

  subtract_normal_expr_fun <- function(row_idx) {
    gene_means <- as.numeric(ref_grp_gene_means[row_idx, , drop=TRUE])
    gene_means_mean <- mean(gene_means)
    x <- as.numeric(object@cnv.data[row_idx, , drop=TRUE])
    row_init = rep(0, length(x))

    row_init <- x - gene_means_mean
    return(row_init)
  }

  object@cnv.data <- do.call(rbind, lapply(seq_len(nrow(object@cnv.data)), subtract_normal_expr_fun))

  rownames(object@cnv.data) <- my.rownames
  colnames(object@cnv.data) <- my.colnames

  return(object)
}


#' clear_noise_via_ref_mean_sd()
#'
#' Define noise based on the standard deviation of the reference spot expression data.
#' The range to remove noise would be mean +- sdev * 1.5.
#' Data points defined as noise are set to average expression level of normal reference.
#'
#' @param object TUSCAN object
#'
#' @return Return TUSCAN object
#'
#' @keywords internal
#' @noRd
#

clear_noise_via_ref_mean_sd <- function(object) {

  sd_amplifier=1.5
  ref_idx = unlist(object@reference_grouped_spot_indices)
  vals = object@cnv.data[,ref_idx]
  mean_ref_vals = mean(vals)
  mean_ref_sd <- mean(apply(vals, 2, function(x) sd(x, na.rm=TRUE))) * sd_amplifier
  upper_bound = mean_ref_vals + mean_ref_sd
  lower_bound = mean_ref_vals - mean_ref_sd
  object@cnv.data[object@cnv.data > lower_bound & object@cnv.data < upper_bound] = mean_ref_vals

  return(object)
}

#' center_spot_expr_across_chromosome()
#'
#' Centers expression data across all genes for each spot, using the spot mean or median expression value as the center.
#'
#' @param object TUSCAN object
#' @param method Method to select the center of the spot expression value. (default: 'mean', options: 'mean,'median')
#'
#' @return Return TUSCAN object
#'
#' @keywords internal
#' @noRd
#'

center_spot_expr_across_chromosome <- function(object, method="median") {
  object@cnv.data <- .center_columns(object@cnv.data, method)
  return(object)
}

#' @keywords internal
#' @noRd
#'

.center_columns <- function(expr_data, method) {
  # Center within columns (spots)
  if (method == "median") {
    median_gene_spot <- apply(expr_data, 2, function(x) { median(x, na.rm=TRUE) } )
    expr_data <- sweep(expr_data, 2, median_gene_spot, "-")
  }
  else {
    # by mean
    mean_gene_spot <- apply(expr_data, 2, function(x) { mean(x, na.rm=TRUE) } )
    expr_data <- sweep(expr_data, 2, mean_gene_spot, "-")
  }
  return(expr_data)
}


#' Subtracting the mean of the reference expr distributions from the observed spots.
#'
#' Remove the average of the genes of the reference observations from all observations' expression.
#'
#' @param object TUSCAN object
#'
#' @return return TUSCAN object containing the reference subtracted values.
#'
#' @keywords internal
#' @noRd
#'

subtract_ref_expr_from_obs <- function(object) {

  ref_groups = object@reference_grouped_spot_indices
  ref_grp_gene_means <- .get_normal_gene_mean_bounds(object@cnv.data, ref_groups)
  object@cnv.data <- .subtract_expr(object@cnv.data, ref_grp_gene_means)

  return(object)

}

#' @keywords internal
#' @noRd
#'

.get_normal_gene_mean_bounds <- function(cnv.data, ref_groups) {

  get_indiv_gene_group_means_bounds_fun <- function(x) {
    grp_means = c()
    grp_means <- vapply(ref_groups, function(ref_group) {
      mean(x[ref_group])
    }, double(1))

    names(grp_means) <- names(ref_groups)
    return(as.data.frame(t(data.frame(grp_means))))
  }

  gene_ref_grp_means <- do.call(rbind, apply(cnv.data, 1, get_indiv_gene_group_means_bounds_fun))
  rownames(gene_ref_grp_means) <- rownames(cnv.data)
  return(gene_ref_grp_means)

}

#' @keywords internal
#' @noRd
#'

.subtract_expr <- function(expr_matrix, ref_grp_gene_means) {

  my.rownames = rownames(expr_matrix)
  my.colnames = colnames(expr_matrix)

  subtract_normal_expr_fun <- function(row_idx) {

    gene_means <- as.numeric(ref_grp_gene_means[row_idx, , drop=TRUE])
    gene_means_mean <- mean(gene_means)
    x <- as.numeric(expr_matrix[row_idx, , drop=TRUE])
    row_init = rep(0, length(x))
    row_init <- x - gene_means_mean

    return(row_init)
  }

  subtr_data <- do.call(rbind, lapply(seq_len(nrow(expr_matrix)), subtract_normal_expr_fun))
  rownames(subtr_data) <- my.rownames
  colnames(subtr_data) <- my.colnames

  return(subtr_data)

}

#' smooth_by_chromosome()
#'
#' @param object TUSCAN object.
#' @param window_length Length of window (number of genes) for the moving average.
#' @param smooth_ends Perform smoothing at the ends of the chromosomes (default: TRUE).
#'
#' @return return TUSCAN object
#'
#' @keywords internal
#' @noRd
#'

smooth_by_chromosome <- function(object, window_length, smooth_ends=TRUE) {

  gene_chr_listing = object@gene_order$chr
  chrs = unlist(unique(gene_chr_listing))

  for (chr in chrs) {
    chr_genes_indices = which(gene_chr_listing == chr)
    chr_data=object@cnv.data[chr_genes_indices, , drop=FALSE]
    if (nrow(chr_data) > 1) {
      smoothed_chr_data = .smooth_window(data=chr_data,
                                         window_length=window_length)
      object@cnv.data[chr_genes_indices, ] <- smoothed_chr_data
    }
  }

  return(object)
}

#' @keywords internal
#' @noRd
#'

.smooth_window <- function(data, window_length) {

  if (window_length < 2){
    flog.warn("window length < 2, returning original unmodified data")
    return(data)
  }

  ## Fix ends that couldn't be smoothed since not spanned by win/2 at ends.
  data_sm <- apply(data,
                   2,
                   .smooth_helper,
                   window_length=window_length)

  ## Set back row and column names
  row.names(data_sm) <- row.names(data)
  colnames(data_sm) <- colnames(data)
  return(data_sm)
}

#' Helper function for smoothing with a moving average
#'
#' @param obs_data Data to smooth.
#' @param window_length Length of the window to smooth.
#'
#' @return Return data smoothed
#'
#' @keywords internal
#' @noRd
#'

.smooth_helper <- function(obs_data, window_length) {
  # strip NAs out and replace after smoothing
  orig_obs_data = obs_data
  nas = is.na(obs_data)
  obs_data = obs_data[!nas]
  obs_length <- length(obs_data)
  end_data <- obs_data
  tail_length = (window_length - 1)/2
  if (obs_length >= window_length) {
    end_data <- .smooth_center_helper(obs_data, window_length)
  }
  # end_data will have the end positions replaced with mean values, smoothing just at the ends.

  obs_count <- length(obs_data)
  numerator_counts_vector = c(seq_len(tail_length), tail_length + 1, c(tail_length:1))

  # defining the iteration range in cases where the window size is larger than the number of genes. In that case we only iterate to the half since the process is applied from both ends.
  iteration_range = ifelse(obs_count > window_length, tail_length, ceiling(obs_count/2))
  for (tail_end in seq_len(iteration_range)) {
    end_tail = obs_count - tail_end + 1
    d_left = tail_end - 1
    d_right = obs_count - tail_end
    d_right = ifelse(d_right > tail_length, tail_length, d_right)
    r_left = tail_length - d_left
    r_right = tail_length - d_right
    denominator = (((window_length - 1)/2)^2 + window_length) - ((r_left * (r_left + 1))/2) - ((r_right * (r_right + 1))/2)
    left_input_vector_chunk = obs_data[seq_len(tail_end + d_right)]
    right_input_vector_chunk = obs_data[(end_tail - d_right):obs_length]
    numerator_range = numerator_counts_vector[(tail_length + 1 - d_left):(tail_length + 1 + d_right)]
    end_data[tail_end] = sum(left_input_vector_chunk * numerator_range)/denominator
    end_data[end_tail] = sum(right_input_vector_chunk * rev(numerator_range))/denominator
  }
  orig_obs_data[! nas] = end_data  # replace original data with end-smoothed data
  return(orig_obs_data)
}


#' Smooth vector of values over the given window length
#'
#' @param obs_data Vector of data to smooth with a moving average.
#' @param window_length Length of the window for smoothing.
#'
#' @return Return vector of values smoothed with a moving average.
#'
#' @keywords internal
#' @noRd
#'

.smooth_center_helper <- function(obs_data, window_length){

  nas = is.na(obs_data)
  vals = obs_data[! nas]
  custom_filter_denominator = ((window_length-1)/2)^2 + window_length
  custom_filter_numerator = c(seq_len((window_length-1)/2), ((window_length-1)/2)+1, c(((window_length-1)/2):1))
  custom_filter = custom_filter_numerator/rep(custom_filter_denominator, window_length)

  smoothed = stats::filter(vals, custom_filter, sides=2)
  ind = which(! is.na(smoothed))
  vals[ind] = smoothed[ind]
  obs_data[! nas] = vals

  return(obs_data)
}

