#' Build a Spatial Weight Matrix
#'
#' Builds the square spatial weight matrix `ww` from an index matrix and a
#' corresponding value matrix.
#'
#' @param ww_index Matrix. Spatial-neighbor index matrix. Each row identifies
#'   the neighbors of the corresponding spatial unit. Positive entries are
#'   interpreted as valid neighbor indices.
#' @param ww_values Matrix. Spatial-weight value matrix. Must have the same
#'   structure as `ww_index`; values are assigned to the corresponding neighbor
#'   positions in the output matrix.
#'
#' @return A square numeric matrix containing the spatial weights. The matrix
#'   has dimension `pp` by `pp`, where `pp = nrow(ww_index)`. Row and column
#'   names are inherited from `rownames(ww_index)`.
#'
#' @details
#' The function checks that the resulting spatial matrix:
#' \itemize{
#'   \item contains no `NA` values;
#'   \item has a zero diagonal;
#'   \item has no rows with all zero entries.
#' }
#'
#' @examples
#' ww_index <- matrix(
#'   c(2, 3,
#'     1, 3,
#'     1, 2),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#'
#' ww_values <- matrix(
#'   c(0.5, 0.5,
#'     0.5, 0.5,
#'     0.5, 0.5),
#'   nrow = 3,
#'   byrow = TRUE
#' )
#'
#' rownames(ww_index) <- c("1", "2", "3")
#' rownames(ww_values) <- c("1", "2", "3")
#'
#' ww <- build_spatial_matrix(
#'   ww_index = ww_index,
#'   ww_values = ww_values
#' )
#'
#' ww
#'
#' @export
build_spatial_matrix <- function(ww_index, ww_values) {

  ## This function builds the square spatial matrix ww.

  pp <- dim(ww_index)[1]

  ww <- matrix(0, nrow = pp, ncol = pp)
  dimnames(ww)[[1]] <- dimnames(ww)[[2]] <- dimnames(ww_index)[[1]]

  for (ii in seq_len(pp)) {
    indices <- as.character(ww_index[ii, ww_index[ii, ] > 0])
    ww[ii, indices] <- ww_values[ii, ww_index[ii, ] > 0]
  }

  # Check the consistency of the spatial matrix.
  if (sum(is.na(ww)) > 0) {
    return(list(error = "The spatial matrix cannot have NA values."))
  }

  if (sum(diag(abs(ww))) > 0) {
    return(list(error = "The spatial matrix must have a zero diagonal."))
  }

  zero_row_index <- apply(abs(ww), 1, sum) == 0

  if (sum(zero_row_index) > 0) {
    return(list(error = "The spatial matrix has one or more zero rows."))
  }

  # Output.
  ww
}


