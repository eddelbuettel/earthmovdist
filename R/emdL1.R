
emdL1 <- function(x, y, dims=NULL, verbose=FALSE)
{
  ## Check if x and y are identical concerning their dimensions
  if ( !identical(dim(x), dim(y)) ) {
    stop("x and y must be of same dimension and size!\n\nExecution aborted!\n")
  }

  if (is.null(dims)) {                  # no dimensions given, so test and determine

      ## if x and y are vectors, check if they are of same length
      if (is.null(dim(x))) {
          if (length(x)!=length(y)) {
              stop("x and y must be of the same size!\n\nExecution aborted!\n")
          } else {
              ## set dims to length if vector
              noDims <- 1
              dim1 <- length(x)
              dim2 <- NA
              dim3 <- NA
          }
      } else {
          ## set dims to dim() if matrix or array
          noDims <- length(dim(x))
          if (noDims > 3) {
              stop("x and y must be a two or three dimensional matrix / array!\n\nExecution aborted!\n")
          }
          dim1 <- dim(x)[1]
          dim2 <- dim(x)[2]
          dim3 <- ifelse(
                         noDims==3,
                         dim(x)[3],
                         NA
                         )
      }

      parms <- list(
                    noDims  = noDims,
                    dim1    = dim1,
                    dim2    = dim2,
                    dim3    = dim3,
                    verbose = verbose
                    )

  } else {
      if ( !is.vector(dims) ) {
          stop("dims supplied, but not a vector. Aborting!\n")
      }
      parms <- list(noDims  = length(dims),
                    dim1    = dims[1],
                    dim2    = ifelse(length(dims) >= 2, dims[2], NA),
                    dim3    = ifelse(length(dims) >= 3, dims[3], NA),
                    verbose = verbose
                    )
  }

  val <- .Call("emdL1", as.vector(x), as.vector(y), parms, PACKAGE="earthmovdist")
  d <- val[["dist"]]

  return(d)
}
