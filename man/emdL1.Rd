\name{emdL1}
\alias{emdL1}
\title{Earth Mover's Distance between histograms}
\description{
  The \code{emdL1} provides a wrapper to the emd-L1 distance norm by
  Haibin Ling and Kazunori Okada.
}
\usage{
  emdL1(x, y, dims=NULL, verbose=FALSE)
}
\arguments{
  \item{x}{vector, matrix or 3D array containing signature 1}
  \item{y}{vector, matrix or 3D array containing signature 2}
  \item{dims}{optional vector of dimensions of input data, see Details}
  \item{verbose}{optional boolean flag indicating verbose operation}
}
\value{A real number with the EMD-L1 distance metric between the two
  histograms}
\details{
  This function provides an implementation of the Earth Movers Distance (EMD),
  using the L1 (taxicab) metric.

  'The earth mover's distance (EMD) is a mathematical measure of the
  distance between two distributions over some region D. Informally, if
  the distributions are interpreted as two different ways of piling up a
  certain amount of dirt over the region D, the EMD is the minimum cost
  of turning one pile into the other; where the cost is assumed to be
  amount of dirt moved times the distance by which its is moved'
  

  The distance 'the dirt' has to be moved is in the implementation of
  emd-L1, the L1 distance: if we have two points (x1, y1) and (x2, y2),
  the L1 distance (or taxicab metric) is abs(x2-x1)+abs(y2-y1).

  The above definition only covers the case, where the 'amount of dirt'
  in both distributionssignatures is the same. In this case, the EMD is identical
  to the 1st Mallows distance or 1st Wasserstein distance.

  For details on how the emd-L1 distance is computed in this algorithm, see
  @ARTICLE{
    author = {Ling, Haibin and Okada, K.},
    title = {An efficient earth mover's distance algorithm for robust histogram
      comparison},
    journal = {IEEE Transactions on Pattern Analysis and Machine Intelligence},
    year = {2007},
    volume = {29},
    pages = {840-53},
    number = {5},
    month = {05},
    doi = {10.1109/TPAMI.2007.1058},
    publication-type = {J},
    publisher = {IEEE Comput. Soc},
    type = {Journal Paper}
  }
  
  The function has to be called by providing two signatures (x and
  y). these can be vectors (1D signatures), matrix (2D signatures) or 3D
  arrays (3D signatures). The dimensions of x and y have to be the same,
  and x and y should not contain any missing or infinite values.

  If the vector dims is supplied, x and y have to be vectors, and dims
  will be used to interpet those as vectors (length(dims)==1), matrix
  (length(dims)==2) or 3D array (length(dims)==3).

  If the vector dims is not supplied and dim(x)==dim(y), dims is set
  internally to dim(x).
}
\author{Dirk Eddelbuettel and Rainer M Krug (Rainer@krugs.de) for the R package; Haibin Ling and Kazunori
  Okada for emd-L1.
}
\examples{
  set.seed(42)
  x <- rnorm(100)
  y <- x + rnorm(100)/10
  emdL1(x, y, verbose=TRUE)
}