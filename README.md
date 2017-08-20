## earthmovdist

[![Build Status](https://travis-ci.org/eddelbuettel/earthmovdist.png)](https://travis-ci.org/eddelbuettel/earthmovdist) 
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html) 


Rcpp bindings to the Emd-L1 distance library for the Earth Mover's Distances

### About

Earth Mover's Distance is a popular distance measure. This package provides the Emd-L1 library.

### Library License

Note the rather stringent _non-commercial_ licensing terms in the C++ source code and
copied into [LICENSE](LICENSE).

### See Also

The [emdist](https://cloud.r-project.org/package=emdist) package provides Earth Mover's
Distance calculations under a more common open source license.

### History

We wrote this package _many_ years ago as a first and rather simple illustration of
[Rcpp](http://dirk.eddelbuettel.com/code/rcpp.html).  At the time, no other Earth Mover's
Distance implementations were available for R. These days, you may wan to use
[emdist](https://cloud.r-project.org/package=emdist) or one of the other, newer
distance-computing packages.

### Authors

Dirk Eddelbuettel and Rainer Krug

### License

GPL (>= 2) for the code added by this package, but note the restrictive library license.
