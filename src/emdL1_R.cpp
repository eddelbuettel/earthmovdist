// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; -*-

// emdL1_R.cc:

#include <Rcpp.h>		

#include "emdL1.h"

// [[Rcpp::export]]
Rcpp::NumericVector emdL1_impl(Rcpp::NumericVector h1, 
                               Rcpp::NumericVector h2,
                               Rcpp::List rparam) {

    bool verbose = Rcpp::as<bool>(rparam["verbose"]);
	
    EmdL1 em;
	
    double d = -1;
    switch( Rcpp::as<int>(rparam["noDims"]) ) {
    case 1:
        d = em.EmdDist(h1.begin(), 
                       h2.begin(), 
                       Rcpp::as<int>(rparam["dim1"]));
        if (verbose) Rprintf("1D - EmdL1(h1,h2)=%lf\n", d);
        break;
    case 2 :
        d = em.EmdDist(h1.begin(), 
                       h2.begin(), 
                       Rcpp::as<int>(rparam["dim1"]),
                       Rcpp::as<int>(rparam["dim2"]));
        if (verbose) Rprintf("2D - EmdL1(h1,h2)=%lf\n", d);
        break;
    case 3 :
        d = em.EmdDist(h1.begin(), 
                       h2.begin(), 
                       Rcpp::as<int>(rparam["dim1"]),
                       Rcpp::as<int>(rparam["dim2"]),
                       Rcpp::as<int>(rparam["dim3"]));
        if (verbose) Rprintf("3D - EmdL1(h1,h2)=%lf\n", d);
        break;
    }

    return Rcpp::NumericVector::create(Rcpp::Named("dist", d));
 
}

