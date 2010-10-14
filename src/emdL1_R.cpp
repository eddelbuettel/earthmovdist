// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; -*-

// emdL1_R.cc:

#include <Rcpp.h>		

#include "emdL1.h"

RcppExport SEXP emdL1(SEXP H1, SEXP H2, SEXP parms) {

    try {

	Rcpp::NumericVector h1(H1);	// double vector based on H1
	Rcpp::NumericVector h2(H2);	// double vector based on H2
	Rcpp::List rparam(parms);       // parameter from R based on parms
	bool verbose = Rcpp::as<bool>(rparam["verbose"]);
	
	EmdL1 em;
	
	double d = -1;
	switch( Rcpp::as<int>(rparam["noDims"]) ) {
	    case 1:
		d = em.EmdDist(h1.begin(), 
			       h2.begin(), 
			       Rcpp::as<int>(rparam["dim1"]));
		if (verbose) 
		    Rprintf("1D - EmdL1(h1,h2)=%lf\n", d);
		break;
	    case 2 :
 		d = em.EmdDist(h1.begin(), 
			       h2.begin(), 
			       Rcpp::as<int>(rparam["dim1"]),
			       Rcpp::as<int>(rparam["dim2"]));
		if (verbose) 
		    Rprintf("2D - EmdL1(h1,h2)=%lf\n", d);
		break;
	    case 3 :
 		d = em.EmdDist(h1.begin(), 
			       h2.begin(), 
			       Rcpp::as<int>(rparam["dim1"]),
			       Rcpp::as<int>(rparam["dim2"]),
			       Rcpp::as<int>(rparam["dim3"]));
		if (verbose) 
		    Rprintf("3D - EmdL1(h1,h2)=%lf\n", d);
		break;
	}

	//return Rcpp::NumericVector::create(Rcpp::Named("dist", Rcpp::wrap(d)));
	return Rcpp::wrap(d);

    } catch(std::exception &ex) { 
        forward_exception_to_r(ex); 
    } catch(...) { 
        ::Rf_error("c++ exception (unknown reason)"); 
    }

    return R_NilValue;
  
}

