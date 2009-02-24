// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; -*-

// emdL1_R.cc:

#include "Rcpp.h"		

#include "emdL1.h"

RcppExport SEXP emdL1(SEXP H1, SEXP H2, SEXP parms) {

    SEXP rl=R_NilValue;
    char* exceptionMesg=NULL;
  
    try {

        RcppVector<double> h1(H1);	// double vector based on H1
	RcppVector<double> h2(H2);	// double vector based on H2
	RcppParams rparam(parms);       // parameter from R based on parms
	bool verbose = rparam.getBoolValue("verbose");
	RcppResultSet rs;
	
	EmdL1 em;
	
	double d = -1;
	switch( rparam.getIntValue("noDims") ) 
	    {
	    case 1:
		d = em.EmdDist( 
			       h1.cVector(), 
			       h2.cVector(), 
			       rparam.getIntValue("dim1") 
				);
		rs.add("dist", d);
		if (verbose) 
		    Rprintf("1D - EmdL1(h1,h2)=%lf\n", d);
		break;
	    case 2 :
 		d = em.EmdDist( 
			       h1.cVector(), 
			       h2.cVector(), 
			       rparam.getIntValue("dim1"), 
			       rparam.getIntValue("dim2") 
				); 
		rs.add("dist", d);
		if (verbose) 
		    Rprintf("2D - EmdL1(h1,h2)=%lf\n", d);
		break;
	    case 3 :
 		d = em.EmdDist( 
			       h1.cVector(), 
			       h2.cVector(), 
			       rparam.getIntValue("dim1"), 
			       rparam.getIntValue("dim2"), 
			       rparam.getIntValue("dim3") 
				);
		rs.add("dist", d);
		if (verbose) 
		    Rprintf("3D - EmdL1(h1,h2)=%lf\n", d);
		break;
	    }


        rl = rs.getReturnList();

    } catch(std::exception& ex) {
        exceptionMesg = copyMessageToR(ex.what());
    } catch(...) {
        exceptionMesg = copyMessageToR("unknown reason");
    }
  
    if(exceptionMesg != NULL)
        error(exceptionMesg);
    
    return rl;
}

