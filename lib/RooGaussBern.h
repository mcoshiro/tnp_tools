#ifndef ROO_GAUSS_BERN
#define ROO_GAUSS_BERN

#include "RooAbsPdf.h"
#include "RooAbsArg.h"
#include "RooAbsBinning.h"
#include "RooArgList.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooListProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "TMath.h"
#include "Riostream.h"
#include <vector>

class RooGaussBern : public RooAbsPdf {
public:
  RooGaussBern() {} ; 
  RooGaussBern(const char *name, const char *title,
		    RooAbsReal& x,
		    RooAbsReal& x0,
		    RooAbsReal& sigma,
		    RooAbsReal& alphaL,
		    RooAbsReal& alphaR,
        RooArgList& bernCoefsL,
        RooArgList& bernCoefsR);

  RooGaussBern (const RooGaussBern& other, const char* name);
  inline virtual TObject* clone(const char* newname) const { return new RooGaussBern(*this,newname);}
  inline ~RooGaussBern(){}
  Double_t evaluate() const ;
  
  ClassDef(RooGaussBern, 2)

protected:

  RooRealProxy x;
  RooRealProxy x0;
	RooRealProxy sigma;
	RooRealProxy alphaL;
	RooRealProxy alphaR;
  RooListProxy bernCoefsL;
  RooListProxy bernCoefsR;
    
};

#endif
