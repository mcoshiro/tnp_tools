#ifndef ROO_MOD_DSCB
#define ROO_MOD_DSCB


#include "RooAbsPdf.h"
#include "RooAbsArg.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "TMath.h"
#include "Riostream.h"

class RooModDSCB : public RooAbsPdf {
public:
  RooModDSCB() {} ; 
  RooModDSCB(const char *name, const char *title,
		    RooAbsReal& x,
		    RooAbsReal& x0,
		    RooAbsReal& sigmaL,
		    RooAbsReal& sigmaR,
		    RooAbsReal& alphaL,
		    RooAbsReal& nL1,
		    RooAbsReal& nL2,
		    RooAbsReal& fL,
		    RooAbsReal& alphaR,
		    RooAbsReal& nR1,
		    RooAbsReal& nR2,
		    RooAbsReal& fR);

  RooModDSCB (const RooModDSCB& other, const char* name);
  inline virtual TObject* clone(const char* newname) const { return new RooModDSCB(*this,newname);}
  inline ~RooModDSCB(){}
  Double_t evaluate() const ;
  
  ClassDef(RooModDSCB, 2)

protected:

  RooRealProxy x;
  RooRealProxy x0;
	RooRealProxy sigmaL;
	RooRealProxy sigmaR;
	RooRealProxy alphaL;
	RooRealProxy nL1;
	RooRealProxy nL2;
	RooRealProxy fL;
	RooRealProxy alphaR;
	RooRealProxy nR1;
	RooRealProxy nR2;
	RooRealProxy fR;
    
};

#endif
