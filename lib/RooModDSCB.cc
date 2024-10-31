#include "RooModDSCB.h" 

ClassImp(RooModDSCB) 

RooModDSCB::RooModDSCB(const char *name, const char *title,
	    RooAbsReal& _x,
	    RooAbsReal& _x0,
	    RooAbsReal& _sigmaL,
	    RooAbsReal& _sigmaR,
	    RooAbsReal& _alphaL,
	    RooAbsReal& _nL1,
	    RooAbsReal& _nL2,
	    RooAbsReal& _fL,
	    RooAbsReal& _alphaR,
	    RooAbsReal& _nR1,
	    RooAbsReal& _nR2,
	    RooAbsReal& _fR) :
  RooAbsPdf(name,title), 
  x("x","x",this,_x),
  x0("x0","x0",this,_x0),
  sigmaL("sigmaL","sigmaL",this,_sigmaL),
  sigmaR("sigmaR","sigmaR",this,_sigmaR),
  alphaL("alphaL","alphaL",this,_alphaL),
  nL1("nL1","nL1",this,_nL1),
  nL2("nL2","nL2",this,_nL2),
  fL("fL","fL",this,_fL),
  alphaR("alphaR","alphaR",this,_alphaR),
  nR1("nR1","nR1",this,_nR1),
  nR2("nR2","nR2",this,_nR2),
  fR("fR","fR",this,_fR)
{}

RooModDSCB::RooModDSCB(const RooModDSCB& other, const char* name):
  RooAbsPdf(other,name), 
  x("x",this,other.x),
  x0("x0",this,other.x0),
  sigmaL("sigmaL",this,other.sigmaL),
  sigmaR("sigmaR",this,other.sigmaR),
  alphaL("alphaL",this,other.alphaL),
  nL1("nL1",this,other.nL1),
  nL2("nL2",this,other.nL2),
  fL("fL",this,other.fL),
  alphaR("alphaR",this,other.alphaR),
  nR1("nR1",this,other.nR1),
  nR2("nR2",this,other.nR2),
  fR("fR",this,other.fR)
{}


Double_t RooModDSCB::evaluate() const
{ 
  Double_t left_sigma = (x-x0)/sigmaL;
  Double_t right_sigma = (x-x0)/sigmaR;
  Double_t rval = 0;

  if (left_sigma < -1.0*alphaL) {
    //left power law region
    //force nL2>=nL1 to avoid multiple minima
    Double_t nL2_eff = nL2;
    if (nL2 < nL1)
      nL2_eff = nL1;
    Double_t abs_alpha = fabs(alphaL);
    Double_t BL1 = nL1/abs_alpha-abs_alpha;
    Double_t BL2 = nL2_eff/abs_alpha-abs_alpha;
    Double_t AL1 = TMath::Power(nL1/abs_alpha,nL1)
                   *exp(-0.5*abs_alpha*abs_alpha);
    Double_t AL2 = TMath::Power(nL2_eff/abs_alpha,nL2_eff)
                   *exp(-0.5*abs_alpha*abs_alpha);
    Double_t power_law1 = AL1*TMath::Power((BL1-left_sigma),-1.0*nL1);
    Double_t power_law2 = AL2*TMath::Power((BL2-left_sigma),-1.0*nL2_eff);
    rval = fL*power_law1+(1.0-fL)*power_law2;
  }
  else if (left_sigma < 0) {
    //left Gaussian region
    Double_t exp_arg = (x-x0)/sigmaL;
    rval = exp(-0.5*exp_arg*exp_arg);
  }
  else if (right_sigma < alphaR) {
    //right Gaussian region
    Double_t exp_arg = (x-x0)/sigmaR;
    rval = exp(-0.5*exp_arg*exp_arg);
  }
  else {
    //right power law region
    //force nR2>=nR1 to avoid multiple minima
    Double_t nR2_eff = nR2;
    if (nR2 < nR1)
      nR2_eff = nR1;
    Double_t abs_alpha = fabs(alphaR);
    Double_t BR1 = nR1/abs_alpha-abs_alpha;
    Double_t BR2 = nR2_eff/abs_alpha-abs_alpha;
    Double_t AR1 = TMath::Power(nR1/abs_alpha,nR1)
                   *exp(-0.5*abs_alpha*abs_alpha);
    Double_t AR2 = TMath::Power(nR2_eff/abs_alpha,nR2_eff)
                   *exp(-0.5*abs_alpha*abs_alpha);
    Double_t power_law1 = AR1*TMath::Power((BR1+right_sigma),-1.0*nR1);
    Double_t power_law2 = AR2*TMath::Power((BR2+right_sigma),-1.0*nR2_eff);
    rval = fR*power_law1+(1.0-fR)*power_law2;
  }

  return rval;
} 
