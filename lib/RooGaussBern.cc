#include "RooGaussBern.h" 

ClassImp(RooGaussBern) 

//Constructor
//note the order of the polynomial is equal to the number of coefficients given
//as the one additional coefficient is fixed by normalization. Note other
//coefficients are relative to fixed one
RooGaussBern::RooGaussBern(const char *name, const char *title,
	    RooAbsReal& _x,
	    RooAbsReal& _x0,
	    RooAbsReal& _sigma,
	    RooAbsReal& _alphaL,
	    RooAbsReal& _alphaR,
      RooArgList& _bernCoefsL,
      RooArgList& _bernCoefsR) :
  RooAbsPdf(name,title), 
  x("x","x",this,_x),
  x0("x0","x0",this,_x0),
  sigma("sigma","sigma",this,_sigma),
  alphaL("alphaL","alphaL",this,_alphaL),
  alphaR("alphaR","alphaR",this,_alphaR),
  bernCoefsL("bernCoefsL","",this),
  bernCoefsR("bernCoefsR","",this)
{
  for (auto *coef : _bernCoefsL) {
    bernCoefsL.add(*coef);
  }
  for (auto *coef : _bernCoefsR) {
    bernCoefsR.add(*coef);
  }
}

RooGaussBern::RooGaussBern(const RooGaussBern& other, const char* name):
  RooAbsPdf(other,name), 
  x("x",this,other.x),
  x0("x0",this,other.x0),
  sigma("sigma",this,other.sigma),
  alphaL("alphaL",this,other.alphaL),
  alphaR("alphaR",this,other.alphaR),
  bernCoefsL("bernCoefsL",this,other.bernCoefsL),
  bernCoefsR("bernCoefsR",this,other.bernCoefsR)
{}


Double_t RooGaussBern::evaluate() const
{ 

  Double_t low_edge = x.min();
  Double_t high_edge = x.max();
  Double_t norm_dist = (x-x0)/sigma;
  Double_t rval = 0.0;

  if (norm_dist < -1.0*alphaL) {
    //left polynomial region
    Double_t boundary = x0-alphaL*sigma;
    Double_t x_scaled = (x-low_edge)/(boundary-low_edge);
    Int_t bern_order = bernCoefsL.getSize();
    Double_t coef0 = exp(-0.5*alphaL*alphaL);
    for (Int_t i = 0; i <= bern_order; i++) {
      Double_t powone = static_cast<double>(i);
      Double_t powtwo = static_cast<double>(bern_order-i);
      Double_t term_coef = coef0;
      if (i < bern_order) {
        term_coef *= static_cast<RooAbsReal&>(bernCoefsL[i]).getVal();
      }
      rval += term_coef
              *TMath::Power(x_scaled, powone)
              *TMath::Power(1.0-x_scaled, powtwo);
    }
  }
  else if (norm_dist < alphaR) {
    //Gaussian region
    Double_t exp_arg = (x-x0)/sigma;
    rval = exp(-0.5*exp_arg*exp_arg);
  }
  else {
    //right polynomial region
    Double_t boundary = x0+alphaR*sigma;
    Double_t x_scaled = (x-boundary)/(high_edge-boundary);
    Int_t bern_order = bernCoefsR.getSize();
    Double_t coef0 = exp(-0.5*alphaR*alphaR);
    for (Int_t i = 0; i <= bern_order; i++) {
      Double_t powone = static_cast<double>(i);
      Double_t powtwo = static_cast<double>(bern_order-i);
      Double_t term_coef = coef0;
      if (i > 0) {
        term_coef *= static_cast<RooAbsReal&>(bernCoefsR[i-1]).getVal();
      }
      rval += term_coef
              *TMath::Power(x_scaled, powone)
              *TMath::Power(1.0-x_scaled, powtwo);
    }
  }

  return rval;
} 
