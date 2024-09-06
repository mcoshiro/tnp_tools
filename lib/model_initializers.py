"""@package docstring
Package containing model initializers that set up fitting models 
(signal+background fit) for tnp_analyzer. See tnp_analyzer.py for more info.
"""

import ROOT

def model_initializer_dscb_p_cms(fit_var):
  '''Model initializer that returns a tnp workspace where the signal model is
  a double sided crystal ball and the background shape is a CMSshape (erf*exp)
  '''

  #set up RooFit stuff
  workspace = ROOT.RooWorkspace()
  getattr(workspace,'import')(fit_var)

  gauss_mu = ROOT.RooRealVar('gauss_mu', 'Z peak Gaussian mean', 85.0, 95.0) 
  gauss_sigma = ROOT.RooRealVar('gauss_sigma', 'Z peak Gaussian width', 0.01, 15.0) 
  cb_alphal = ROOT.RooRealVar('cb_alphal', 'Z peak CB left switchover', 0.1, 10.0) 
  cb_nl = ROOT.RooRealVar('cb_nl', 'Z peak CB left power', 0.1, 10.0) 
  cb_alphar = ROOT.RooRealVar('cb_alphar', 'Z peak CB right switchover', 0.1, 10.0) 
  cb_nr = ROOT.RooRealVar('cb_nr', 'Z peak CB right power', 0.1, 10.0) 
  erf_mu = ROOT.RooRealVar('erf_mu', 'Nonresonant (erf) turn-on midpoint', 30.0, 85.0)
  erf_sigma = ROOT.RooRealVar('erf_sigma', 'Nonresonant (erf) turn-on width', 0.001, 2.0) 
  exp_lambda = ROOT.RooRealVar('exp_lambda', 'Nonresonant exponential parameter', 0.0001, 10.0) 
  nSig = ROOT.RooRealVar('nSig', 'Z peak normalization', 0.0, 100000.0) 
  nBkg = ROOT.RooRealVar('nBkg', 'Nonresonant normalization', 0.0, 100000.0) 
  getattr(workspace,'import')(gauss_mu)
  getattr(workspace,'import')(gauss_sigma)
  getattr(workspace,'import')(cb_alphal)
  getattr(workspace,'import')(cb_nl)
  getattr(workspace,'import')(cb_alphar)
  getattr(workspace,'import')(cb_nr)
  getattr(workspace,'import')(erf_mu)
  getattr(workspace,'import')(erf_sigma)
  getattr(workspace,'import')(exp_lambda)
  getattr(workspace,'import')(nSig)
  getattr(workspace,'import')(nBkg)

  pdf_s  = ROOT.RooCrystalBall('pdf_s','pdf_s', fit_var, gauss_mu, gauss_sigma, cb_alphal, cb_nl, cb_alphar, cb_nr)
  pdf_b = ROOT.RooGenericPdf('pdf_b','pdf_b',
      '(TMath::Erf((@0-@1)*@2)+1.0)/2.0*exp(-1.0*@3*(@0-60.0)/40.0)',
         ROOT.RooArgList(fit_var, erf_mu, erf_sigma, exp_lambda))
  pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(nSig, nBkg))
  getattr(workspace,'import')(pdf_sb)
  return workspace
