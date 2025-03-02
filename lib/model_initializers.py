"""@package docstring
Package containing model initializers that set up fitting models 
(signal+background fit) for tnp_analyzer. See tnp_analyzer.py for more info.
"""

import ROOT
import json

MAX_SIGNAL = 100000000.0
MAX_BACKGROUND = 100000000.0

ROOT.gInterpreter.ProcessLine('.L lib/RooCBExGaussShapeTNP.cc+')
ROOT.gInterpreter.ProcessLine('.L lib/RooModDSCB.cc+')
ROOT.gInterpreter.ProcessLine('.L lib/RooGaussBern.cc+')

def model_initializer_dscb_p_cms(fit_var, ibin, is_pass):
  '''Model initializer that returns a tnp workspace where the signal model is
  a double sided crystal ball and the background shape is a CMSshape (erf*exp)

  fit_var  RooRealVar representing variable to fit
  ibin     int, bin number
  is_pass  bool, if is passing leg
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
  nSig = ROOT.RooRealVar('nSig', 'Z peak normalization', 0.0, MAX_SIGNAL) 
  nBkg = ROOT.RooRealVar('nBkg', 'Nonresonant normalization', 0.0, MAX_BACKGROUND) 

  gauss_sigma.setVal(3.0)
  cb_nl.setVal(3.0)
  cb_nr.setVal(3.0)
  cb_alphal.setVal(3.0)
  cb_alphar.setVal(3.0)
  nBkg.setVal(5000.0)

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

def model_initializer_dscb(fit_var, ibin, is_pass):
  '''Model initializer that reutrns a tnp workspace where the signal model
  is a DSCB shape and there is no background

  fit_var      RooRealVar representing variable to fit
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  '''
  workspace = ROOT.RooWorkspace()
  getattr(workspace,'import')(fit_var)

  gauss_mu = ROOT.RooRealVar('mean', 'Gaussian mean', 85.0, 95.0) 
  sigmal = ROOT.RooRealVar('sigmal', 'Gaussian sigma', 0.01, 15.0) 
  sigmar = ROOT.RooRealVar('sigmar', 'Gaussian sigma', 0.01, 15.0) 
  cb_alphal = ROOT.RooRealVar('alphal', 'CB left switchover', 0.1, 10.0) 
  cb_nl = ROOT.RooRealVar('nl', 'CB left power', 0.1, 10.0) 
  cb_alphar = ROOT.RooRealVar('alphar', 'CB right switchover', 0.1, 10.0) 
  cb_nr = ROOT.RooRealVar('nr', 'CB right power', 0.1, 10.0) 
  sigmal.setVal(3.0)
  sigmar.setVal(3.0)
  cb_alphal.setVal(2.0)
  cb_alphar.setVal(2.0)
  cb_nl.setVal(2.0)
  cb_nr.setVal(2.0)

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, MAX_SIGNAL) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, 0.0) 

  getattr(workspace,'import')(gauss_mu)
  getattr(workspace,'import')(sigmal)
  getattr(workspace,'import')(sigmar)
  getattr(workspace,'import')(cb_alphal)
  getattr(workspace,'import')(cb_nl)
  getattr(workspace,'import')(cb_alphar)
  getattr(workspace,'import')(cb_nr)
  getattr(workspace,'import')(nSig)
  getattr(workspace,'import')(nBkg)

  pdf_sb  = ROOT.RooCrystalBall('pdf_sb','pdf_sb', fit_var, gauss_mu, sigmal, 
                                sigmar, cb_alphal, cb_nl, cb_alphar, cb_nr)
  getattr(workspace,'import')(pdf_sb)
  return workspace


def model_initializer_dscbgaus(fit_var, ibin, is_pass):
  '''Model initializer that reutrns a tnp workspace where the signal model
  is a DSCB shape plus a Gaussian and there is no background

  fit_var      RooRealVar representing variable to fit
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  '''
  workspace = ROOT.RooWorkspace()
  getattr(workspace,'import')(fit_var)

  mu = ROOT.RooRealVar('mu', 'Gaussian mean', 85.0, 95.0) 
  sigma = ROOT.RooRealVar('sigma', 'Gaussian sigma', 0.5, 10.0) 
  cb_alphal = ROOT.RooRealVar('alphal', 'CB left switchover', 0.1, 10.0) 
  cb_nl = ROOT.RooRealVar('nl', 'CB left power', 0.1, 10.0) 
  cb_alphar = ROOT.RooRealVar('alphar', 'CB right switchover', 0.1, 10.0) 
  cb_nr = ROOT.RooRealVar('nr', 'CB right power', 0.1, 10.0) 

  gauss_mu = ROOT.RooRealVar('gauss_mu', 'Gaussian mean', 60.0, 120.0) 
  gauss_sigma = ROOT.RooRealVar('gauss_sigma', 'Gaussian sigma', 0.5, 30.0) 
  gauss_frac = ROOT.RooRealVar('gauss_frac', 'Gaussian fraction', 0.0, 1.0)

  sigma.setVal(3.0)
  cb_alphal.setVal(2.0)
  cb_alphar.setVal(2.0)
  cb_nl.setVal(2.0)
  cb_nr.setVal(2.0)
  gauss_mu.setVal(75.0)
  gauss_sigma.setVal(5.0)
  gauss_frac.setVal(0.0)

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, MAX_SIGNAL) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, 0.0) 

  getattr(workspace,'import')(mu)
  getattr(workspace,'import')(sigma)
  getattr(workspace,'import')(cb_alphal)
  getattr(workspace,'import')(cb_nl)
  getattr(workspace,'import')(cb_alphar)
  getattr(workspace,'import')(cb_nr)
  getattr(workspace,'import')(gauss_mu)
  getattr(workspace,'import')(gauss_sigma)
  getattr(workspace,'import')(gauss_frac)
  getattr(workspace,'import')(nSig)
  getattr(workspace,'import')(nBkg)

  pdf_dscb = ROOT.RooCrystalBall('pdf_dscb','pdf_dscb', fit_var, mu, sigma, 
                                 cb_alphal, cb_nl, cb_alphar, cb_nr)
  pdf_gaus = ROOT.RooGaussian('pdf_gaus','pdf_gaus', fit_var, gauss_mu,
                              gauss_sigma)
  pdf_sb = ROOT.RooAddPdf('pdf_sb','pdf_sb',
                          ROOT.RooArgList(pdf_gaus, pdf_dscb), 
                          ROOT.RooArgList(gauss_frac))
  getattr(workspace,'import')(pdf_sb)
  return workspace

def model_initializer_cbconvgen(fit_var, ibin, is_pass):
  '''Adds CBExGaussShapeTNP convoluted with Z lineshape as signal model and no
  background

  fit_var      RooRealVar representing variable to fit
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  '''
  workspace = ROOT.RooWorkspace()
  getattr(workspace,'import')(fit_var)

  m0 = ROOT.RooRealVar('m0', 'm0', -5.0, 5.0) 
  sigma = ROOT.RooRealVar('sigma', 'sigma', 0.7, 15.0) 
  alpha = ROOT.RooRealVar('alpha', 'alpha', 0.5, 10.0) 
  n = ROOT.RooRealVar('n', 'n', -5.0, 5.0) 
  sigma_2 = ROOT.RooRealVar('sigma_2', 'sigma_2', 0.5, 6.0) 
  tailLeft = ROOT.RooRealVar('tailLeft', 'tailLeft', 0.5, 5.0) 
  sigma.setVal(2.0)
  alpha.setVal(2.0)
  n.setVal(3.0)
  sigma_2.setVal(2.0)
  tailLeft.setVal(1.0)
  getattr(workspace,'import')(m0)
  getattr(workspace,'import')(sigma)
  getattr(workspace,'import')(alpha)
  getattr(workspace,'import')(n)
  getattr(workspace,'import')(sigma_2)
  getattr(workspace,'import')(tailLeft)

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, MAX_SIGNAL) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, 0.0) 
  getattr(workspace,'import')(nSig)
  getattr(workspace,'import')(nBkg)

  input_file = ROOT.TFile('lib/ZeeGenLevel.root','READ')
  hist = input_file.Get('Mass')
  hist.SetDirectory(ROOT.nullptr)
  input_file.Close()
  core_hist = ROOT.RooDataHist('pdf_s_core_hist',
                               'pdf_s_core_hist',
                               ROOT.RooArgList(fit_var), hist)
  getattr(workspace,'import')(core_hist)
  pdf_s_core = ROOT.RooHistPdf('pdf_s_core','pdf_s_core',
                               ROOT.RooArgSet(fit_var),core_hist)
  pdf_s_res = ROOT.RooCBExGaussShapeTNP('pdf_s_res','pdf_s_res',fit_var,
                                        m0,sigma,
                                        alpha,n,sigma_2,tailLeft)
  pdf_s = ROOT.RooFFTConvPdf('pdf_sb','pdf_sb',fit_var,pdf_s_core,pdf_s_res)
  getattr(workspace,'import')(pdf_s)
  return workspace

def model_initializer_moddscb(fit_var, ibin, is_pass):
  '''Adds modified double-sided crystal ball distribution as signal model
  with no background

  fit_var      RooRealVar representing variable to fit
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  '''
  workspace = ROOT.RooWorkspace()
  getattr(workspace,'import')(fit_var)

  gauss_mu = ROOT.RooRealVar('mean', 'Gaussian mean', 85.0, 95.0) 
  sigmal = ROOT.RooRealVar('sigmal', 'Gaussian sigma', 0.01, 15.0) 
  sigmar = ROOT.RooRealVar('sigmar', 'Gaussian sigma', 0.01, 15.0) 
  cb_alphal = ROOT.RooRealVar('alphal', 'CB left switchover', 0.1, 10.0) 
  cb_nl1 = ROOT.RooRealVar('nl1', 'CB left power', 0.1, 10.0) 
  cb_nl2 = ROOT.RooRealVar('nl2', 'CB left power', 0.1, 10.0) 
  cb_fl = ROOT.RooRealVar('fl', 'left power law fraction', 0.0, 1.0) 
  cb_alphar = ROOT.RooRealVar('alphar', 'CB right switchover', 0.1, 10.0) 
  cb_nr1 = ROOT.RooRealVar('nr1', 'CB right power', 0.1, 10.0) 
  cb_nr2 = ROOT.RooRealVar('nr2', 'CB right power', 0.1, 10.0) 
  cb_fr = ROOT.RooRealVar('fr', 'right power law fraction', 0.0, 1.0) 

  sigmal.setVal(3.0)
  sigmar.setVal(3.0)
  cb_alphal.setVal(2.0)
  cb_alphar.setVal(2.0)
  cb_nl1.setVal(1.0)
  cb_nl2.setVal(2.0)
  cb_nr1.setVal(1.0)
  cb_nr2.setVal(2.0)
  getattr(workspace,'import')(gauss_mu)
  getattr(workspace,'import')(sigmal)
  getattr(workspace,'import')(sigmar)
  getattr(workspace,'import')(cb_alphal)
  getattr(workspace,'import')(cb_nl1)
  getattr(workspace,'import')(cb_nl2)
  getattr(workspace,'import')(cb_fl)
  getattr(workspace,'import')(cb_alphar)
  getattr(workspace,'import')(cb_nr1)
  getattr(workspace,'import')(cb_nr2)
  getattr(workspace,'import')(cb_fr)

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, MAX_SIGNAL) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, 0.0) 
  getattr(workspace,'import')(nSig)
  getattr(workspace,'import')(nBkg)
  
  pdf_s  = ROOT.RooModDSCB('pdf_sb','pdf_sb', fit_var, gauss_mu, sigmal, 
                           sigmar, cb_alphal, cb_nl1, cb_nl2, cb_fl, cb_alphar,
                           cb_nr1, cb_nr2, cb_fr)
  getattr(workspace,'import')(pdf_s)
  return workspace

def model_initializer_gaussbern(fit_var, ibin, is_pass):
  '''Adds a piecewise gaussian-polynomial distributino as the signal model and
  no background

  fit_var      RooRealVar representing variable to fit
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  '''
  workspace = ROOT.RooWorkspace()
  getattr(workspace,'import')(fit_var)

  mu = ROOT.RooRealVar('mean', 'Gaussian mean', 85.0, 95.0) 
  sigma = ROOT.RooRealVar('sigma', 'Gaussian sigma', 0.01, 15.0) 
  cb_alphal = ROOT.RooRealVar('alphal', 'CB left switchover', 0.1, 10.0) 
  cb_alphar = ROOT.RooRealVar('alphar', 'CB right switchover', 0.1, 10.0) 
  al1 = ROOT.RooRealVar('al1', 'Bernstein left coef 0', 0.0, 1.0) 
  al2 = ROOT.RooRealVar('al2', 'Bernstein left coef 1', 0.0, 1.0) 
  al3 = ROOT.RooRealVar('al3', 'Bernstein left coef 2', 0.0, 1.0) 
  al4 = ROOT.RooRealVar('al4', 'Bernstein left coef 3', 0.0, 1.0) 
  ar1 = ROOT.RooRealVar('ar1', 'Bernstein right coef 0', 0.0, 1.0) 
  ar2 = ROOT.RooRealVar('ar2', 'Bernstein right coef 1', 0.0, 1.0) 
  ar3 = ROOT.RooRealVar('ar3', 'Bernstein right coef 2', 0.0, 1.0) 
  ar4 = ROOT.RooRealVar('ar4', 'Bernstein right coef 3', 0.0, 1.0) 

  sigma.setVal(3.0)
  cb_alphal.setVal(2.0)
  cb_alphar.setVal(2.0)
  getattr(workspace,'import')(mu)
  getattr(workspace,'import')(sigma)
  getattr(workspace,'import')(cb_alphal)
  getattr(workspace,'import')(cb_alphar)
  getattr(workspace,'import')(al1)
  getattr(workspace,'import')(al2)
  getattr(workspace,'import')(al3)
  getattr(workspace,'import')(al4)
  getattr(workspace,'import')(ar1)
  getattr(workspace,'import')(ar2)
  getattr(workspace,'import')(ar3)
  getattr(workspace,'import')(ar4)

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, MAX_SIGNAL) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, 0.0) 
  getattr(workspace,'import')(nSig)
  getattr(workspace,'import')(nBkg)
  
  pdf_s  = ROOT.RooGaussBern('pdf_sb','pdf_sb', fit_var, mu, sigma, cb_alphal,
                             cb_alphar, ROOT.RooArgList(al1,al2,al3,al4),
                             ROOT.RooArgList(ar1,ar2,ar3,ar4))
  getattr(workspace,'import')(pdf_s)
  return workspace

def make_signal_background_model(fit_var, ibin, is_pass, add_signal_model, 
                                 add_background_model):
  '''Model initializer that reutrns a tnp workspace with given signal and background 
  models

  fit_var              RooRealVar representing variable to fit
  ibin                 int, bin number
  is_pass              bool, if is passing leg
  add_signal_model     function that adds signal pdf (pdf_s) to workspace
  add_background_model function that adds background pdf (pdf_b) to workspace
  '''
  workspace = ROOT.RooWorkspace()
  getattr(workspace,'import')(fit_var)

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, MAX_SIGNAL) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, MAX_BACKGROUND) 
  getattr(workspace,'import')(nSig)
  getattr(workspace,'import')(nBkg)
  add_signal_model(workspace, ibin, is_pass)
  add_background_model(workspace, ibin, is_pass)

  pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', 
      ROOT.RooArgList(workspace.pdf('pdf_s'), workspace.pdf('pdf_b')), 
      ROOT.RooArgList(nSig, nBkg))
  getattr(workspace,'import')(pdf_sb)
  return workspace

def add_signal_model_mcsmear(workspace, ibin, is_pass, get_histogram):
  '''Adds signal model that is template smeared by a Gaussian

  workspace      workspace to add model to; must have variable fit_var
  ibin           int, bin number
  is_pass        bool, if is passing leg
  get_histogram  function that returns TH1D given ibin and is_pass
  '''

  fit_var = workspace.var('fit_var')
  mean = ROOT.RooRealVar('mean', 'Smearing gaussian peak', -10.0, 10.0) 
  sigma = ROOT.RooRealVar('sigma', 'Smearing gaussian sigma', 0.5, 5.0) 
  getattr(workspace,'import')(mean)
  getattr(workspace,'import')(sigma)

  hist = get_histogram(ibin, is_pass)
  core_hist = ROOT.RooDataHist('pdf_s_core_hist',
                               'pdf_s_core_hist',
                               ROOT.RooArgList(fit_var), hist)
  getattr(workspace,'import')(core_hist)
  pdf_s_core = ROOT.RooHistPdf('pdf_s_core','pdf_s_core',
                               ROOT.RooArgSet(fit_var),core_hist)
  pdf_s_res = ROOT.RooGaussian('pdf_s_res','pdf_s_res',fit_var,mean,sigma)
  pdf_s = ROOT.RooFFTConvPdf('pdf_s','pdf_s',fit_var,pdf_s_core,pdf_s_res)
  getattr(workspace,'import')(pdf_s)

def add_signal_model_mcsumsmear(workspace, ibin, is_pass, get_histogram):
  '''Adds signal model that is template formed from sum of passing and
  failing MC smeared by a Gaussian

  workspace      workspace to add model to; must have variable fit_var
  ibin           int, bin number
  is_pass        bool, if is passing leg
  get_histogram  function that returns TH1D given ibin and is_pass
  '''

  fit_var = workspace.var('fit_var')
  mean = ROOT.RooRealVar('mean', 'Smearing gaussian peak', -10.0, 10.0) 
  sigma = ROOT.RooRealVar('sigma', 'Smearing gaussian sigma', 0.5, 5.0) 
  pass_frac = ROOT.RooRealVar('pass_frac', 'pass template coef', 0.0, 1.0)
  if (is_pass):
    pass_frac.setVal(1.0)
  else:
    pass_frac.setVal(0.0)
  getattr(workspace,'import')(mean)
  getattr(workspace,'import')(sigma)
  getattr(workspace,'import')(pass_frac)

  pass_hist = get_histogram(ibin, True)
  fail_hist = get_histogram(ibin, False)
  core_hist_pass = ROOT.RooDataHist('core_hist_pass',
                                    'core_hist_pass',
                                    ROOT.RooArgList(fit_var), pass_hist)
  core_hist_fail = ROOT.RooDataHist('core_hist_fail',
                                    'core_hist_fail',
                                    ROOT.RooArgList(fit_var), fail_hist)
  getattr(workspace,'import')(core_hist_pass)
  getattr(workspace,'import')(core_hist_fail)
  pdf_s_core_pass = ROOT.RooHistPdf('pdf_s_core_pass','pdf_s_core_pass',
                                    ROOT.RooArgSet(fit_var),core_hist_pass)
  pdf_s_core_fail = ROOT.RooHistPdf('pdf_s_core_fail','pdf_s_core_fail',
                                    ROOT.RooArgSet(fit_var),core_hist_fail)
  pdf_s_core = ROOT.RooAddPdf('pdf_s_core','pdf_s_core',pdf_s_core_pass,
                              pdf_s_core_fail,pass_frac)
  pdf_s_res = ROOT.RooGaussian('pdf_s_res','pdf_s_res',fit_var,mean,sigma)
  pdf_s = ROOT.RooFFTConvPdf('pdf_s','pdf_s',fit_var,pdf_s_core,pdf_s_res)
  getattr(workspace,'import')(pdf_s)

def add_signal_model_mcdscbsmear(workspace, ibin, is_pass, get_histogram):
  '''Adds signal model that is template convoluted with a crystal ball

  workspace      workspace to add model to; must have variable fit_var
  ibin           int, bin number
  is_pass        bool, if is passing leg
  get_histogram  function that returns TH1D given ibin and is_pass
  '''

  fit_var = workspace.var('fit_var')
  mu = ROOT.RooRealVar('mu', 'Smearing gaussian peak', -5.0, 5.0) 
  sigma = ROOT.RooRealVar('sigma', 'Smearing gaussian sigma', 0.5, 5.0) 
  cb_alphal = ROOT.RooRealVar('alphal', 'CB left switchover', 0.1, 10.0) 
  cb_nl = ROOT.RooRealVar('nl', 'CB left power', 0.1, 10.0) 
  cb_alphar = ROOT.RooRealVar('alphar', 'CB right switchover', 0.1, 10.0) 
  cb_nr = ROOT.RooRealVar('nr', 'CB right power', 0.1, 10.0) 
  sigma.setVal(3.0)
  cb_alphal.setVal(2.0)
  cb_alphar.setVal(2.0)
  cb_nl.setVal(2.0)
  cb_nr.setVal(2.0)
  getattr(workspace,'import')(mu)
  getattr(workspace,'import')(sigma)
  getattr(workspace,'import')(cb_alphal)
  getattr(workspace,'import')(cb_nl)
  getattr(workspace,'import')(cb_alphar)
  getattr(workspace,'import')(cb_nr)

  hist = get_histogram(ibin, is_pass)
  core_hist = ROOT.RooDataHist('pdf_s_core_hist',
                               'pdf_s_core_hist',
                               ROOT.RooArgList(fit_var), hist)
  getattr(workspace,'import')(core_hist)
  pdf_s_core = ROOT.RooHistPdf('pdf_s_core','pdf_s_core',
                               ROOT.RooArgSet(fit_var),core_hist)
  pdf_s_res = ROOT.RooCrystalBall('pdf_s_res','pdf_s_res',fit_var,mu,sigma,
                                  cb_alphal, cb_nl, cb_alphar, cb_nr)
  pdf_s = ROOT.RooFFTConvPdf('pdf_s','pdf_s',fit_var,pdf_s_core,pdf_s_res)
  getattr(workspace,'import')(pdf_s)

def add_signal_model_dscb(workspace, ibin, is_pass):
  '''Adds double-sided crystal ball distribution as signal model

  workspace    workspace to add model to; must have variable fit_var
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  '''
  fit_var = workspace.var('fit_var')
  gauss_mu = ROOT.RooRealVar('mean', 'Gaussian mean', 85.0, 95.0) 
  sigmal = ROOT.RooRealVar('sigmal', 'Gaussian sigma', 0.01, 15.0) 
  sigmar = ROOT.RooRealVar('sigmar', 'Gaussian sigma', 0.01, 15.0) 
  cb_alphal = ROOT.RooRealVar('alphal', 'CB left switchover', 0.1, 10.0) 
  cb_nl = ROOT.RooRealVar('nl', 'CB left power', 0.1, 10.0) 
  cb_alphar = ROOT.RooRealVar('alphar', 'CB right switchover', 0.1, 10.0) 
  cb_nr = ROOT.RooRealVar('nr', 'CB right power', 0.1, 10.0) 
  sigmal.setVal(3.0)
  sigmar.setVal(3.0)
  cb_alphal.setVal(2.0)
  cb_alphar.setVal(2.0)
  cb_nl.setVal(2.0)
  cb_nr.setVal(2.0)
  getattr(workspace,'import')(gauss_mu)
  getattr(workspace,'import')(sigmal)
  getattr(workspace,'import')(sigmar)
  getattr(workspace,'import')(cb_alphal)
  getattr(workspace,'import')(cb_nl)
  getattr(workspace,'import')(cb_alphar)
  getattr(workspace,'import')(cb_nr)
  pdf_s  = ROOT.RooCrystalBall('pdf_s','pdf_s', fit_var, gauss_mu, sigmal, 
                               sigmar, cb_alphal, cb_nl, cb_alphar, cb_nr)
  getattr(workspace,'import')(pdf_s)

def add_signal_model_dscbgaus(workspace, ibin, is_pass):
  '''Adds double-sided crystal ball distribution plus Gaussian as signal model

  workspace    workspace to add model to; must have variable fit_var
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  '''
  fit_var = workspace.var('fit_var')
  mu = ROOT.RooRealVar('mu', 'Gaussian mean', 85.0, 95.0) 
  sigma = ROOT.RooRealVar('sigma', 'Gaussian sigma', 0.1, 10.0) 
  cb_alphal = ROOT.RooRealVar('alphal', 'CB left switchover', 0.1, 10.0) 
  cb_nl = ROOT.RooRealVar('nl', 'CB left power', 0.1, 10.0) 
  cb_alphar = ROOT.RooRealVar('alphar', 'CB right switchover', 0.1, 10.0) 
  cb_nr = ROOT.RooRealVar('nr', 'CB right power', 0.1, 10.0) 

  gauss_mu = ROOT.RooRealVar('gauss_mu', 'Gaussian mean', 60.0, 120.0) 
  gauss_sigma = ROOT.RooRealVar('gauss_sigma', 'Gaussian sigma', 0.01, 30.0) 
  gauss_frac = ROOT.RooRealVar('gauss_frac', 'Gaussian fraction', 0.0, 1.0)

  sigma.setVal(3.0)
  cb_alphal.setVal(2.0)
  cb_alphar.setVal(2.0)
  cb_nl.setVal(2.0)
  cb_nr.setVal(2.0)
  gauss_mu.setVal(75.0)
  gauss_sigma.setVal(5.0)
  gauss_frac.setVal(0.0)

  getattr(workspace,'import')(mu)
  getattr(workspace,'import')(sigma)
  getattr(workspace,'import')(cb_alphal)
  getattr(workspace,'import')(cb_nl)
  getattr(workspace,'import')(cb_alphar)
  getattr(workspace,'import')(cb_nr)
  getattr(workspace,'import')(gauss_mu)
  getattr(workspace,'import')(gauss_sigma)
  getattr(workspace,'import')(gauss_frac)

  pdf_dscb = ROOT.RooCrystalBall('pdf_dscb','pdf_dscb', fit_var, mu, sigma, 
                                 cb_alphal, cb_nl, cb_alphar, cb_nr)
  pdf_gaus = ROOT.RooGaussian('pdf_gaus','pdf_gaus', fit_var, gauss_mu,
                              gauss_sigma)
  pdf_s = ROOT.RooAddPdf('pdf_s','pdf_s',
                          ROOT.RooArgList(pdf_gaus, pdf_dscb), 
                          ROOT.RooArgList(gauss_frac))
  getattr(workspace,'import')(pdf_s)

def add_signal_model_moddscb(workspace, ibin, is_pass):
  '''Adds modified double-sided crystal ball distribution as signal model

  workspace    workspace to add model to; must have variable fit_var
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  '''
  fit_var = workspace.var('fit_var')

  gauss_mu = ROOT.RooRealVar('mean', 'Gaussian mean', 85.0, 95.0) 
  sigmal = ROOT.RooRealVar('sigmal', 'Gaussian sigma', 0.01, 15.0) 
  sigmar = ROOT.RooRealVar('sigmar', 'Gaussian sigma', 0.01, 15.0) 
  cb_alphal = ROOT.RooRealVar('alphal', 'CB left switchover', 0.1, 10.0) 
  cb_nl1 = ROOT.RooRealVar('nl1', 'CB left power', 0.1, 10.0) 
  cb_nl2 = ROOT.RooRealVar('nl2', 'CB left power', 0.1, 10.0) 
  cb_fl = ROOT.RooRealVar('fl', 'left power law fraction', 0.0, 1.0) 
  cb_alphar = ROOT.RooRealVar('alphar', 'CB right switchover', 0.1, 10.0) 
  cb_nr1 = ROOT.RooRealVar('nr1', 'CB right power', 0.1, 10.0) 
  cb_nr2 = ROOT.RooRealVar('nr2', 'CB right power', 0.1, 10.0) 
  cb_fr = ROOT.RooRealVar('fr', 'right power law fraction', 0.0, 1.0) 

  sigmal.setVal(3.0)
  sigmar.setVal(3.0)
  cb_alphal.setVal(2.0)
  cb_alphar.setVal(2.0)
  cb_nl1.setVal(1.0)
  cb_nl2.setVal(2.0)
  cb_nr1.setVal(1.0)
  cb_nr2.setVal(2.0)
  getattr(workspace,'import')(gauss_mu)
  getattr(workspace,'import')(sigmal)
  getattr(workspace,'import')(sigmar)
  getattr(workspace,'import')(cb_alphal)
  getattr(workspace,'import')(cb_nl1)
  getattr(workspace,'import')(cb_nl2)
  getattr(workspace,'import')(cb_fl)
  getattr(workspace,'import')(cb_alphar)
  getattr(workspace,'import')(cb_nr1)
  getattr(workspace,'import')(cb_nr2)
  getattr(workspace,'import')(cb_fr)
  pdf_s  = ROOT.RooModDSCB('pdf_s','pdf_s', fit_var, gauss_mu, sigmal, sigmar, 
                           cb_alphal, cb_nl1, cb_nl2, cb_fl, cb_alphar, cb_nr1,
                           cb_nr2, cb_fr)
  getattr(workspace,'import')(pdf_s)

def add_signal_model_cbconvgen(workspace, ibin, is_pass):
  '''Adds CBExGaussShapeTNP convoluted with Z lineshape

  workspace    workspace to add model to; must have variable fit_var
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  '''
  fit_var = workspace.var('fit_var')
  m0 = ROOT.RooRealVar('m0', 'm0', -5.0, 5.0) 
  sigma = ROOT.RooRealVar('sigma', 'sigma', 0.7, 15.0) 
  alpha = ROOT.RooRealVar('alpha', 'alpha', 0.5, 10.0) 
  n = ROOT.RooRealVar('n', 'n', -5.0, 5.0) 
  sigma_2 = ROOT.RooRealVar('sigma_2', 'sigma_2', 0.5, 6.0) 
  tailLeft = ROOT.RooRealVar('tailLeft', 'tailLeft', 0.5, 5.0) 
  sigma.setVal(2.0)
  alpha.setVal(2.0)
  n.setVal(3.0)
  sigma_2.setVal(2.0)
  tailLeft.setVal(1.0)
  getattr(workspace,'import')(m0)
  getattr(workspace,'import')(sigma)
  getattr(workspace,'import')(alpha)
  getattr(workspace,'import')(n)
  getattr(workspace,'import')(sigma_2)
  getattr(workspace,'import')(tailLeft)

  input_file = ROOT.TFile('lib/ZeeGenLevel.root','READ')
  hist = input_file.Get('Mass')
  hist.SetDirectory(ROOT.nullptr)
  input_file.Close()
  core_hist = ROOT.RooDataHist('pdf_s_core_hist',
                               'pdf_s_core_hist',
                               ROOT.RooArgList(fit_var), hist)
  getattr(workspace,'import')(core_hist)
  pdf_s_core = ROOT.RooHistPdf('pdf_s_core','pdf_s_core',
                               ROOT.RooArgSet(fit_var),core_hist)
  pdf_s_res = ROOT.RooCBExGaussShapeTNP('pdf_s_res','pdf_s_res',fit_var,m0,
                                        sigma,
                                        alpha,n,sigma_2,tailLeft)
  pdf_s = ROOT.RooFFTConvPdf('pdf_s','pdf_s',fit_var,pdf_s_core,pdf_s_res)
  getattr(workspace,'import')(pdf_s)

def add_background_model_cmsshape(workspace, ibin, is_pass):
  '''Adds background model that is CMS shape (erf*exp)

  workspace      workspace to add model to; must have variable fit_var
  ibin           int, bin number
  is_pass        bool, if is passing leg
  '''
  fit_var = workspace.var('fit_var')
  acms = ROOT.RooRealVar('acms', 'erf turn-on point', 50.0, 80.0)
  beta = ROOT.RooRealVar('beta', 'erf width', 0.01, 0.08)
  gamma = ROOT.RooRealVar('gamma', 'exp parameter', -2.0, 2.0)
  peak = ROOT.RooRealVar('peak', 'peak', 90.0, 90.0)
  acms.setVal(60.0)
  beta.setVal(0.01)
  gamma.setVal(0.005)
  getattr(workspace,'import')(acms)
  getattr(workspace,'import')(beta)
  getattr(workspace,'import')(gamma)
  getattr(workspace,'import')(peak)
  pdf_b = ROOT.RooGenericPdf('pdf_b','pdf_b',
      'TMath::Erfc((@1-@0)*@2)*exp((@4-@0)*@3)',
      ROOT.RooArgList(fit_var, acms, beta, gamma, peak))
  getattr(workspace,'import')(pdf_b)

def add_background_model_chebyshev(workspace, ibin, is_pass):
  '''Adds Chebyshev polynomial of degree 3 as background model

  workspace    workspace to add model to; must have variable fit_var
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  '''
  fit_var = workspace.var('fit_var')
  a0 = ROOT.RooRealVar('a0', '0th Chebyshev coefficient', -1.0, 1.0)
  a1 = ROOT.RooRealVar('a1', '1st Chebyshev coefficient', -1.0, 1.0)
  a2 = ROOT.RooRealVar('a2', '2nd Chebyshev coefficient', -1.0, 1.0)
  getattr(workspace,'import')(a0)
  getattr(workspace,'import')(a1)
  getattr(workspace,'import')(a2)
  #pdf_b = ROOT.RooChebyshev('pdf_b', 'pdf_b', fit_var, ROOT.RooArgList(a0, a1, a2))
  pdf_b = ROOT.RooGenericPdf('pdf_b','pdf_b',
                             '1+@1*@0+@2*(2*@0*@0-1)+@3*(4*@0*@0*@0-3*@0)',
                             ROOT.RooArgList(fit_var, a0, a1, a2))
  getattr(workspace,'import')(pdf_b)

def add_background_model_bernstein(workspace, ibin, is_pass):
  '''Adds Bernstein polynomial of degree 4 as background model

  workspace    workspace to add model to; must have variable fit_var
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  '''
  fit_var = workspace.var('fit_var')
  a0 = ROOT.RooRealVar('a0', '0th Bernstein coefficient', 0.0, 1.0)
  a1 = ROOT.RooRealVar('a1', '1st Bernstein coefficient', 0.0, 1.0)
  a2 = ROOT.RooRealVar('a2', '2nd Bernstein coefficient', 0.0, 1.0)
  a3 = ROOT.RooRealVar('a3', '2nd Bernstein coefficient', 0.0, 1.0)
  a4 = ROOT.RooRealVar('a4', '3rd Bernstein coefficient', 0.0, 1.0)
  a0.setVal(0.8)
  a1.setVal(0.7)
  a2.setVal(0.3)
  a3.setVal(0.1)
  a4.setVal(0.1)
  getattr(workspace,'import')(a0)
  getattr(workspace,'import')(a1)
  getattr(workspace,'import')(a2)
  getattr(workspace,'import')(a3)
  getattr(workspace,'import')(a4)
  pdf_b = ROOT.RooBernstein('pdf_b', 'pdf_b', fit_var, ROOT.RooArgList(a0, a1, a2, a3, a4))
  getattr(workspace,'import')(pdf_b)

def add_background_model_bernstein8(workspace, ibin, is_pass):
  '''Adds Bernstein polynomial of degree 8 as background model

  workspace    workspace to add model to; must have variable fit_var
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  '''
  fit_var = workspace.var('fit_var')
  a0 = ROOT.RooRealVar('a0', '0th Bernstein coefficient', 0.0, 1.0)
  a1 = ROOT.RooRealVar('a1', '1st Bernstein coefficient', 0.0, 1.0)
  a2 = ROOT.RooRealVar('a2', '2nd Bernstein coefficient', 0.0, 1.0)
  a3 = ROOT.RooRealVar('a3', '2nd Bernstein coefficient', 0.0, 1.0)
  a4 = ROOT.RooRealVar('a4', '3rd Bernstein coefficient', 0.0, 1.0)
  a5 = ROOT.RooRealVar('a5', '4th Bernstein coefficient', 0.0, 1.0)
  a6 = ROOT.RooRealVar('a6', '5th Bernstein coefficient', 0.0, 1.0)
  a7 = ROOT.RooRealVar('a7', '6th Bernstein coefficient', 0.0, 1.0)
  a8 = ROOT.RooRealVar('a8', '7th Bernstein coefficient', 0.0, 1.0)
  a0.setVal(0.8)
  a1.setVal(0.7)
  a2.setVal(0.3)
  a3.setVal(0.1)
  a4.setVal(0.1)
  a5.setVal(0.1)
  a6.setVal(0.1)
  a7.setVal(0.1)
  a8.setVal(0.1)
  getattr(workspace,'import')(a0)
  getattr(workspace,'import')(a1)
  getattr(workspace,'import')(a2)
  getattr(workspace,'import')(a3)
  getattr(workspace,'import')(a4)
  getattr(workspace,'import')(a5)
  getattr(workspace,'import')(a6)
  getattr(workspace,'import')(a7)
  getattr(workspace,'import')(a8)
  pdf_b = ROOT.RooBernstein('pdf_b', 'pdf_b', fit_var, ROOT.RooArgList(a0, a1, a2, a3, a4, a5, a6, a7, a8))
  getattr(workspace,'import')(pdf_b)

def add_background_model_exponential(workspace, ibin, is_pass):
  '''Adds exponential as background model

  workspace    workspace to add model to; must have variable fit_var
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  '''
  fit_var = workspace.var('fit_var')
  alpha = ROOT.RooRealVar('alpha', 'exponential parameter', -5.0, 5.0)
  getattr(workspace,'import')(alpha)
  pdf_b = ROOT.RooExponential('pdf_b', 'pdf_b', fit_var, alpha)
  getattr(workspace,'import')(pdf_b)

def add_background_model_gamma(workspace, ibin, is_pass):
  '''Adds gamma distribution as background model

  workspace    workspace to add model to; must have variable fit_var
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  '''
  fit_var = workspace.var('fit_var')
  gamma = ROOT.RooRealVar('gamma', 'RooGamma gamma', 0.01, 20.0)
  beta = ROOT.RooRealVar('beta', 'RooGamma beta', 0.5, 100.0)
  gamma_mu = ROOT.RooRealVar('gamma_mu', 'RooGamma mu', 0.0, 60.0)
  gamma.setVal(18.0)
  beta.setVal(1.5)
  gamma_mu.setVal(59.0)
  gamma_mu.setConstant(True)
  getattr(workspace,'import')(gamma)
  getattr(workspace,'import')(beta)
  getattr(workspace,'import')(gamma_mu)
  pdf_b = ROOT.RooGamma('pdf_b', 'pdf_b', fit_var, gamma, beta, gamma_mu)
  getattr(workspace,'import')(pdf_b)

def add_background_model_gammagauss(workspace, ibin, is_pass):
  '''Adds gamma+gauss distribution as background model

  workspace    workspace to add model to; must have variable fit_var
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  '''
  fit_var = workspace.var('fit_var')
  gamma = ROOT.RooRealVar('gamma', 'RooGamma gamma', 0.01, 20.0)
  beta = ROOT.RooRealVar('beta', 'RooGamma beta', 0.5, 100.0)
  gamma_mu = ROOT.RooRealVar('gamma_mu', 'RooGamma mu', 0.0, 60.0)
  gamadd_mu = ROOT.RooRealVar('gamadd_mu', 'RooGaussian mu', 40.0, 140.0)
  gamadd_sigma = ROOT.RooRealVar('gamadd_sigma', 'RooGaussian sigma', 0.5, 60.0)
  gamadd_frac = ROOT.RooRealVar('gamadd_frac', 'Gaussian frac', 0.0, 1.0)
  gamma.setVal(18.0)
  beta.setVal(4.0)
  gamma_mu.setVal(55.0)
  gamadd_mu.setVal(65.0)
  gamadd_sigma.setVal(7.0)
  gamadd_frac.setVal(0.34)
  getattr(workspace,'import')(gamma)
  getattr(workspace,'import')(beta)
  getattr(workspace,'import')(gamma_mu)
  getattr(workspace,'import')(gamadd_mu)
  getattr(workspace,'import')(gamadd_sigma)
  getattr(workspace,'import')(gamadd_frac)
  pdf_b_gamma = ROOT.RooGamma('pdf_b_gamma', 'pdf_b_gamma', fit_var, gamma, 
                              beta, gamma_mu)
  pdf_b_gauss = ROOT.RooGaussian('pdf_b_gauss', 'pdf_b_gauss', fit_var, 
                              gamadd_mu, gamadd_sigma)
  pdf_b = ROOT.RooAddPdf('pdf_b','pdf_b',
                          ROOT.RooArgList(pdf_b_gauss, pdf_b_gamma), 
                          ROOT.RooArgList(gamadd_frac))
  getattr(workspace,'import')(pdf_b)
