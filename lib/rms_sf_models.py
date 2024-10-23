#!/usr/bin/env python3
"""@package docstring
T&P meta analyzer that follows standard EGM procedures to generate scale factors. In particular:
"""

from array import array
import ROOT
import json

MAX_SIGNAL = 100000000.0
MAX_BACKGROUND = 100000000.0

def model_initializer_nom_egm_meta(fit_var, ibin, is_pass, mc_analyzer, highpt_bins):
  '''Model initializer that reutrns a tnp workspace where the signal model
  is an MC shape convolved with a Gaussian and the background model is a 
  CMSshape (erf*exp)

  fit_var      RooRealVar representing variable to fit
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  mc_analyzer  TnpAnalyzer for MC samples
  highpt_bins  list of ints, bins where passing histogram is used always
  '''

  #set up RooFit stuff
  workspace = ROOT.RooWorkspace()
  getattr(workspace,'import')(fit_var)

  mean = ROOT.RooRealVar('mean', 'Smearing gaussian peak', -5.0, 5.0) 
  sigma = ROOT.RooRealVar('sigma', 'Smearing gaussian sigma', 0.5, 5.0) 

  acms = ROOT.RooRealVar('acms', 'erf turn-on point', 50.0, 80.0)
  beta = ROOT.RooRealVar('beta', 'erf width', 0.01, 0.08)
  gamma = ROOT.RooRealVar('gamma', 'exp parameter', -2.0, 2.0)
  peak = ROOT.RooRealVar('peak', 'peak', 90.0, 90.0)

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, MAX_SIGNAL) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, MAX_BACKGROUND) 

  getattr(workspace,'import')(mean)
  getattr(workspace,'import')(sigma)
  getattr(workspace,'import')(acms)
  getattr(workspace,'import')(beta)
  getattr(workspace,'import')(gamma)
  getattr(workspace,'import')(peak)
  getattr(workspace,'import')(nSig)
  getattr(workspace,'import')(nBkg)

  pass_fail = 'pass'
  if not is_pass and not (ibin in highpt_bins):
    pass_fail = 'fail'
  mc_filename = 'out/'+mc_analyzer.temp_name+'/'+mc_analyzer.temp_name+'.root'
  hist_name = 'hist_'+pass_fail+'_bin{}'.format(ibin)
  input_file = ROOT.TFile(mc_filename,'READ')
  hist = input_file.Get(hist_name)
  hist.SetDirectory(ROOT.nullptr)
  input_file.Close()

  core_hist = ROOT.RooDataHist('pdf_s_core_hist',
                               'pdf_s_core_hist',
                               ROOT.RooArgList(fit_var), hist)
  getattr(workspace,'import')(core_hist)
  pdf_s_core = ROOT.RooHistPdf('pdf_s_core','pdf_s_core',
                               ROOT.RooArgSet(fit_var),core_hist)
  pdf_s_res = ROOT.RooGaussian('pdf_s_res','pdf_s_res',fit_var,mean,sigma)
  pdf_s = ROOT.RooFFTConvPdf('pdf_s','pdf_s',fit_var,pdf_s_core,pdf_s_res)
  pdf_b = ROOT.RooCMSShape('pdf_b', 'pdf_b', fit_var, acms, beta, gamma, peak)
  pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(nSig, nBkg))
  getattr(workspace,'import')(pdf_sb)
  return workspace

def model_initializer_cheby_egm_meta(fit_var, ibin, is_pass, mc_analyzer, highpt_bins):
  '''Model initializer that reutrns a tnp workspace where the signal model
  is an MC shape convolved with a Gaussian and the background model is a 
  Chebyshev polynomial

  fit_var      RooRealVar representing variable to fit
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  mc_analyzer  TnpAnalyzer for MC samples
  highpt_bins  list of ints, bins where passing histogram is used always
  '''

  #set up RooFit stuff
  workspace = ROOT.RooWorkspace()
  getattr(workspace,'import')(fit_var)

  mean = ROOT.RooRealVar('mean', 'Smearing gaussian peak', -5.0, 5.0) 
  sigma = ROOT.RooRealVar('sigma', 'Smearing gaussian sigma', 0.5, 5.0) 

  a0 = ROOT.RooRealVar('a0', '0th Chebyshev coefficient', -1.0, 1.0)
  a1 = ROOT.RooRealVar('a1', '1st Chebyshev coefficient', -1.0, 1.0)
  a2 = ROOT.RooRealVar('a2', '2nd Chebyshev coefficient', -1.0, 1.0)

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, MAX_SIGNAL) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, MAX_BACKGROUND) 

  getattr(workspace,'import')(mean)
  getattr(workspace,'import')(sigma)
  getattr(workspace,'import')(a0)
  getattr(workspace,'import')(a1)
  getattr(workspace,'import')(a2)
  getattr(workspace,'import')(nSig)
  getattr(workspace,'import')(nBkg)

  pass_fail = 'pass'
  if not is_pass and not (ibin in highpt_bins):
    pass_fail = 'fail'
  mc_filename = 'out/'+mc_analyzer.temp_name+'/'+mc_analyzer.temp_name+'.root'
  hist_name = 'hist_'+pass_fail+'_bin{}'.format(ibin)
  input_file = ROOT.TFile(mc_filename,'READ')
  hist = input_file.Get(hist_name)
  hist.SetDirectory(ROOT.nullptr)
  input_file.Close()

  core_hist = ROOT.RooDataHist('pdf_s_core_hist',
                               'pdf_s_core_hist',
                               ROOT.RooArgList(fit_var), hist)
  getattr(workspace,'import')(core_hist)
  pdf_s_core = ROOT.RooHistPdf('pdf_s_core','pdf_s_core',
                               ROOT.RooArgSet(fit_var),core_hist)
  pdf_s_res = ROOT.RooGaussian('pdf_s_res','pdf_s_res',fit_var,mean,sigma)
  pdf_s = ROOT.RooFFTConvPdf('pdf_s','pdf_s',fit_var,pdf_s_core,pdf_s_res)
  #pdf_b = ROOT.RooChebyshev('pdf_b', 'pdf_b', fit_var, ROOT.RooArgList(a0, a1, a2))
  pdf_b = ROOT.RooGenericPdf('pdf_b','pdf_b',
                             '1+@1*@0+@2*(2*@0*@0-1)+@3*(4*@0*@0*@0-3*@0)',
                             ROOT.RooArgList(fit_var, a0, a1, a2))
  pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(nSig, nBkg))
  getattr(workspace,'import')(pdf_sb)
  return workspace

def model_initializer_bern_egm_meta(fit_var, ibin, is_pass, mc_analyzer, highpt_bins):
  '''Model initializer that reutrns a tnp workspace where the signal model
  is an MC shape convolved with a Gaussian and the background model is a 
  Bernstein polynomial

  fit_var      RooRealVar representing variable to fit
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  mc_analyzer  TnpAnalyzer for MC samples
  highpt_bins  list of ints, bins where passing histogram is used always
  '''

  #set up RooFit stuff
  workspace = ROOT.RooWorkspace()
  getattr(workspace,'import')(fit_var)

  mean = ROOT.RooRealVar('mean', 'Smearing gaussian peak', -5.0, 5.0) 
  sigma = ROOT.RooRealVar('sigma', 'Smearing gaussian sigma', 0.5, 5.0) 

  a0 = ROOT.RooRealVar('a0', '0th Bernstein coefficient', 0.0, 1.0)
  a1 = ROOT.RooRealVar('a1', '1st Bernstein coefficient', 0.0, 1.0)
  a2 = ROOT.RooRealVar('a2', '2nd Bernstein coefficient', 0.0, 1.0)
  a3 = ROOT.RooRealVar('a3', '2nd Bernstein coefficient', 0.0, 1.0)

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, MAX_SIGNAL) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, MAX_BACKGROUND) 

  getattr(workspace,'import')(mean)
  getattr(workspace,'import')(sigma)
  getattr(workspace,'import')(a0)
  getattr(workspace,'import')(a1)
  getattr(workspace,'import')(a2)
  getattr(workspace,'import')(a3)
  getattr(workspace,'import')(nSig)
  getattr(workspace,'import')(nBkg)

  pass_fail = 'pass'
  if not is_pass and not (ibin in highpt_bins):
    pass_fail = 'fail'
  mc_filename = 'out/'+mc_analyzer.temp_name+'/'+mc_analyzer.temp_name+'.root'
  hist_name = 'hist_'+pass_fail+'_bin{}'.format(ibin)
  input_file = ROOT.TFile(mc_filename,'READ')
  hist = input_file.Get(hist_name)
  hist.SetDirectory(ROOT.nullptr)
  input_file.Close()

  core_hist = ROOT.RooDataHist('pdf_s_core_hist',
                               'pdf_s_core_hist',
                               ROOT.RooArgList(fit_var), hist)
  getattr(workspace,'import')(core_hist)
  pdf_s_core = ROOT.RooHistPdf('pdf_s_core','pdf_s_core',
                               ROOT.RooArgSet(fit_var),core_hist)
  pdf_s_res = ROOT.RooGaussian('pdf_s_res','pdf_s_res',fit_var,mean,sigma)
  pdf_s = ROOT.RooFFTConvPdf('pdf_s','pdf_s',fit_var,pdf_s_core,pdf_s_res)
  pdf_b = ROOT.RooBernstein('pdf_b', 'pdf_b', fit_var, ROOT.RooArgList(a0, a1, a2, a3))
  pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(nSig, nBkg))
  getattr(workspace,'import')(pdf_sb)
  return workspace

def model_initializer_cheby_dscb(fit_var, ibin, is_pass):
  '''Model initializer that reutrns a tnp workspace where the signal model
  is a DSCB shape convolved with a Gaussian and the background model is a 
  Chebyshev polynomial

  fit_var      RooRealVar representing variable to fit
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  '''

  #set up RooFit stuff
  workspace = ROOT.RooWorkspace()
  getattr(workspace,'import')(fit_var)

  gauss_mu = ROOT.RooRealVar('mean', 'Gaussian mean', 85.0, 95.0) 
  gauss_sigma = ROOT.RooRealVar('sigma', 'Gaussian sigma', 0.01, 15.0) 
  cb_alphal = ROOT.RooRealVar('alphal', 'CB left switchover', 0.1, 10.0) 
  cb_nl = ROOT.RooRealVar('nl', 'CB left power', 0.1, 10.0) 
  cb_alphar = ROOT.RooRealVar('alphar', 'CB right switchover', 0.1, 10.0) 
  cb_nr = ROOT.RooRealVar('nr', 'CB right power', 0.1, 10.0) 

  mean = ROOT.RooRealVar('mean', 'Smearing gaussian peak', -5.0, 5.0) 
  sigma = ROOT.RooRealVar('sigma', 'Smearing gaussian sigma', 0.5, 5.0) 

  a0 = ROOT.RooRealVar('a0', '0th Chebyshev coefficient', -1.0, 1.0)
  a1 = ROOT.RooRealVar('a1', '1st Chebyshev coefficient', -1.0, 1.0)
  a2 = ROOT.RooRealVar('a2', '2nd Chebyshev coefficient', -1.0, 1.0)

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, MAX_SIGNAL) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, MAX_BACKGROUND) 

  getattr(workspace,'import')(gauss_mu)
  getattr(workspace,'import')(gauss_sigma)
  getattr(workspace,'import')(cb_alphal)
  getattr(workspace,'import')(cb_nl)
  getattr(workspace,'import')(cb_alphar)
  getattr(workspace,'import')(cb_nr)
  getattr(workspace,'import')(a0)
  getattr(workspace,'import')(a1)
  getattr(workspace,'import')(a2)
  getattr(workspace,'import')(nSig)
  getattr(workspace,'import')(nBkg)

  pdf_s  = ROOT.RooCrystalBall('pdf_s','pdf_s', fit_var, gauss_mu, gauss_sigma, cb_alphal, cb_nl, cb_alphar, cb_nr)
  #pdf_b = ROOT.RooChebyshev('pdf_b', 'pdf_b', fit_var, ROOT.RooArgList(a0, a1, a2))
  pdf_b = ROOT.RooGenericPdf('pdf_b','pdf_b',
                             '1+@1*@0+@2*(2*@0*@0-1)+@3*(4*@0*@0*@0-3*@0)',
                             ROOT.RooArgList(fit_var, a0, a1, a2))
  pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(nSig, nBkg))
  getattr(workspace,'import')(pdf_sb)
  return workspace

def model_initializer_dscb(fit_var, ibin, is_pass):
  '''Model initializer that reutrns a tnp workspace where the signal model
  is a DSCB shape 

  fit_var      RooRealVar representing variable to fit
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  '''

  #set up RooFit stuff
  workspace = ROOT.RooWorkspace()
  getattr(workspace,'import')(fit_var)

  gauss_mu = ROOT.RooRealVar('mean', 'Gaussian mean', 85.0, 95.0) 
  gauss_sigma = ROOT.RooRealVar('sigma', 'Gaussian sigma', 0.01, 15.0) 
  cb_alphal = ROOT.RooRealVar('alphal', 'CB left switchover', 0.1, 10.0) 
  cb_nl = ROOT.RooRealVar('nl', 'CB left power', 0.1, 10.0) 
  cb_alphar = ROOT.RooRealVar('alphar', 'CB right switchover', 0.1, 10.0) 
  cb_nr = ROOT.RooRealVar('nr', 'CB right power', 0.1, 10.0) 

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, MAX_SIGNAL) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, 0.0) 

  getattr(workspace,'import')(gauss_mu)
  getattr(workspace,'import')(gauss_sigma)
  getattr(workspace,'import')(cb_alphal)
  getattr(workspace,'import')(cb_nl)
  getattr(workspace,'import')(cb_alphar)
  getattr(workspace,'import')(cb_nr)
  getattr(workspace,'import')(nSig)
  getattr(workspace,'import')(nBkg)

  pdf_sb  = ROOT.RooCrystalBall('pdf_sb','pdf_sb', fit_var, gauss_mu, gauss_sigma, cb_alphal, cb_nl, cb_alphar, cb_nr)
  getattr(workspace,'import')(pdf_sb)
  return workspace

def model_initializer_bern_dscb(fit_var, ibin, is_pass):
  '''Model initializer that reutrns a tnp workspace where the signal model
  is a DSCB shape convolved with a Gaussian and the background model is a 
  Bernstein polynomial

  fit_var      RooRealVar representing variable to fit
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  '''

  #set up RooFit stuff
  workspace = ROOT.RooWorkspace()
  getattr(workspace,'import')(fit_var)

  gauss_mu = ROOT.RooRealVar('mean', 'Gaussian mean', 85.0, 95.0) 
  gauss_sigma = ROOT.RooRealVar('sigma', 'Gaussian sigma', 0.01, 15.0) 
  cb_alphal = ROOT.RooRealVar('alphal', 'CB left switchover', 0.1, 10.0) 
  cb_nl = ROOT.RooRealVar('nl', 'CB left power', 0.1, 10.0) 
  cb_alphar = ROOT.RooRealVar('alphar', 'CB right switchover', 0.1, 10.0) 
  cb_nr = ROOT.RooRealVar('nr', 'CB right power', 0.1, 10.0) 

  mean = ROOT.RooRealVar('mean', 'Smearing gaussian peak', -5.0, 5.0) 
  sigma = ROOT.RooRealVar('sigma', 'Smearing gaussian sigma', 0.5, 5.0) 

  a0 = ROOT.RooRealVar('a0', '0th Bernstein coefficient', 0.0, 1.0)
  a1 = ROOT.RooRealVar('a1', '1st Bernstein coefficient', 0.0, 1.0)
  a2 = ROOT.RooRealVar('a2', '2nd Bernstein coefficient', 0.0, 1.0)
  a3 = ROOT.RooRealVar('a3', '2nd Bernstein coefficient', 0.0, 1.0)

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, MAX_SIGNAL) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, MAX_BACKGROUND) 

  getattr(workspace,'import')(gauss_mu)
  getattr(workspace,'import')(gauss_sigma)
  getattr(workspace,'import')(cb_alphal)
  getattr(workspace,'import')(cb_nl)
  getattr(workspace,'import')(cb_alphar)
  getattr(workspace,'import')(cb_nr)
  getattr(workspace,'import')(a0)
  getattr(workspace,'import')(a1)
  getattr(workspace,'import')(a2)
  getattr(workspace,'import')(a3)
  getattr(workspace,'import')(nSig)
  getattr(workspace,'import')(nBkg)

  pdf_s  = ROOT.RooCrystalBall('pdf_s','pdf_s', fit_var, gauss_mu, gauss_sigma, cb_alphal, cb_nl, cb_alphar, cb_nr)
  pdf_b = ROOT.RooBernstein('pdf_b', 'pdf_b', fit_var, ROOT.RooArgList(a0, a1, a2, a3))
  pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(nSig, nBkg))
  getattr(workspace,'import')(pdf_sb)
  return workspace

def param_initializer_zero_bkg_cheby(ibin, is_pass, workspace):
  '''
  Parameter initializer for cheby_dscb model that fixes the background to zero
  
  ibin       int, bin number
  is_pass    bool, indicates if passing leg
  workspace  RooWorkspace for this bin
  '''
  workspace.var('nBkg').setVal(0.0)
  workspace.var('a0').setConstant()
  workspace.var('a1').setConstant()
  workspace.var('a2').setConstant()
  workspace.var('nBkg').setConstant()

def param_initializer_zero_bkg_bern(ibin, is_pass, workspace):
  '''
  Parameter initializer for bern_dscb model that fixes the background to zero
  
  ibin       int, bin number
  is_pass    bool, indicates if passing leg
  workspace  RooWorkspace for this bin
  '''
  workspace.var('nBkg').setVal(0.0)
  workspace.var('a0').setConstant()
  workspace.var('a1').setConstant()
  workspace.var('a2').setConstant()
  workspace.var('a3').setConstant()
  workspace.var('nBkg').setConstant()

def param_initializer_dscb_from_mc_meta(ibin, is_pass, workspace, mc_analyzer):
  '''
  Parameter initializer for cheby_dscb model that fixes DSCB parameters except
  mean and sigma to MC result
  
  ibin         int, bin number
  is_pass      bool, indicates if passing leg
  workspace    RooWorkspace for this bin
  mc_analyzer  TnpAnalyzer for MC samples
  '''
  pass_fail = 'pass'
  if not is_pass:
    pass_fail = 'fail'
  mc_json_filename = 'out/'+mc_analyzer.temp_name+'/fitinfo_bin'+str(ibin)+'_'+pass_fail+'.json'
  with open(mc_json_filename,'r') as mc_file:
    param_dict = json.loads(mc_file.read())
    for var in ['mean','sigma','alphal','nl','alphar','nr']:
      workspace.var(var).setVal(param_dict[var])
    for var in ['alphal','nl','alphar','nr']:
      workspace.var('alphal').setConstant()
      workspace.var('nl').setConstant()
      workspace.var('alphar').setConstant()
      workspace.var('nr').setConstant()

def make_param_initializer_dscb_from_mc(mc_analyzer):
  '''
  make a param_initializer from the meta version

  mc_analyzer  TnpAnalyzer for MC samples
  '''
  return (lambda ibin, is_pass, workspace : 
          param_initializer_dscb_from_mc_meta(ibin, is_pass, workspace, 
                                              mc_analyzer))

def model_initializer_exp_egm_meta(fit_var, ibin, is_pass, mc_analyzer, highpt_bins):
  '''Model initializer that reutrns a tnp workspace where the signal model
  is an MC shape convolved with a Gaussian and the background model is an exponential

  fit_var      RooRealVar representing variable to fit
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  mc_analyzer  TnpAnalyzer for MC samples
  highpt_bins  list of ints, bins where passing histogram is used always
  '''

  #set up RooFit stuff
  workspace = ROOT.RooWorkspace()
  getattr(workspace,'import')(fit_var)

  mean = ROOT.RooRealVar('mean', 'Smearing gaussian peak', -5.0, 5.0) 
  sigma = ROOT.RooRealVar('sigma', 'Smearing gaussian sigma', 0.5, 5.0) 

  alpha = ROOT.RooRealVar('alpha', 'exponential parameter', -5.0, 5.0)

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, MAX_SIGNAL) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, MAX_BACKGROUND) 

  getattr(workspace,'import')(mean)
  getattr(workspace,'import')(sigma)
  getattr(workspace,'import')(alpha)
  getattr(workspace,'import')(nSig)
  getattr(workspace,'import')(nBkg)

  pass_fail = 'pass'
  if not is_pass and not (ibin in highpt_bins):
    pass_fail = 'fail'
  mc_filename = 'out/'+mc_analyzer.temp_name+'/'+mc_analyzer.temp_name+'.root'
  hist_name = 'hist_'+pass_fail+'_bin{}'.format(ibin)
  input_file = ROOT.TFile(mc_filename,'READ')
  hist = input_file.Get(hist_name)
  hist.SetDirectory(ROOT.nullptr)
  input_file.Close()

  core_hist = ROOT.RooDataHist('pdf_s_core_hist',
                               'pdf_s_core_hist',
                               ROOT.RooArgList(fit_var), hist)
  getattr(workspace,'import')(core_hist)
  pdf_s_core = ROOT.RooHistPdf('pdf_s_core','pdf_s_core',
                               ROOT.RooArgSet(fit_var),core_hist)
  pdf_s_res = ROOT.RooGaussian('pdf_s_res','pdf_s_res',fit_var,mean,sigma)
  pdf_s = ROOT.RooFFTConvPdf('pdf_s','pdf_s',fit_var,pdf_s_core,pdf_s_res)
  pdf_b = ROOT.RooExponential('pdf_b', 'pdf_b', fit_var, alpha)
  pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(nSig, nBkg))
  getattr(workspace,'import')(pdf_sb)
  return workspace

def model_initializer_gamma_egm_meta(fit_var, ibin, is_pass, mc_analyzer, highpt_bins):
  '''Model initializer that reutrns a tnp workspace where the signal model
  is an MC shape convolved with a Gaussian and the background model is a gamma
  distribution

  fit_var      RooRealVar representing variable to fit
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  mc_analyzer  TnpAnalyzer for MC samples
  highpt_bins  list of ints, bins where passing histogram is used always
  '''

  #set up RooFit stuff
  workspace = ROOT.RooWorkspace()
  getattr(workspace,'import')(fit_var)

  mean = ROOT.RooRealVar('mean', 'Smearing gaussian peak', -5.0, 5.0) 
  sigma = ROOT.RooRealVar('sigma', 'Smearing gaussian sigma', 0.5, 5.0) 

  a0 = ROOT.RooRealVar('gamma', 'RooGamma gamma', 0.01, 20.0)
  a1 = ROOT.RooRealVar('beta', 'RooGamma beta', 0.1, 20.0)
  a2 = ROOT.RooRealVar('mu', 'RooGamma mu', 0.0, 60.0)

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, MAX_SIGNAL) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, MAX_BACKGROUND) 

  getattr(workspace,'import')(mean)
  getattr(workspace,'import')(sigma)
  getattr(workspace,'import')(a0)
  getattr(workspace,'import')(a1)
  getattr(workspace,'import')(a2)
  getattr(workspace,'import')(nSig)
  getattr(workspace,'import')(nBkg)

  pass_fail = 'pass'
  if not is_pass and not (ibin in highpt_bins):
    pass_fail = 'fail'
  mc_filename = 'out/'+mc_analyzer.temp_name+'/'+mc_analyzer.temp_name+'.root'
  hist_name = 'hist_'+pass_fail+'_bin{}'.format(ibin)
  input_file = ROOT.TFile(mc_filename,'READ')
  hist = input_file.Get(hist_name)
  hist.SetDirectory(ROOT.nullptr)
  input_file.Close()

  core_hist = ROOT.RooDataHist('pdf_s_core_hist',
                               'pdf_s_core_hist',
                               ROOT.RooArgList(fit_var), hist)
  getattr(workspace,'import')(core_hist)
  pdf_s_core = ROOT.RooHistPdf('pdf_s_core','pdf_s_core',
                               ROOT.RooArgSet(fit_var),core_hist)
  pdf_s_res = ROOT.RooGaussian('pdf_s_res','pdf_s_res',fit_var,mean,sigma)
  pdf_s = ROOT.RooFFTConvPdf('pdf_s','pdf_s',fit_var,pdf_s_core,pdf_s_res)
  pdf_b = ROOT.RooGamma('pdf_b', 'pdf_b', fit_var, a0, a1, a2)
  pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(nSig, nBkg))
  getattr(workspace,'import')(pdf_sb)
  return workspace

def model_initializer_gamma_dscb(fit_var, ibin, is_pass):
  '''Model initializer that reutrns a tnp workspace where the signal model
  is a DSCB shape convolved with a Gaussian and the background model is a 
  gamma distribution

  fit_var      RooRealVar representing variable to fit
  ibin         int bin number
  is_pass      bool, indicates if passing or failing leg
  '''

  #set up RooFit stuff
  workspace = ROOT.RooWorkspace()
  getattr(workspace,'import')(fit_var)

  gauss_mu = ROOT.RooRealVar('mean', 'Gaussian mean', 85.0, 95.0) 
  gauss_sigma = ROOT.RooRealVar('sigma', 'Gaussian sigma', 0.01, 15.0) 
  cb_alphal = ROOT.RooRealVar('alphal', 'CB left switchover', 0.1, 10.0) 
  cb_nl = ROOT.RooRealVar('nl', 'CB left power', 0.1, 10.0) 
  cb_alphar = ROOT.RooRealVar('alphar', 'CB right switchover', 0.1, 10.0) 
  cb_nr = ROOT.RooRealVar('nr', 'CB right power', 0.1, 10.0) 

  mean = ROOT.RooRealVar('mean', 'Smearing gaussian peak', -5.0, 5.0) 
  sigma = ROOT.RooRealVar('sigma', 'Smearing gaussian sigma', 0.5, 5.0) 

  a0 = ROOT.RooRealVar('gamma', 'RooGamma gamma', 0.01, 20.0)
  a1 = ROOT.RooRealVar('beta', 'RooGamma beta', 0.1, 20.0)
  a2 = ROOT.RooRealVar('mu', 'RooGamma mu', 0.0, 60.0)

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, MAX_SIGNAL) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, MAX_BACKGROUND) 

  getattr(workspace,'import')(gauss_mu)
  getattr(workspace,'import')(gauss_sigma)
  getattr(workspace,'import')(cb_alphal)
  getattr(workspace,'import')(cb_nl)
  getattr(workspace,'import')(cb_alphar)
  getattr(workspace,'import')(cb_nr)
  getattr(workspace,'import')(a0)
  getattr(workspace,'import')(a1)
  getattr(workspace,'import')(a2)
  getattr(workspace,'import')(nSig)
  getattr(workspace,'import')(nBkg)

  pdf_s  = ROOT.RooCrystalBall('pdf_s','pdf_s', fit_var, gauss_mu, gauss_sigma, cb_alphal, cb_nl, cb_alphar, cb_nr)
  pdf_b = ROOT.RooGamma('pdf_b', 'pdf_b', fit_var, a0, a1, a2)
  pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(nSig, nBkg))
  getattr(workspace,'import')(pdf_sb)
  return workspace

def make_model_initializer(meta_initializer, mc_analyzer, highpt_bins):
  '''
  creates model_initializer function lambda from a meta initializer

  meta_initializer function returning model with extra parameters
  mc_analyzer      TnpAnalyzer for MC samples
  highpt_bins      list of ints, bins where passing histogram is always used
  '''
  return (lambda fit_var, ibin, is_pass : meta_initializer(
      fit_var, ibin, is_pass, mc_analyzer, highpt_bins))
