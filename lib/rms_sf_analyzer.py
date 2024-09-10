#!/usr/bin/env python3
"""@package docstring
T&P meta analyzer that follows standard EGM procedures to generate scale factors. In particular:
"""

import os
import ROOT

from tnp_analyzer import *
from model_initializers import *

def model_initializer_nom_egm_meta(fit_var, ibin=0, is_pass=True, mc_analyzer, highpt_bins):
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

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, 1000000.0) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, 1000000.0) 

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
  hist_name = 'hist_'+pass_fail+'_'+mc_analyzer.get_binname(ibin)
  with ROOT.TFile(mc_filename,'READ') as input_file:
    hist = input_file.Get(hist_name)
    hist.SetDirectory(ROOT.nullptr)

  pdf_s_core = ROOT.RooHistPdf('pdf_s_core','pdf_s_core',
                               ROOT.RooArgSet(fit_var),
                               ROOT.RooDataHist('pdf_s_core_hist',
                                                'pdf_s_core_hist',
                                                ROOT.RooArgList(fit_var), hist))
  pdf_s_res = ROOT.RooGaussian('pdf_s_res','pdf_s_res',fit_var,mean,sigma)
  pdf_s = ROOT.RooFFTConvPdf('pdf_s','pdf_s',fit_var,pdf_s_core,pdf_s_res)
  pdf_b = ROOT.RooCMSShape('pdf_b', 'pdf_b', fit_var, acms, beta, gamma, peak)
  pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(nSig, nBkg))
  getattr(workspace,'import')(pdf_sb)
  return workspace

def model_initializer_cheby_egm_meta(fit_var, ibin=0, is_pass=True, mc_analyzer, highpt_bins):
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
  a1 = ROOT.RooRealVar('a0', '1st Chebyshev coefficient', -1.0, 1.0)
  a2 = ROOT.RooRealVar('a0', '2nd Chebyshev coefficient', -1.0, 1.0)

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, 1000000.0) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, 1000000.0) 

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
  hist_name = 'hist_'+pass_fail+'_'+mc_analyzer.get_binname(ibin)
  with ROOT.TFile(mc_filename,'READ') as input_file:
    hist = input_file.Get(hist_name)
    hist.SetDirectory(ROOT.nullptr)

  pdf_s_core = ROOT.RooHistPdf('pdf_s_core','pdf_s_core',
                               ROOT.RooArgSet(fit_var),
                               ROOT.RooDataHist('pdf_s_core_hist',
                                                'pdf_s_core_hist',
                                                ROOT.RooArgList(fit_var), hist))
  pdf_s_res = ROOT.RooGaussian('pdf_s_res','pdf_s_res',fit_var,mean,sigma)
  pdf_s = ROOT.RooFFTConvPdf('pdf_s','pdf_s',fit_var,pdf_s_core,pdf_s_res)
  pdf_b = ROOT.RooChebyshev('pdf_b', 'pdf_b', fit_var, ROOT.RooArgList(a0, a1, a2))
  pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(nSig, nBkg))
  getattr(workspace,'import')(pdf_sb)
  return workspace

def model_initializer_cheby_dscb(fit_var, ibin=0, is_pass=True):
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

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, 1000000.0) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, 1000000.0) 

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
  pdf_b = ROOT.RooChebyshev('pdf_b', 'pdf_b', fit_var, ROOT.RooArgList(a0, a1, a2))
  pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(nSig, nBkg))
  getattr(workspace,'import')(pdf_sb)
  return workspace

def param_initializer_zero_bkg(ibin, is_pass, workspace):
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

def model_initializer_exp_egm_meta(fit_var, ibin=0, is_pass=True, mc_analyzer, highpt_bins):
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

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, 1000000.0) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, 1000000.0) 

  getattr(workspace,'import')(mean)
  getattr(workspace,'import')(sigma)
  getattr(workspace,'import')(alpha)
  getattr(workspace,'import')(nSig)
  getattr(workspace,'import')(nBkg)

  pass_fail = 'pass'
  if not is_pass and not (ibin in highpt_bins):
    pass_fail = 'fail'
  mc_filename = 'out/'+mc_analyzer.temp_name+'/'+mc_analyzer.temp_name+'.root'
  hist_name = 'hist_'+pass_fail+'_'+mc_analyzer.get_binname(ibin)
  with ROOT.TFile(mc_filename,'READ') as input_file:
    hist = input_file.Get(hist_name)
    hist.SetDirectory(ROOT.nullptr)

  pdf_s_core = ROOT.RooHistPdf('pdf_s_core','pdf_s_core',
                               ROOT.RooArgSet(fit_var),
                               ROOT.RooDataHist('pdf_s_core_hist',
                                                'pdf_s_core_hist',
                                                ROOT.RooArgList(fit_var), hist))
  pdf_s_res = ROOT.RooGaussian('pdf_s_res','pdf_s_res',fit_var,mean,sigma)
  pdf_s = ROOT.RooFFTConvPdf('pdf_s','pdf_s',fit_var,pdf_s_core,pdf_s_res)
  pdf_b = ROOT.RooExponential('pdf_b', 'pdf_b', fit_var, alpha)
  pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(nSig, nBkg))
  getattr(workspace,'import')(pdf_sb)
  return workspace

def model_initializer_gamma_egm_meta(fit_var, ibin=0, is_pass=True, mc_analyzer, highpt_bins):
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

  a0 = ROOT.RooRealVar('gamma', 'RooGamma gamma', 0.001, 0.1)
  a1 = ROOT.RooRealVar('beta', 'RooGamma beta', 10.0, 200.0)
  a2 = ROOT.RooRealVar('mu', 'RooGamma mu', 10.0, 5.0)

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, 1000000.0) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, 1000000.0) 

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
  hist_name = 'hist_'+pass_fail+'_'+mc_analyzer.get_binname(ibin)
  with ROOT.TFile(mc_filename,'READ') as input_file:
    hist = input_file.Get(hist_name)
    hist.SetDirectory(ROOT.nullptr)

  pdf_s_core = ROOT.RooHistPdf('pdf_s_core','pdf_s_core',
                               ROOT.RooArgSet(fit_var),
                               ROOT.RooDataHist('pdf_s_core_hist',
                                                'pdf_s_core_hist',
                                                ROOT.RooArgList(fit_var), hist))
  pdf_s_res = ROOT.RooGaussian('pdf_s_res','pdf_s_res',fit_var,mean,sigma)
  pdf_s = ROOT.RooFFTConvPdf('pdf_s','pdf_s',fit_var,pdf_s_core,pdf_s_res)
  pdf_b = ROOT.RooGamma('pdf_b', 'pdf_b', fit_var, a0, a1, a2)
  pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(nSig, nBkg))
  getattr(workspace,'import')(pdf_sb)
  return workspace

def model_initializer_gamma_dscb(fit_var, ibin=0, is_pass=True):
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

  a0 = ROOT.RooRealVar('gamma', 'RooGamma gamma', 0.001, 0.1)
  a1 = ROOT.RooRealVar('beta', 'RooGamma beta', 10.0, 200.0)
  a2 = ROOT.RooRealVar('mu', 'RooGamma mu', 10.0, 5.0)

  nSig = ROOT.RooRealVar('nSig', 'Signal normalization', 0.0, 1000000.0) 
  nBkg = ROOT.RooRealVar('nBkg', 'Background normalization', 0.0, 1000000.0) 

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

def make_model_initializer(mc_analyzer, highpt_bins, meta_initializer)
  '''
  creates model_initializer function lambda from a meta initializer

  mc_analyzer      TnpAnalyzer for MC samples
  highpt_bins      list of ints, bins where passing histogram is always used
  meta_initializer function returning model with extra parameters
  '''
  return (lambda fit_var, ibin, is_pas : meta_initializer(
      fit_var, ibin, is_pass, mc_analyzer, highpt_bins))

class RmsSFAnalyzer:

  def __init__(self, name):
    '''
    Constructor

    name  string, name of analyzer
    '''
    ROOT.EnableImplicitMT()
    self.name = name
    self.data_nom_tnp_analyzer = TnpAnalyzer(name+'_data_nom')
    self.data_altsig_tnp_analyzer = TnpAnalyzer(name+'_data_altsig')
    self.data_altbkg_tnp_analyzer = TnpAnalyzer(name+'_data_altbkg')
    self.data_altsigbkg_tnp_analyzer = TnpAnalyzer(name+'_data_altsigbkg')
    self.mc_tnp_analyzer = TnpAnalyzer(name+'_mc_nom')
    self.mc_alt_tnp_analyzer = TnpAnalyzer(name+'_mc_alt')

  def set_input_files(self, data_files, mc_files, mc_alt_files, data_tree, 
                      mc_tree='', mc_alt_tree=''):
    '''
    Sets input files

    data_files    list of strings, data files to process
    mc_files      list of strings, mc files to process
    mc_alt_files  list of strings, alternative MC files to process
    data_tree     string, name of TTree in data files
    mc_tree       string, name of TTree in MC files, defaults to data_tree if ''
    mc_alt_tree   string, name of TTree in MC alt file, defaults to mc_tree if ''
    '''
    if mc_tree=='':
      mc_tree = data_tree
    if mc_alt_tree=='':
      mc_alt_tree = data_tree
    self.data_nom_tnp_analyzer.set_input_files(data_files, data_tree)
    self.data_altsig_tnp_analyzer.set_input_files(data_files, data_tree)
    self.data_altbkg_tnp_analyzer.set_input_files(data_files, data_tree)
    self.data_altsigbkg_tnp_analyzer.set_input_files(data_files, data_tree)
    self.mc_nom_tnp_analyzer.set_input_files(mc_files, mc_tree)
    self.mc_alt_tnp_analyzer.set_input_files(mc_alt_files, mc_alt_tree)
    
  def set_fitting_variable(self, name, description, nbins, 
                           var_range, weight='1'):
    '''
    Adds information about fitting variable

    name         string, name of branch in TTree or C++ expression
    description  string, name used in plots (with TLaTeX)
    nbins        int, number of bins to use for fit variable
    var_range    tuple of two floats, start and end of fit range
    weight       string, expression for weight to use
    '''
    self.data_nom_tnp_analyzer.set_fitting_variable(name, description, nbins, 
                                                    var_range, weight)
    self.data_altsig_tnp_analyzer.set_fitting_variable(name, description, 
                                                       nbins, var_range, weight)
    self.data_altbkg_tnp_analyzer.set_fitting_variable(name, description, 
                                                       nbins, var_range, weight)
    self.data_altsigbkg_tnp_analyzer.set_fitting_variable(name, description, 
                                                          nbins, var_range, 
                                                          weight)
    self.mc_nom_tnp_analyzer.set_fitting_variable(name, description, nbins, 
                                                  var_range, weight)
    self.mc_alt_tnp_analyzer.set_fitting_variable(name, description, nbins, 
                                                  var_range, weight)

  def set_measurement_variable(self, var):
    '''
    Sets selection efficiency to measure with tag & probe

    var  string, name of branch in TTree or C++ expression
    '''
    self.data_nom_tnp_analyzer.set_measurement_variable(var)
    self.data_altsig_tnp_analyzer.set_measurement_variable(var)
    self.data_altbkg_tnp_analyzer.set_measurement_variable(var)
    self.data_altsigbkg_tnp_analyzer.set_measurement_variable(var)
    self.mc_nom_tnp_analyzer.set_measurement_variable(var)
    self.mc_alt_tnp_analyzer.set_measurement_variable(var)

  def set_preselection(self, preselection, desc):
    '''
    Sets basic preselection applied to all bins

    preselection  string, selection as a C++ expression
    desc          string, description of selection in TLaTeX
    '''
    self.data_nom_tnp_analyzer.set_preselection(preselection, desc)
    self.data_altsig_tnp_analyzer.set_preselection(preselection, desc)
    self.data_altbkg_tnp_analyzer.set_preselection(preselection, desc)
    self.data_altsigbkg_tnp_analyzer.set_preselection(preselection, desc)
    self.mc_nom_tnp_analyzer.set_preselection(preselection, desc)
    self.mc_alt_tnp_analyzer.set_preselection(preselection, desc)

  #def add_nd_binning(self,dimensions):
  #  '''
  #  Sets even n-dimensional binning for TnP analysis

  #  dimensions  list of tuples generated by make_nd_bin_dimension
  #  '''
  #  self.data_nom_tnp_analyzer.add_nd_binning(dimensions)
  #  self.data_altsig_tnp_analyzer.add_nd_binning(dimensions)
  #  self.data_altbkg_tnp_analyzer.add_nd_binning(dimensions)
  #  self.data_altsigbkg_tnp_analyzer.add_nd_binning(dimensions)
  #  self.mc_nom_tnp_analyzer.add_nd_binning(dimensions)
  #  self.mc_alt_tnp_analyzer.add_nd_binning(dimensions)

  def add_custom_binning(self, bin_selections, bin_names, is_high_pt):
    '''
    Creates custom bins for TnP analysis

    bin_selections  list of strings, describes selection for each bin
    bin_names       list of strings, names of each bin that appear in plots
    is_high_pt      list of bools, indicates whether each bin is high pT
    '''
    self.data_nom_tnp_analyzer.add_custom_binning(bin_selections, bin_names)
    self.data_altsig_tnp_analyzer.add_custom_binning(bin_selections, bin_names)
    self.data_altbkg_tnp_analyzer.add_custom_binning(bin_selections, bin_names)
    self.data_altsigbkg_tnp_analyzer.add_custom_binning(bin_selections, bin_names)
    self.mc_nom_tnp_analyzer.add_custom_binning(bin_selections, bin_names)
    self.mc_alt_tnp_analyzer.add_custom_binning(bin_selections, bin_names)
    for ibin in range(len(is_high_pt)):
      if is_high_pt[ibin]:
        self.highpt_bins.append(ibin)

  def add_models(self):
    '''
    Adds standard models and parameter initializers for fitting
    '''
    #self.data_nom_tnp_analyzer.add_model('mc_cms',make_model_initializer(
    #    model_initializer_nom_egm_meta, self.mc_nom_tnp_analyzer, 
    #    self.highpt_bins))
    self.data_nom_tnp_analyzer.add_model('mc_cheby',make_model_initializer(
        model_initializer_cheby_egm_meta, self.mc_nom_tnp_analyzer, 
        self.highpt_bins))
    self.data_altbkg_tnp_analyzer.add_model('mc_gamma',make_model_initializer(
        model_initializer_gamma_egm_meta, self.mc_nom_tnp_analyzer, 
        self.highpt_bins))
    self.data_altsig_tnp_analyzer.add_model('dscb_cheby',model_initializer_cheby_dscb)
    self.data_altsigbkg_tnp_analyzer.add_model('dscb_gamma',model_initializer_cheby_dscb)
    self.mc_nom_tnp_analyzer.add_model('dscb_cheby',model_initializer_cheby_dscb)
    #MC fit is only to constrain DSCB in data, so no fit to alt MC needed
    self.mc_nom_tnp_analyzer.add_param_initializer('zero_bkg',param_initializer_zero_bkg)
    self.data_altsig_tnp.add_param_initializer('dscb_init',
        make_param_initializer_dscb_from_mc(self.mc_nom_tnp_analyzer))
    self.data_altsigbkg_tnp.add_param_initializer('dscb_init',
        make_param_initializer_dscb_from_mc(self.mc_nom_tnp_analyzer))

  def produce_histograms(self):
    '''
    Produce histograms. Only performs one loop over data files for efficiency
    '''
    if not os.path.isdir('out'):
      os.mkdir('out')
    nomdat_name = self.name+'_data_nom'
    altsig_name = self.name+'_data_altsig'
    altbkg_name = self.name+'_data_altbkg'
    altsnb_name = self.name+'_data_altsigbkg'
    nomsim_name = self.name+'_mc_nom'
    altsim_name = self.name+'_mc_alt'
    if os.path.isdir('out/'+altsig_name) or
       os.path.isdir('out/'+altbkg_name) or
       os.path.isdir('out/'+nomdat_name) or
       os.path.isdir('out/'+nomsim_name) or
       os.path.isdir('out/'+altsim_name) or
       os.path.isdir('out/'+altsnb_name):
      raise RuntimeError('Output directories already exist, aborting to avoid overwrites.')
    self.data_nom_tnp_analyzer.produce_histograms()
    os.mkdir('out/'+altsig_name)
    os.mkdir('out/'+altbkg_name)
    os.mkdir('out/'+altsnb_name)
    os.system('cp out/'+nomdat_name+'/'+nomdat_name+'.root out/'+altsig_name+'/')
    os.system('cp out/'+nomdat_name+'/'+nomdat_name+'.root out/'+altbkg_name+'/')
    os.system('cp out/'+nomdat_name+'/'+nomdat_name+'.root out/'+alttnp_name+'/')
    self.mc_nom_tnp_analyzer.produce_histograms()
    self.mc_alt_tnp_analyzer.produce_histograms()

  def run_interactive(self):
    '''
    Run an interactive T&P analysis
    '''
    #TODO left off, need to add tools for aggregating results, make some instructions for using the interactive fitters, and test
    pass

