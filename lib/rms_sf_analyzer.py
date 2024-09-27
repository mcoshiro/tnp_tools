#!/usr/bin/env python3
"""@package docstring
T&P meta analyzer that follows standard EGM procedures to generate scale factors. In particular:
"""

from array import array
from correctionlib import schemav2 
import os
import ROOT
import statistics

from tnp_analyzer import *
from tnp_utils import CMS_COLORS
from model_initializers import *
from root_plot_lib import RplPlot

LUMI_TAGS = {'2016APV' : [(20,13)],
             '2016' : [(17,13)],
             '2017' : [(41,13)],
             '2018' : [(60,13)],
             '2022' : [(8,13.6)],
             '2022EE' : [(27,13.6)],
             '2023' : [(18,13.6)],
             '2023BPix' : [(10,13.6)]}

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
  hist_name = 'hist_'+pass_fail+'_bin{}'.format(ibin)
  with ROOT.TFile(mc_filename,'READ') as input_file:
    hist = input_file.Get(hist_name).Clone()
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
  hist_name = 'hist_'+pass_fail+'_bin{}'.format(ibin)
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
  hist_name = 'hist_'+pass_fail+'_bin{}'.format(ibin)
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
  hist_name = 'hist_'+pass_fail+'_bin{}'.format(ibin)
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

def make_model_initializer(mc_analyzer, highpt_bins, meta_initializer):
  '''
  creates model_initializer function lambda from a meta initializer

  mc_analyzer      TnpAnalyzer for MC samples
  highpt_bins      list of ints, bins where passing histogram is always used
  meta_initializer function returning model with extra parameters
  '''
  return (lambda fit_var, ibin, is_pas : meta_initializer(
      fit_var, ibin, is_pass, mc_analyzer, highpt_bins))

def add_gap_eta_bins(original_bins):
  '''
  Modifies eta binning to include EB-EE gap bins (-1.566,-1.4442) and 
  (1.4442,1.566). Returns a tuple (new_bins, -gap_index, +gap_index)

  original_bins   sorted list of floats, must have an entry in the ranges
                  specified above
  '''
  new_bins = original_bins
  neg_gap_location = -1
  pos_gap_location = -1
  for i in range(len(original_bins)-1):
    if new_bins[i]>-1.566 and new_bins[i+1]<-1.4442:
      neg_gap_location = i
    if new_bins[i]>1.4442 and new_bins[i+1]<1.566:
      pos_gap_location = i
  if neg_gap_location==-1 or pos_gap_location==-1:
    raise ValueError('Input binning must have borders in gap region.')
  new_bins.insert(neg_gap_location,-1.566)
  new_bins[neg_gap_location+1] = -1.4442
  new_bins.insert(pos_gap_location+1,1.4442)
  new_bins[pos_gap_location+2] = 1.566
  return (new_bins, neg_gap_location, pos_gap_location+1)


def make_data_mc_graph(x, ex, data_y, data_ey, sim_y, sim_ey, name, data_names, 
                       mc_names, x_title, y_title, lumi):
  '''
  Makes a nice plot with multiple graphs overlayed. Each data graph will be
  drawn in solid line while each simulation graph in dotted lines. The number
  of y/ey points should match

  x           list of floats, x values for points
  ex          list of floats, x error bars for points
  data_y      list of list of floats, y values for data points
  data_ey     list of list of floats, y error bars for data points
  sim_y       list of list of floats, y values for simulation points
  sim_ey      list of list of floats, y error bars for simulation points
  name        string filename
  data_names  list of string names for data graphs
  mc_names    list of string names for mc graphs
  x_title     X-axis label
  y_title     y-axis label
  lumi        list of tuple of two floats representing lumi and CM energy
  '''
  ROOT.gStyle.SetOptStat(0)
  x_vals = array('d',x)
  ex_vals = array('d',ex)
  data_y_vals = []
  data_ey_vals = []
  sim_y_vals = []
  sim_ey_vals = []
  graphs = []
  for idata in range(len(data_y)):
    data_y_vals.append(array.array('d',data_y[idata]))
    data_ey_vals.append(array.array('d',data_ey[idata]))
    graphs.append(ROOT.TGraphErrors(len(x),x_vals,data_y_vals[idata],ex_vals,
                                    data_ey_vals[idata]))
    graphs[-1].SetTitle(data_names[idata])
    graphs[-1].SetLineStyle(ROOT.kSolid)
    graphs[-1].SetLineColor(CMS_COLORS[idata])
  for isim in range(len(sim_y)):
    sim_y_vals.append(array.array('d',sim_y[isim]))
    sim_ey_vals.append(array.array('d',sim_ey[isim]))
    graphs.append(ROOT.TGraphErrors(len(x),x_vals,sim_y_vals[isim],ex_vals,
                                    sim_ey_vals[isim]))
    graphs[-1].SetTitle(mc_names[imc])
    graphs[-1].SetLineStyle(ROOT.kDashed)
    graphs[-1].SetLineColor(CMS_COLORS[isim])
  sf_plot = RplPlot()
  sf_plot.lumi_data = lumi
  for graph in graphs
    sf_plot.plot_graph(graph, graph.GetColor())
  sf_plot.x_title = x_title
  sf_plot.y_title = y_title
  sf_plot.draw(name)

def make_sf_graph(x, ex, y, ey, name, graph_names, x_title, y_title, lumi):
  '''
  Makes a nice plot with multiple graphs overlayed. 

  x           list of floats, x values for points
  ex          list of floats, x error bars for points
  y           list of list of floats, y values  for points
  ey          list of list of floats, y error bars for points
  name        string filename
  graph_names list of string names for graphs
  x_title     X-axis label
  y_title     y-axis label
  lumi        list of tuple of two floats representing lumi and CM energy
  '''
  ROOT.gStyle.SetOptStat(0)
  x_vals = array('d',x)
  ex_vals = array('d',ex)
  y_vals = []
  ey_vals = []
  graphs = []
  for idata in range(len(y)):
    y_vals.append(array('d',y[idata]))
    ey_vals.append(array('d',ey[idata]))
    graphs.append(ROOT.TGraphErrors(len(x),x_vals,y_vals[idata],
                                    ex_vals,ey_vals[idata]))
    graphs[-1].SetTitle(graph_names[idata])
    graphs[-1].SetLineColor(CMS_COLORS[idata])
  sf_plot = RplPlot()
  sf_plot.lumi_data = lumi
  for graph in graphs
    sf_plot.plot_graph(graph, graph.GetColor())
  sf_plot.x_title = x_title
  sf_plot.y_title = y_title
  sf_plot.draw(name)

def make_heatmap(x, y, z, name, x_title, y_title, z_title, lumi):
  '''
  Makes a heatmap (2D histogram/colz)

  x           list of floats, x axis bin divisions
  y           list of floats, y axis bin divisions
  z           list of floats, heatmap values
  name        string filename
  x_title     X-axis label
  y_title     y-axis label
  z_title     z-axis label
  lumi        list of tuple of two floats representing lumi and CM energy
  '''
  x_bins = array('d',x)
  y_bins = array('d',y)
  hist = ROOT.TH2D('heatmap',';'+x_title+';'+y_title+';'+z_title,
                   len(x)-1,x_bins,len(y)-1,y_bins)
  for ix in range(len(x)-1):
    for iy in range(len(y)-1):
      hist.SetBinContent(ix+1, iy+1, z[ix][iy])
  sf_plot = RplPlot()
  sf_plot.lumi_data = lumi
  sf_plot.plot_colormap(hist)
  sf_plot.draw(name)

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
    self.binning_type = 'custom'

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
                           var_range=(60.0,120.0), custom_mc_range=(50.0,130.0), weight='1'):
    '''
    Adds information about fitting variable

    name             string, name of branch in TTree or C++ expression
    description      string, name used in plots (with TLaTeX)
    nbins            int, number of bins to use for fit variable
    var_range        tuple of two floats, start and end of fit range
    custom_mc_range  tuple of two floats, start and end of MC histogram range
                     (may be different from var_range for convolutions)
    weight           string, expression for weight to use
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
                                                  custom_mc_range, weight)
    self.mc_alt_tnp_analyzer.set_fitting_variable(name, description, nbins, 
                                                  custom_mc_range, weight)
    self.mc_nom_tnp_analyzer.set_custom_fit_range(var_range)
    self.mc_alt_tnp_analyzer.set_custom_fit_range(var_range)

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

  def set_preselection(self, preselection_data, preselection_mc, desc):
    '''
    Sets basic preselection applied to all bins

    preselection_data  string, selection for data as a C++ expression
    preselection_mc    string, selection for MC as a C++ expression
    desc               string, description of selection in TLaTeX
    '''
    self.data_nom_tnp_analyzer.set_preselection(preselection_data, desc)
    self.data_altsig_tnp_analyzer.set_preselection(preselection_data, desc)
    self.data_altbkg_tnp_analyzer.set_preselection(preselection_data, desc)
    self.data_altsigbkg_tnp_analyzer.set_preselection(preselection_data, desc)
    self.mc_nom_tnp_analyzer.set_preselection(preselection_mc, desc)
    self.mc_alt_tnp_analyzer.set_preselection(preselection_mc, desc)

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

  def add_standard_gap_binning(self, pt_bins, eta_bins, gap_pt_bins, pt_var_name, eta_var_name):
    '''
    Creates standard binning including specialized pt bins for gap region

    pt_bins       list of floats, pt bin edges
    eta_bins      list of floats, eta bin edges
    gap_pt_bins   list of floats, gap region pt bin edges
    pt_var_name   string, name of pt variable
    eta_var_name  string, name of eta variable
    '''
    bin_selections = []
    bin_names = []
    is_high_pt = []
    for ipt in range(len(pt_bins)-1):
      for ieta in range(len(pt_bins)-1):
        bin_selections.append('{}<{}&&{}<{}&&{}<{}&&{}<{}'.format(
            pt_bins[ipt],pt_var_name,pt_var_name,pt_bins[ipt+1],
            eta_bins[ieta],eta_var_name,eta_var_name,eta_bins[ieta+1]))
        bin_names.append('{}<p_{T}<{} GeV, {}<#eta<{}'.format(
            pt_bins[ipt],pt_bins[ipt+1],eta_bins[ieta],eta_bins[ieta+1]))
        if (pt_bins[ipt]>70.0):
          is_high_pt.append(True)
        else:
          is_high_pt.append(False)
    for ipt in range(len(gap_pt_bins)-1):
      bin_selections.append('{}<{}&&{}<{}&&-1.566<{}&&{}<-1.4442'.format(
          pt_bins[ipt],pt_var_name,pt_var_name,pt_bins[ipt+1],
          eta_var_name,eta_var_name))
      bin_names.append('{}<p_{T}<{} GeV, -1.566<#eta<-1.4442'.format(
          pt_bins[ipt],pt_bins[ipt+1],eta_bins[ieta],eta_bins[ieta+1]))
      bin_selections.append('{}<{}&&{}<{}&&1.4442<{}&&{}<1.566'.format(
          pt_bins[ipt],pt_var_name,pt_var_name,pt_bins[ipt+1],
          eta_var_name,eta_var_name))
      bin_names.append('{}<p_{T}<{} GeV, 1.4442<#eta<1.566'.format(
          pt_bins[ipt],pt_bins[ipt+1],eta_bins[ieta],eta_bins[ieta+1]))
      if (pt_bins[ipt]>70.0):
        is_high_pt.append(True)
        is_high_pt.append(True)
      else:
        is_high_pt.append(False)
        is_high_pt.append(False)
    self.binning_type = 'std_gap'
    self.add_custom_binning(bin_selections, bin_names, gap_pt_bins)
    self.pt_bins = pt_bins
    self.eta_bins = eta_bins
    self.gap_pt_bins = gap_pt_bins

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

  def generate_output(self):
    '''
    Generate final SFs and histograms
    '''

    #TODO left off, need to add tools for aggregating results, make some instructions for using the interactive fitters, and test

    #do checks, generate individual outputs if not already generated
    if (self.binning_type != 'std_gap'):
      print('ERROR: only standard gap binning currently supported for ',end='')
      print('automatic output formatting')
      return
    nomdat_name = self.name+'_data_nom'
    altsig_name = self.name+'_data_altsig'
    altbkg_name = self.name+'_data_altbkg'
    altsnb_name = self.name+'_data_altsigbkg'
    nomsim_name = self.name+'_mc_nom'
    altsim_name = self.name+'_mc_alt'
    if not os.path.isfile('out/'+nomdat_name+'/efficiencies.json'):
      self.data_nom_tnp_analyzer.generate_final_output()
    if not os.path.isfile('out/'+altsig_name+'/efficiencies.json'):
      self.data_altsig_tnp_analyzer.generate_final_output()
    if not os.path.isfile('out/'+altbkg_name+'/efficiencies.json'):
      self.data_altbkg_tnp_analyzer.generate_final_output()
    if not os.path.isfile('out/'+altsnb_name+'/efficiencies.json'):
      self.data_altsigbkg_tnp_analyzer.generate_final_output()
    if not os.path.isfile('out/'+nomsim_name+'/cnc_efficiencies.json'):
      self.mc_nom_tnp_analyzer.generate_cut_and_count_output()
    if not os.path.isfile('out/'+altsim_name+'/cnc_efficiencies.json'):
      self.mc_alt_tnp_analyzer.generate_cut_and_count_output()
    if ((not os.path.isfile('out/'+nomdat_name+'/efficiencies.json')) or:
        (not os.path.isfile('out/'+altsig_name+'/efficiencies.json')) or:
        (not os.path.isfile('out/'+altbkg_name+'/efficiencies.json')) or:
        (not os.path.isfile('out/'+altsnb_name+'/efficiencies.json')) or:
        (not os.path.isfile('out/'+nomsim_name+'/cnc_efficiencies.json')) or:
        (not os.path.isfile('out/'+altsim_name+'/cnc_efficiencies.json'))):
      print('ERROR: Could not generate individual outputs, please ensure '+
            'all fits have been performed.')
      return

    #get efficiencies from JSON files
    eff_dat_nom = []
    eff_alt_sig = []
    eff_alt_bkg = []
    eff_alt_snb = []
    eff_sim_nom = []
    eff_sim_alt = []
    with open('out/'+nomdat_name+'/efficiencies.json','r') as input_file:
      eff_dat_nom = json.loads(input_file.read())
    with open('out/'+altsig_name+'/efficiencies.json','r') as input_file:
      eff_alt_sig = json.loads(input_file.read())
    with open('out/'+altbkg_name+'/efficiencies.json','r') as input_file:
      eff_alt_bkg = json.loads(input_file.read())
    with open('out/'+altsnb_name+'/efficiencies.json','r') as input_file:
      eff_alt_snb = json.loads(input_file.read())
    with open('out/'+nomsim_name+'/efficiencies.json','r') as input_file:
      eff_sim_nom = json.loads(input_file.read())
    with open('out/'+altsim_name+'/efficiencies.json','r') as input_file:
      eff_sim_alt = json.loads(input_file.read())

    #calculate interesting quantities
    data_eff = []
    data_unc = []
    mc_eff = []
    mc_unc = []
    pass_sf = []
    pass_unc = []
    fail_sf = []
    fail_unc = []
    for ibin in range(len(eff_dat_nom)):
      eff_dat1 = eff_dat_nom[ibin][0]
      eff_dat2 = eff_alt_sig[ibin][0]
      eff_dat3 = eff_alt_bkg[ibin][0]
      eff_dat4 = eff_alt_snb[ibin][0]
      eff_sim1 = eff_sim_nom[ibin][0]
      eff_sim2 = eff_sim_nom[ibin][0]
      unc_dat1 = eff_dat_nom[ibin][1]
      unc_sim1 = eff_sim_nom[ibin][1]

      eff_dat_mean = statistics.mean((eff_dat1, eff_dat2, eff_dat3, eff_dat4))
      eff_dat_unc = math.hypot(unc_dat[1], math.hypot((eff_dat1-eff_dat_mean),
          (eff_dat2-eff_dat_mean),(eff_dat3-eff_dat_mean),
          (eff_dat4-eff_dat_mean))/math.sqrt(12.0))
      data_eff.append(eff_dat_mean)
      data_unc.append(eff_dat_unc)
      
      mc_eff.append(eff_sim1)
      mc_unc.append(max(unc_sim1,abs(eff_sim1-eff_sim2)))

      if (eff_sim1<=0.0 or eff_sim2<=0.0):
        print('WARNING: Nonpositive efficiency')
        pass_sf.append(1.0)
        pass_unc.append(1.0)
      else:
        sfp1 = eff_dat1/eff_sim1
        sfp2 = eff_dat2/eff_sim1
        sfp3 = eff_dat3/eff_sim1
        sfp4 = eff_dat4/eff_sim1
        sfpm = statistics.mean([sfp1,sfp2,sfp3,sfp4])
        sfprms = math.hypot(sfp1-sfpm,sfp2-sfpm,sfp3-sfpm,sfp4-sfpm)/math.sqrt(3.0)
        sfpdstat = sfp1*unc_dat1/eff_dat1
        sfpmstat = sfp1*unc_sim1/eff_sim1
        sfpmcalt = abs(sfp1-eff_dat1/eff_sim2)
        pass_sf.append(sfpm)
        pass_unc.append(math.hypot(sfprms/math.sqrt(4.0),sfpdstat,max(sfpmstat,sfpmcalt)))

      if (eff_sim1>=1.0 or eff_sim2>=1.0):
        print('WARNING: Unity efficiency')
        fail_sf.append(1.0)
        fail_unc.append(1.0)
      else:
        sff1 = (1.0-eff_dat1)/(1.0-eff_sim1)
        sff2 = (1.0-eff_dat2)/(1.0-eff_sim1)
        sff3 = (1.0-eff_dat3)/(1.0-eff_sim1)
        sff4 = (1.0-eff_dat4)/(1.0-eff_sim1)
        sffm = statistics.mean([sfp1,sfp2,sfp3,sfp4])
        sffrms = math.hypot(sff1-sffm,sff2-sffm,sff3-sffm,sff4-sffm)/math.sqrt(3.0)
        sffdstat = sff1*unc_dat1/(1.0-eff_dat1)
        sffmstat = sff1*unc_sim1/(1.0-eff_sim1)
        sffmcalt = abs(sff1-(1.0-eff_dat1)/(1.0-eff_sim2))
        fail_sf.append(sffm)
        fail_unc.append(math.hypot(sffrms/math.sqrt(4.0),sffdstat,max(sffmstat,sffmcalt)))

    #organize SFs as they will be saved in the JSON
    gapincl_eta_bins, neg_gap_idx, pos_gap_idx = add_gap_eta_bins(self.eta_bins)
    pass_json_sfs = []
    pass_json_uns = []
    fail_json_sfs = []
    fail_json_uns = []
    for ipt in range(len(self.pt_bins)-1):
      mean_bin_pt = (self.pt_bins[ipt]+self.pt_bins[ipt+1])/2.0
      for ieta in range(len(gapincl_eta_bins)-1):
        tnp_bin = -1
        if ieta < neg_gap_idx:
          tnp_bin = ipt*len(self.eta_bins)+ieta
        elif ieta == neg_gap_idx:
          tnp_bin = (len(self.pt_bins)*len(self.eta_bins)
                     +get_bin(mean_bin_pt, self.gap_pt_bins)*2)
        elif ieta > neg_gap_idx and ieta < pos_gap_idx:
          tnp_bin = ipt*len(self.eta_bins)+(ieta-1)
        elif ieta == pos_gap_idx:
          tnp_bin = (len(self.pt_bins)*len(self.eta_bins)
                     +get_bin(mean_bin_pt, self.gap_pt_bins)*2+1)
        elif ieta > pos_gap_idx:
          tnp_bin = ipt*len(self.eta_bins)+(ieta-2)
        pass_json_sfs.append(pass_sf[tnp_bin])
        pass_json_uns.append(pass_unc[tnp_bin])
        fail_json_sfs.append(fail_sf[tnp_bin])
        fail_json_uns.append(fail_unc[tnp_bin])

    #write JSON
    clib_sfs_pass = schemav2.correction(
        name='sf_pass',
        version=1,
        inputs=[schemav2.variable(name='pt', type='real', description='pt'),
                schemav2.variable(name='eta', type='real', description='eta')],
        output=schemav2.variable(name='sf', type='real', description='data-mc sf'),
        data=schemav2.multibinning(
            nodetype='multibinning',
            inputs=['pt','eta'],
            edges=[self.pt_bins,gapincl_eta_bins],
            content=pass_json_sfs,
            flow='clamp',
            ),
        )
    clib_uns_pass = schemav2.correction(
        name='unc_pass',
        version=1,
        inputs=[schemav2.variable(name='pt', type='real', description='pt'),
                schemav2.variable(name='eta', type='real', description='eta')],
        output=schemav2.variable(name='sf', type='real', description='data-mc sf'),
        data=schemav2.multibinning(
            nodetype='multibinning',
            inputs=['pt','eta'],
            edges=[self.pt_bins,gapincl_eta_bins],
            content=pass_json_uns,
            flow='clamp',
            ),
        )
    clib_sfs_fail = schemav2.correction(
        name='sf_fail',
        version=1,
        inputs=[schemav2.variable(name='pt', type='real', description='pt'),
                schemav2.variable(name='eta', type='real', description='eta')],
        output=schemav2.variable(name='sf', type='real', description='data-mc sf'),
        data=schemav2.multibinning(
            nodetype='multibinning',
            inputs=['pt','eta'],
            edges=[self.pt_bins,gapincl_eta_bins],
            content=fail_json_sfs,
            flow='clamp',
            ),
        )
    clib_uns_fail = schemav2.correction(
        name='unc_fail',
        version=1,
        inputs=[schemav2.variable(name='pt', type='real', description='pt'),
                schemav2.variable(name='eta', type='real', description='eta')],
        output=schemav2.variable(name='sf', type='real', description='data-mc sf'),
        data=schemav2.multibinning(
            nodetype='multibinning',
            inputs=['pt','eta'],
            edges=[self.pt_bins,gapincl_eta_bins],
            content=fail_json_uns,
            flow='clamp',
            ),
        )
    with open('out/'+self.name+'_scalefactors.json','w') as output_file:
      output_file.write(fix_correctionlib_json(
        [clib_sfs_pass.json(exclude_unset=True),
         clib_uns_pass.json(exclude_unset=True),
         clib_sfs_fail.json(exclude_unset=True),
         clib_uns_fail.json(exclude_unset=True)]))
    self.generate_summary_plots(data_eff, data_unc, mc_eff, mc_unc, pass_sf, 
                                pass_unc, fail_sf, fail_unc)

  def generate_summary_plots(data_eff, data_unc, mc_eff, mc_unc, pass_sf, 
                             pass_unc, fail_sf, fail_unc):
    '''generate following plots: 1D eta efficiency plot with MC & data
                              1D pt efficiency plot with MC & data
                              gap 1D pt efficiency plot with MC & data
                              2D efficiency plot with MC
                              2D efficiency plot with data
                              1D eta SF plot pass
                              1D eta SF plot fail
                              1D pt SF plot pass
                              1D pt SF plot fail
                              gap 1D pt SF plot pass
                              gap 1D pt SF plot fail
                              2D SF plot pass
                              2D SF plot fail
                              2D SF uncertainty plot pass
                              2D SF uncertainty plot fail

    data_eff  list of data efficiencies
    data_unc  list of data uncertainties
    mc_eff    list of mc efficiencies
    mc_unc    list of mc uncertainties
    pass_sf   list of scale factors (SFs) for passing selection
    pass_unc  list of uncertainties on passing SFs
    fail_sf   list of scale factors for failing selection
    fail_unc  list of uncertainties on failing SFs
    '''
    eta_plot_x = []
    eta_plot_ex = []
    pt_plot_x = []
    pt_plot_ex = []
    gap_pt_plot_x = []
    gap_pt_plot_ex = []
    eff_eta_plot_data_y = []
    eff_eta_plot_data_ey = []
    eff_eta_plot_mc_y = []
    eff_eta_plot_mc_ey = []
    eff_pt_plot_data_y = []
    eff_pt_plot_data_ey = []
    eff_pt_plot_mc_y = []
    eff_pt_plot_mc_ey = []
    eff_gappt_plot_data_y = []
    eff_gappt_plot_data_ey = []
    eff_gappt_plot_mc_y = []
    eff_gappt_plot_mc_ey = []
    sf_eta_plot_pass_y = []
    sf_eta_plot_pass_ey = []
    sf_eta_plot_fail_y = []
    sf_eta_plot_fail_ey = []
    sf_pt_plot_pass_y = []
    sf_pt_plot_pass_ey = []
    sf_pt_plot_fail_y = []
    sf_pt_plot_fail_ey = []
    sf_gappt_plot_pass_y = []
    sf_gappt_plot_pass_ey = []
    sf_gappt_plot_fail_y = []
    sf_gappt_plot_fail_ey = []
    eta_plot_names = []
    eta_plot_data_names = []
    eta_plot_mc_names = []
    pt_plot_names = []
    pt_plot_data_names = []
    pt_plot_mc_names = []
    gappt_plot_names = []
    gappt_plot_data_names = []
    gappt_plot_mc_names = []

    for ieta in range(len(self.eta_bins)-1):
      eta_plot_x.append((self.eta_bins[ieta]+self.eta_bins[ieta+1])/2.0)
      eta_plot_ex.append((self.eta_bins[ieta]-self.eta_bins[ieta+1])/2.0)
      pt_plot_names.append('{}<#eta<{}'.format(self.eta_bins[ieta],self.eta_bins[ieta+1]))
      pt_plot_data_names.append('Data {}<#eta<{}'.format(self.eta_bins[ieta],self.eta_bins[ieta+1]))
      pt_plot_mc_names.append('MC {}<#eta<{}'.format(self.eta_bins[ieta],self.eta_bins[ieta+1]))
      eff_pt_plot_data_y.append([])
      eff_pt_plot_data_ey.append([])
      eff_pt_plot_mc_y.append([])
      eff_pt_plot_mc_ey.append([])
      sf_pt_plot_pass_y.append([])
      sf_pt_plot_pass_ey.append([])
      sf_pt_plot_fail_y.append([])
      sf_pt_plot_fail_ey.append([])
    for ipt in range(len(self.pt_bins)-1):
      pt_plot_x.append((self.pt_bins[ipt]+self.pt_bins[ipt+1])/2.0)
      pt_plot_ex.append((self.pt_bins[ipt]-self.pt_bins[ipt+1])/2.0)
      eta_plot_names.append('{}<p_{{T}}<{} GeV'.format(self.pt_bins[ipt],self.pt_bins[ipt+1]))
      eta_plot_data_names.append('Data {}<p_{{T}}<{} GeV'.format(self.pt_bins[ipt],self.pt_bins[ipt+1]))
      eta_plot_mc_names.append('MC {}<p_{{T}}<{} GeV'.format(self.pt_bins[ipt],self.pt_bins[ipt+1]))
      eff_eta_plot_data_y.append([])
      eff_eta_plot_data_ey.append([])
      eff_eta_plot_mc_y.append([])
      eff_eta_plot_mc_ey.append([])
      sf_eta_plot_pass_y.append([])
      sf_eta_plot_pass_ey.append([])
      sf_eta_plot_fail_y.append([])
      sf_eta_plot_fail_ey.append([])
    for ipt in range(len(self.gap_pt_bins)-1):
      gap_pt_plot_x.append((self.gap_pt_bins[ipt]+self.gap_pt_bins[ipt+1])/2.0)
      gap_pt_plot_ex.append((self.gap_pt_bins[ipt]-self.gap_pt_bins[ipt+1])/2.0)
    gappt_plot_names.append('-1.566<#eta<-1.4442'))
    gappt_plot_names.append('1.4442<#eta<1.566'))
    gappt_plot_data_names.append('Data -1.566<#eta<-1.4442'))
    gappt_plot_data_names.append('Data 1.4442<#eta<1.566'))
    gappt_plot_mc_names.append('MC -1.566<#eta<-1.4442'))
    gappt_plot_mc_names.append('MC 1.4442<#eta<1.566'))
    for ieta in range(2):
      eff_gappt_plot_data_y.append([])
      eff_gappt_plot_data_ey.append([])
      eff_gappt_plot_mc_y.append([])
      eff_gappt_plot_mc_ey.append([])
      sf_gappt_plot_pass_y.append([])
      sf_gappt_plot_pass_ey.append([])
      sf_gappt_plot_fail_y.append([])
      sf_gappt_plot_fail_ey.append([])

    for ipt in range(len(self.pt_bins)-1):
      for ieta in range(len(self.eta_bins)-1):
        tnp_bin = ieta+ipt*len(self.eta_bins)
        eff_eta_plot_data_y[ipt].append(data_eff[tnp_bin])
        eff_eta_plot_data_ey[ipt].append(data_unc[tnp_bin])
        eff_eta_plot_mc_y[ipt].append(mc_eff[tnp_bin])
        eff_eta_plot_mc_ey[ipt].append(mc_unc[tnp_bin])
        eff_pt_plot_data_y[ieta].append(data_eff[tnp_bin])
        eff_pt_plot_data_ey[ieta].append(data_unc[tnp_bin])
        eff_pt_plot_mc_y[ieta].append(mc_eff[tnp_bin])
        eff_pt_plot_mc_ey[ieta].append(mc_unc[tnp_bin])
        sf_eta_plot_pass_y[ipt].append(pass_sf[tnp_bin])
        sf_eta_plot_pass_ey[ipt].append(pass_unc[tnp_bin])
        sf_eta_plot_fail_y[ipt].append(fail_sf[tnp_bin])
        sf_eta_plot_fail_ey[ipt].append(fail_unc[tnp_bin])
        sf_pt_plot_pass_y[ieta].append(pass_sf[tnp_bin])
        sf_pt_plot_pass_ey[ieta].append(pass_unc[tnp_bin])
        sf_pt_plot_fail_y[ieta].append(fail_sf[tnp_bin])
        sf_pt_plot_fail_ey[ieta].append(fail_unc[tnp_bin])
    for ipt in range(len(self.gap_pt_bins)-1):
      for ieta in range(2):
        tnp_bin = len(self.pt_bins)*len(self.eta_bins)+ieta+ipt*2
        eff_gappt_plot_data_y[ieta].append(data_eff[tnp_bin])
        eff_gappt_plot_data_ey[ieta].append(data_unc[tnp_bin])
        eff_gappt_plot_mc_y[ieta].append(mc_eff[tnp_bin])
        eff_gappt_plot_mc_ey[ieta].append(mc_unc[tnp_bin])
        sf_gappt_plot_pass_y[ieta].append(pass_sf[tnp_bin])
        sf_gappt_plot_pass_ey[ieta].append(pass_unc[tnp_bin])
        sf_gappt_plot_fail_y[ieta].append(fail_sf[tnp_bin])
        sf_gappt_plot_fail_ey[ieta].append(fail_unc[tnp_bin])

    eff_string = 'Efficiency '+self.data_nom_tnp_analyzer.measurement_variable
    make_data_mc_graph(eta_plot_x, eta_plot_ex, eff_eta_plot_data_y, 
                       eff_eta_plot_data_ey, eff_eta_plot_mc_y, 
                       eff_eta_plot_mc_ey, 
                       'out/'+self.name+'_eff_etabinned.pdf', 
                       eta_plot_data_names, eta_plot_mc_names,
                       '#eta', eff_string, LUMI_TAGS[year])
                       )
    make_data_mc_graph(pt_plot_x, pt_plot_ex, eff_pt_plot_data_y, 
                       eff_pt_plot_data_ey, eff_pt_plot_mc_y, 
                       eff_pt_plot_mc_ey, 
                       'out/'+self.name+'_eff_ptbinned.pdf',
                       pt_plot_data_names, pt_plot_mc_names,
                       'p_{T} [GeV]', eff_string, LUMI_TAGS[year])
    make_data_mc_graph(gappt_plot_x, gappt_plot_ex, eff_gappt_plot_data_y, 
                       eff_gappt_plot_data_ey, eff_gappt_plot_mc_y, 
                       eff_gappt_plot_mc_ey, 
                       'out/'+self.name+'_eff_gapptbinned.pdf',
                       gappt_plot_data_names, gappt_plot_mc_names,
                       'p_{T} [GeV]', eff_string, LUMI_TAGS[year])
    make_heatmap(self.eta_bins, self.pt_bins, eff_pt_plot_data_y, 
                 'out/'+self.name+'_eff_data.pdf', '#eta', 'p_{T} [GeV]', 
                 eff_string, LUMI_TAGS[year])
    make_heatmap(self.eta_bins, self.pt_bins, eff_pt_plot_mc_y, 
                 'out/'+self.name+'_eff_mc.pdf', '#eta', 'p_{T} [GeV]',
                 eff_string, LUMI_TAGS[year])
    make_sf_graph(eta_plot_x, eta_plot_ex, sf_eta_plot_pass_y, 
                  sf_eta_plot_pass_ey, 
                  'out/'+self.name+'_sfpass_etabinned.pdf',
                  eta_plot_names, '#eta', 'Pass SF', LUMI_TAGS[year])
    make_sf_graph(eta_plot_x, eta_plot_ex, sf_eta_plot_fail_y, 
                  sf_eta_plot_fail_ey, 
                  'out/'+self.name+'_sffail_etabinned.pdf',
                  eta_plot_names, '#eta', 'Fail SF', LUMI_TAGS[year])
    make_sf_graph(pt_plot_x, pt_plot_ex, sf_pt_plot_pass_y, sf_pt_plot_pass_ey,
                  'out/'+self.name+'_sfpass_ptbinned.pdf',
                  pt_plot_names, 'p_{T} [GeV]', 'Pass SF', LUMI_TAGS[year])
    make_sf_graph(pt_plot_x, pt_plot_ex, sf_pt_plot_fail_y, sf_pt_plot_fail_ey, 
                  'out/'+self.name+'_sffail_ptbinned.pdf',
                  pt_plot_names, 'p_{T} [GeV]', 'Fail SF', LUMI_TAGS[year])
    make_sf_graph(gappt_plot_x, gappt_plot_ex, sf_gappt_plot_pass_y, 
                  sf_gappt_plot_pass_ey, 
                  'out/'+self.name+'_sfpass_gapptbinned.pdf',
                  gappt_plot_names, 'p_{T} [GeV]', 'Pass SF', LUMI_TAGS[year])
    make_sf_graph(gappt_plot_x, gappt_plot_ex, sf_gappt_plot_fail_y, 
                  sf_gappt_plot_fail_ey, 
                  'out/'+self.name+'_sffail_gapptbinned.pdf',
                  gappt_plot_names, 'p_{T} [GeV]', 'Fail SF', LUMI_TAGS[year])
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_pass_y, 
                 'out/'+self.name+'_sfpass.pdf', '#eta', 'p_{T} [GeV]', 
                 eff_string, LUMI_TAGS[year])
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_fail_y, 
                 'out/'+self.name+'_sffail.pdf', '#eta', 'p_{T} [GeV]',
                 eff_string, LUMI_TAGS[year])
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_pass_ey, 
                 'out/'+self.name+'_sfpass_unc.pdf', '#eta', 'p_{T} [GeV]',
                 eff_string, LUMI_TAGS[year])
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_fail_ey, 
                 'out/'+self.name+'_sffail_unc.pdf', '#eta', 'p_{T} [GeV]',
                 eff_string, LUMI_TAGS[year])

  def run_interactive(self):
    '''
    Run an interactive T&P analysis
    '''

    #run main interactive loop
    exit_loop = False
    print('Welcome to tnp_tools interactive analysis. Type [h]elp for more information.')
    while not exit_loop:
      user_input = input('>:')
      user_input = user_input.split()
      if len(user_input)<1:
        continue
      elif (user_input[0] == 'h' or user_input[0] == 'help'):
        print('Available commands:')
        print('h(elp)                         print help information')
        print('p(roduce)                      produce histograms to fit; run this first')
        print('f(it) <nom/alts/altb/altsb/mc> run interactive session (nominal/alt signal/')
        print('                               alt background/alt signal+alt background/')
        print('                               mc alt signal constraint')
        print('o(utput)                       generate final outputs')
      elif (user_input[0] == 'p' or user_input[0] == 'produce'):
        self.produce_histograms()
      elif (user_input[0] == 'f' or user_input[0] == 'fit'):
        if (len(user_input)<2):
          print('ERROR: f(it) requires an argument')
        elif (user_input[1] == 'nom'):
          self.data_nom_tnp_analyzer.run_interactive()
        elif (user_input[1] == 'alts' or user_input[1] == 'altsignal'):
          self.data_altsig_tnp_analyzer.run_interactive()
        elif (user_input[1] == 'altb' or user_input[1] == 'altbackground'):
          self.data_altbkg_tnp_analyzer.run_interactive()
        elif (user_input[1] == 'altsb' or user_input[1] == 'altsignalbackground'):
          self.data_altsigbkg_tnp_analyzer.run_interactive()
        elif (user_input[1] == 'mc' or user_input[1] == 'mcalt'):
          self.mc_nom_tnp_analyzer.run_interactive()
        else:
          print('ERROR: unrecognized argument to f(it)')
      elif (user_input[0] == 'o' or user_input[0] == 'output'):
        #TODO add calls to outputs for each individual piece
        self.generate_output()
      else:
        print('ERROR: unrecognized command')

