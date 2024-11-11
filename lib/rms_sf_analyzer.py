#!/usr/bin/env python3
"""@package docstring
T&P meta analyzer that follows standard EGM procedures to generate scale factors. In particular:
"""

from array import array
from correctionlib import schemav2 
from functools import partial
import os
import ROOT
import statistics
import json

from tnp_analyzer import *
from tnp_utils import CMS_COLORS, LUMI_TAGS
from model_initializers import *
from root_plot_lib import RplPlot

def param_initializer_dscb_from_mc(ibin, is_pass, workspace, mc_analyzer):
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
  mc_json_filename = 'out/{}/fitinfo_bin{}_{}.json'.format(
          mc_analyzer.temp_name,str(ibin),pass_fail)
  with open(mc_json_filename,'r') as mc_file:
    param_dict = json.loads(mc_file.read())
    for var in ['mean','sigmal','sigmar','alphal','nl','alphar','nr']:
      workspace.var(var).setVal(param_dict[var])
    for var in ['alphal','nl','alphar','nr']:
      workspace.var(var).setConstant()

def param_initializer_moddscb_from_mc(ibin, is_pass, workspace, mc_analyzer):
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
  mc_json_filename = 'out/{}/fitinfo_bin{}_{}.json'.format(
          mc_analyzer.temp_name,str(ibin),pass_fail)
  with open(mc_json_filename,'r') as mc_file:
    param_dict = json.loads(mc_file.read())
    for var in ['mean','sigmal','sigmar','alphal','nl1','nl2','fl','alphar',
                'nr1','nr2','fr']:
      workspace.var(var).setVal(param_dict[var])
    for var in ['alphal','nl1','nl2','fl','alphar','nr1','nr2','fr']:
      workspace.var(var).setConstant()

def param_initializer_dscbgaus_from_mc(ibin, is_pass, workspace, mc_analyzer):
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
  mc_json_filename = 'out/{}/fitinfo_bin{}_{}.json'.format(
          mc_analyzer.temp_name,str(ibin),pass_fail)
  with open(mc_json_filename,'r') as mc_file:
    param_dict = json.loads(mc_file.read())
    for var in ['mu','sigma','alphal','nl','alphar','nr','gauss_mu',
                'gauss_sigma','gauss_frac']:
      workspace.var(var).setVal(param_dict[var])
    for var in ['alphal','nl','alphar','nr','gauss_mu','gauss_sigma',
                'gauss_frac']:
      workspace.var(var).setConstant()

def param_initializer_cbconvgen_from_mc(ibin, is_pass, workspace, mc_analyzer):
  '''
  Parameter initializer for cbconvgen model that fixes CB parameters to MC
  result
  
  ibin         int, bin number
  is_pass      bool, indicates if passing leg
  workspace    RooWorkspace for this bin
  mc_analyzer  TnpAnalyzer for MC samples
  '''
  pass_fail = 'pass'
  if not is_pass:
    pass_fail = 'fail'
  mc_json_filename = 'out/{}/fitinfo_bin{}_{}.json'.format(
          mc_analyzer.temp_name,str(ibin),pass_fail)
  with open(mc_json_filename,'r') as mc_file:
    param_dict = json.loads(mc_file.read())
    for var in ['m0','sigma','alpha','n','sigma_2','tailLeft']:
      workspace.var(var).setVal(param_dict[var])
    for var in ['alpha','n','sigma_2','tailLeft']:
      workspace.var(var).setConstant()

def get_mc_histogram(ibin, is_pass, mc_analyzer, highpt_bins):
  '''Helper function used to get appropriate TH1D from analyzer
  '''
  pass_fail = 'pass'
  if not is_pass and not (ibin in highpt_bins):
    pass_fail = 'fail'
  mc_filename = 'out/'+mc_analyzer.temp_name+'/'+mc_analyzer.temp_name+'.root'
  hist_name = 'hist_{}_bin{}'.format(pass_fail,ibin)
  input_file = ROOT.TFile(mc_filename,'READ')
  hist = input_file.Get(hist_name)
  #deal with rare case of empty histogram
  if hist.Integral()<=0.0:
    if pass_fail=='pass':
      hist_name = 'hist_fail_bin{}'.format(ibin)
      hist = input_file.Get(hist_name)
    else:
      hist_name = 'hist_pass_bin{}'.format(ibin)
      hist = input_file.Get(hist_name)
  hist.SetDirectory(ROOT.nullptr)
  input_file.Close()
  return hist

def add_gap_eta_bins(original_bins):
  '''
  Modifies eta binning to include EB-EE gap bins (-1.566,-1.4442) and 
  (1.4442,1.566). Returns a tuple (new_bins, -gap_index, +gap_index)

  original_bins   sorted list of floats, must have an entry in the ranges
                  specified above
  '''
  new_bins = original_bins.copy()
  neg_gap_location = -1
  pos_gap_location = -1
  for i in range(len(original_bins)-1):
    if new_bins[i]>-1.566 and new_bins[i]<-1.4442:
      neg_gap_location = i
    if new_bins[i]>1.4442 and new_bins[i]<1.566:
      pos_gap_location = i
  if neg_gap_location==-1 or pos_gap_location==-1:
    raise ValueError('Input binning must have borders in gap region.')
  new_bins.insert(neg_gap_location,-1.566)
  new_bins[neg_gap_location+1] = -1.4442
  new_bins.insert(pos_gap_location+1,1.4442)
  new_bins[pos_gap_location+2] = 1.566
  return (new_bins, neg_gap_location, pos_gap_location+1)

def calculate_sfs(eff_dat1, eff_dat2, eff_dat3, eff_dat4, 
                  eff_sim1, eff_sim2, unc_dat1, unc_sim1,
                  unc_sim2):
  '''Calculates scale factors from efficiencies and associated uncertainties 
  using RMS method

  eff_dat1  data efficiency measurement 1
  eff_dat2  data efficiency measurement 2
  eff_dat3  data efficiency measurement 3
  eff_dat4  data efficiency measurement 4
  eff_sim1  simulation efficiency measurement 1
  eff_sim2  simulation efficiency measurement 2
  unc_dat1  data efficiency 1 uncertainty
  unc_sim1  simulation efficiency 1 uncertainty
  unc_sim2  simulation efficiency 2 uncertainty

  returns (scale factor for passing events,
           associated uncertainty,
           scale factor for failing events,
           associated uncertainty)
  '''
  pass_sf = 1.0
  pass_unc = 1.0
  fail_sf = 1.0
  fail_unc = 1.0
  if not (eff_sim1<=0.0 or eff_sim2<=0.0):
    sfp1 = eff_dat1/eff_sim1
    sfp2 = eff_dat2/eff_sim1
    sfp3 = eff_dat3/eff_sim1
    sfp4 = eff_dat4/eff_sim1
    sfpm = statistics.mean([sfp1,sfp2,sfp3,sfp4])
    sfprms = math.hypot(sfp1-sfpm,sfp2-sfpm,sfp3-sfpm,sfp4-sfpm)/math.sqrt(3.0)
    sfpdstat = sfp1*unc_dat1/eff_dat1
    sfpmstat = sfp1*unc_sim1/eff_sim1
    sfpmcalt = abs(sfp1-eff_dat1/eff_sim2)
    pass_sf = sfpm
    pass_unc = math.hypot(sfprms/math.sqrt(4.0),sfpdstat,max(sfpmstat,sfpmcalt))
  elif (eff_sim1 > 0.0 or eff_sim2 > 0.0):
    nonzero_eff_sim = eff_sim1
    nonzero_unc_sim = unc_sim1
    if eff_sim2 > 0.0:
      nonzero_eff_sim = eff_sim2
      nonzero_unc_sim = unc_sim2
    sfp1 = eff_dat1/nonzero_eff_sim
    sfp2 = eff_dat2/nonzero_eff_sim
    sfp3 = eff_dat3/nonzero_eff_sim
    sfp4 = eff_dat4/nonzero_eff_sim
    sfpm = statistics.mean([sfp1,sfp2,sfp3,sfp4])
    sfprms = math.hypot(sfp1-sfpm,sfp2-sfpm,sfp3-sfpm,sfp4-sfpm)/math.sqrt(3.0)
    sfpdstat = sfp1*unc_dat1/eff_dat1
    sfpmstat = sfp1*nonzero_unc_sim/nonzero_eff_sim
    pass_sf = sfpm
    pass_unc = math.hypot(sfprms/math.sqrt(4.0),sfpdstat,sfpmstat)
  else:
    print('WARNING: zero efficiency found')
    print(eff_sim1)
    print(eff_sim2)

  if not (eff_sim1>=1.0 or eff_sim2>=1.0):
    sff1 = (1.0-eff_dat1)/(1.0-eff_sim1)
    sff2 = (1.0-eff_dat2)/(1.0-eff_sim1)
    sff3 = (1.0-eff_dat3)/(1.0-eff_sim1)
    sff4 = (1.0-eff_dat4)/(1.0-eff_sim1)
    sffm = statistics.mean([sff1,sff2,sff3,sff4])
    sffrms = math.hypot(sff1-sffm,sff2-sffm,sff3-sffm,sff4-sffm)/math.sqrt(3.0)
    sffdstat = sff1*unc_dat1/(1.0-eff_dat1)
    sffmstat = sff1*unc_sim1/(1.0-eff_sim1)
    sffmcalt = abs(sff1-(1.0-eff_dat1)/(1.0-eff_sim2))
    fail_sf = sffm
    fail_unc = math.hypot(sffrms/math.sqrt(4.0),sffdstat,max(sffmstat,sffmcalt))
  elif (eff_sim1 < 1.0 or eff_sim2 < 1.0):
    nonunity_eff_sim = eff_sim1
    nonunity_unc_sim = unc_sim1
    if eff_sim2 < 1.0:
      nonunity_eff_sim = eff_sim2
      nonunity_unc_sim = unc_sim2
    sff1 = (1.0-eff_dat1)/(1.0-nonunity_eff_sim)
    sff2 = (1.0-eff_dat2)/(1.0-nonunity_eff_sim)
    sff3 = (1.0-eff_dat3)/(1.0-nonunity_eff_sim)
    sff4 = (1.0-eff_dat4)/(1.0-nonunity_eff_sim)
    sffm = statistics.mean([sff1,sff2,sff3,sff4])
    sffrms = math.hypot(sff1-sffm,sff2-sffm,sff3-sffm,sff4-sffm)/math.sqrt(3.0)
    sffdstat = sff1*unc_dat1/(1.0-eff_dat1)
    sffmstat = sff1*nonunity_unc_sim/(1.0-nonunity_eff_sim)
    fail_sf = sffm
    fail_unc = math.hypot(sffrms/math.sqrt(4.0),sffdstat,sffmstat)
  else:
    print('WARNING: unity efficiency found')
  return pass_sf, pass_unc, fail_sf, fail_unc

def make_data_mc_graph(x, ex, data_y, data_ey, sim_y, sim_ey, name, data_names, 
                       mc_names, x_title, y_title, lumi, log_x=False):
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
  log_x       bool, if true makes x-axis logarithmic
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
    data_y_vals.append(array('d',data_y[idata]))
    data_ey_vals.append(array('d',data_ey[idata]))
    graphs.append(ROOT.TGraphErrors(len(x),x_vals,data_y_vals[idata],ex_vals,
                                    data_ey_vals[idata]))
    graphs[-1].SetTitle(data_names[idata])
    graphs[-1].SetLineStyle(ROOT.kSolid)
    graphs[-1].SetLineColor(CMS_COLORS[idata])
  for isim in range(len(sim_y)):
    sim_y_vals.append(array('d',sim_y[isim]))
    sim_ey_vals.append(array('d',sim_ey[isim]))
    graphs.append(ROOT.TGraphErrors(len(x),x_vals,sim_y_vals[isim],ex_vals,
                                    sim_ey_vals[isim]))
    graphs[-1].SetTitle(mc_names[isim])
    graphs[-1].SetLineStyle(ROOT.kDashed)
    graphs[-1].SetLineColor(CMS_COLORS[isim])
  sf_plot = RplPlot()
  sf_plot.lumi_data = lumi
  for graph in graphs:
    sf_plot.plot_graph(graph, graph.GetLineColor())
  sf_plot.x_title = x_title
  sf_plot.y_title = y_title
  sf_plot.log_x = log_x
  sf_plot.draw(name)

def make_correction(name, desc, pt_bins, eta_bins, content):
  '''generates correctionlib correction object

  name     string correction name
  desc     string correction description
  pt_bins  list of floats, bin edges
  eta_bins list of floats, bin edges
  content  list of floats, content
  '''
  return schemav2.Correction(
      name=name,
      version=1,
      inputs=[schemav2.Variable(name='pt', type='real', description='pt'),
              schemav2.Variable(name='eta', type='real', description='eta')],
      output=schemav2.Variable(name='sf', type='real', description=desc),
      data=schemav2.MultiBinning(
          nodetype='multibinning',
          inputs=['pt','eta'],
          edges=[pt_bins, eta_bins],
          content=content,
          flow='clamp',
          ),
      )

def make_sf_graph(x, ex, y, ey, name, graph_names, x_title, y_title, lumi,
                  log_x=False):
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
  log_x       boolean, if true sets x-axis to be logarithmic
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
  for graph in graphs:
    sf_plot.plot_graph(graph, graph.GetLineColor())
  sf_plot.x_title = x_title
  sf_plot.y_title = y_title
  sf_plot.log_x = log_x
  sf_plot.draw(name)

def make_heatmap(x, y, z, name, x_title, y_title, z_title, lumi, log_x=False, log_y=False):
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
  log_x       boolean, if true sets x-axis to be logarithmic
  log_y       boolean, if true sets y-axis to be logarithmic
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
  sf_plot.log_x = log_x
  sf_plot.log_y = log_y
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
    self.mc_nom_tnp_analyzer = TnpAnalyzer(name+'_mc_nom')
    self.mc_alt_tnp_analyzer = TnpAnalyzer(name+'_mc_alt')
    self.binning_type = 'custom'
    self.year = '2016APV'

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
    
  def set_fitting_variable(self, name, description, nbins=60, nbins_mc=80,
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
    self.mc_nom_tnp_analyzer.set_fitting_variable(name, description, nbins_mc, 
                                                  custom_mc_range, weight)
    self.mc_alt_tnp_analyzer.set_fitting_variable(name, description, nbins_mc, 
                                                  custom_mc_range, weight)
    self.mc_nom_tnp_analyzer.set_custom_fit_range(var_range)
    self.mc_alt_tnp_analyzer.set_custom_fit_range(var_range)

  def set_measurement_variable(self, var, desc=''):
    '''
    Sets selection efficiency to measure with tag & probe

    var   string, name of branch in TTree or C++ expression
    desc  string, description o fmeasurement variable
    '''
    self.data_nom_tnp_analyzer.set_measurement_variable(var,desc)
    self.data_altsig_tnp_analyzer.set_measurement_variable(var,desc)
    self.data_altbkg_tnp_analyzer.set_measurement_variable(var,desc)
    self.data_altsigbkg_tnp_analyzer.set_measurement_variable(var,desc)
    self.mc_nom_tnp_analyzer.set_measurement_variable(var,desc)
    self.mc_alt_tnp_analyzer.set_measurement_variable(var,desc)

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
    self.highpt_bins = []
    for ibin in range(len(is_high_pt)):
      if is_high_pt[ibin]:
        self.highpt_bins.append(ibin)

  def add_standard_binning(self, pt_bins, eta_bins, pt_var_name, eta_var_name):
    '''
    Creates standard pt-eta binning 

    pt_bins       list of floats, pt bin edges
    eta_bins      list of floats, eta bin edges
    pt_var_name   string, name of pt variable
    eta_var_name  string, name of eta variable
    '''
    bin_selections = []
    bin_names = []
    is_high_pt = []
    for ipt in range(len(pt_bins)-1):
      for ieta in range(len(eta_bins)-1):
        bin_selections.append('{}<{}&&{}<{}&&{}<{}&&{}<{}'.format(
            pt_bins[ipt],pt_var_name,pt_var_name,pt_bins[ipt+1],
            eta_bins[ieta],eta_var_name,eta_var_name,eta_bins[ieta+1]))
        bin_names.append('{}<p_{{T}}<{} GeV, {}<#eta<{}'.format(
            pt_bins[ipt],pt_bins[ipt+1],eta_bins[ieta],eta_bins[ieta+1]))
        if (pt_bins[ipt]>70.0):
          is_high_pt.append(True)
        else:
          is_high_pt.append(False)
    self.binning_type = 'std'
    self.add_custom_binning(bin_selections, bin_names, is_high_pt)
    self.pt_bins = pt_bins
    self.eta_bins = eta_bins
    self.gap_pt_bins = []

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
      for ieta in range(len(eta_bins)-1):
        bin_selections.append('{}<{}&&{}<{}&&{}<{}&&{}<{}'.format(
            pt_bins[ipt],pt_var_name,pt_var_name,pt_bins[ipt+1],
            eta_bins[ieta],eta_var_name,eta_var_name,eta_bins[ieta+1]))
        bin_names.append('{}<p_{{T}}<{} GeV, {}<#eta<{}'.format(
            pt_bins[ipt],pt_bins[ipt+1],eta_bins[ieta],eta_bins[ieta+1]))
        if (pt_bins[ipt]>70.0):
          is_high_pt.append(True)
        else:
          is_high_pt.append(False)
    for ipt in range(len(gap_pt_bins)-1):
      bin_selections.append('{}<{}&&{}<{}&&-1.566<{}&&{}<-1.4442'.format(
          pt_bins[ipt],pt_var_name,pt_var_name,pt_bins[ipt+1],
          eta_var_name,eta_var_name))
      bin_names.append('{}<p_{{T}}<{} GeV, -1.566<#eta<-1.4442'.format(
          pt_bins[ipt],pt_bins[ipt+1],eta_bins[ieta],eta_bins[ieta+1]))
      bin_selections.append('{}<{}&&{}<{}&&1.4442<{}&&{}<1.566'.format(
          pt_bins[ipt],pt_var_name,pt_var_name,pt_bins[ipt+1],
          eta_var_name,eta_var_name))
      bin_names.append('{}<p_{{T}}<{} GeV, 1.4442<#eta<1.566'.format(
          pt_bins[ipt],pt_bins[ipt+1],eta_bins[ieta],eta_bins[ieta+1]))
      if (pt_bins[ipt]>70.0):
        is_high_pt.append(True)
        is_high_pt.append(True)
      else:
        is_high_pt.append(False)
        is_high_pt.append(False)
    self.binning_type = 'std_gap'
    self.add_custom_binning(bin_selections, bin_names, is_high_pt)
    self.pt_bins = pt_bins
    self.eta_bins = eta_bins
    self.gap_pt_bins = gap_pt_bins

  def add_models(self):
    '''
    Adds standard models and parameter initializers for fitting
    '''
    #this makes liberal use of functools.partial in order to 
    # 1. merge signal and background models
    # 2. pass meta parameters such as location of MC template histograms
    nom_signal_model = partial(add_signal_model_mcsumsmear, 
        get_histogram = partial(get_mc_histogram, 
            mc_analyzer = self.mc_nom_tnp_analyzer, 
            highpt_bins = self.highpt_bins))
    nom_signal_model_name = 'mc'
    nom_background_model = add_background_model_bernstein
    nom_background_model_name = 'bern'
    contingency_background_model = add_background_model_bernstein8
    contingency_background_model_name = 'bern8'
    alt_signal_model = add_signal_model_dscbgaus
    alt_signal_model_standalone = model_initializer_dscbgaus
    alt_signal_model_initializer = partial(param_initializer_dscbgaus_from_mc,
        mc_analyzer = self.mc_nom_tnp_analyzer)
    alt_signal_model_name = 'dscbgaus'
    alt_background_model = add_background_model_gamma
    alt_background_model_name = 'gamma'

    self.nom_fn_name = '{}_{}'.format(nom_signal_model_name, 
                                      nom_background_model_name)
    self.contingency_fn_name = '{}_{}'.format(nom_signal_model_name, 
                                              contingency_background_model_name)
    self.alts_fn_name = '{}_{}'.format(alt_signal_model_name, 
                                       nom_background_model_name)
    self.contingencyalts_fn_name = '{}_{}'.format(alt_signal_model_name, 
        contingency_background_model_name)
    self.altb_fn_name = '{}_{}'.format(nom_signal_model_name, 
                                       alt_background_model_name)
    self.altsb_fn_name = '{}_{}'.format(alt_signal_model_name, 
                                        alt_background_model_name)
    self.alts_fn_init = '{}_init'.format(alt_signal_model_name)
    self.alts_fn_name_sa = alt_signal_model_name

    self.data_nom_tnp_analyzer.add_model(self.nom_fn_name,
        partial(make_signal_background_model, 
            add_signal_model = nom_signal_model,
            add_background_model = nom_background_model))
    self.data_nom_tnp_analyzer.add_model(self.contingency_fn_name,
        partial(make_signal_background_model, 
            add_signal_model = nom_signal_model,
            add_background_model = contingency_background_model))
    self.data_altsig_tnp_analyzer.add_model(self.alts_fn_name,
        partial(make_signal_background_model, 
            add_signal_model = alt_signal_model,
            add_background_model = nom_background_model))
    self.data_altsig_tnp_analyzer.add_model(self.contingencyalts_fn_name,
        partial(make_signal_background_model, 
            add_signal_model = alt_signal_model,
            add_background_model = contingency_background_model))
    self.data_altbkg_tnp_analyzer.add_model(self.altb_fn_name,
        partial(make_signal_background_model, 
            add_signal_model = nom_signal_model,
            add_background_model = alt_background_model))
    self.data_altsigbkg_tnp_analyzer.add_model(self.altsb_fn_name,
        partial(make_signal_background_model, 
            add_signal_model = alt_signal_model,
            add_background_model = alt_background_model))
    self.mc_nom_tnp_analyzer.add_model(alt_signal_model_name,
        alt_signal_model_standalone)
    self.data_altsig_tnp_analyzer.add_param_initializer(self.alts_fn_init,
        alt_signal_model_initializer)
    self.data_altsigbkg_tnp_analyzer.add_param_initializer(self.alts_fn_init,
        alt_signal_model_initializer)

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
    if (os.path.isdir('out/'+altsig_name) or
        os.path.isdir('out/'+altbkg_name) or
        os.path.isdir('out/'+nomdat_name) or
        os.path.isdir('out/'+nomsim_name) or
        os.path.isdir('out/'+altsim_name) or
        os.path.isdir('out/'+altsnb_name)):
      raise RuntimeError('Output folders exist, aborting to avoid overwrite.')
    self.data_nom_tnp_analyzer.produce_histograms()
    os.mkdir('out/'+altsig_name)
    os.mkdir('out/'+altbkg_name)
    os.mkdir('out/'+altsnb_name)
    self.data_nom_tnp_analyzer.close_file()
    os.system('cp out/'+nomdat_name+'/'+nomdat_name+'.root '
              +'out/'+altsig_name+'/'+altsig_name+'.root')
    os.system('cp out/'+nomdat_name+'/'+nomdat_name+'.root '
              +'out/'+altbkg_name+'/'+altbkg_name+'.root')
    os.system('cp out/'+nomdat_name+'/'+nomdat_name+'.root '
              +'out/'+altsnb_name+'/'+altsnb_name+'.root')
    self.mc_nom_tnp_analyzer.produce_histograms()
    self.mc_alt_tnp_analyzer.produce_histograms()
    self.mc_nom_tnp_analyzer.close_file()
    self.mc_alt_tnp_analyzer.close_file()

  def clean_output(self):
    '''Cleans the output so efficiencies will be regenerated
    '''
    self.data_nom_tnp_analyzer.clean_output()
    self.data_altsig_tnp_analyzer.clean_output()
    self.data_altbkg_tnp_analyzer.clean_output()
    self.data_altsigbkg_tnp_analyzer.clean_output()

  def generate_individual_outputs(self):
    '''Generates individual efficiency measurements if they have not already 
    been generated and checks output files are generated correctly. Returns
    true if the outputs are generated correctly
    '''
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
    if not os.path.isfile('out/'+nomsim_name+'/efficiencies.json'):
      self.mc_nom_tnp_analyzer.generate_final_output() #just for fit plots
    if not os.path.isfile('out/'+nomsim_name+'/cnc_efficiencies.json'):
      self.mc_nom_tnp_analyzer.generate_cut_and_count_output()
    if not os.path.isfile('out/'+altsim_name+'/cnc_efficiencies.json'):
      self.mc_alt_tnp_analyzer.generate_cut_and_count_output()
    if ((not os.path.isfile('out/'+nomdat_name+'/efficiencies.json')) or
        (not os.path.isfile('out/'+altsig_name+'/efficiencies.json')) or
        (not os.path.isfile('out/'+altbkg_name+'/efficiencies.json')) or
        (not os.path.isfile('out/'+altsnb_name+'/efficiencies.json')) or
        (not os.path.isfile('out/'+nomsim_name+'/cnc_efficiencies.json')) or
        (not os.path.isfile('out/'+altsim_name+'/cnc_efficiencies.json'))):
      print('ERROR: Could not generate individual outputs, please ensure '+
            'all fits have been performed.')
      return False
    return True

  def generate_output(self):
    '''Generate final SFs and histograms
    '''

    if not self.generate_individual_outputs():
      return

    #get efficiencies from JSON files
    nomdat_name = self.name+'_data_nom'
    altsig_name = self.name+'_data_altsig'
    altbkg_name = self.name+'_data_altbkg'
    altsnb_name = self.name+'_data_altsigbkg'
    nomsim_name = self.name+'_mc_nom'
    altsim_name = self.name+'_mc_alt'
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
    with open('out/'+nomsim_name+'/cnc_efficiencies.json','r') as input_file:
      eff_sim_nom = json.loads(input_file.read())
    with open('out/'+altsim_name+'/cnc_efficiencies.json','r') as input_file:
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
      eff_sim2 = eff_sim_alt[ibin][0]
      unc_dat1 = eff_dat_nom[ibin][1]
      unc_sim1 = eff_sim_nom[ibin][1]
      unc_sim2 = eff_sim_alt[ibin][1]

      eff_dat_mean = statistics.mean((eff_dat1, eff_dat2, eff_dat3, eff_dat4))
      eff_dat_unc = math.hypot(unc_dat1, math.hypot((eff_dat1-eff_dat_mean),
          (eff_dat2-eff_dat_mean),(eff_dat3-eff_dat_mean),
          (eff_dat4-eff_dat_mean))/math.sqrt(12.0))
      data_eff.append(eff_dat_mean)
      data_unc.append(eff_dat_unc)
      
      if (eff_sim1>0.0 and eff_sim1<1.0):
        mc_eff.append(eff_sim1)
        mc_unc.append(max(unc_sim1,abs(eff_sim1-eff_sim2)))
      else:
        mc_eff.append(eff_sim2)
        mc_unc.append(max(unc_sim2,abs(eff_sim1-eff_sim2)))
      sfs = calculate_sfs(eff_dat1, eff_dat2, eff_dat3, eff_dat4, 
                          eff_sim1, eff_sim2, unc_dat1, unc_sim1,
                          unc_sim2)
      pass_sf.append(sfs[0])
      pass_unc.append(sfs[1])
      fail_sf.append(sfs[2])
      fail_unc.append(sfs[3])

    if self.binning_type=='std':
      self.generate_jsons_nogap(data_eff, data_unc, mc_eff, mc_unc, 
                                pass_sf, pass_unc, fail_sf, fail_unc)
      self.generate_summary_plots_nogap(data_eff, data_unc, mc_eff, mc_unc,
                                        pass_sf, pass_unc, fail_sf, fail_unc)
    elif self.binning_type=='std_gap':
      self.generate_jsons_gap(data_eff, data_unc, mc_eff, mc_unc, pass_sf, 
                              pass_unc, fail_sf, fail_unc)
      self.generate_summary_plots_gap(data_eff, data_unc, mc_eff, mc_unc,
                                      pass_sf, pass_unc, fail_sf, fail_unc)
    else:
      raise RuntimeError('Unsupported binning')

  def generate_jsons_nogap(self, data_eff, data_unc, mc_eff, mc_unc, pass_sf,
                           pass_unc, fail_sf, fail_unc):
    '''Generate output assuming standard binning with no gap

    data_eff  list of data efficiencies
    data_unc  list of data uncertainties
    mc_eff    list of mc efficiencies
    mc_unc    list of mc uncertainties
    pass_sf   list of scale factors (SFs) for passing selection
    pass_unc  list of uncertainties on passing SFs
    fail_sf   list of scale factors for failing selection
    fail_unc  list of uncertainties on failing SFs
    '''

    if not os.path.isdir('out/'+self.name):
      print('Output directory not found, making new output directory')
      os.mkdir('out/'+self.name)

    #write JSON
    clib_sfs_pass = make_correction('sf_pass', 'data-MC SF', self.pt_bins, 
                                    self.eta_bins, pass_sf)
    clib_uns_pass = make_correction('unc_pass', 'data-MC unc', self.pt_bins, 
                                    self.eta_bins, pass_unc)
    clib_sfs_fail = make_correction('sf_fail', 'data-MC SF', self.pt_bins, 
                                    self.eta_bins, fail_sf)
    clib_uns_fail = make_correction('unc_fail', 'data-MC unc', self.pt_bins, 
                                    self.eta_bins, fail_unc)
    clib_dat_eff = make_correction('effdata', 'data eff', self.pt_bins, 
                                   self.eta_bins, data_eff)
    clib_dat_unc = make_correction('systdata', 'data unc', self.pt_bins, 
                                   self.eta_bins, data_unc)
    clib_sim_eff = make_correction('effmc', 'MC eff', self.pt_bins, 
                                   self.eta_bins, mc_eff)
    clib_sim_unc = make_correction('systmc', 'MC unc', self.pt_bins, 
                                   self.eta_bins, mc_unc)

    sf_filename = 'out/{0}/{0}_scalefactors.json'.format(self.name)
    eff_filename = 'out/{0}/{0}_efficiencies.json'.format(self.name)
    with open(eff_filename,'w') as output_file:
      output_file.write(fix_correctionlib_json(
        [clib_dat_eff.json(exclude_unset=True),
         clib_dat_unc.json(exclude_unset=True),
         clib_sim_eff.json(exclude_unset=True),
         clib_sim_unc.json(exclude_unset=True)]))
    with open(sf_filename,'w') as output_file:
      output_file.write(fix_correctionlib_json(
        [clib_sfs_pass.json(exclude_unset=True),
         clib_uns_pass.json(exclude_unset=True),
         clib_sfs_fail.json(exclude_unset=True),
         clib_uns_fail.json(exclude_unset=True)]))

  def generate_jsons_gap(self, data_eff, data_unc, mc_eff, mc_unc, pass_sf, 
                         pass_unc, fail_sf, fail_unc):
    '''Generate output assuming standard gap binning

    data_eff  list of data efficiencies
    data_unc  list of data uncertainties
    mc_eff    list of mc efficiencies
    mc_unc    list of mc uncertainties
    pass_sf   list of scale factors (SFs) for passing selection
    pass_unc  list of uncertainties on passing SFs
    fail_sf   list of scale factors for failing selection
    fail_unc  list of uncertainties on failing SFs
    '''
    #organize SFs as they will be saved in the JSON
    gapincl_eta_bins, neg_gap_idx, pos_gap_idx = add_gap_eta_bins(self.eta_bins)
    pass_json_sfs = []
    pass_json_uns = []
    fail_json_sfs = []
    fail_json_uns = []
    json_dat_eff = []
    json_dat_unc = []
    json_sim_eff = []
    json_sim_unc = []
    num_bins_pt = len(self.pt_bins)-1
    num_bins_eta = len(self.eta_bins)-1
    for ipt in range(num_bins_pt):
      mean_bin_pt = (self.pt_bins[ipt]+self.pt_bins[ipt+1])/2.0
      for ieta in range(num_bins_eta+2):
        tnp_bin = -1
        if ieta < neg_gap_idx:
          tnp_bin = ipt*num_bins_eta+ieta
        elif ieta == neg_gap_idx:
          tnp_bin = (num_bins_pt*num_bins_eta
                     +get_bin(mean_bin_pt, self.gap_pt_bins)*2)
        elif ieta > neg_gap_idx and ieta < pos_gap_idx:
          tnp_bin = ipt*num_bins_eta+(ieta-1)
        elif ieta == pos_gap_idx:
          tnp_bin = (num_bins_pt*num_bins_eta
                     +get_bin(mean_bin_pt, self.gap_pt_bins)*2+1)
        elif ieta > pos_gap_idx:
          tnp_bin = ipt*num_bins_eta+(ieta-2)
        pass_json_sfs.append(pass_sf[tnp_bin])
        pass_json_uns.append(pass_unc[tnp_bin])
        fail_json_sfs.append(fail_sf[tnp_bin])
        fail_json_uns.append(fail_unc[tnp_bin])
        json_dat_eff.append(data_eff[tnp_bin])
        json_dat_unc.append(data_unc[tnp_bin])
        json_sim_eff.append(mc_eff[tnp_bin])
        json_sim_unc.append(mc_unc[tnp_bin])

    if not os.path.isdir('out/'+self.name):
      print('Output directory not found, making new output directory')
      os.mkdir('out/'+self.name)

    #write JSON
    clib_sfs_pass = make_correction('sf_pass', 'data-MC SF', self.pt_bins, 
                                    gapincl_eta_bins, pass_json_sfs)
    clib_uns_pass = make_correction('unc_pass', 'data-MC unc', self.pt_bins, 
                                    gapincl_eta_bins, pass_json_uns)
    clib_sfs_fail = make_correction('sf_fail', 'data-MC SF', self.pt_bins, 
                                    gapincl_eta_bins, fail_json_sfs)
    clib_uns_fail = make_correction('unc_fail', 'data-MC unc', self.pt_bins, 
                                    gapincl_eta_bins, fail_json_uns)
    clib_dat_eff = make_correction('effdata', 'data eff', self.pt_bins, 
                                   gapincl_eta_bins, json_dat_eff)
    clib_dat_unc = make_correction('systdata', 'data unc', self.pt_bins, 
                                   gapincl_eta_bins, json_dat_unc)
    clib_sim_eff = make_correction('effmc', 'MC eff', self.pt_bins, 
                                   gapincl_eta_bins, json_sim_eff)
    clib_sim_unc = make_correction('systmc', 'MC unc', self.pt_bins, 
                                   gapincl_eta_bins, json_sim_unc)

    sf_filename = 'out/{0}/{0}_scalefactors.json'.format(self.name)
    eff_filename = 'out/{0}/{0}_efficiencies.json'.format(self.name)
    with open(eff_filename,'w') as output_file:
      output_file.write(fix_correctionlib_json(
        [clib_dat_eff.json(exclude_unset=True),
         clib_dat_unc.json(exclude_unset=True),
         clib_sim_eff.json(exclude_unset=True),
         clib_sim_unc.json(exclude_unset=True)]))
    with open(sf_filename,'w') as output_file:
      output_file.write(fix_correctionlib_json(
        [clib_sfs_pass.json(exclude_unset=True),
         clib_uns_pass.json(exclude_unset=True),
         clib_sfs_fail.json(exclude_unset=True),
         clib_uns_fail.json(exclude_unset=True)]))

  def generate_summary_plots_nogap(self, data_eff, data_unc, mc_eff, mc_unc, 
                                 pass_sf, pass_unc, fail_sf, fail_unc):
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
    eff_eta_plot_data_y = []
    eff_eta_plot_data_ey = []
    eff_eta_plot_mc_y = []
    eff_eta_plot_mc_ey = []
    eff_pt_plot_data_y = []
    eff_pt_plot_data_ey = []
    eff_pt_plot_mc_y = []
    eff_pt_plot_mc_ey = []
    sf_eta_plot_pass_y = []
    sf_eta_plot_pass_ey = []
    sf_eta_plot_fail_y = []
    sf_eta_plot_fail_ey = []
    sf_pt_plot_pass_y = []
    sf_pt_plot_pass_ey = []
    sf_pt_plot_fail_y = []
    sf_pt_plot_fail_ey = []
    eta_plot_names = []
    eta_plot_data_names = []
    eta_plot_mc_names = []
    pt_plot_names = []
    pt_plot_data_names = []
    pt_plot_mc_names = []

    for ieta in range(len(self.eta_bins)-1):
      eta_plot_x.append((self.eta_bins[ieta+1]+self.eta_bins[ieta])/2.0)
      eta_plot_ex.append((self.eta_bins[ieta+1]-self.eta_bins[ieta])/2.0)
      pt_plot_names.append('{}<|#eta|<{}'.format(self.eta_bins[ieta],self.eta_bins[ieta+1]))
      pt_plot_data_names.append('Data {}<|#eta|<{}'.format(self.eta_bins[ieta],self.eta_bins[ieta+1]))
      pt_plot_mc_names.append('MC {}<|#eta|<{}'.format(self.eta_bins[ieta],self.eta_bins[ieta+1]))
      eff_pt_plot_data_y.append([])
      eff_pt_plot_data_ey.append([])
      eff_pt_plot_mc_y.append([])
      eff_pt_plot_mc_ey.append([])
      sf_pt_plot_pass_y.append([])
      sf_pt_plot_pass_ey.append([])
      sf_pt_plot_fail_y.append([])
      sf_pt_plot_fail_ey.append([])
    for ipt in range(len(self.pt_bins)-1):
      pt_plot_x.append((self.pt_bins[ipt+1]+self.pt_bins[ipt])/2.0)
      pt_plot_ex.append((self.pt_bins[ipt+1]-self.pt_bins[ipt])/2.0)
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

    for ipt in range(len(self.pt_bins)-1):
      for ieta in range(len(self.eta_bins)-1):
        tnp_bin = ieta+ipt*(len(self.eta_bins)-1)
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

    eff_string = 'Efficiency '+self.data_nom_tnp_analyzer.measurement_desc
    eff_string = 'Efficiency '+self.data_nom_tnp_analyzer.measurement_desc
    unc_string = 'Eff. Unc. '+self.data_nom_tnp_analyzer.measurement_desc
    passsf_string = 'Pass SF '+self.data_nom_tnp_analyzer.measurement_desc
    failsf_string = 'Fail SF '+self.data_nom_tnp_analyzer.measurement_desc
    passunc_string = 'Pass SF Unc. '+self.data_nom_tnp_analyzer.measurement_desc
    failunc_string = 'Fail SF Unc. '+self.data_nom_tnp_analyzer.measurement_desc
    for fit_syst in ['data_nom','data_altsig','data_altbkg','data_altsigbkg','mc_nom']:
      os.system('cp out/{0}_{1}/allfits.pdf out/{0}/{1}_allfits.pdf'.format(self.name,fit_syst))
    make_data_mc_graph(eta_plot_x, eta_plot_ex, eff_eta_plot_data_y, 
                       eff_eta_plot_data_ey, eff_eta_plot_mc_y, 
                       eff_eta_plot_mc_ey, 
                       'out/{0}/{0}_eff_etabinned.pdf'.format(self.name), 
                       eta_plot_data_names, eta_plot_mc_names,
                       '|#eta|', eff_string, LUMI_TAGS[self.year])
    make_data_mc_graph(pt_plot_x, pt_plot_ex, eff_pt_plot_data_y, 
                       eff_pt_plot_data_ey, eff_pt_plot_mc_y, 
                       eff_pt_plot_mc_ey, 
                       'out/{0}/{0}_eff_ptbinned.pdf'.format(self.name),
                       pt_plot_data_names, pt_plot_mc_names,
                       'p_{T} [GeV]', eff_string, LUMI_TAGS[self.year],True)
    make_heatmap(self.eta_bins, self.pt_bins, eff_pt_plot_data_y, 
                 'out/{0}/{0}_eff_data.pdf'.format(self.name), '|#eta|', 
                 'p_{T} [GeV]', 
                 eff_string, LUMI_TAGS[self.year],False,True)
    make_heatmap(self.eta_bins, self.pt_bins, eff_pt_plot_mc_y, 
                 'out/{0}/{0}_eff_mc.pdf'.format(self.name), '|#eta|', 
                 'p_{T} [GeV]',
                 eff_string, LUMI_TAGS[self.year],False,True)
    make_sf_graph(eta_plot_x, eta_plot_ex, sf_eta_plot_pass_y, 
                  sf_eta_plot_pass_ey, 
                  'out/{0}/{0}_sfpass_etabinned.pdf'.format(self.name),
                  eta_plot_names, '|#eta|', 'Pass SF', LUMI_TAGS[self.year])
    make_sf_graph(eta_plot_x, eta_plot_ex, sf_eta_plot_fail_y, 
                  sf_eta_plot_fail_ey, 
                  'out/{0}/{0}_sffail_etabinned.pdf'.format(self.name),
                  eta_plot_names, '|#eta|', 'Fail SF', LUMI_TAGS[self.year])
    make_sf_graph(pt_plot_x, pt_plot_ex, sf_pt_plot_pass_y, sf_pt_plot_pass_ey,
                  'out/{0}/{0}_sfpass_ptbinned.pdf'.format(self.name),
                  pt_plot_names, 'p_{T} [GeV]', 'Pass SF', LUMI_TAGS[self.year],
                  True)
    make_sf_graph(pt_plot_x, pt_plot_ex, sf_pt_plot_fail_y, sf_pt_plot_fail_ey, 
                  'out/{0}/{0}_sffail_ptbinned.pdf'.format(self.name),
                  pt_plot_names, 'p_{T} [GeV]', 'Fail SF', LUMI_TAGS[self.year],
                  True)
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_pass_y, 
                 'out/{0}/{0}_sfpass.pdf'.format(self.name), '|#eta|', 
                 'p_{T} [GeV]', passsf_string, LUMI_TAGS[self.year], False, 
                 True)
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_fail_y, 
                 'out/{0}/{0}_sffail.pdf'.format(self.name), '|#eta|', 
                 'p_{T} [GeV]', failsf_string, LUMI_TAGS[self.year], False, 
                 True)
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_pass_ey, 
                 'out/{0}/{0}_sfpass_unc.pdf'.format(self.name), '|#eta|', 
                 'p_{T} [GeV]', passunc_string, LUMI_TAGS[self.year], False, 
                 True)
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_fail_ey, 
                 'out/{0}/{0}_sffail_unc.pdf'.format(self.name), '|#eta|', 
                 'p_{T} [GeV]', failunc_string, LUMI_TAGS[self.year], False, 
                 True)

  def generate_summary_plots_gap(self, data_eff, data_unc, mc_eff, mc_unc, 
                                 pass_sf, pass_unc, fail_sf, fail_unc):
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
    gappt_plot_x = []
    gappt_plot_ex = []
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
      eta_plot_x.append((self.eta_bins[ieta+1]+self.eta_bins[ieta])/2.0)
      eta_plot_ex.append((self.eta_bins[ieta+1]-self.eta_bins[ieta])/2.0)
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
      pt_plot_x.append((self.pt_bins[ipt+1]+self.pt_bins[ipt])/2.0)
      pt_plot_ex.append((self.pt_bins[ipt+1]-self.pt_bins[ipt])/2.0)
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
      gappt_plot_x.append((self.gap_pt_bins[ipt+1]+self.gap_pt_bins[ipt])/2.0)
      gappt_plot_ex.append((self.gap_pt_bins[ipt+1]-self.gap_pt_bins[ipt])/2.0)
    gappt_plot_names.append('-1.566<#eta<-1.4442')
    gappt_plot_names.append('1.4442<#eta<1.566')
    gappt_plot_data_names.append('Data -1.566<#eta<-1.4442')
    gappt_plot_data_names.append('Data 1.4442<#eta<1.566')
    gappt_plot_mc_names.append('MC -1.566<#eta<-1.4442')
    gappt_plot_mc_names.append('MC 1.4442<#eta<1.566')
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
        tnp_bin = ieta+ipt*(len(self.eta_bins)-1)
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
        tnp_bin = (len(self.pt_bins)-1)*(len(self.eta_bins)-1)+ieta+ipt*2
        eff_gappt_plot_data_y[ieta].append(data_eff[tnp_bin])
        eff_gappt_plot_data_ey[ieta].append(data_unc[tnp_bin])
        eff_gappt_plot_mc_y[ieta].append(mc_eff[tnp_bin])
        eff_gappt_plot_mc_ey[ieta].append(mc_unc[tnp_bin])
        sf_gappt_plot_pass_y[ieta].append(pass_sf[tnp_bin])
        sf_gappt_plot_pass_ey[ieta].append(pass_unc[tnp_bin])
        sf_gappt_plot_fail_y[ieta].append(fail_sf[tnp_bin])
        sf_gappt_plot_fail_ey[ieta].append(fail_unc[tnp_bin])

    eff_string = 'Efficiency '+self.data_nom_tnp_analyzer.measurement_desc
    unc_string = 'Eff. Unc. '+self.data_nom_tnp_analyzer.measurement_desc
    passsf_string = 'Pass SF '+self.data_nom_tnp_analyzer.measurement_desc
    failsf_string = 'Fail SF '+self.data_nom_tnp_analyzer.measurement_desc
    passunc_string = 'Pass SF Unc. '+self.data_nom_tnp_analyzer.measurement_desc
    failunc_string = 'Fail SF Unc. '+self.data_nom_tnp_analyzer.measurement_desc
    for fit_syst in ['data_nom','data_altsig','data_altbkg','data_altsigbkg','mc_nom']:
      os.system('cp out/{0}_{1}/allfits.pdf out/{0}/{1}_allfits.pdf'.format(self.name,fit_syst))
    make_data_mc_graph(eta_plot_x, eta_plot_ex, eff_eta_plot_data_y, 
                       eff_eta_plot_data_ey, eff_eta_plot_mc_y, 
                       eff_eta_plot_mc_ey, 
                       'out/{0}/{0}_eff_etabinned.pdf'.format(self.name), 
                       eta_plot_data_names, eta_plot_mc_names,
                       '#eta', eff_string, LUMI_TAGS[self.year])
    make_data_mc_graph(pt_plot_x, pt_plot_ex, eff_pt_plot_data_y, 
                       eff_pt_plot_data_ey, eff_pt_plot_mc_y, 
                       eff_pt_plot_mc_ey, 
                       'out/{0}/{0}_eff_ptbinned.pdf'.format(self.name),
                       pt_plot_data_names, pt_plot_mc_names,
                       'p_{T} [GeV]', eff_string, LUMI_TAGS[self.year],True)
    make_data_mc_graph(gappt_plot_x, gappt_plot_ex, eff_gappt_plot_data_y, 
                       eff_gappt_plot_data_ey, eff_gappt_plot_mc_y, 
                       eff_gappt_plot_mc_ey, 
                       'out/{0}/{0}_eff_gapptbinned.pdf'.format(self.name),
                       gappt_plot_data_names, gappt_plot_mc_names,
                       'p_{T} [GeV]', eff_string, LUMI_TAGS[self.year],True)
    make_heatmap(self.eta_bins, self.pt_bins, eff_pt_plot_data_y, 
                 'out/{0}/{0}_eff_data.pdf'.format(self.name), '#eta', 'p_{T} [GeV]', 
                 eff_string, LUMI_TAGS[self.year],False,True)
    make_heatmap(self.eta_bins, self.pt_bins, eff_pt_plot_mc_y, 
                 'out/{0}/{0}_eff_mc.pdf'.format(self.name), '#eta', 'p_{T} [GeV]',
                 eff_string, LUMI_TAGS[self.year],False,True)
    make_sf_graph(eta_plot_x, eta_plot_ex, sf_eta_plot_pass_y, 
                  sf_eta_plot_pass_ey, 
                  'out/{0}/{0}_sfpass_etabinned.pdf'.format(self.name),
                  eta_plot_names, '#eta', 'Pass SF', LUMI_TAGS[self.year])
    make_sf_graph(eta_plot_x, eta_plot_ex, sf_eta_plot_fail_y, 
                  sf_eta_plot_fail_ey, 
                  'out/{0}/{0}_sffail_etabinned.pdf'.format(self.name),
                  eta_plot_names, '#eta', 'Fail SF', LUMI_TAGS[self.year])
    make_sf_graph(pt_plot_x, pt_plot_ex, sf_pt_plot_pass_y, sf_pt_plot_pass_ey,
                  'out/{0}/{0}_sfpass_ptbinned.pdf'.format(self.name),
                  pt_plot_names, 'p_{T} [GeV]', 'Pass SF', LUMI_TAGS[self.year],
                  True)
    make_sf_graph(pt_plot_x, pt_plot_ex, sf_pt_plot_fail_y, sf_pt_plot_fail_ey, 
                  'out/{0}/{0}_sffail_ptbinned.pdf'.format(self.name),
                  pt_plot_names, 'p_{T} [GeV]', 'Fail SF', LUMI_TAGS[self.year],
                  True)
    make_sf_graph(gappt_plot_x, gappt_plot_ex, sf_gappt_plot_pass_y, 
                  sf_gappt_plot_pass_ey, 
                  'out/{0}/{0}_sfpass_gapptbinned.pdf'.format(self.name),
                  gappt_plot_names, 'p_{T} [GeV]', 'Pass SF', 
                  LUMI_TAGS[self.year], True)
    make_sf_graph(gappt_plot_x, gappt_plot_ex, sf_gappt_plot_fail_y, 
                  sf_gappt_plot_fail_ey, 
                  'out/{0}/{0}_sffail_gapptbinned.pdf'.format(self.name),
                  gappt_plot_names, 'p_{T} [GeV]', 'Fail SF', 
                  LUMI_TAGS[self.year], True)
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_pass_y, 
                 'out/{0}/{0}_sfpass.pdf'.format(self.name), '#eta', 
                 'p_{T} [GeV]', passsf_string, LUMI_TAGS[self.year], False, 
                 True)
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_fail_y, 
                 'out/{0}/{0}_sffail.pdf'.format(self.name), '#eta', 
                 'p_{T} [GeV]', failsf_string, LUMI_TAGS[self.year], False, 
                 True)
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_pass_ey, 
                 'out/{0}/{0}_sfpass_unc.pdf'.format(self.name), '#eta', 
                 'p_{T} [GeV]', passunc_string, LUMI_TAGS[self.year], False, 
                 True)
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_fail_ey, 
                 'out/{0}/{0}_sffail_unc.pdf'.format(self.name), '#eta', 
                 'p_{T} [GeV]', failunc_string, LUMI_TAGS[self.year], False, 
                 True)

  def run_interactive(self):
    '''
    Run an interactive T&P analysis
    '''
    self.add_models()

    #run main interactive loop
    exit_loop = False
    print('Welcome to rms_analyzer interactive analysis.')
    print('Type [h]elp for more information.')
    past_commands = []
    while not exit_loop:
      user_input = input('>:')
      past_commands.append(user_input)
      user_input = user_input.split()
      if len(user_input)<1:
        continue
      elif (user_input[0] == 'h' or user_input[0] == 'help'):
        print('Available commands:')
        print('h(elp)                         print help information')
        print('p(roduce)                      produce histograms to fit; run this first')
        print('f(it) <sample> [<bin> <pass>]  run interactive session (sample=nom,alts,')
        print('                               altb,altsb,mc) use bin and pass(=p/f) to ')
        print('                               start from a particular bin and category')
        print('o(utput)                       generate final outputs')
        print('c(lean)                        cleans previous output')
        print('q(uit)                         exit')
        print('(pre)v(ious)                   list previous commands entered')
      elif (user_input[0] == 'p' or user_input[0] == 'produce'):
        self.produce_histograms()
      elif (user_input[0] == 'f' or user_input[0] == 'fit'):
        if len(user_input)<2:
          print('ERROR: f(it) requires an argument')
        else:
          starting_bin = '0'
          starting_cat = 'p'
          if (len(user_input)>=3):
            starting_bin = user_input[2]
          if (len(user_input)>=4):
            starting_cat = user_input[3]
          if (user_input[1] == 'nom'):
            #avoid accessing the same ROOT file in multiple tnp_analyzers
            self.mc_nom_tnp_analyzer.close_file()
            self.data_nom_tnp_analyzer.close_file()
            self.data_nom_tnp_analyzer.fit_histogram(starting_bin,starting_cat,
                                                     self.nom_fn_name)
          elif (user_input[1] == 'nomcont'):
            #avoid accessing the same ROOT file in multiple tnp_analyzers
            self.mc_nom_tnp_analyzer.close_file()
            self.data_nom_tnp_analyzer.close_file()
            self.data_nom_tnp_analyzer.fit_histogram(starting_bin,starting_cat,
                                                     self.contingency_fn_name)
          elif (user_input[1] == 'alts' or user_input[1] == 'altsignal'):
            self.data_altsig_tnp_analyzer.close_file()
            self.data_altsig_tnp_analyzer.fit_histogram(starting_bin,
                                                        starting_cat,
                                                        self.alts_fn_name,
                                                        self.alts_fn_init)
          elif (user_input[1] == 'altscont'):
            self.data_altsig_tnp_analyzer.close_file()
            self.data_altsig_tnp_analyzer.fit_histogram(starting_bin,
                starting_cat,self.contingencyalts_fn_name,self.alts_fn_init)
          elif (user_input[1] == 'altb' or user_input[1] == 'altbackground'):
            #avoid accessing the same ROOT file in multiple tnp_analyzers
            self.mc_nom_tnp_analyzer.close_file()
            self.data_altbkg_tnp_analyzer.close_file()
            self.data_altbkg_tnp_analyzer.fit_histogram(starting_bin,
                                                        starting_cat,
                                                        self.altb_fn_name)
          elif (user_input[1] == 'altsb' or 
                user_input[1] == 'altsignalbackground'):
            self.data_altsigbkg_tnp_analyzer.close_file()
            self.data_altsigbkg_tnp_analyzer.fit_histogram(starting_bin,
                                                           starting_cat,
                                                           self.altsb_fn_name,
                                                           self.alts_fn_init)
          elif (user_input[1] == 'mc' or user_input[1] == 'mcalt'):
            self.mc_nom_tnp_analyzer.close_file()
            self.mc_nom_tnp_analyzer.fit_histogram(starting_bin,
                                                   starting_cat,
                                                   self.alts_fn_name_sa)
          else:
            print('ERROR: unrecognized argument to f(it)')
      elif (user_input[0] == 'o' or user_input[0] == 'output'):
        #TODO add calls to outputs for each individual piece
        self.generate_output()
      elif (user_input[0] == 'q' or user_input[0] == 'quit'):
        exit_loop = True
      elif (user_input[0] == 'c' or user_input[0] == 'clean'):
        self.clean_output()
      elif (user_input[0] == 'v' or user_input[0] == 'previous'):
        for past_command in past_commands:
          print(past_command)
      else:
        print('ERROR: unrecognized command')

