#!/usr/bin/env python3
"""@package docstring
T&P meta analyzer that follows standard EGM procedures to generate scale factors. In particular:
"""

from array import array
from correctionlib import schemav2 
from functools import partial
import gc
import os
import ROOT
import statistics
import subprocess
import json

from tnp_analyzer import *
from tnp_utils import CMS_COLORS, LUMI_TAGS
from model_initializers import *
from root_plot_lib import RplPlot

def param_initializer_dscb_from_mc(ibin: int, is_pass: bool, 
                                   workspace: ROOT.RooWorkspace, 
                                   mc_analyzer: TnpAnalyzer):
  '''Parameter initializer for DSCB based on MC fit

  Parameter initializer for cheby_dscb model that fixes DSCB parameters 
  except mean and sigma to MC result
  
  Args:
    ibin: bin number
    is_pass: indicates if passing leg
    workspace: RooWorkspace for this bin
    mc_analyzer: TnpAnalyzer for MC samples
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

def param_initializer_moddscb_from_mc(ibin: int, is_pass: bool, 
                                      workspace: ROOT.RooWorkspace, 
                                      mc_analyzer: TnpAnalyzer):
  '''Parameter initializer for Cheby+DSCB based on MC fit

  Parameter initializer for cheby_dscb model that fixes DSCB parameters except
  mean and sigma to MC result
  
  Args:
    ibin: bin number
    is_pass: indicates if passing leg
    workspace: RooWorkspace for this bin
    mc_analyzer: TnpAnalyzer for MC samples
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

def param_initializer_dscbgaus_from_mc(ibin: int, is_pass: bool, 
                                       workspace: ROOT.RooWorkspace, 
                                       mc_analyzer: TnpAnalyzer):
  '''Parameter initializer for Cheby and DSCB+Gauss based on MC fit

  Parameter initializer for cheby_dscb model that fixes DSCB parameters except
  mean and sigma to MC result
  
  Args:
    ibin: bin number
    is_pass: indicates if passing leg
    workspace: RooWorkspace for this bin
    mc_analyzer: TnpAnalyzer for MC samples
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

def param_initializer_cbconvgen_from_mc(ibin: int, is_pass: bool, 
                                        workspace: ROOT.RooWorkspace, 
                                        mc_analyzer: TnpAnalyzer):
  '''Parameter initializer for CB * template parameters to MC fit

  Parameter initializer for cbconvgen model that fixes CB parameters to MC
  result
  
  Args:
    ibin: bin number
    is_pass: indicates if passing leg
    workspace: RooWorkspace for this bin
    mc_analyzer: TnpAnalyzer for MC samples
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

def get_mc_histogram(ibin: int, is_pass: bool, mc_analyzer: TnpAnalyzer, 
                     highpt_bins: list[int]) -> ROOT.TH1D:
  '''Helper function used to get appropriate histogram from analyzer

  Args:
    ibin: bin index
    is_pass: indicates if passing leg
    mc_analyzer: TnpAnalyzer to extract histogram from
    highpt_bins: list of high pt bins where the failing leg is ignored

  Returns:
    specified histogram
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

def add_gap_eta_bins(original_bins: list[float]) -> tuple[list[float],int,int]:
  '''Modifies binning to include gap bins

  Modifies eta binning to include EB-EE gap bins (-1.566,-1.4442) and 
  (1.4442,1.566). Original_bins MUST include -1.5 and 1.5

  Args:
    original_bins: ordered bin edges

  Returns: 
    (new_bins, -gap_index, +gap_index)
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

def calculate_sfs(eff_dat1: float, eff_dat2: float, eff_dat3: float, 
                  eff_dat4: float, eff_sim1: float, eff_sim2: float, 
                  unc_dat1: float, unc_sim1: float, 
                  unc_sim2: float) -> tuple[float,float,float,float]:
  '''Computes SFs factors from efficiencies and uncertainties using RMS method

  Args:
    eff_dat1: data efficiency measurement 1
    eff_dat2: data efficiency measurement 2
    eff_dat3: data efficiency measurement 3
    eff_dat4: data efficiency measurement 4
    eff_sim1: simulation efficiency measurement 1
    eff_sim2: simulation efficiency measurement 2
    unc_dat1: data efficiency 1 uncertainty
    unc_sim1: simulation efficiency 1 uncertainty
    unc_sim2: simulation efficiency 2 uncertainty

  Returns:
    (scale factor for passing events,
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
    nonzero_eff_sim = unc_sim1
    nonzero_unc_sim = unc_sim1
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

def make_data_mc_graph(x: list[float], ex: list[float], 
                       data_y: list[list[float]], data_ey: list[list[float]], 
                       sim_y: list[list[float]], sim_ey: list[list[float]], 
                       name: str, data_names: list[float], 
                       mc_names: list[str], x_title: str, 
                       y_title: str, lumi: tuple[float, float], 
                       log_x: bool=False):
  '''Generate Data-MC comparison plot

  Makes a nice plot with multiple graphs overlayed. Each data graph will be
  drawn in solid line while each simulation graph in dotted lines. The number
  of y/ey points should match

  Args:
    x: x values for points
    ex: x error bars for points
    data_y: y values for data points
    data_ey: y error bars for data points
    sim_y: y values for simulation points
    sim_ey: y error bars for simulation points
    name: filename
    data_names: names for data graphs
    mc_names: names for mc graphs
    x_title: X-axis label
    y_title: y-axis label
    lumi: lumi and CM energy
    log_x: if true, makes x-axis logarithmic
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

def make_correction(name: str, desc: str, pt_bins: list[float], 
                    eta_bins: list[float], 
                    content: list[float]) -> schemav2.Correction:
  '''Generates correctionlib correction object

  Args:
    name: correction name
    desc: description
    pt_bins: bin edges
    eta_bins: bin edges
    content: content
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

def make_sf_graph(x: list[float], ex: list[float], y: list[list[float]], 
                  ey: list[list[float]], name: str, graph_names: list[str], 
                  x_title: str, y_title: str, 
                  lumi: list[tuple[float, float]], log_x: bool=False):
  '''Makes a nice plot with multiple graphs overlayed. 

  Args:
    x: x values for points
    ex: x error bars for points
    y: y values  for points
    ey: y error bars for points
    name: filename
    graph_names: names for graphs
    x_title: X-axis label
    y_title: y-axis label
    lumi: list of tuple of two floats representing lumi and CM energy
    log_x: if true sets x-axis to be logarithmic
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

def make_heatmap(x: list[float], y: list[float], z: list[float], name: str, 
                 x_title: str, y_title: str, z_title: str, 
                 lumi: list[tuple[float,float]], log_x: bool=False, 
                 log_y: bool=False):
  '''Makes a heatmap (2D histogram/colz)

  Args:
    x: x axis bin divisions
    y: y axis bin divisions
    z: heatmap values
    name: filename
    x_title: x-axis label
    y_title: y-axis label
    z_title: z-axis label
    lumi: list of tuple of two floats representing lumi and CM energy
    log_x: if true sets x-axis to be logarithmic
    log_y: if true sets y-axis to be logarithmic
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
  '''Class that manages RMS SF/efficiency derivation

  Attributes:
    name: analyzer name
    data_nom_tnp_analyzer: analyzer for nominal data efficiencies
    data_altsig_tnp_analyzer: analyzer for alt sig data efficiencies
    data_altbkg_tnp_analyzer: analyzer for alt bkg data efficiencies
    data_altsigbkg_tnp_analyzer: analyzer for alt sig+bkg data efficiencies
    mc_nom_tnp_analyzer: analyzer for nominal MC efficiencies
    mc_alt_tnp_analyzer: analyzer for alternative MC efficiencies
    binning_type: 'std', 'std_gap', or 'custom' to indicate binning
    year: data taking era
    pt_bins: pt bin edges
    eta_bins: eta bin edges
    gap_pt_bins: pt bin edges for the gap region
    highpt_bins: list of bins considered high-pt
    nom_fn_name: name of nominal fit function
    contingency_fn_name: name of nominal fit function used in contingency
    alts_fn_name: name of alternate signal fit function
    contingencyalts_fn_name: name of alternate signal fit function for cont.
    altb_fn_name: name of alternate background fit function
    altsb_fn_name: name of alternate signal+background fit function
    alts_fn_init: initializer for alternate signal fit function
    alts_fn_name_sa: name of alternate signal fit function without background
  '''

  def __init__(self, name: str):
    '''Constructor

    Args:
      name: name of analyzer
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

    self.pt_bins = []
    self.eta_bins = []
    self.gap_pt_bins = []
    self.highpt_bins = []
    self.nom_fn_name = ''
    self.contingency_fn_name = ''
    self.alts_fn_name = ''
    self.contingencyalts_fn_name = ''
    self.altb_fn_name = ''
    self.altsb_fn_name = ''
    self.alts_fn_init = ''
    self.alts_fn_name_sa = ''

  def set_input_files(self, data_files: list[str], mc_files: list[str], 
                      mc_alt_files: list[str], data_tree: str, 
                      mc_tree: str='', mc_alt_tree: str=''):
    '''Sets input files

    Args:
      data_files: data files to process
      mc_files: mc files to process
      mc_alt_files: alternative MC files to process
      data_tree: name of TTree in data files
      mc_tree: name of TTree in MC files, defaults to data_tree if ''
      mc_alt_tree: name of TTree in MC alt file, defaults to mc_tree if ''
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
    
  def set_fitting_variable(self, name: str, description: str, nbins: int=60, 
                           nbins_mc: int=80,
                           var_range: tuple[float,float]=(60.0,120.0), 
                           custom_mc_range: tuple[float,float]=(50.0,130.0), 
                           weight: str='1', 
                           weight_mc: str='1'):
    '''Adds information about fitting variable

    Args:
      name: name of branch in TTree or C++ expression
      description: name used in plots (with TLaTeX)
      nbins: number of bins to use for fit variable
      var_range: start and end of fit range
      custom_mc_range: start and end of MC histogram range
      weight: expression for weight to use
      weight_mc: expression for weight to use for MC
    '''
    self.data_nom_tnp_analyzer.set_fitting_variable(name, description, nbins, 
                                                    var_range, weight)
    self.data_nom_tnp_analyzer.set_template_range(custom_mc_range)
    self.data_altsig_tnp_analyzer.set_fitting_variable(name, description, 
                                                       nbins, var_range, weight)
    self.data_altbkg_tnp_analyzer.set_fitting_variable(name, description, 
                                                       nbins, var_range, weight)
    self.data_altbkg_tnp_analyzer.set_template_range(custom_mc_range)
    self.data_altsigbkg_tnp_analyzer.set_fitting_variable(name, description, 
                                                          nbins, var_range, 
                                                          weight)
    self.mc_nom_tnp_analyzer.set_fitting_variable(name, description, nbins_mc, 
                                                  custom_mc_range, weight_mc)
    self.mc_alt_tnp_analyzer.set_fitting_variable(name, description, nbins_mc, 
                                                  custom_mc_range, weight_mc)
    self.mc_nom_tnp_analyzer.set_custom_fit_range(var_range)
    self.mc_alt_tnp_analyzer.set_custom_fit_range(var_range)

  def set_measurement_variable(self, var: str, desc: str=''):
    '''Sets selection efficiency to measure with tag & probe

    Args:
      var: name of branch in TTree or C++ expression
      desc: description o fmeasurement variable
    '''
    self.data_nom_tnp_analyzer.set_measurement_variable(var,desc)
    self.data_altsig_tnp_analyzer.set_measurement_variable(var,desc)
    self.data_altbkg_tnp_analyzer.set_measurement_variable(var,desc)
    self.data_altsigbkg_tnp_analyzer.set_measurement_variable(var,desc)
    self.mc_nom_tnp_analyzer.set_measurement_variable(var,desc)
    self.mc_alt_tnp_analyzer.set_measurement_variable(var,desc)

  def set_preselection(self, preselection_data: str, preselection_mc: str, 
                       desc: str):
    '''Sets basic preselection applied to all bins

    Args:
      preselection_data: selection for data as a C++ expression
      preselection_mc: selection for MC as a C++ expression
      desc: description of selection in TLaTeX
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

  def add_custom_binning(self, bin_selections: list[str], 
                         bin_names: list[str], is_high_pt: list[bool]):
    '''Creates custom bins for TnP analysis

    Args:
      bin_selections: describes selection for each bin
      bin_names: names of each bin that appear in plots
      is_high_pt: indicates whether each bin is high pT
    '''
    self.data_nom_tnp_analyzer.add_custom_binning(bin_selections, bin_names)
    self.data_altsig_tnp_analyzer.add_custom_binning(bin_selections, bin_names)
    self.data_altbkg_tnp_analyzer.add_custom_binning(bin_selections, bin_names)
    self.data_altsigbkg_tnp_analyzer.add_custom_binning(bin_selections, 
                                                        bin_names)
    self.mc_nom_tnp_analyzer.add_custom_binning(bin_selections, bin_names)
    self.mc_alt_tnp_analyzer.add_custom_binning(bin_selections, bin_names)
    self.highpt_bins = []
    for ibin in range(len(is_high_pt)):
      if is_high_pt[ibin]:
        self.highpt_bins.append(ibin)

  def add_standard_binning(self, pt_bins: list[float], eta_bins: list[float], 
                           pt_var_name: str, eta_var_name: str):
    '''Creates standard pt-eta binning 

    Args:
      pt_bins: pt bin edges
      eta_bins: eta bin edges
      pt_var_name: name of pt variable
      eta_var_name: name of eta variable
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

  def add_standard_gap_binning(self, pt_bins: list[float], 
                               eta_bins: list[float], 
                               gap_pt_bins: list[float], pt_var_name: str, 
                               eta_var_name: str):
    '''Creates standard binning including specialized pt bins for gap region

    Args:
      pt_bins: pt bin edges
      eta_bins: eta bin edges
      gap_pt_bins: gap region pt bin edges
      pt_var_name: name of pt variable
      eta_var_name: name of eta variable
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

  def add_models(self, gamma_add_gauss: bool=False):
    '''Adds standard models and parameter initializers for fitting

    Args:
      gamma_add_gauss: use extra gaussian for background model
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
    if gamma_add_gauss:
      alt_background_model = add_background_model_gammagauss
      alt_background_model_name = 'gammagauss'

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
    '''Produce histograms. Performs one loop over data files for efficiency
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

  def generate_individual_outputs(self, cnc: bool=False):
    '''Generates individual efficiency measurements if they have not already 
    been generated and checks output files are generated correctly

    Args:
      cnc: use only cut and count rather than fits

    Returns:
      true if the outputs are generated correctly
    '''
    nomdat_name = self.name+'_data_nom'
    altsig_name = self.name+'_data_altsig'
    altbkg_name = self.name+'_data_altbkg'
    altsnb_name = self.name+'_data_altsigbkg'
    nomsim_name = self.name+'_mc_nom'
    altsim_name = self.name+'_mc_alt'
    eff_name = 'efficiencies.json'
    if cnc:
      eff_name = 'cnc_efficiencies.json'
    if not cnc:
      if not os.path.isfile('out/'+nomdat_name+'/efficiencies.json'):
        self.data_nom_tnp_analyzer.generate_final_output()
      if not os.path.isfile('out/'+altsig_name+'/efficiencies.json'):
        self.data_altsig_tnp_analyzer.generate_final_output()
      if not os.path.isfile('out/'+altbkg_name+'/efficiencies.json'):
        self.data_altbkg_tnp_analyzer.generate_final_output()
      if not os.path.isfile('out/'+altsnb_name+'/efficiencies.json'):
        self.data_altsigbkg_tnp_analyzer.generate_final_output()
      if not os.path.isfile('out/'+nomsim_name+'/efficiencies.json'):
        #just for fit plots
        self.mc_nom_tnp_analyzer.generate_final_output() 
    else:
      if not os.path.isfile('out/'+nomdat_name+'/cnc_efficiencies.json'):
        self.data_nom_tnp_analyzer.generate_cut_and_count_output()
      if not os.path.isfile('out/'+altsig_name+'/cnc_efficiencies.json'):
        self.data_altsig_tnp_analyzer.generate_cut_and_count_output()
      if not os.path.isfile('out/'+altbkg_name+'/cnc_efficiencies.json'):
        self.data_altbkg_tnp_analyzer.generate_cut_and_count_output()
      if not os.path.isfile('out/'+altsnb_name+'/cnc_efficiencies.json'):
        self.data_altsigbkg_tnp_analyzer.generate_cut_and_count_output()
      if not os.path.isfile('out/'+nomsim_name+'/cnc_efficiencies.json'):
        #just for fit plots
        self.mc_nom_tnp_analyzer.generate_cut_and_count_output() 
    if not os.path.isfile('out/'+nomsim_name+'/cnc_efficiencies.json'):
      self.mc_nom_tnp_analyzer.generate_cut_and_count_output()
    if not os.path.isfile('out/'+altsim_name+'/cnc_efficiencies.json'):
      self.mc_alt_tnp_analyzer.generate_cut_and_count_output()
    if ((not os.path.isfile('out/'+nomdat_name+'/'+eff_name)) or
        (not os.path.isfile('out/'+altsig_name+'/'+eff_name)) or
        (not os.path.isfile('out/'+altbkg_name+'/'+eff_name)) or
        (not os.path.isfile('out/'+altsnb_name+'/'+eff_name)) or
        (not os.path.isfile('out/'+nomsim_name+'/cnc_efficiencies.json')) or
        (not os.path.isfile('out/'+altsim_name+'/cnc_efficiencies.json'))):
      print('ERROR: Could not generate individual outputs, please ensure '+
            'all fits have been performed.')
      return False
    return True

  def generate_web_output(self):
    '''Generates webpage that summarizes output
    '''
    #do checks
    sf_filename = 'out/{0}/{0}_scalefactors.json'.format(self.name)
    if not os.path.exists(sf_filename):
      print('ERROR: general output must be created before web output')
      return
    if os.path.isdir('out/web_{0}'.format(self.name)):
      print('ERROR: web output already exists, aborting')
      return

    #generate plots
    webdir = 'out/web_{0}'.format(self.name)
    os.mkdir(webdir)
    self.data_nom_tnp_analyzer.generate_web_output(webdir)
    self.data_altsig_tnp_analyzer.generate_web_output(webdir)
    self.data_altbkg_tnp_analyzer.generate_web_output(webdir)
    self.data_altsigbkg_tnp_analyzer.generate_web_output(webdir)
    self.mc_nom_tnp_analyzer.generate_web_output(webdir)
    self.generate_output(True)
    plot_types = ['eff_data','eff_mc','eff_ptbinned','eff_etabinned','sfpass',
                  'sfpass_unc','sfpass_ptbinned','sfpass_etabinned','sffail',
                  'sffail_unc','sffail_ptbinned','sffail_etabinned']
    for plot_type in plot_types:
      subprocess.run(('mv out/{0}/{0}_{1}.png out/web_{0}/'.format(self.name,
                      plot_type)).split())

    nbins = self.data_nom_tnp_analyzer.nbins

    #setup web page/directory
    with open('data/index_template.html','r') as template_file:
      index_content = template_file.read()
    index_content = index_content.replace('{0}',self.name)
    with open('{0}/index.html'.format(webdir),'w') as webindex_file:
      webindex_file.write(index_content)
    with open('{0}/fits.html'.format(webdir),'w') as fits_file:
      fits_file.write('<!DOCTYPE html>\n<html>\n<head>\n')
      fits_file.write('<link rel="stylesheet" href="../style.css">\n')
      fits_file.write('</head>\n<body>\n')
      fits_file.write('<p>This page was auto-generated.</p>')
      fits_file.write('<p>Data in black, background in blue,')
      fits_file.write(' signal+background in red. ')
      fits_file.write('All fit models included (scroll to see).</p>')
      fits_file.write('<p><a href="index.html">Back to summary</a></p>\n')
      for analyzer, desc in [(self.data_nom_tnp_analyzer, 'Nominal fits'), 
                       (self.data_altsig_tnp_analyzer, 'Altsig fits'), 
                       (self.data_altbkg_tnp_analyzer, 'Altbkg fits'), 
                       (self.data_altsigbkg_tnp_analyzer, 'Altsigbkg fits'),
                       (self.mc_nom_tnp_analyzer, 'MC constraint fits')]:
        fits_file.write('<p>{0}</p>'.format(desc))
        for ibin in range(nbins):
          fits_file.write(
            '<img src="{0}_fit{1}.png" width="765" height="256">'.format(
            analyzer.temp_name,ibin))
      fits_file.write('</body>\n</html>\n')

  def generate_output(self, is_webversion: bool=False, cnc: bool=False):
    '''Generate final SFs and histograms

    Args:
      is_webvserion: flag indicating to make png versions of plots
      cnc: flag to use cut-and-count efficiencies rather than fits
    '''

    if not is_webversion:
      if not self.generate_individual_outputs(cnc):
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
    eff_name = 'efficiencies.json'
    if cnc:
      eff_name = 'cnc_efficiencies.json'
    with open('out/'+nomdat_name+'/'+eff_name,'r') as input_file:
      eff_dat_nom = json.loads(input_file.read())
    with open('out/'+altsig_name+'/'+eff_name,'r') as input_file:
      eff_alt_sig = json.loads(input_file.read())
    with open('out/'+altbkg_name+'/'+eff_name,'r') as input_file:
      eff_alt_bkg = json.loads(input_file.read())
    with open('out/'+altsnb_name+'/'+eff_name,'r') as input_file:
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
      if not is_webversion:
        self.generate_jsons_nogap(data_eff, data_unc, mc_eff, mc_unc, 
                                  pass_sf, pass_unc, fail_sf, fail_unc)
      self.generate_summary_plots_nogap(data_eff, data_unc, mc_eff, mc_unc,
                                        pass_sf, pass_unc, fail_sf, fail_unc,
                                        is_webversion)
    elif self.binning_type=='std_gap':
      if not is_webversion:
        self.generate_jsons_gap(data_eff, data_unc, mc_eff, mc_unc, pass_sf, 
                                pass_unc, fail_sf, fail_unc)
      self.generate_summary_plots_gap(data_eff, data_unc, mc_eff, mc_unc,
                                      pass_sf, pass_unc, fail_sf, fail_unc,
                                      is_webversion)
    else:
      raise RuntimeError('Unsupported binning')

  def generate_jsons_nogap(self, data_eff: list[float], data_unc: list[float], 
                           mc_eff: list[float], mc_unc: list[float], 
                           pass_sf: list[float], pass_unc: list[float], 
                           fail_sf: list[float], fail_unc: list[float]):
    '''Generate output assuming standard binning with no gap

    Args:
      data_eff: list of data efficiencies
      data_unc: list of data uncertainties
      mc_eff: list of mc efficiencies
      mc_unc: list of mc uncertainties
      pass_sf: list of scale factors (SFs) for passing selection
      pass_unc: list of uncertainties on passing SFs
      fail_sf: list of scale factors for failing selection
      fail_unc: list of uncertainties on failing SFs
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

  def generate_jsons_gap(self, data_eff: list[float], data_unc: list[float], 
                         mc_eff: list[float], mc_unc: list[float], 
                         pass_sf: list[float], pass_unc: list[float], 
                         fail_sf: list[float], fail_unc: list[float]):
    '''Generate output assuming standard gap binning

    Args:
      data_eff: list of data efficiencies
      data_unc: list of data uncertainties
      mc_eff: list of mc efficiencies
      mc_unc: list of mc uncertainties
      pass_sf: list of scale factors (SFs) for passing selection
      pass_unc: list of uncertainties on passing SFs
      fail_sf: list of scale factors for failing selection
      fail_unc: list of uncertainties on failing SFs
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

  def generate_summary_plots_nogap(self, data_eff: list[float], 
                                   data_unc: list[float], mc_eff: list[float], 
                                   mc_unc: list[float], pass_sf: list[float], 
                                   pass_unc: list[float], fail_sf: list[float],
                                   fail_unc: list[float], 
                                   is_webversion: bool=False):
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

    Args:
      data_eff: list of data efficiencies
      data_unc: list of data uncertainties
      mc_eff: list of mc efficiencies
      mc_unc: list of mc uncertainties
      pass_sf: list of scale factors (SFs) for passing selection
      pass_unc: list of uncertainties on passing SFs
      fail_sf: list of scale factors for failing selection
      fail_unc: list of uncertainties on failing SFs
      is_webvserion: flag to indicate if output should be web-friendly
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
    gc.collect()
    gc.disable()
    file_extension = 'pdf'
    if is_webversion:
      file_extension = 'png'
    make_data_mc_graph(eta_plot_x, eta_plot_ex, eff_eta_plot_data_y, 
                       eff_eta_plot_data_ey, eff_eta_plot_mc_y, 
                       eff_eta_plot_mc_ey, 
                       'out/{0}/{0}_eff_etabinned.{1}'.format(self.name, 
                                                              file_extension), 
                       eta_plot_data_names, eta_plot_mc_names,
                       '|#eta|', eff_string, LUMI_TAGS[self.year])
    make_data_mc_graph(pt_plot_x, pt_plot_ex, eff_pt_plot_data_y, 
                       eff_pt_plot_data_ey, eff_pt_plot_mc_y, 
                       eff_pt_plot_mc_ey, 
                       'out/{0}/{0}_eff_ptbinned.{1}'.format(self.name, 
                                                             file_extension),
                       pt_plot_data_names, pt_plot_mc_names,
                       'p_{T} [GeV]', eff_string, LUMI_TAGS[self.year],True)
    make_heatmap(self.eta_bins, self.pt_bins, eff_pt_plot_data_y, 
                 'out/{0}/{0}_eff_data.{1}'.format(self.name, file_extension), 
                 '|#eta|', 
                 'p_{T} [GeV]', 
                 eff_string, LUMI_TAGS[self.year],False,True)
    make_heatmap(self.eta_bins, self.pt_bins, eff_pt_plot_mc_y, 
                 'out/{0}/{0}_eff_mc.{1}'.format(self.name, file_extension), 
                 '|#eta|', 
                 'p_{T} [GeV]',
                 eff_string, LUMI_TAGS[self.year],False,True)
    make_sf_graph(eta_plot_x, eta_plot_ex, sf_eta_plot_pass_y, 
                  sf_eta_plot_pass_ey, 
                  'out/{0}/{0}_sfpass_etabinned.{1}'.format(self.name, 
                                                            file_extension),
                  eta_plot_names, '|#eta|', 'Pass SF', LUMI_TAGS[self.year])
    make_sf_graph(eta_plot_x, eta_plot_ex, sf_eta_plot_fail_y, 
                  sf_eta_plot_fail_ey, 
                  'out/{0}/{0}_sffail_etabinned.{1}'.format(self.name, 
                                                            file_extension),
                  eta_plot_names, '|#eta|', 'Fail SF', LUMI_TAGS[self.year])
    make_sf_graph(pt_plot_x, pt_plot_ex, sf_pt_plot_pass_y, sf_pt_plot_pass_ey,
                  'out/{0}/{0}_sfpass_ptbinned.{1}'.format(self.name, 
                                                           file_extension),
                  pt_plot_names, 'p_{T} [GeV]', 'Pass SF', LUMI_TAGS[self.year],
                  True)
    make_sf_graph(pt_plot_x, pt_plot_ex, sf_pt_plot_fail_y, sf_pt_plot_fail_ey, 
                  'out/{0}/{0}_sffail_ptbinned.{1}'.format(self.name, 
                                                           file_extension),
                  pt_plot_names, 'p_{T} [GeV]', 'Fail SF', LUMI_TAGS[self.year],
                  True)
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_pass_y, 
                 'out/{0}/{0}_sfpass.{1}'.format(self.name, file_extension), 
                 '|#eta|', 
                 'p_{T} [GeV]', passsf_string, LUMI_TAGS[self.year], False, 
                 True)
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_fail_y, 
                 'out/{0}/{0}_sffail.{1}'.format(self.name, file_extension), 
                 '|#eta|', 
                 'p_{T} [GeV]', failsf_string, LUMI_TAGS[self.year], False, 
                 True)
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_pass_ey, 
                 'out/{0}/{0}_sfpass_unc.{1}'.format(self.name, 
                                                     file_extension), 
                 '|#eta|', 
                 'p_{T} [GeV]', passunc_string, LUMI_TAGS[self.year], False, 
                 True)
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_fail_ey, 
                 'out/{0}/{0}_sffail_unc.{1}'.format(self.name, 
                                                     file_extension), 
                 '|#eta|', 
                 'p_{T} [GeV]', failunc_string, LUMI_TAGS[self.year], False, 
                 True)
    gc.enable()

  def generate_summary_plots_gap(self, data_eff: list[float], 
                                 data_unc: list[float], mc_eff: list[float], 
                                 mc_unc: list[float], pass_sf: list[float], 
                                 pass_unc: list[float], fail_sf: list[float], 
                                 fail_unc: list[float], 
                                 is_webversion: bool=False):
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

    Args:
      data_eff: list of data efficiencies
      data_unc: list of data uncertainties
      mc_eff: list of mc efficiencies
      mc_unc: list of mc uncertainties
      pass_sf: list of scale factors (SFs) for passing selection
      pass_unc: list of uncertainties on passing SFs
      fail_sf: list of scale factors for failing selection
      fail_unc: list of uncertainties on failing SFs
      is_webvserion: flag to indicate if output should be web-friendly
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
    gc.collect()
    gc.disable()
    file_extension = 'pdf'
    if is_webversion:
      file_extension = 'png'
    make_data_mc_graph(eta_plot_x, eta_plot_ex, eff_eta_plot_data_y, 
                       eff_eta_plot_data_ey, eff_eta_plot_mc_y, 
                       eff_eta_plot_mc_ey, 
                       'out/{0}/{0}_eff_etabinned.{1}'.format(self.name,
                                                              file_extension), 
                       eta_plot_data_names, eta_plot_mc_names,
                       '#eta', eff_string, LUMI_TAGS[self.year])
    make_data_mc_graph(pt_plot_x, pt_plot_ex, eff_pt_plot_data_y, 
                       eff_pt_plot_data_ey, eff_pt_plot_mc_y, 
                       eff_pt_plot_mc_ey, 
                       'out/{0}/{0}_eff_ptbinned.{1}'.format(self.name,
                                                             file_extension),
                       pt_plot_data_names, pt_plot_mc_names,
                       'p_{T} [GeV]', eff_string, LUMI_TAGS[self.year],True)
    make_data_mc_graph(gappt_plot_x, gappt_plot_ex, eff_gappt_plot_data_y, 
                       eff_gappt_plot_data_ey, eff_gappt_plot_mc_y, 
                       eff_gappt_plot_mc_ey, 
                       'out/{0}/{0}_eff_gapptbinned.{1}'.format(self.name,
                           file_extension),
                       gappt_plot_data_names, gappt_plot_mc_names,
                       'p_{T} [GeV]', eff_string, LUMI_TAGS[self.year],True)
    make_heatmap(self.eta_bins, self.pt_bins, eff_pt_plot_data_y, 
                 'out/{0}/{0}_eff_data.{1}'.format(self.name,file_extension), 
                 '#eta', 'p_{T} [GeV]', 
                 eff_string, LUMI_TAGS[self.year],False,True)
    make_heatmap(self.eta_bins, self.pt_bins, eff_pt_plot_mc_y, 
                 'out/{0}/{0}_eff_mc.{1}'.format(self.name, file_extension), 
                 '#eta', 'p_{T} [GeV]',
                 eff_string, LUMI_TAGS[self.year],False,True)
    make_sf_graph(eta_plot_x, eta_plot_ex, sf_eta_plot_pass_y, 
                  sf_eta_plot_pass_ey, 
                  'out/{0}/{0}_sfpass_etabinned.{1}'.format(self.name, 
                      file_extension),
                  eta_plot_names, '#eta', 'Pass SF', LUMI_TAGS[self.year])
    make_sf_graph(eta_plot_x, eta_plot_ex, sf_eta_plot_fail_y, 
                  sf_eta_plot_fail_ey, 
                  'out/{0}/{0}_sffail_etabinned.{1}'.format(self.name, 
                      file_extension),
                  eta_plot_names, '#eta', 'Fail SF', LUMI_TAGS[self.year])
    make_sf_graph(pt_plot_x, pt_plot_ex, sf_pt_plot_pass_y, sf_pt_plot_pass_ey,
                  'out/{0}/{0}_sfpass_ptbinned.{1}'.format(self.name,
                      file_extension),
                  pt_plot_names, 'p_{T} [GeV]', 'Pass SF', LUMI_TAGS[self.year],
                  True)
    make_sf_graph(pt_plot_x, pt_plot_ex, sf_pt_plot_fail_y, sf_pt_plot_fail_ey, 
                  'out/{0}/{0}_sffail_ptbinned.{1}'.format(self.name,
                      file_extension),
                  pt_plot_names, 'p_{T} [GeV]', 'Fail SF', LUMI_TAGS[self.year],
                  True)
    make_sf_graph(gappt_plot_x, gappt_plot_ex, sf_gappt_plot_pass_y, 
                  sf_gappt_plot_pass_ey, 
                  'out/{0}/{0}_sfpass_gapptbinned.{1}'.format(self.name,
                      file_extension),
                  gappt_plot_names, 'p_{T} [GeV]', 'Pass SF', 
                  LUMI_TAGS[self.year], True)
    make_sf_graph(gappt_plot_x, gappt_plot_ex, sf_gappt_plot_fail_y, 
                  sf_gappt_plot_fail_ey, 
                  'out/{0}/{0}_sffail_gapptbinned.{1}'.format(self.name,
                      file_extension),
                  gappt_plot_names, 'p_{T} [GeV]', 'Fail SF', 
                  LUMI_TAGS[self.year], True)
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_pass_y, 
                 'out/{0}/{0}_sfpass.{1}'.format(self.name, file_extension), 
                 '#eta', 
                 'p_{T} [GeV]', passsf_string, LUMI_TAGS[self.year], False, 
                 True)
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_fail_y, 
                 'out/{0}/{0}_sffail.{1}'.format(self.name, file_extension), 
                 '#eta', 
                 'p_{T} [GeV]', failsf_string, LUMI_TAGS[self.year], False, 
                 True)
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_pass_ey, 
                 'out/{0}/{0}_sfpass_unc.{1}'.format(self.name,
                                                     file_extension), 
                 '#eta', 
                 'p_{T} [GeV]', passunc_string, LUMI_TAGS[self.year], False, 
                 True)
    make_heatmap(self.eta_bins, self.pt_bins, sf_pt_plot_fail_ey, 
                 'out/{0}/{0}_sffail_unc.{1}'.format(self.name, 
                                                     file_extension), 
                 '#eta', 
                 'p_{T} [GeV]', failunc_string, LUMI_TAGS[self.year], False, 
                 True)
    gc.enable()

  def run_interactive(self, gamma_add_gauss: bool=False):
    '''Run an interactive T&P analysis

    Args:
      gamma_add_gauss: use gamma+gaussian background model
    '''
    self.add_models(gamma_add_gauss)

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
        print('o(utput) [cnc]                 generate final outputs optionally with cut&count')
        print('w(eboutput)                    generate web output')
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
            self.data_nom_tnp_analyzer.fit_histogram_wrapper(starting_bin,
                starting_cat,self.nom_fn_name)
          elif (user_input[1] == 'nomcont'):
            #avoid accessing the same ROOT file in multiple tnp_analyzers
            self.mc_nom_tnp_analyzer.close_file()
            self.data_nom_tnp_analyzer.close_file()
            self.data_nom_tnp_analyzer.fit_histogram_wrapper(starting_bin,
                starting_cat,self.contingency_fn_name)
          elif (user_input[1] == 'alts' or user_input[1] == 'altsignal'):
            self.data_altsig_tnp_analyzer.close_file()
            self.data_altsig_tnp_analyzer.fit_histogram_wrapper(starting_bin,
                starting_cat,self.alts_fn_name,self.alts_fn_init)
          elif (user_input[1] == 'altscont'):
            self.data_altsig_tnp_analyzer.close_file()
            self.data_altsig_tnp_analyzer.fit_histogram_wrapper(starting_bin,
                starting_cat,self.contingencyalts_fn_name,self.alts_fn_init)
          elif (user_input[1] == 'altb' or user_input[1] == 'altbackground'):
            #avoid accessing the same ROOT file in multiple tnp_analyzers
            self.mc_nom_tnp_analyzer.close_file()
            self.data_altbkg_tnp_analyzer.close_file()
            self.data_altbkg_tnp_analyzer.fit_histogram_wrapper(starting_bin,
                starting_cat,self.altb_fn_name)
          elif (user_input[1] == 'altsb' or 
                user_input[1] == 'altsignalbackground'):
            self.data_altsigbkg_tnp_analyzer.close_file()
            self.data_altsigbkg_tnp_analyzer.fit_histogram_wrapper(
                starting_bin,starting_cat,self.altsb_fn_name,self.alts_fn_init)
          elif (user_input[1] == 'mc' or user_input[1] == 'mcalt'):
            self.mc_nom_tnp_analyzer.close_file()
            self.mc_nom_tnp_analyzer.fit_histogram_wrapper(starting_bin,
                starting_cat, self.alts_fn_name_sa)
          else:
            print('ERROR: unrecognized argument to f(it)')
      elif (user_input[0] == 'o' or user_input[0] == 'output'):
        do_cnc = False
        if len(user_input)>=2:
          if (user_input[1] == 'cnc'):
            do_cnc = True
        #TODO add calls to outputs for each individual piece
        self.generate_output(False, do_cnc)
      elif (user_input[0] == 'q' or user_input[0] == 'quit'):
        exit_loop = True
      elif (user_input[0] == 'c' or user_input[0] == 'clean'):
        self.clean_output()
      elif (user_input[0] == 'w' or user_input[0] == 'weboutput'):
        self.generate_web_output()
      elif (user_input[0] == 'v' or user_input[0] == 'previous'):
        for past_command in past_commands:
          print(past_command)
      else:
        print('ERROR: unrecognized command')

