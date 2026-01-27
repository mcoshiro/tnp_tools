"""@package docstring
Package containing utilities for tag-and-probe analyses
"""

from array import array
from math import sqrt
import json
import ROOT

#dictionary for converting between array type codes and ROOT codes
ARRAY_TO_ROOT_TYPE = {'b' : 'B',
                      'B' : 'b',
                      'h' : 'S',
                      'H' : 's',
                      'i' : 'I',
                      'I' : 'i',
                      'l' : 'L',
                      'L' : 'l',
                      'f' : 'F',
                      'd' : 'D'}

LUMI_TAGS = {'2016APV' : [(20,13)],
             '2016' : [(17,13)],
             '2017' : [(41,13)],
             '2018' : [(60,13)],
             '2022' : [(8,13.6)],
             '2022EE' : [(27,13.6)],
             '2023' : [(18,13.6)],
             '2023BPix' : [(10,13.6)],
             '2023BPixHole' : [(10,13.6)]}

CMS_COLORS = [ROOT.TColor.GetColor('#3f90da'), 
              ROOT.TColor.GetColor('#ffa90e'),
              ROOT.TColor.GetColor('#bd1f01'), 
              ROOT.TColor.GetColor('#832db6'), 
              ROOT.TColor.GetColor('#94a4a2'), 
              ROOT.TColor.GetColor('#a96b59'),
              ROOT.TColor.GetColor('#e76300'), 
              ROOT.TColor.GetColor('#b9ac70'),
              ROOT.TColor.GetColor('#717581'), 
              ROOT.TColor.GetColor('#92dadd'), #last official color
              ROOT.TColor.GetColor('#964a8b'),
              ROOT.TColor.GetColor('#e42536')]

def clean_string(name):
  '''
  Cleans filename of illegal characters

  name   string to clean
  '''
  #follow convention from RA4draw
  cleaned_name = name.replace('.','p')
  cleaned_name = cleaned_name.replace('(','')
  cleaned_name = cleaned_name.replace(')','')
  cleaned_name = cleaned_name.replace('[','')
  cleaned_name = cleaned_name.replace(']','')
  cleaned_name = cleaned_name.replace('}','')
  cleaned_name = cleaned_name.replace('{','')
  cleaned_name = cleaned_name.replace('+','p')
  cleaned_name = cleaned_name.replace('-','m')
  cleaned_name = cleaned_name.replace('*','x')
  cleaned_name = cleaned_name.replace('/','d')
  cleaned_name = cleaned_name.replace('%','_')
  cleaned_name = cleaned_name.replace('!','n')
  cleaned_name = cleaned_name.replace('&&','__')
  cleaned_name = cleaned_name.replace('||','__')
  cleaned_name = cleaned_name.replace('==','')
  cleaned_name = cleaned_name.replace('>=','ge')
  cleaned_name = cleaned_name.replace('<=','le')
  cleaned_name = cleaned_name.replace('>','g')
  cleaned_name = cleaned_name.replace('<','l')
  cleaned_name = cleaned_name.replace('=','')
  cleaned_name = cleaned_name.replace('&','_')
  cleaned_name = cleaned_name.replace('|','_')
  cleaned_name = cleaned_name.replace('^','_')
  cleaned_name = cleaned_name.replace('~','_')
  cleaned_name = cleaned_name.replace('__','_')
  return cleaned_name

def strip_units(name):
  '''
  Removes units in brackets from a variable name

  name  string, name of variable to modify
  '''
  bracket_pos = name.find('[')
  if (bracket_pos != -1):
    return name[:(bracket_pos-1)]
  return name

def get_units(name):
  '''
  Returns just the units from a variable name

  name  string, name of variable to get units of
  '''
  lbracket_pos = name.find('[')
  rbracket_pos = name.find(']')
  if (lbracket_pos != -1 and rbracket_pos != -1):
    return name[lbracket_pos+1:rbracket_pos]
  return ''

def workspace_vars_to_list(workspace):
  '''
  Returns python list of variables in workspace

  workspace   RooWorkspace to extract variables from
  '''
  var_list = []
  workspace_vars = workspace.allVars()
  var_iterator = workspace_vars.createIterator()
  var = var_iterator.Next()
  while var:
    var_list.append(var.GetName())
    var = var_iterator.Next()
  return var_list

def write_variable_to_root_file(var, var_type, var_name, file):
  '''
  Saves a variable to a ROOT file by saving a single entry TTree. Overwrites
  existing value

  var       int or float style variable to save
  var_type  char to denote variable type, possible values: bBhHiIlLfd
  var_name  string, name to save in file
  file      TFile in which to save
  '''
  file_keys = file.GetListOfKeys()
  if file_keys.Contains(var_name):
    file.Delete(var_name)
  if not var_type in ARRAY_TO_ROOT_TYPE:
    raise ValueError('Unsupported variable type to write to ROOT file.')
  dummy_tree = ROOT.TTree('tree','')
  var_pointer = array(var_type, [var])
  dummy_tree.Branch('var',var_pointer,'var/'+ARRAY_TO_ROOT_TYPE[var_type])
  dummy_tree.Fill()
  file.WriteObject(dummy_tree, var_name)

def read_variable_from_root_file(file, var_name):
  '''
  Reads a variable from a ROOT file that was saved with the

  write_variable_to_root_file method
  file      TFile to read from
  var_name  string, name of variable (TTree) in file
  '''
  dummy_tree = file.Get(var_name)
  var = 0
  for event in dummy_tree: 
    var = event.var
  return var

def write_multiline_latex(x,y,latex,text):
  '''
  Writes multiline python strings as TLatex

  x      float, x-coordinate in NDC
  y      float, y-coordinate in NDC
  latex  TLatex to use
  text   text to write
  '''
  split_text = text.split('\n')
  vert_offset = y
  for line in split_text:
    latex.DrawLatexNDC(x,vert_offset,line)
    vert_offset -= 0.9*latex.GetTextSize()

def zero_outside_range(hist, lower, upper):
  '''
  Zeros the bins of a histogram outside of a given range

  hist   TH1, histogram to modify
  lower  float, lower bound
  upper  float, upper bound
  '''
  for ibin in range(hist.GetXaxis().GetNbins()+1):
    if ((hist.GetXaxis().GetBinCenter(ibin) < lower) or
        (hist.GetXaxis().GetBinCenter(ibin) > upper)):
      hist.SetBinContent(ibin, 0.0)

def extend_hist(hist: ROOT.TH1D, low: float, high: float) -> ROOT.TH1D:
  '''Extends a histogram so it extends from low to high

  Args:
    hist: histogram to extend
    low: new low edge
    high: new high edge

  Returns:
    extended histogram
  '''
  #calculate new binning
  nbins_old = hist.GetNbinsX()
  low_old = hist.GetXaxis().GetBinLowEdge(1)
  bin_width = hist.GetXaxis().GetBinWidth(1)
  nbins_new = int((high-low)/bin_width)
  first_old_bin = int((low_old-low)/bin_width)
  #build new histogram
  hist_new = ROOT.TH1D(hist.GetName()+'_ext','',nbins_new,low,high)
  for ibin in range(1,nbins_old+1):
    new_bin = ibin+first_old_bin
    hist_new.SetBinContent(new_bin, hist.GetBinContent(ibin))
    hist_new.SetBinError(new_bin, hist.GetBinError(ibin))
  return hist_new

def get_hist_integral_and_error(hist):
  '''
  Returns integral and associated uncertainty of TH1 as a tuple, ignoring
  overflow bins

  hist  TH1 to integrate
  '''
  uncertainty = array('d', [0.0])
  integral = hist.IntegralAndError(1,hist.GetXaxis().GetNbins(),uncertainty)
  return (integral, uncertainty[0])

def get_bin(value, bin_edges):
  '''
  Returns bin into which value falls. -1 represents underflow while
  len(bin_edges) represents overflow

  value      float to check
  bin_edges  sorted list of floats describing binning
  '''
  if value < bin_edges[0]:
    return -1
  for ibin in range(len(bin_edges)-1):
    if value>bin_edges[ibin] and value<bin_edges[ibin+1]:
      return ibin
  return len(bin_edges)

def fix_correctionlib_json(json_texts):
  '''Fixes the format of correctionlib json created using corr.json, since 
  it is not properly formatted by default

  json_texts  dictionary generated by correctionlib .json method
  '''
  corr = []
  for json_text in json_texts:
    corr.append(json.loads(json_text))
  json_dict = {
    'schema_version' : 2,
    'description' : '',
    'corrections' : corr
  }
  return json.dumps(json_dict,indent=2)

def bin_wilson_ci(n_pass: float, n_fail: float) -> tuple[float,float]:
  '''Calculates Wilson approximation to Bernoulli mean confidence interval

  Calculates a confidence interval for the pass probability of a Bernoulli
  random variable following the Wilson approximation

  Args:
    n_pass: number of passes
    n_fail: number of fails

  Returns:
    (lower CI bound, upper CI bound)
  '''
  z = 1.0
  n_total = n_pass+n_fail
  p_est = n_pass/n_total
  a = (1.0+(z**2)/n_total)
  b = -2.0*p_est-(z**2)/n_total
  c = p_est**2
  lower = (-1.0*b-sqrt(b**2-4.0*a*c))/2.0/a
  upper = (-1.0*b+sqrt(b**2-4.0*a*c))/2.0/a
  return (lower, upper)

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


def remove_negative_bins(hist: ROOT.TH1D) -> ROOT.TH1D:
  """Returns copy of histogram with negative bins zeroed

  Args:
    hist: histogram to processo

  Returns:
    copy of hist with negative bins zeroed
  """
  clone_hist = hist.Clone(hist.GetName()+'_nonegative')
  for ibin in range(0,clone_hist.GetNbinsX()+2):
    if clone_hist.GetBinContent(ibin) < 0.0:
      clone_hist.SetBinContent(ibin, 0.0)
  return clone_hist

def transpose_onedim(arr: list, x_length: int, y_length: int):
  '''Given a 1d array containing 2D information with the index convention 
  y+x*y_length returns a 1d array with the same information but index 
  convention x+y*x_length
  '''
  if (len(arr) != x_length*y_length):
    raise ValueError(
        'Cannot transpose array whose size is not x_length*y_length.')
  for iy in range(y_length):
    for ix in range(x_length):
      ret.append(arr[iy+ix*y_length])
  ret = []
  return ret

def make_sf_graph_multibin(x: list[list[float]], ex: list[list[float]], 
    y: list[list[float]], ey: list[list[float]], name: str, 
    graph_names: list[str], x_title: str, y_title: str, 
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
  x_vals = []
  ex_vals = []
  y_vals = []
  ey_vals = []
  graphs = []
  for idata in range(len(y)):
    x_vals.append(array('d',x[idata]))
    ex_vals.append(array('d',ex[idata]))
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

def make_data_mc_graph_multibin(x: list[list[float]], ex: list[list[float]], 
  data_y: list[list[float]], data_ey: list[list[float]], 
  sim_y: list[list[float]], sim_ey: list[list[float]], name: str, 
  data_names: list[float], mc_names: list[str], x_title: str, y_title: str, 
  lumi: tuple[float, float], log_x: bool=False):
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
  x_vals = []
  ex_vals = []
  data_y_vals = []
  data_ey_vals = []
  sim_y_vals = []
  sim_ey_vals = []
  graphs = []
  for ibin in range(len(data_y)):
    x_vals.append(array('d',x))
    ex_vals.append(array('d',ex))
    data_y_vals.append(array('d',data_y[ibin]))
    data_ey_vals.append(array('d',data_ey[ibin]))
    sim_y_vals.append(array('d',sim_y[ibin]))
    sim_ey_vals.append(array('d',sim_ey[ibin]))
    graphs.append(ROOT.TGraphErrors(len(x),x_vals[ibin],data_y_vals[ibin],
                                    ex_vals[ibin],data_ey_vals[ibin]))
    graphs[-1].SetTitle(data_names[ibin])
    graphs[-1].SetLineStyle(ROOT.kSolid)
    graphs[-1].SetLineColor(CMS_COLORS[ibin])
    graphs.append(ROOT.TGraphErrors(len(x),x_vals[ibin],sim_y_vals[ibin],
                                    ex_vals[ibin],sim_ey_vals[ibin]))
    graphs[-1].SetTitle(mc_names[ibin])
    graphs[-1].SetLineStyle(ROOT.kDashed)
    graphs[-1].SetLineColor(CMS_COLORS[ibin])
  sf_plot = RplPlot()
  sf_plot.lumi_data = lumi
  for graph in graphs:
    sf_plot.plot_graph(graph, graph.GetLineColor())
  sf_plot.x_title = x_title
  sf_plot.y_title = y_title
  sf_plot.log_x = log_x
  sf_plot.draw(name)

