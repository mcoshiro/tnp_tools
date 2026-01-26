"""@package docstring
Functions implementing various custom binning schemes for T&P analyses
"""

from array import array
from correctionlib import schemav2 
from functools import partial
import ROOT

from tnp_utils import CMS_COLORS, LUMI_TAGS, make_heatmap
from rms_sf_analyzer import RmsSFAnalyzer

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

def add_standard_gap_highpt_binning(analyzer: RmsSFAnalyzer, 
                                    pt_bins: list[float], 
                                    eta_bins: list[float], 
                                    gap_pt_bins: list[float], 
                                    pt_var_name: str, 
                                    eta_var_name: str):
    '''Creates standard binning including specialized pt bins for gap region
    and combined bins for high pT

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
    for eta_lo, eta_hi in [(-1.566, -1.4442),(1.4442, 1.566)]:
      for ipt in range(len(gap_pt_bins)-1):
        bin_selections.append(f'{gap_pt_bins[ipt]}<{pt_var_name}'
            +f'&&{pt_var_name}<{gap_pt_bins[ipt+1]}&&{eta_lo}<{eta_var_name}'
            +f'&&{eta_var_name}<{eta_hi}')
        bin_selections.append(f'{gap_pt_bins[ipt]}<p_{{T}}'
            +f'<{gap_pt_bins[ipt+1]} GeV, {eta_lo}<#eta<{eta_hi}')
        if (pt_bins[ipt]>70.0):
          is_high_pt.append(True)
        else:
          is_high_pt.append(False)
    bin_selections.append(f'100<{pt_var_name}&&{pt_var_name}<500'
        +f'&&0.0<fabs({eta_var_name})&&fabs({eta_var_name})<1.5')
    bin_names.append('100<p_{T}<500 GeV, 0<|#eta|<1.5')
    bin_selections.append(f'100<{pt_var_name}&&{pt_var_name}<500'
        +f'&&1.5<fabs({eta_var_name})&&fabs({eta_var_name})<2.5')
    bin_names.append('100<p_{T}<500 GeV, 1.5<|#eta|<2.5')
    is_high_pt.append(True)
    is_high_pt.append(True)
    analyzer.add_custom_binning(bin_selections, bin_names, is_high_pt,
        'data/index_template_gap_highpt.html', 
        partial(generate_jsons_plots_gap_highpt, pt_bins=pt_bins, 
                eta_bins=eta_bins, gap_pt_bins=gap_pt_bins))

def generate_jsons_plots_gap_highpt(self, data_eff: list[float], 
    data_unc: list[float], mc_eff: list[float], mc_unc: list[float], 
    pass_sf: list[float], pass_unc: list[float], fail_sf: list[float], 
    fail_unc: list[float], name: str, desc: str, year: str, 
    pt_bins: list[float], eta_bins: list[float], gap_pt_bins: list[float]):
  '''Generate output assuming gap + only barrel/endcap for pT>100

  Generates json output along with the following plots:
    1D eta efficiency plot with MC & data
    1D pt efficiency plot with MC & data
    2D efficiency plot with MC
    2D efficiency plot with data
    1D eta SF plot pass
    1D eta SF plot fail
    1D pt SF plot pass
    1D pt SF plot fail
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
    name: analyzer name
    desc: analyzer description
    year: year
    pt_bins: pt bin boundaries for pt<100 excluding gap
    eta_bins: eta bin boundaries for pt<100 and excluding gap
    gap_pt_bins: pt bin boundaries for -1.566<|eta|<-1.4442
  '''
  #organize SFs as they will be saved in the JSON
  gapincl_eta_bins, neg_gap_idx, pos_gap_idx = add_gap_eta_bins(eta_bins)
  pass_json_sfs = []
  pass_json_uns = []
  fail_json_sfs = []
  fail_json_uns = []
  json_dat_eff = []
  json_dat_unc = []
  json_sim_eff = []
  json_sim_unc = []

  pteta_simueff = []
  pteta_dataeff = []
  pteta_sfpass = []
  pteta_sffail = []
  ptplot_x = [[], []]
  ptplot_dateff = [[], []]
  ptplot_simeff = [[], []]
  ptplot_sfpass = [[], []]
  ptplot_sffail = [[], []]
  etaplot_x = [[], []]
  etaplot_dateff = [[], []]
  etaplot_simeff = [[], []]
  etaplot_sfpass = [[], []]
  etaplot_sffail = [[], []]
  ptplot_names = []
  etaplot_names = []

  num_bins_pt = len(pt_bins)-1
  num_bins_eta = len(eta_bins)-1
  num_bins_ptgap = len(gap_pt_bins)

  for ipt in range(num_bins_pt+1):
    mean_bin_pt = (pt_bins[ipt]+pt_bins[ipt+1])/2.0
    for ieta in range(num_bins_eta+2):
      tnp_bin = -1
      if ieta < neg_gap_idx:
        if ipt == num_bins_pt:
          tnp_bin = num_bins_pt*num_bins_eta+2*num_bins_ptgap
        else:
          tnp_bin = ipt*num_bins_eta+ieta
      elif ieta == neg_gap_idx:
        tnp_bin = (num_bins_pt*num_bins_eta
                   +get_bin(mean_bin_pt, gap_pt_bins)*2)
      elif ieta > neg_gap_idx and ieta < pos_gap_idx:
        if ipt == num_bins_pt:
          tnp_bin = num_bins_pt*num_bins_eta+2*num_bins_ptgap+1
        else:
          tnp_bin = ipt*num_bins_eta+(ieta-1)
      elif ieta == pos_gap_idx:
        tnp_bin = (num_bins_pt*num_bins_eta
                   +get_bin(mean_bin_pt, gap_pt_bins)*2+1)
      elif ieta > pos_gap_idx:
        if ipt == num_bins_pt:
          tnp_bin = num_bins_pt*num_bins_eta+2*num_bins_ptgap
        else:
          tnp_bin = ipt*num_bins_eta+(ieta-2)
      pass_json_sfs.append(pass_sf[tnp_bin])
      pass_json_uns.append(pass_unc[tnp_bin])
      fail_json_sfs.append(fail_sf[tnp_bin])
      fail_json_uns.append(fail_unc[tnp_bin])
      json_dat_eff.append(data_eff[tnp_bin])
      json_dat_unc.append(data_unc[tnp_bin])
      json_sim_eff.append(mc_eff[tnp_bin])
      json_sim_unc.append(mc_unc[tnp_bin])
  pteta_dateff = transpose_onedim(json_dat_eff, num_bins_pt+1, 
                                  num_bins_eta+2)
  pteta_simeff = transpose_onedim(json_sim_eff, num_bins_pt+1, 
                                  num_bins_eta+2)
  pteta_sfpass = transpose_onedim(pass_json_sfs, num_bins_pt+1, 
                                  num_bins_eta+2)
  pteta_sffail = transpose_onedim(fail_json_sfs, num_bins_pt+1, 
                                  num_bins_eta+2)
  pteta_unpass = transpose_onedim(pass_json_uns, num_bins_pt+1, 
                                  num_bins_eta+2)
  pteta_unfail = transpose_onedim(fail_json_uns, num_bins_pt+1, 
                                  num_bins_eta+2)

  #get info for pt-binned plots
  for ieta in range(num_bins_eta):
    eta_lo = eta_bins[ieta]
    eta_hi = eta_bins[ieta+1]
    if eta_lo == -1.5:
      eta_lo = -1.4442
    elif eta_lo == 1.5:
      eta_lo = 1.566
    if eta_hi == -1.5:
      eta_hi = -1.566
    elif eta_hi == 1.5:
      eta_hi = 1.4442
    ptplot_names.append(f'{eta_lo}<#eta<{eta_hi}')
    ptplot_x[0].append([])
    ptplot_x[1].append([])
    ptplot_dateff[0].append([])
    ptplot_dateff[1].append([])
    ptplot_simeff[0].append([])
    ptplot_simeff[1].append([])
    ptplot_passsf[0].append([])
    ptplot_passsf[1].append([])
    ptplot_failsf[0].append([])
    ptplot_failsf[1].append([])
    for ipt in range(num_bins_pt):
      tnp_bin = ipt*num_bins_eta+ieta
      pt_mean = (pt_bins[ipt+1]+pt_bins[ipt])/2.0
      pt_diff = (pt_bins[ipt+1]-pt_bins[ipt])/2.0
      ptplot_x[0][ieta].append(pt_mean)
      ptplot_x[1][ieta].append(pt_diff)
      ptplot_dateff[0][ieta].append(data_eff[tnp_bin])
      ptplot_dateff[1][ieta].append(data_unc[tnp_bin])
      ptplot_simeff[0][ieta].append(mc_eff[tnp_bin])
      ptplot_simeff[1][ieta].append(mc_unc[tnp_bin])
      ptplot_passsf[0][ieta].append(pass_sf[tnp_bin])
      ptplot_passsf[1][ieta].append(pass_unc[tnp_bin])
      ptplot_failsf[0][ieta].append(fail_sf[tnp_bin])
      ptplot_failsf[1][ieta].append(fail_unc[tnp_bin])
  eta_gap_lo = [-1.566, 1.4442]
  eta_gap_hi = [-1.4442, 1.566]
  for ieta in range(2):
    ptplot_names.append(f'{eta_gap_lo[ieta]}<#eta<{eta_gap_hi[ieta]}')
    ptplot_x[0].append([])
    ptplot_x[1].append([])
    ptplot_dateff[0].append([])
    ptplot_dateff[1].append([])
    ptplot_simeff[0].append([])
    ptplot_simeff[1].append([])
    ptplot_passsf[0].append([])
    ptplot_passsf[1].append([])
    ptplot_failsf[0].append([])
    ptplot_failsf[1].append([])
    for ipt in range(num_bins_ptgap):
      tnp_bin = num_bins_pt*num_bins_eta+ipt*2+ieta
      pt_mean = (gap_pt_bins[ipt+1]+gap_pt_bins[ipt])/2.0
      pt_diff = (gap_pt_bins[ipt+1]-gap_pt_bins[ipt])/2.0
      ptplot_x[0][ieta].append(pt_mean)
      ptplot_x[1][ieta].append(pt_diff)
      ptplot_dateff[0][ieta].append(data_eff[tnp_bin])
      ptplot_dateff[1][ieta].append(data_unc[tnp_bin])
      ptplot_simeff[0][ieta].append(mc_eff[tnp_bin])
      ptplot_simeff[1][ieta].append(mc_unc[tnp_bin])
      ptplot_passsf[0][ieta].append(pass_sf[tnp_bin])
      ptplot_passsf[1][ieta].append(pass_unc[tnp_bin])
      ptplot_failsf[0][ieta].append(fail_sf[tnp_bin])
      ptplot_failsf[1][ieta].append(fail_unc[tnp_bin])
  highpt_eta_names = ['|#eta|<1.4442', '|#eta|>1.566']
  for ieta in range(2):
    ptplot_names.append(highpt_eta_names[ieta])
    tnp_bin = num_bins_pt*num_bins_eta+num_bins_ptgap*2+ieta
    pt_mean = (500.0+100.0)/2.0
    pt_diff = (500.0-100.0)/2.0
    ptplot_x[0].append([pt_mean])
    ptplot_x[1].append([pt_diff])
    ptplot_dateff[0].append([data_eff[tnp_bin]])
    ptplot_dateff[1].append([data_unc[tnp_bin]])
    ptplot_simeff[0].append([mc_eff[tnp_bin]])
    ptplot_simeff[1].append([mc_unc[tnp_bin]])
    ptplot_passsf[0].append([pass_sf[tnp_bin]])
    ptplot_passsf[1].append([pass_unc[tnp_bin]])
    ptplot_failsf[0].append([fail_sf[tnp_bin]])
    ptplot_failsf[1].append([fail_unc[tnp_bin]])

  #get info for eta-binned plots
  for ipt in range(num_bins_pt):
    pt_lo = pt_bins[ipt]
    pt_hi = pt_bins[ipt+1]
    etaplot_names.append(f'{pt_lo}<#pt<{pt_hi} GeV')
    etaplot_x[0].append([])
    etaplot_x[1].append([])
    etaplot_dateff[0].append([])
    etaplot_dateff[1].append([])
    etaplot_simeff[0].append([])
    etaplot_simeff[1].append([])
    etaplot_passsf[0].append([])
    etaplot_passsf[1].append([])
    etaplot_failsf[0].append([])
    etaplot_failsf[1].append([])
    for ieta in range(num_bins_eta):
      tnp_bin = ipt*num_bins_eta+ieta
      eta_lo = eta_bins[ieta]
      eta_hi = eta_bins[ieta+1]
      if eta_lo == -1.5:
        eta_lo = -1.4442
      elif eta_lo == 1.5:
        eta_lo = 1.566
      if eta_hi == -1.5:
        eta_hi = -1.566
      elif eta_hi == 1.5:
        eta_hi = 1.4442
      eta_mean = (eta_hi+eta_lo)/2.0
      eta_diff = (eta_hi-eta_lo)/2.0
      etaplot_x[0][ipt].append(eta_mean)
      etaplot_x[1][ipt].append(eta_diff)
      etaplot_dateff[0][ipt].append(data_eff[tnp_bin])
      etaplot_dateff[1][ipt].append(data_unc[tnp_bin])
      etaplot_simeff[0][ipt].append(mc_eff[tnp_bin])
      etaplot_simeff[1][ipt].append(mc_unc[tnp_bin])
      etaplot_passsf[0][ipt].append(pass_sf[tnp_bin])
      etaplot_passsf[1][ipt].append(pass_unc[tnp_bin])
      etaplot_failsf[0][ipt].append(fail_sf[tnp_bin])
      etaplot_failsf[1][ipt].append(fail_unc[tnp_bin])
  for ipt in range(num_bins_ptgap):
    etaplot_names.append(f'{gap_pt_bins[ipt]}<p_{T}<{gap_pt_bins[ipt+1]} GeV')
    etaplot_x[0].append([])
    etaplot_x[1].append([])
    etaplot_dateff[0].append([])
    etaplot_dateff[1].append([])
    etaplot_simeff[0].append([])
    etaplot_simeff[1].append([])
    etaplot_passsf[0].append([])
    etaplot_passsf[1].append([])
    etaplot_failsf[0].append([])
    etaplot_failsf[1].append([])
    for ieta in range(2):
      tnp_bin = num_bins_eta*num_bins_pt+ipt*2+ieta
      eta_mean = (gap_gap_hi[ieta]+eta_gap_lo[ieta])/2.0
      eta_diff = (gap_gap_hi[ieta]-eta_gap_lo[ieta])/2.0
      etaplot_x[0][ipt].append(eta_mean)
      etaplot_x[1][ipt].append(eta_diff)
      etaplot_dateff[0][ipt].append(data_eff[tnp_bin])
      etaplot_dateff[1][ipt].append(data_unc[tnp_bin])
      etaplot_simeff[0][ipt].append(mc_eff[tnp_bin])
      etaplot_simeff[1][ipt].append(mc_unc[tnp_bin])
      etaplot_passsf[0][ipt].append(pass_sf[tnp_bin])
      etaplot_passsf[1][ipt].append(pass_unc[tnp_bin])
      etaplot_failsf[0][ipt].append(fail_sf[tnp_bin])
      etaplot_failsf[1][ipt].append(fail_unc[tnp_bin])
  etaplot_names.append('100<p_{T}<500 GeV')
  tnp_bin = num_bins_eta*num_bins_pt+num_bins_ptgap*2
  etaplot_x[0].append([-2.033, 0.0, 2.033])
  etaplot_x[1].append([0.467, 1.4442, 0.467])
  etaplot_dateff[0].append([data_eff[tnp_bin+1], data_eff[tnp_bin], 
                            data_eff[tnp_bin+1]])
  etaplot_dateff[1].append([data_unc[tnp_bin+1], data_unc[tnp_bin], 
                            data_unc[tnp_bin+1]])
  etaplot_simeff[0].append([mc_eff[tnp_bin+1], mc_eff[tnp_bin], 
                            mc_eff[tnp_bin+1]])
  etaplot_simeff[1].append([mc_unc[tnp_bin+1], mc_unc[tnp_bin], 
                            mc_unc[tnp_bin+1]])
  etaplot_passsf[0].append([pass_sf[tnp_bin+1], pass_sf[tnp_bin], 
                            pass_sf[tnp_bin+1]])
  etaplot_passsf[1].append([pass_unc[tnp_bin+1], pass_unc[tnp_bin], 
                            pass_unc[tnp_bin+1]])
  etaplot_failsf[0].append([fail_sf[tnp_bin+1], fail_sf[tnp_bin], 
                            fail_sf[tnp_bin+1]])
  etaplot_failsf[1].append([fail_unc[tnp_bin+1], fail_unc[tnp_bin], 
                            fail_unc[tnp_bin+1]])

  if not os.path.isdir('out/'+name):
    print('Output directory not found, making new output directory')
    os.mkdir('out/'+name)

  #write JSON
  clib_sfs_pass = make_correction('sf_pass', 'data-MC SF', pt_bins, 
                                  gapincl_eta_bins, pass_json_sfs)
  clib_uns_pass = make_correction('unc_pass', 'data-MC unc', pt_bins, 
                                  gapincl_eta_bins, pass_json_uns)
  clib_sfs_fail = make_correction('sf_fail', 'data-MC SF', pt_bins, 
                                  gapincl_eta_bins, fail_json_sfs)
  clib_uns_fail = make_correction('unc_fail', 'data-MC unc', pt_bins, 
                                  gapincl_eta_bins, fail_json_uns)
  clib_dat_eff = make_correction('effdata', 'data eff', pt_bins, 
                                 gapincl_eta_bins, json_dat_eff)
  clib_dat_unc = make_correction('systdata', 'data unc', pt_bins, 
                                 gapincl_eta_bins, json_dat_unc)
  clib_sim_eff = make_correction('effmc', 'MC eff', pt_bins, 
                                 gapincl_eta_bins, json_sim_eff)
  clib_sim_unc = make_correction('systmc', 'MC unc', pt_bins, 
                                 gapincl_eta_bins, json_sim_unc)

  sf_filename = 'out/{0}/{0}_scalefactors.json'.format(name)
  eff_filename = 'out/{0}/{0}_efficiencies.json'.format(name)
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

  eff_string = 'Efficiency '+desc
  dateff_string = 'Data Efficiency '+desc
  simeff_string = 'MC Efficiency '+desc
  passsf_string = 'Pass SF '+desc
  failsf_string = 'Fail SF '+desc
  passun_string = 'Pass SF Unc. '+desc
  failun_string = 'Fail SF Unc. '+desc
  ptplot_dat_names = [f'Data {name}' for name in ptplot_names]
  ptplot_sim_names = [f'MC {name}' for name in ptplot_names]
  etaplot_dat_names = [f'Data {name}' for name in etaplot_names]
  etaplot_sim_names = [f'MC {name}' for name in etaplot_names]
  gc.collect()
  gc.disable()
  file_extension = 'pdf'
  if is_webversion:
    file_extension = 'png'

  for file_extension, web_tag in [('pdf',''),('png','web_')]:
    #2D plots
    for plot_data, plot_name, plot_desc in [
        (pteta_dateff, '_eff_data', dateff_string),
        (pteta_simeff, '_eff_mc',   simeff_string),
        (pteta_sfpass, '_sfpass', passsf_string),
        (pteta_sffail, '_sffail', failsf_string),
        (pteta_unpass, '_sfpass_unc', passun_string),
        (pteta_unfail, '_sffail_unc', failun_string)]:
      make_heatmap(gapincl_eta_bins, pt_bins, pteta_dataeff, 
                   f'out/{web_tag}{name}/{name}_eff_data.{file_extension}'
                   '#eta', 'p_{T} [GeV]', 
                   eff_string, LUMI_TAGS[year],False,True)
    #1D plots
    make_data_mc_graph_multibin(ptplot_x[0], ptplot_x[1], 
        ptplot_dateff_y[0], ptplot_dateff_y[1], ptplot_simeff_y[0],
        ptplot_simeff_y[1], 
        f'out/{web_tag}{name}/{name}_eff_ptbinned.{file_extension}',
        ptplot_dat_names, ptplot_sim_names, 'p_{T} [GeV]', eff_string,
        LUMI_TAGS[years], True)
    make_data_mc_graph_multibin(etaplot_x[0], etaplot_x[1], 
        etaplot_dateff_y[0], etaplot_dateff_y[1], etaplot_simeff_y[0],
        etaplot_simeff_y[1], 
        f'out/{web_tag}{name}/{name}_eff_etabinned.{file_extension}',
        etaplot_dat_names, etaplot_sim_names, '#eta', eff_string,
        LUMI_TAGS[years], False)
    make_sf_graph_multibin(ptplot_x[0], ptplot_x[1], 
        ptplot_sfpass_y[0], ptplot_sfpass_y[1], 
        f'out/{web_tag}{name}/{name}_sfpass_ptbinned.{file_extension}',
        ptplot_names, 'p_{T} [GeV]', passsf_string, LUMI_TAGS[years], True)
    make_sf_graph_multibin(ptplot_x[0], ptplot_x[1], 
        ptplot_sffail_y[0], ptplot_sffail_y[1], 
        f'out/{web_tag}{name}/{name}_sffail_ptbinned.{file_extension}',
        ptplot_names, 'p_{T} [GeV]', failsf_string, LUMI_TAGS[years], True)
    make_sf_graph_multibin(etaplot_x[0], etaplot_x[1], 
        etaplot_sfpass_y[0], etaplot_sfpass_y[1], 
        f'out/{web_tag}{name}/{name}_sfpass_etabinned.{file_extension}',
        etaplot_names, 'p_{T} [GeV]', passsf_string, LUMI_TAGS[years], False)
    make_sf_graph_multibin(etaplot_x[0], etaplot_x[1], 
        etaplot_sffail_y[0], etaplot_sffail_y[1], 
        f'out/{web_tag}{name}/{name}_sffail_etabinned.{file_extension}',
        etaplot_names, 'p_{T} [GeV]', failsf_string, LUMI_TAGS[years], False)

  gc.enable()
