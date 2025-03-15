"""@package docstring
Package containing analyzer that can be used for tag-and-probe analyses

For CMS members, see https://twiki.cern.ch/twiki/bin/view/CMS/ElectronScaleFactorsRun2
                 and https://indico.cern.ch/event/1288547/contributions/5414376/attachments/2654622/4597206/26052023_RMS_EGamma.pdf
                 and https://indico.cern.ch/event/1360957/contributions/5994827/attachments/2873039/5030766/SF2022update_fnal.pdf
"""

import json
import math
import os
import ROOT
from functools import partial
from tnp_utils import *
from merge_pdfs import merge_pdfs
from collections.abc import Callable
from interactive_fit import InteractiveFit

ROOT.gROOT.LoadMacro('./lib/getpusf.cpp')

class TnpAnalyzer:
  '''Class used to do a T&P analysis. See scripts/example.py for example usage.

  Attributes:
    temp_name: name for analyzer
    fit_var_name: name of fitting variable
    fit_var_nbins: number of bins for fitting variable
    fit_var_range: range for fit variable
    fit_var_weight: weight for fit variable histogram
    preselection: preselection for fit histogram
    preselection_desc: description of preslection
    measurement_variable: selection of which to measure efficiency
    measurement_desc: description of selection
    bin_selections: selections defining bins (categories) of measurement
    bin_names: names of measurement bins
    nbins: number of bins for measurement
    nbins_x: number of x-bins for measurement, when bins displayed in 2D
    input_filenames: data n-tuple filenames
    tree_name: ROOT TTree name in file
    temp_file: file used to hold generated histograms
    model_initializers: collection of model initializers
    param_initializers: collection of model parameter initializers
    flag_set_input: indicates input initialized
    flag_set_fitting_variable: indicates fitting initialized
    flag_set_measurement_variable: indicated measurement initialized
    flag_set_binning: indicates binning initialized
    custom_fit_range: used to set fit range separately from hist range
    template_range: used to set range for any template fit functions
  '''

  def __init__(self, name: str):
    '''Default constructor

    Args:
      name: name for this tnp_analysis, used in temp file names
    '''
    self.temp_name = name
    self.fit_var_name = ''
    self.fit_var_nbins = 0
    self.fit_var_range = (0.0,0.0)
    self.fit_var_weight = '1'
    self.preselection = '1'
    self.preselection_desc = ''
    self.measurement_variable = '1'
    self.measurement_desc = ''
    self.bin_selections = []
    self.bin_names = []
    self.nbins = 0
    self.nbins_x = 0
    self.input_filenames = []
    self.tree_name = ''
    self.temp_file = None
    self.model_initializers = dict()
    self.param_initializers = dict()
    self.flag_set_input = False
    self.flag_set_fitting_variable = False
    self.flag_set_measurement_variable = False
    self.flag_set_binning = False
    self.custom_fit_range = None
    self.template_range = None

  def set_input_files(self, filenames: list[str], tree_name: str):
    '''Sets input file

    Args:
      filenames: list of strings, name of input files
      tree_name: string, name of TTree in file
    '''
    self.input_filenames = filenames
    self.tree_name = tree_name
    self.flag_set_input = True

  def set_fitting_variable(self, name: str, description: str, 
                           nbins: int, var_range: tuple[float,float], 
                           weight: str='1'):
    '''Adds information about fitting variable

    Args:
      name: name of branch in TTree or C++ expression
      description: name used in plots (with TLaTeX)
      nbins: number of bins to use for fit variable
      var_range: start and end of fit range
      weight: expression for weight to use
    '''
    self.fit_var_name = name
    self.fit_var_desc = description
    self.fit_var_nbins = nbins
    self.fit_var_range = var_range
    self.fit_var_weight = weight
    self.flag_set_fitting_variable = True

  def set_custom_fit_range(self, fit_range: tuple[float,float]):
    '''Sets a custom fit range separate from the measurement variable range

    Args:
      fit_range: tuple of two floats indicating range
    '''
    self.custom_fit_range = fit_range

  def set_template_range(self, template_range: tuple[float,float]):
    '''Sets a custom range for MC templates to avoid edge effects

    Args:
      template_range: tuple of two floats indicating range
    '''
    self.template_range = template_range


  def set_measurement_variable(self, var: str, desc: str=''):
    '''Sets selection efficiency to measure with tag & probe

    Args:
      var: name of branch in TTree or C++ expression
      desc: description of measurement variable
    '''
    self.measurement_variable = var 
    self.flag_set_measurement_variable = True
    if desc=='':
      desc = var
    self.measurement_desc = desc

  def set_preselection(self, preselection: str, desc: str):
    '''Sets basic preselection applied to all bins

    Args:
      preselection: selection as a C++ expression
      desc: description of selection in TLaTeX
    '''
    self.preselection = preselection
    self.preselection_desc = desc

  @staticmethod
  def make_nd_bin_dimension(name: str, description: str, edges: list[float]):
    '''Generates a dimension for an even ND binning

    Args:
      name: name of branch in TTree or C++ expression
      description: name used in plots (with TLaTeX)
      edges: bin edges

    Returns:
      tuple (name, description, edges)
    '''
    return (name, description, edges)

  def add_nd_binning(self, dimensions: list[tuple[str,str,list[float]]]):
    '''Sets even n-dimensional binning for TnP analysis

    Args:
      dimensions: list of tuples generated by make_nd_bin_dimension
    '''
    #calculate number of bins
    ndims = len(dimensions)
    self.nbins = 1
    self.nbins_x = len(dimensions[0][2])-1
    for dim in dimensions:
      self.nbins *= (len(dim[2])-1)

    #initialize each bin (which is just a selection and a name)
    for ibin in range(0,self.nbins):
      bin_selection = ''
      bin_name = ''
      remaining_dim_index = ibin
      for idim in range(len(dimensions)):
        this_dim = dimensions[ndims-idim-1]
        this_dim_len = len(this_dim[2])-1
        dim_index = remaining_dim_index % this_dim_len
        if idim != 0:
          bin_selection += '&&'
          bin_name += ', '
        bin_selection += (str(this_dim[2][dim_index])+'<'+this_dim[0])
        bin_selection += ('&&'+this_dim[0]+'<'+str(this_dim[2][dim_index+1]))
        bin_name += (str(this_dim[2][dim_index])+' < '+strip_units(this_dim[1]))
        bin_name += (' < '+str(this_dim[2][dim_index+1])+' '+get_units(this_dim[1]))
        remaining_dim_index //= this_dim_len
      self.bin_selections.append(bin_selection)
      self.bin_names.append(bin_name)
    self.flag_set_binning = True

  def add_custom_binning(self, bin_selections: list[str], 
                         bin_names: list[str]):
    '''Creates custom bins for TnP analysis

    Args:
      bin_selections: describes selection for each bin
      bin_names: names of each bin that appear in plots
    '''
    self.nbins = len(bin_selections)
    self.nbins_x = 0
    self.bin_selections = bin_selections
    self.bin_names = bin_names
    self.flag_set_binning = True

  def add_model(self, model_name: str, model_initializer: Callable):
    '''Add fitting model

    Adds an initializer for a RooWorkspace, i.e. signal and background shapes
    with appropriate Parameters. See below for examples of model initializers
    The model must take as an argument the fitting variable (as a RooRealVar)
    as well as the bin number and include RooRealVars nSig and nBkg 
    representing the norm of signal and background as well as a RooRealVar 
    fit_var representing the fit variable and a RooAbsPdf pdf_sb representing 
    the S+B model

    Args:
      model_name: string, name for this model
      model_initializer: function that returns RooWorkspace, see examples
    '''
    self.model_initializers[model_name] = model_initializer

  def add_param_initializer(self, name: str, param_initializer: Callable):
    '''Add fit model parameter initializer

    Adds a parameter initializer for a RooWorkspace, i.e. a function that
    takes in the bin number, a bool representing probe pass/fail, and a 
    workspace that sets the parameter values for the workspace

    Args:
      name: name for the parameter initializer
      param_initializer: function that sets Workspace parameters
    '''
    self.param_initializers[name] = param_initializer

  def get_binname(self, ibin: int) -> str:
    '''Get name from bin indexo

    Args:
      ibin: bin indexo

    Returns:
      bin name
    '''
    return clean_string(self.bin_selections[ibin])

  def get_fit_range(self) -> tuple[float,float]:
    '''Provides range for fit
    '''
    fit_range_lower = self.fit_var_range[0]
    fit_range_upper = self.fit_var_range[1]
    if not (self.custom_fit_range == None):
      fit_range_lower = self.custom_fit_range[0]
      fit_range_upper = self.custom_fit_range[1]
    return (fit_range_lower,fit_range_upper)

  def get_var_range(self) -> tuple[float,float]:
    '''Provides range for fitting variable
    '''
    var_range_lower = self.fit_var_range[0]
    var_range_upper = self.fit_var_range[1]
    if not (self.template_range == None):
      var_range_lower = self.template_range[0]
      var_range_upper = self.template_range[1]
    return (var_range_lower,var_range_upper)

  def make_simple_tnp_plot(self, workspace: ROOT.RooWorkspace, 
                           canvas: ROOT.TCanvas):
    '''Draws a simple plot for debugging fits

    Args:
      workspace: T&P workspace
      canvas: canvas to draw plot on
    '''
    fit_range_lower, fit_range_upper = self.get_fit_range()
    plot = workspace.var('fit_var').frame(fit_range_lower, fit_range_upper)
    workspace.data('data').plotOn(plot)
    pdf_sb = workspace.pdf('pdf_sb')
    pdf_sb.plotOn(plot,ROOT.RooFit.Name('fit_sb'))
    sig_norm_temp = workspace.var('nSig').getValV()
    bak_norm_temp = workspace.var('nBkg').getValV()
    workspace.var('nSig').setVal(0.0)
    pdf_sb.plotOn(plot,ROOT.RooFit.Normalization(bak_norm_temp/(sig_norm_temp+bak_norm_temp),ROOT.RooAbsReal.Relative),ROOT.RooFit.Name('fit_b'))
    workspace.var('nSig').setVal(sig_norm_temp)
    canvas.cd()
    plot.Draw()
    canvas.Update()

  def check_initialization(self):
    #check things that are needed are already set
    if ((not self.flag_set_input) or
        (not self.flag_set_fitting_variable) or
        (not self.flag_set_measurement_variable) or
        (not self.flag_set_binning)):
      raise ValueError('Must initialize binning, fit variable, measurement, '+
                       'input files, and model before calling run methods.')

  def initialize_files_directories(self):
    '''Generates output file and directory if they do not already exist
    '''
    #make output directory if it doesn't already exist
    if not os.path.isdir('out'):
      os.mkdir('out')
    if not os.path.isdir('out/'+self.temp_name):
      print('Output directory not found, making new output directory')
      os.mkdir('out/'+self.temp_name)

    #open working (swap/temp) file
    temp_file_name = 'out/'+self.temp_name+'/'+self.temp_name+'.root'
    if not os.path.isfile(temp_file_name):
      print('Temp ROOT file not found, making new temp file.')
      self.temp_file = ROOT.TFile(temp_file_name,'CREATE')
    else:
      if self.temp_file == None:
        self.temp_file = ROOT.TFile(temp_file_name,'UPDATE')

  def close_file(self):
    '''Closes temp file
    '''
    if (self.temp_file != None):
      self.temp_file.Close()
      self.temp_file = None

  def produce_histograms(self):
    '''Processes the input file(s) and generates histograms for use in fitting
    '''
    self.check_initialization()
    self.initialize_files_directories()
    print('Preparing file(s) for processing.')
    filenames_vec = ROOT.std.vector('string')()
    for filename in self.input_filenames:
      filenames_vec.push_back(filename)
    df = ROOT.RDataFrame(self.tree_name,filenames_vec)
    df = df.Filter(self.preselection)
    df = df.Define('tnpanalysis_fit_var', self.fit_var_name)
    df = df.Define('tnpanalysis_fit_var_weight', self.fit_var_weight)
    pass_hist_ptrs = []
    fail_hist_ptrs = []
    for ibin in range(0,self.nbins):
      df_bin = df.Filter(self.bin_selections[ibin])
      df_bin_pass = df_bin.Filter(self.measurement_variable)
      df_bin_fail = df_bin.Filter('!('+self.measurement_variable+')')
      pass_hist_ptrs.append(df_bin_pass.Histo1D((
          'hist_pass_bin{}'.format(ibin),
          ';'+self.fit_var_desc+';Events/bin',
          self.fit_var_nbins,
          self.fit_var_range[0],
          self.fit_var_range[1]),
          'tnpanalysis_fit_var',
          'tnpanalysis_fit_var_weight'))
      fail_hist_ptrs.append(df_bin_fail.Histo1D((
          'hist_fail_bin{}'.format(ibin),
          ';'+self.fit_var_desc+';Events/bin',
          self.fit_var_nbins,
          self.fit_var_range[0],
          self.fit_var_range[1]),
          'tnpanalysis_fit_var',
          'tnpanalysis_fit_var_weight'))
    print('Performing event loop. This may take a while.')
    for ibin in range(0,self.nbins):
      pass_hist_ptrs[ibin].Write()
      fail_hist_ptrs[ibin].Write()

  def fit_histogram_wrapper(self, ibin_str: str, pass_probe: str, model: str, 
                            param_initializer: str=''):
    """Performs interactive fit to data

    Args:
      ibin_str: bin to fit
      pass_probe: pass or fail
      model: model name
      param_initializer: parameter initializer name
    """
    #do initial checks
    ibin = 0
    try:
      ibin = int(ibin_str)
    except ValueError:
      print('ERROR: Unable to cast bin, skipping fit.')
      return
    if ibin<0 or ibin>=self.nbins:
      print('ERROR: Bin out of range.')
      return
    if pass_probe not in ['pass','p','P','fail','f','F']:
      print('ERROR: <pass> parameter must be (p)ass or (f)ail.')
      return
    if not model in self.model_initializers:
      print('ERROR: unknown fit model, aborting fit.')
      return
    if not param_initializer in self.param_initializers:
      if not param_initializer=='':
        print('ERROR: unknown param initializer, aborting fit.')
        return
    exit_code = 1
    while exit_code==1:
      exit_code = self.fit_histogram(ibin, pass_probe, model, 
                                     param_initializer)
      if exit_code == 1:
        if (ibin != (self.nbins-1)):
          ibin += 1
        else:
          if (pass_probe in ['pass','p','P']):
            ibin = 0
            pass_probe = 'f'
          else:
            exit_code = 0
      if exit_code==0:
        print('Exiting interactive fitter bin {} category {}.'
            .format(ibin, pass_probe))

  def fit_histogram(self, ibin: int, pass_probe: str, model: str, 
                    param_initializer: str='') -> int:
    """Performs interactive fit to data

    Args:
      ibin: bin to fit
      pass_probe: pass or fail
      model: model name
      param_initializer: parameter initializer name

    Returns:
      0 to exit, 1 to go to next fit
    """

    self.check_initialization()
    self.initialize_files_directories()
    #do initial checks
    pass_bool = False
    pass_fail = 'fail'
    if pass_probe in ['pass','p','P']:
      pass_bool = True
      pass_fail = 'pass'
    hist_name = 'hist_'+pass_fail+'_bin{}'.format(ibin)
    if not self.temp_file.GetListOfKeys().Contains(hist_name):
      print('ERROR: Please produce histograms before calling fit.')
      return 0

    #initialize workspace
    fit_range_lower, fit_range_upper = self.get_fit_range()
    var_range_lower, var_range_upper = self.get_var_range()
    fit_var = ROOT.RooRealVar('fit_var', self.fit_var_name, 
                              var_range_lower, var_range_upper)
    #fit_var.setBins(10000,"cache")
    #fit_var.setMin("cache",(self.fit_var_range[0]+fit_range_lower)/2.0)
    #fit_var.setMax("cache",(self.fit_var_range[1]+fit_range_upper)/2.0)
    #fit_var.setRange(self.fit_var_range[0], self.fit_var_range[1])
    fit_var.setRange('fitMassRange', fit_range_lower, fit_range_upper)
    workspace = self.model_initializers[model](fit_var, ibin, pass_bool)
    fit_hist = extend_hist(self.temp_file.Get(hist_name), var_range_lower,
                           var_range_upper)
    #zero_outside_range(fit_hist, fit_range_lower, fit_range_upper)
    data = ROOT.RooDataHist('data','fit variable', ROOT.RooArgList(fit_var), 
                            fit_hist)
    getattr(workspace,'import')(data)
    if not param_initializer=='':
      self.param_initializers[param_initializer](ibin, pass_bool, workspace)

    #run interactive fit
    canvas = ROOT.TCanvas()
    print('Fitting '+pass_fail+' bin {}'.format(ibin))
    plot_callback = partial(self.make_simple_tnp_plot,workspace,canvas)
    output_filename = 'out/{}/fitinfo_bin{}_{}.json'.format(self.temp_name,
                                                            str(ibin),
                                                            pass_fail)
    help_str = (('This is an interactive fitting session for bin {} '
                 +'category {}.').format(ibin, pass_probe))
    fit_session = InteractiveFit(workspace, plot_callback, canvas.SaveAs,
                                 output_filename, model, help_str)
    fit_session.interactive_fit(['f'])
    return fit_session.run_interactive()

  def draw_fit_set(self, ibin: int, filename: str, 
                   canvas_size: int=320) -> tuple[float,float]:
    '''Draws fits and write fit parameters

    Draws fits and writes fit parameters. Writes to the current gPad (global 
    ROOT pad pointer) and returns a tuple (effificiency, uncertainty) for the
    current bin

    Args:
      ibin: bin to analyze
      filename: output filename
      canvas_size: canvas height (canvas width/3) in pixels

    Returns:
       efficiency and uncertainty
    '''

    #load fit parameters and calculate efficiencies
    pass_param_dict = dict()
    fail_param_dict = dict()
    pass_param_filename = 'out/{}/fitinfo_bin{}_pass.json'.format(
        self.temp_name,ibin)
    fail_param_filename = 'out/{}/fitinfo_bin{}_fail.json'.format(
        self.temp_name,ibin)
    with open(pass_param_filename,'r') as input_file:
      pass_param_dict = json.loads(input_file.read())
    with open(fail_param_filename,'r') as input_file:
      fail_param_dict = json.loads(input_file.read())
    nsig_pass = pass_param_dict['nSig']
    nsig_fail = fail_param_dict['nSig']
    nsig_pass_unc = pass_param_dict['nSig_unc']
    nsig_fail_unc = fail_param_dict['nSig_unc']
    nsig_total = nsig_pass+nsig_fail
    #nsig_total_unc = math.hypot(nsig_pass_unc, nsig_fail_unc)
    eff = 0.0
    unc = 1.0
    if (nsig_pass > 0.0):
      eff = nsig_pass/nsig_total
      #unc = eff*math.hypot(nsig_pass_unc/nsig_pass, nsig_total_unc/nsig_total)
      unc = ((eff/nsig_pass-eff**2/nsig_pass)*nsig_pass_unc
             -(eff**2/nsig_pass)*nsig_fail_unc)
    else:
      print('WARNING: no fit passing signal in bin '+str(ibin))

    #if fragment already exists, just return efficiency
    if os.path.exists(filename):
      return (eff, unc)

    #draw fit to passing leg on middle subpad
    canvas = ROOT.TCanvas('can','can',3*canvas_size,canvas_size)
    canvas.Divide(3,1,0.0,0.0)
    canvas.cd(2)
    ROOT.gPad.SetMargin(0.1,0.1,0.1,0.1)

    #initialize workspace
    fit_range_lower, fit_range_upper = self.get_fit_range()
    var_range_lower, var_range_upper = self.get_var_range()
    fit_var_pass = ROOT.RooRealVar('fit_var', self.fit_var_name, 
                                   var_range_lower, var_range_upper) 
    fit_var_pass.setRange('fitMassRange', fit_range_lower, fit_range_upper)
    pass_workspace = self.model_initializers[pass_param_dict['fit_model']](
        fit_var_pass, ibin, True)
    hist_name = 'hist_pass_bin{}'.format(ibin)
    fit_hist_pass = extend_hist(self.temp_file.Get(hist_name), var_range_lower,
                                var_range_upper)
    data_pass = ROOT.RooDataHist('data','fit variable', ROOT.RooArgList(
        fit_var_pass),fit_hist_pass)
    getattr(pass_workspace,'import')(data_pass)
    for param_name in pass_param_dict:
      if ((not (param_name == 'fit_model' or param_name == 'fit_status')) and
          (not '_unc' in param_name)):
        pass_workspace.var(param_name).setVal(pass_param_dict[param_name])
    pass_plot = pass_workspace.var('fit_var').frame(ROOT.RooFit.Range(
        fit_range_lower,fit_range_upper),ROOT.RooFit.Title('Passing leg'))
    pass_workspace.data('data').plotOn(pass_plot, ROOT.RooFit.MarkerStyle(1))
    sig_norm_temp = pass_workspace.var('nSig').getValV()
    bak_norm_temp = pass_workspace.var('nBkg').getValV()
    pass_workspace.var('nSig').setVal(0.0)
    pass_workspace.pdf('pdf_sb').plotOn(pass_plot,
        ROOT.RooFit.Normalization(bak_norm_temp/(sig_norm_temp+bak_norm_temp),
        ROOT.RooAbsReal.Relative),ROOT.RooFit.Name('fit_b'),
        ROOT.RooFit.NormRange('fitMassRange'),
        ROOT.RooFit.LineWidth(1),ROOT.RooFit.LineColor(ROOT.kBlue))
    pass_workspace.var('nSig').setVal(sig_norm_temp)
    pass_workspace.pdf('pdf_sb').plotOn(pass_plot,ROOT.RooFit.Name('fit_sb'),
                                        ROOT.RooFit.LineWidth(1),
                                        ROOT.RooFit.NormRange('fitMassRange'),
                                        ROOT.RooFit.LineColor(ROOT.kRed))
    pass_plot.Draw()
    ROOT.gPad.Update()

    #draw fit to failing leg on right subpad
    canvas.cd(3)
    ROOT.gPad.SetMargin(0.1,0.1,0.1,0.1)
    fit_var_fail = ROOT.RooRealVar('fit_var', self.fit_var_name, 
                              var_range_lower, var_range_upper) 
    fit_var_fail.setRange('fitMassRange', fit_range_lower, fit_range_upper)
    fail_workspace = self.model_initializers[fail_param_dict['fit_model']](
        fit_var_fail, ibin, False)
    hist_name = 'hist_fail_bin{}'.format(ibin)
    fit_hist_fail = extend_hist(self.temp_file.Get(hist_name), var_range_lower,
                                var_range_upper)
    data_fail = ROOT.RooDataHist('data','fit variable', ROOT.RooArgList(
        fit_var_fail),fit_hist_fail)
    getattr(fail_workspace,'import')(data_fail)
    for param_name in fail_param_dict:
      if ((not (param_name == 'fit_model' or param_name == 'fit_status')) and
          (not '_unc' in param_name)):
        fail_workspace.var(param_name).setVal(fail_param_dict[param_name])
    fail_plot = fail_workspace.var('fit_var').frame(ROOT.RooFit.Range(
        fit_range_lower,fit_range_upper),ROOT.RooFit.Title('Failing leg'))
    fail_workspace.data('data').plotOn(fail_plot,ROOT.RooFit.MarkerStyle(1))
    sig_norm_temp = fail_workspace.var('nSig').getValV()
    bak_norm_temp = fail_workspace.var('nBkg').getValV()
    fail_workspace.var('nSig').setVal(0.0)
    fail_workspace.pdf('pdf_sb').plotOn(fail_plot,
        ROOT.RooFit.Normalization(bak_norm_temp/(sig_norm_temp+bak_norm_temp),
        ROOT.RooAbsReal.Relative),ROOT.RooFit.Name('fit_b'),
        ROOT.RooFit.NormRange('fitMassRange'),
        ROOT.RooFit.LineWidth(1),ROOT.RooFit.LineColor(ROOT.kBlue))
    fail_workspace.var('nSig').setVal(sig_norm_temp)
    fail_workspace.pdf('pdf_sb').plotOn(fail_plot,ROOT.RooFit.Name('fit_sb'),
                                        ROOT.RooFit.LineWidth(1),
                                        ROOT.RooFit.NormRange('fitMassRange'),
                                        ROOT.RooFit.LineColor(ROOT.kRed))
    fail_plot.Draw()
    ROOT.gPad.Update()

    #write text on left subpad
    canvas.cd(1)
    pad_text = self.bin_names[ibin]+'\n'
    pad_text += ('Fit status pass: {}\n'.format(pass_param_dict['fit_status']))
    pad_text += ('Fit status fail: {}\n'.format(fail_param_dict['fit_status']))
    pad_text += ('Efficiency: {:.4f} #pm {:.4f}\n'.format(eff,unc))
    pad_text += 'Fit parameters:\n'
    for param_name in pass_param_dict:
      if ((not (param_name == 'fit_model' or param_name == 'fit_status')) and
          (not '_unc' in param_name)):
        pad_text += ('{}P = {:.4f} #pm {:.4f}\n'.format(
            param_name, pass_param_dict[param_name], 
            pass_param_dict[param_name+'_unc']))
    for param_name in fail_param_dict:
      if ((not (param_name == 'fit_model' or param_name == 'fit_status')) and
          (not '_unc' in param_name)):
        pad_text += ('{}F = {:.4f} #pm {:.4f}\n'.format(
            param_name, fail_param_dict[param_name], 
            fail_param_dict[param_name+'_unc']))
    latex = ROOT.TLatex()
    latex.SetTextSize(0.025)
    write_multiline_latex(0.1,0.9,latex,pad_text)
    canvas.SaveAs(filename)

    return (eff, unc)

  def clean_output(self):
    '''Cleans output directories
    '''
    plot_filename = 'out/'+self.temp_name+'/allfits.pdf'
    effi_filename = 'out/'+self.temp_name+'/efficiencies.json'
    cnce_filename = 'out/'+self.temp_name+'/cnc_efficiencies.json'
    for remv_filename in [plot_filename, effi_filename, cnce_filename]:
      if os.path.exists(remv_filename):
        os.remove(remv_filename)
    for ibin in range(0,self.nbins):
      fragment_name = 'out/'+self.temp_name+'/allfits_fragment'+str(ibin)+'.pdf'
      if os.path.exists(fragment_name):
        os.remove(fragment_name)

  def generate_web_output(self, webdir: str):
    '''Generates web-friendly fit plots and efficiencies

    Args:
      webdir: output directory
    '''
    #do checks
    self.initialize_files_directories()
    file_keys = self.temp_file.GetListOfKeys()
    for ibin in range(0,self.nbins):
      for pass_fail in ('pass','fail'):
        param_filename = ('out/'+self.temp_name+'/fitinfo_bin'+str(ibin)+
                          '_'+pass_fail+'.json')
        if not os.path.exists(param_filename):
          print('ERROR:'+param_filename+' not found.')
          return

    #make png plots
    file_extension = '.png'
    for ibin in range(0,self.nbins):
      fragment_name = '{0}/{1}_fit{2}.png'.format(webdir,self.temp_name,ibin)
      self.draw_fit_set(ibin,fragment_name,640)

  def generate_final_output(self):
    '''Generates efficiencies and fit plots
    '''
    #do checks
    self.initialize_files_directories()
    plot_filename = 'out/'+self.temp_name+'/allfits.pdf'
    effi_filename = 'out/'+self.temp_name+'/efficiencies.json'
    if os.path.exists(plot_filename) or os.path.exists(effi_filename):
      print('ERROR: output files already exist')
      return
    file_keys = self.temp_file.GetListOfKeys()
    for ibin in range(0,self.nbins):
      for pass_fail in ('pass','fail'):
        param_filename = ('out/'+self.temp_name+'/fitinfo_bin'+str(ibin)+
                          '_'+pass_fail+'.json')
        if not os.path.exists(param_filename):
          print('ERROR:'+param_filename+' not found.')
          return

    #draw final plots and calculate final efficiencies
    #each pass/fail/parameters pad is 600x200
    nbins_x = self.nbins_x
    if nbins_x == 0:
      nbins_x = 6
    effs = []
    fragment_names = []
    file_extension = '.pdf'
    for ibin in range(0,self.nbins):
      fragment_name = ('out/'+self.temp_name+'/allfits_fragment'+str(ibin)
                       +file_extension)
      fragment_names.append(fragment_name)
      effs.append(self.draw_fit_set(ibin,fragment_name))
    fit_plot_name = 'out/'+self.temp_name+'/allfits.pdf'
    merge_pdfs(fragment_names,nbins_x,fit_plot_name)
    with open(effi_filename,'w') as output_file:
      output_file.write(json.dumps(effs))
    print('Wrote '+effi_filename)

  def generate_cut_and_count_output(self):
    '''Generates simple output file via cut-and-count method i.e. assuming 
    100% of events are signal, as is the case in signal MC
    '''

    #do checks
    self.initialize_files_directories()
    if not self.temp_file.GetListOfKeys().Contains('hist_pass_bin0'):
      print('ERROR: Please produce histograms before calling fit.')
      return
    effi_filename = 'out/'+self.temp_name+'/cnc_efficiencies.json'
    if os.path.exists(effi_filename):
      print('ERROR: output files already exist')
      return

    #calculate efficiencies
    effs = []
    for ibin in range(0,self.nbins):
      eff = 0.0
      unc = 1.0
      npass, npass_unc = get_hist_integral_and_error(
           self.temp_file.Get('hist_pass_bin{}'.format(ibin)))
      nfail, nfail_unc = get_hist_integral_and_error(
           self.temp_file.Get('hist_fail_bin{}'.format(ibin)))
      ntotal = npass+nfail
      #ntotal_unc = math.hypot(npass_unc, nfail_unc)
      if (npass > 0.0):
        eff = npass/ntotal
        #unc = eff*math.hypot(npass_unc/npass, ntotal_unc/ntotal)
        unc = ((eff/npass-eff**2/npass)*npass_unc-(eff**2/npass)*nfail_unc)
      else:
        print('WARNING: no passing signal in bin '+str(ibin))
        unc = 1.5/ntotal
      effs.append([eff,unc])
    with open(effi_filename,'w') as output_file:
      output_file.write(json.dumps(effs))
    print('Wrote '+effi_filename)

  def print_info(self):
    '''Prints basic information about T&P Analyzer
    '''
    print('Number of bin: '+str(self.nbins))
    print('Available models: ')
    for model in self.model_initializers:
      print('  '+model)
    print('Available parameter initializers: ')
    for param in self.param_initializers:
      print('  '+param)

  def run_interactive(self):
    '''Begin interactive T&P analysis run
    '''

    self.check_initialization()
    self.initialize_files_directories()

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
        print('h(elp)                                     print help information')
        print('p(roduce)                                  produce histograms to fit')
        print('f(it) <bin> <pass> <model> [<initializer>] perform interactive fit to histograms')
        print('i(nfo)                                     list info about available bins/models/initializers')
        print('c(utncount)                                generate output based on cut & count')
        print('o(utput)                                   generate final output from fits (not yet implemented)')
        print('q(uit)                                     exit the interactive session')
      elif (user_input[0] == 'p' or user_input[0] == 'produce'):
        self.produce_histograms()
      elif (user_input[0] == 'i' or user_input[0] == 'info'):
        self.print_info()
      elif (user_input[0] == 'f' or user_input[0] == 'fit'):
        print('Beginning interactive fitting session.')
        if len(user_input)>=5:
          self.fit_histogram_wrapper(user_input[1],user_input[2],user_input[3],
                                     user_input[4])
        elif len(user_input)==4:
          self.fit_histogram_wrapper(user_input[1],user_input[2],user_input[3])
        else:
          print('ERROR: (f)it requires at least 3 arguments')
      elif (user_input[0] == 'o' or user_input[0] == 'output'):
        self.generate_final_output()
      elif (user_input[0] == 'c' or user_input[0] == 'cutncount'):
        self.generate_cut_and_count_output()
      elif (user_input[0] == 'q' or user_input[0] == 'quit'):
        print('Exiting interactive session and saving temp file')
        exit_loop = True
        self.temp_file.Close()




