#!/usr/bin/env python3
"""@package docstring
Package containing utilities for tag-and-probe analyses

For CMS members, see https://twiki.cern.ch/twiki/bin/view/CMS/ElectronScaleFactorsRun2
                 and https://indico.cern.ch/event/1288547/contributions/5414376/attachments/2654622/4597206/26052023_RMS_EGamma.pdf
                 and https://indico.cern.ch/event/1360957/contributions/5994827/attachments/2873039/5030766/SF2022update_fnal.pdf
"""

import ctypes
import json
import math
import ROOT
import os

#rng = ROOT.TRandom3()

#auxiliary functions

def clean_string(name):
  '''
  Cleans filename of illegal characters

  name   string to clean
  '''
  #follow convention from RA4draw
  cleaned_named = name.replace('.','p')
  cleaned_named = cleaned_name.replace('(','')
  cleaned_named = cleaned_name.replace(')','')
  cleaned_named = cleaned_name.replace('[','')
  cleaned_named = cleaned_name.replace(']','')
  cleaned_named = cleaned_name.replace('}','')
  cleaned_named = cleaned_name.replace('{','')
  cleaned_named = cleaned_name.replace('+','p')
  cleaned_named = cleaned_name.replace('-','m')
  cleaned_named = cleaned_name.replace('*','x')
  cleaned_named = cleaned_name.replace('/','d')
  cleaned_named = cleaned_name.replace('%','_')
  cleaned_named = cleaned_name.replace('!','n')
  cleaned_named = cleaned_name.replace('&&','__')
  cleaned_named = cleaned_name.replace('||','__')
  cleaned_named = cleaned_name.replace('==','')
  cleaned_named = cleaned_name.replace('>=','ge')
  cleaned_named = cleaned_name.replace('<=','le')
  cleaned_named = cleaned_name.replace('>','g')
  cleaned_named = cleaned_name.replace('<','l')
  cleaned_named = cleaned_name.replace('=','')
  cleaned_named = cleaned_name.replace('&','_')
  cleaned_named = cleaned_name.replace('|','_')
  cleaned_named = cleaned_name.replace('^','_')
  cleaned_named = cleaned_name.replace('~','_')
  cleaned_named = cleaned_name.replace('__','_')
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
  var = var_iterator.next()
  while var:
    var_list.append(var.GetName())
    var = var_iterator.next()
  return var_list

#run all variations including constraints from other variations(?)
#  nom, altsig, altbkg, altsigbkg, altmc
#  with outputs, compute SFs, and make pretty plots

#basic:
#input -> differential bins, fit var, fit shapes
#         allows one to run fits in some controlled way
class tnp_analyzer:
  '''
  Class used to do a T&P analysis. Example usage as follows:
  
  #initialize
  my_tnp_analysis = tnp_analyzer()

  #set input file
  my_tnp_analysis.set_input_filename('my_file.root')
  
  #set fitting variable
  my_tnp_analysis.set_fitting_variable('mll','m_{ll} [GeV]')
  
  #add binning
  my_tnp_analysis.add_nd_binning(
      [tnp_analyzer.make_nd_bin_dimension('pt','p_{T} [GeV]',[10.0,20.0,35.0,50.0,100.0]),
      tnp_analyzer.make_nd_bin_dimension('fabs(eta)','|#eta|',[0.0,0.8,1.5,2.0,2.5])]
      )

  #do interactive run
  my_tnp_analysis.run_interactive()
  '''

  def __init__(self, name):
    '''
    Default constructor

    name  string, name for this tnp_analysis, used in temp file names
    '''
    self.temp_name = name
    self.fit_var_name = ''
    self.fit_var_nbins = 0
    self.fit_var_range = (0.0,0.0)
    self.fit_var_weight = '1'
    self.preselection = '1'
    self.preselection_desc = ''
    self.measurement_variable = '1'
    self.bin_selections = []
    self.bin_names = []
    self.nbins = 0
    self.nbins_x = 0
    self.input_filenames = []
    self.tree_name = ''
    self.temp_file = None
    self.model_initializers = dict()
    self.param_initializers = dict()

  def set_input_files(self, filenames, tree_name):
    '''
    Sets input file

    filenames  list of strings, name of input files
    tree_name  string, name of TTree in file
    '''
    self.input_filenames = filenames
    self.tree_name = tree_name

  def set_fitting_variable(self, name, definition, description, nbins, 
                           var_range, weight='1'):
    '''
    Adds information about fitting variable

    name         string, name of branch in TTree or C++ expression
    description  string, name used in plots (with TLaTeX)
    nbins        int, number of bins to use for fit variable
    var_range    tuple of two floats, start and end of fit range
    weight       string, expression for weight to use
    '''
    self.fit_var_name = name
    self.fit_var_desc = description
    self.fit_var_nbins = nbins
    self.fit_var_range = var_range
    self.fit_var_weight = weight

  def set_measurement_variable(self, var):
    '''
    Sets selection efficiency to measure with tag & probe

    var  string, name of branch in TTree or C++ expression
    '''
    self.measurement_variable = var 

  def set_preselection(self, preselection, desc):
    '''
    Sets basic preselection applied to all bins

    preselection  string, selection as a C++ expression
    desc          string, description of selection in TLaTeX
    '''
    self.preselection = preselection
    self.preselection_desc = desc

  @staticmethod
  def make_nd_bin_dimension(name, description, edges):
    '''
    Generates a dimension for an even ND binning

    name         string, name of branch in TTree or C++ expression
    description  string, name used in plots (with TLaTeX)
    edges        list of floats, bin edges
    '''
    return (name, description, edges)

  def add_nd_binning(self,dimensions):
    '''
    Sets even n-dimensional binning for TnP analysis

    dimensions  list of tuples generated by make_nd_bin_dimension
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

  def add_model(self, model_name, model_initializer):
    '''
    Adds an initializer for a RooWorkspace, i.e. signal and background shapes
    with appropriate Parameters. See below for examples of model initializers
    The model must take as an argument the fitting variable (as a RooRealVar)
    and include RooRealVars nSig and nBkg representing the norm of signal and
    background as well as a RooRealVar fit_var representing the fit variable
    and a RooAbsPdf pdf_sb representing the S+B model

    model_name         string, name for this model
    model_initializer  function that returns RooWorkspace, see examples
    '''
    self.model_initializers[model_name] = model_initializer

  def add_param_initializer(self, name, param_initializer):
    '''
    Adds a parameter initializer for a RooWorkspace, i.e. a function that
    takes in the bin number, a bool representing probe pass/fail, and a 
    workspace that sets the parameter values for the workspace

    name               string, name for the parameter initializer
    param_initializer  function that sets Workspace parameters
    '''
    self.param_initializers[name] = param_initializer

  def get_binname(self, ibin):
    return clean_string(self.bin_selection[ibin])

  def make_simple_tnp_plot(self, workspace, canvas):
    '''
    Draws a simple plot for debugging fits

    workspace  T&P workspace
    canvas     TCanvas to draw plot on
    '''
    plot = workspace.var('fit_var').frame()
    workspace.data('data').plotOn(plot)
    pdf_sb.plotOn(plot,ROOT.RooFit.Name('fit_sb'))
    sig_norm_temp = workspace.var('nSig').getValV()
    bak_norm_temp = workspace.var('nBkg').getValV()
    workspace.var('nSig').setVal(0.0)
    pdf_sb.plotOn(plot,ROOT.RooFit.Normalization(bak_norm_temp/(sig_norm_temp+bak_norm_temp),ROOT.RooAbsReal.Relative),ROOT.RooFit.Name('fit_b'))
    workspace.var('nSig').setVal(sig_norm_temp)
    canvas.cd()
    plot.Draw()
    canvas.Update()

  def produce_histograms(self):
    '''
    Processes the input file(s) and generates histograms for use in fitting
    '''
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
      df_bin = df.Filter(self.bin_selection[ibin])
      df_bin_pass = df_bin.Filter(self.measurement_var)
      df_bin_fail = df_bin.Filter('!('+self.measurement_var+')')
      pass_hist_ptrs.append(df_bin_pass.Histo1D((
          'hist_pass_'+self.get_binname(ibin),
          ';'+fit_var_desc+';Events/bin',
          self.fit_var_nbins,
          self.fit_var_range[0],
          self.fit_var_range[1]),
          'tnpanalysis_fit_var',
          'tnpanalysis_fit_var_weight'))
      fail_hist_ptrs.append(df_bin_fail.Histo1D((
          'hist_fail_'+self.get_binname(ibin),
          ';'+fit_var_desc+';Events/bin',
          self.fit_var_nbins,
          self.fit_var_range[0],
          self.fit_var_range[1]),
          'tnpanalysis_fit_var',
          'tnpanalysis_fit_var_weight'))
    for ibin in range(0,self.nbins):
      pass_hist_ptrs[ibin].Write()
      fail_hist_ptrs[ibin].Write()

  def fit_histogram(self, ibin_str, pass_probe, model, param_initializer=''):
    '''
    Performs interactive fit to data

    ibin_str           string, bin to fit
    pass_probe         string, pass or fail
    model              string, model name
    param_initializer  string, parameter initializer name
    '''

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
    pass_bool = False
    pass_fail = 'fail'
    if pass_probe in ['pass','p','P']:
      pass_bool = True
      pass_fail = 'pass'
    if not model in self.model_initializers:
      print('ERROR: unknown fit model, aborting fit.')
      return
    if not param_initializer in self.param_initializers:
      if not param_initializer=='':
        print('ERROR: unknown param initializer, aborting fit.')
        return
    hist_name = 'hist_'+pass_fail+'_'+self.get_binname(ibin)
    if not self.temp_file.GetListOfKeys().Containts(hist_name):
      print('ERROR: Please produce histograms before calling fit.')
      return

    #initialize workspace
    fit_var = ROOT.RooRealVar('fit_var', self.fit_var_name, 
                              self.fit_var_range[0], self.fit_var_range[1]) 
    workspace = self.model_initializers[model]()
    data = ROOT.RooDataHist('data','fit variable', ROOT.RooArgList(fit_var), 
                            self.temp_file.Get(hist_name))
    getattr(workspace,'import')(data)
    if not param_initializer=='':
      self.param_initializers[param_initializer](ibin, pass_bool, workspace)

    #run interactive fit
    exit_loop = False
    canvas = ROOT.TCanvas()
    data = workspace.data('data')
    pdf_sb = workspace.pdf('pdf_sb')
    param_names = workspace_vars_to_list(workspace)
    param_names.remove('fit_var')
    while not exit_loop:
      user_input = input()
      user_input = user_input.split()
      if len(user_input)<1:
        continue
      elif user_input[0]=='list' or user_input[0]=='l':
        workspace_vars = workspace_vars_to_list(workspace)
        for param_name in workspace_vars:
          print(param_name+': '+str(workspace.var(param_name).getValV()))
        #workspace.Print('v') 
      elif user_input[0]=='help' or user_input[0]=='h':
        print('This is an interactive fitting session. Commands include:')
        print('f(it)                    attempt a fit')
        print('l(ist)                   display values of variables')
        print('c(onstant) <var> <value> set a value to constant or not')
        print('n(ext)                   finish fit and proceed to next bin')
        print('q(uit)                   exit interactive fitting session')
        print('q(uit)!                  exit session without saving result')
        print('j(son) <s/l> <fname>     saves or loads parameters to JSON file')
        print('r(evert)                 revert to prefit parameter values')
        print('s(et) <var> <value>      set variable <var> to <value>')
        print('w(rite) <fname>          write current canvas to a file')
      elif user_input[0]=='constant' or user_input[0]=='c':
        if len(user_input)<3:
          print('ERROR: (c)onstant takes two arguments: set <var> <value>')
          continue
        if not (user_input[1] in param_names):
          print('ERROR: unknown parameter '+user_input[1])
          continue
        constant_value = True
        if (user_input[2] == 'False' or user_input[2] == 'false'
            or user_input[2] == 'F' or user_input[2] == 'f'):
          constant_value = False
        workspace.var(user_input[1]).setConstant()
      elif user_input[0]=='set' or user_input[0]=='s':
        if len(user_input)<3:
          print('ERROR: (s)et takes two arguments: set <var> <value>')
          continue
        if not (user_input[1] in param_names):
          print('ERROR: unknown parameter '+user_input[1])
          continue
        try:
          float(user_input[2])
          workspace.var(user_input[1]).setVal(float(user_input[2]))
          self.make_simple_tnp_plot(workspace, canvas)
        except ValueError:
          print('ERROR: Unable to cast value, skipping input.')
      elif user_input[0]=='json' or user_input[0]=='j':
        if len(user_input)<3:
          print('ERROR: (j)son takes two arguments: json <s/l> <fname>')
          continue
        if not user_input[1] in ['s','save','l','load','S','L']:
          print('ERROR: (j)son second argument must be (s)ave or (l)oad')
          continue
        is_save = (user_input[1] in ['s','save','S'])
        filename = user_input[2] 
        if (not is_save) and (not os.path.exists(filename)):
          print('ERROR: File not found.')
          continue
        if is_save and os.path.exists(filename):
          print('WARNING: file already exists, enter "y" to overwrite')
          user_input_2 = input()
          if user_input_2 != 'y':
            print('Aborting write.')
            continue
        if (is_save):
          param_dict = dict()
          for param_name in param_names:
            param_dict[param_name] = workspace.var(param_name).getValV()
          with open(filename,'w') as output_file:
            output_file.write(json.dumps(param_dict))
        if (not is_save):
          with open(filename,'r') as input_file:
            param_dict = json.loads(input_file.read())
            for param_name in param_dict:
              workspace.var(param_name).setVal(param_dict[param_name])
          self.make_simple_tnp_plot(workspace, canvas)
      elif user_input[0]=='fit' or user_input[0]=='f':
        workspace.saveSnapshot('prefit',','.join(param_names))
        pdf_sb.fitTo(data)
        self.make_simple_tnp_plot(workspace, canvas)
      elif user_input[0]=='revert' or user_input[0]=='r':
        workspace.loadSnapshot('prefit')
        self.make_simple_tnp_plot(workspace, canvas)
      elif user_input[0]=='write' or user_input[0]=='w':
        if len(user_input)<3:
          print('ERROR: (w)rite takes one argument: write <fname>')
          continue
        canvas.SaveAs(user_input[1])
      elif user_input[0]=='quit!' or user_input[0]=='q!':
        exit_loop = True
      elif (user_input[0]=='quit' or user_input[0]=='q' or user_input[0]=='n' or
           user_input[0]=='next'):
        output_filename = ('out/'+self.temp_name+'/'+'fitparams_bin'+str(ibin)
                           +'_'+pass_fail+'.json')
        if os.path.exists(output_filename):
          print('WARNING: output file already exists, enter "y" to overwrite')
          user_input_2 = input()
          if user_input_2 != 'y':
            print('Aborting quit.')
            continue
        with open(output_filename,'w') as output_file:
          param_dict = dict()
          for param_name in param_names:
            param_dict[param_name] = workspace.var(param_name).getValV()
          output_file = open(filename,'w')
          output_file.write(json.dumps(param_dict))
        if user_input[0]=='n' or user_input[0]=='next':
          if (ibin != self.nbins-1):
            self.fit_histogram(self, str(ibin+1), pass_probe, model, param_initializer)
        exit_loop = True

  def generate_final_output(self):
    '''
    Generates plots of all the fits and a text file with the measured efficiencies and uncertainties
    '''
    nbins_y = self.nbins//self.nbins_x
    #left off here, calculate pixels for canvas based on number of plots and margins, set margins, etc.
    canvas = ROOT.TCanvas()
    canvas.Divide(self.nbins_x, self.nbins_y, x_margin, y_margin)
    #then combine outputs of all the jsons into a final(?) json with efficiencies and uncertainties

  def run_interactive(self):
    '''
    Begin interactive run
    '''

    #check things that are needed are already set
    if (self.nbins == 0 or self.fit_var_name == '' or len(self.input_filename) == 0):
      raise ValueError('Must initialize binning, fit variable, and input files before calling run_interactive()')

    #make output directory if it doesn't already exist
    if not os.path.isdir('out'):
      os.mkdir('out')
    elif not os.path.isdir('out/'+self.temp_name):
      print('Output directory not found, making new output directory')
      os.mkdir('out/'+self.temp_name)

    #open working (swap/temp) file
    if not os.path.isfile(self.temp_name+'.root'):
      print('Temp ROOT file not found, making new temp file.')
      self.temp_file = ROOT.TFile(self.temp_name+'.root','CREATE')
    else:
      self.temp_file = ROOT.TFile(self.temp_name+'.root','UPDATE')

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
        print('p(roduce)                                  produce histograms to fit')
        print('f(it) <bin> <pass> <model> [<initializer>] perform interactive fit to histograms')
        print('h(elp)                                     print help information')
        print('o(utput)                                   generate final output from fits (not yet implemented)')
        print('q(uit)                                     exit the interactive session')
      elif (user_input[0] == 'p' or user_input[0] == 'produce'):
        self.produce_histograms()
      elif (user_input[0] == 'f' or user_input[0] == 'fit'):
        if len(user_input)>=5:
          self.fit_histogram(user_input[1],user_input[2],user_input[3],user_input[4])
        elif len(user_input)==4:
          self.fit_histogram(user_input[1],user_input[2],user_input[3])
        else:
          print('ERROR: (f)it requires at least 3 arguments')
      elif (user_input[0] == 'q' or user_input[0] == 'quit'):
        print('Exiting interactive session and saving temp file')
        exit_loop = True
        temp_file.Close()

#initialize
my_tnp_analysis = tnp_analyzer('my_tnp_analysis')

my_tnp_analysis.set_input_filename('test.root')

#set fitting variable
my_tnp_analysis.set_fitting_variable('mll','m_{ll} [GeV]',50,(60,110))

#add binning
my_tnp_analysis.add_nd_binning(
    [tnp_analyzer.make_nd_bin_dimension('pt','p_{T} [GeV]',[10.0,20.0,35.0,50.0,100.0]),
    tnp_analyzer.make_nd_bin_dimension('fabs(eta)','|#eta|',[0.0,0.8,1.5,2.0,2.5])]
    )

#run interactive analysis
my_tnp_analysis.run_interactive()

#model initializers

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
  getattr(workspace,'import')(gauss_norm)
  getattr(workspace,'import')(cb_alphal)
  getattr(workspace,'import')(cb_nl)
  getattr(workspace,'import')(cb_alphar)
  getattr(workspace,'import')(cb_nr)
  getattr(workspace,'import')(erf_mu)
  getattr(workspace,'import')(erf_sigma)
  getattr(workspace,'import')(exp_lambda)
  getattr(workspace,'import')(bak_norm)
  param_names = ['gauss_mu','gauss_sigma','sig_norm','cb_alphal','cb_nl',
                 'cb_alphar','cb_nr','erf_mu','erf_sigma','exp_lambda',
                 'bak_norm']

  pdf_s  = ROOT.RooCrystalBall('pdf_s','pdf_s', mll, gauss_mu, gauss_sigma, cb_alphal, cb_nl, cb_alphar, cb_nr)
  pdf_b = ROOT.RooGenericPdf('pdf_b','pdf_b',
      '(TMath::Erf((@0-@1)*@2)+1.0)/2.0*exp(-1.0*@3*(@0-50.0)/80.0)',
         ROOT.RooArgList(mll, erf_mu, erf_sigma, exp_lambda))
  pdf_sb = ROOT.RooAddPdf('pdf_sb', 'pdf_sb', ROOT.RooArgList(pdf_s, pdf_b), ROOT.RooArgList(nSig, nBkg))
  getattr(workspace,'import')(pdf_sb)
  return workspace


