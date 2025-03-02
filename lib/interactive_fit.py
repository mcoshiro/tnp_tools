"""@package docstring
Package containing class to perform interactive fitting, split from 
tnp_analyzer
"""

import json
import os
import random
import ROOT
from tnp_utils import *
from collections.abc import Callable

ROOT.gROOT.LoadMacro('./lib/tnp_utils.cpp')

class InteractiveFit():
  """Class that manages interactive fitting

  Attributes:
    workspace: workspace used in fit
    param_names: parameter names in workspace
    plot_callback: function to plot fit
    saveplot_callback: function to save image of fit
    output_filename: name of json file to save fit output
    model_name: name of fit model used in json output
    help_info: extra information to be printed with help command
    fit_range: name of fit range to use
    fit_var: fitting variable
    data: fitting dataset
    pdf: fitting probability distribution function
    fit_status: tracks fit status
    exit_code: used to send information to caller
  """
  def __init__(self, workspace: ROOT.RooWorkspace, plot_callback: Callable, 
               saveplot_callback: Callable, output_filename: str,
               model_name: str, help_info: str='',
               fit_var_name: str='fit_var', data_name: str='data', 
               pdf_name: str='pdf_sb', fit_range: str='fitMassRange'):
    """Constructor

    Args:
      workspace: workspace for fit
      plot_callback: function that generates fit plot for user
      saveplot_callback: function that saves plot
      output_filename: name of json to save output to
      model_name: name of fit model (used for output only)
      help_info: extra info to print in help menu
      fit_var_name: name of fit variable in workspace
      data_name: name of data in workspace
      pdf_name: name of pdf in workspace
      fit_range: name of range in fitting variable to use
    """
    self.workspace = workspace
    self.param_names = workspace_vars_to_list(self.workspace)
    self.param_names.remove(fit_var_name)
    self.plot_callback = plot_callback
    self.saveplot_callback = saveplot_callback
    self.output_filename = output_filename
    self.model_name = model_name
    self.help_info = help_info
    self.fit_range = fit_range
    self.fit_var = workspace.var(fit_var_name)
    self.data = workspace.data(data_name)
    self.pdf = workspace.pdf(pdf_name)
    self.fit_status = 0
    self.exit_code = 0

  def save_output(self) -> bool:
    """Attempts to save fit output

    Returns:
      true if successful
    """
    #originally was going to save Workspace, but ROOT is too unstable, so
    #just saving fit parameter info to recreate later
    if os.path.exists(self.output_filename):
      print('WARNING: output file already exists, enter "y" to overwrite')
      user_input_2 = input()
      if user_input_2 != 'y':
        print('Aborting write.')
        return False
    param_dict = dict()
    for param_name in self.param_names:
      param_dict[param_name] = self.workspace.var(param_name).getValV()
      param_dict[param_name+'_unc'] = self.workspace.var(param_name).getError()
    param_dict['fit_status'] = self.fit_status
    param_dict['fit_model'] = self.model_name
    with open(self.output_filename,'w') as output_file:
      output_file.write(json.dumps(param_dict))
    return True

  def interactive_constant(self, user_input: list[str]) -> bool:
    """Sets a variable to be constant or nonconstant

    Args:
      user_input: user input 

    Returns:
      true if user should exit interactive session
    """
    if len(user_input)<3:
      print('ERROR: (c)onstant takes two arguments: set <var> <value>')
      return False
    if not (user_input[1] in self.param_names):
      print('ERROR: unknown parameter '+user_input[1])
      return False
    constant_value = True
    if (user_input[2] in ['True','true','T','t']):
      constant_value = True
    elif (user_input[2] in ['False','false','F','f']):
      constant_value = False
    else:
      print('ERROR: unable to cast '+user_input[2])
      return False
    self.workspace.var(user_input[1]).setConstant(constant_value)
    if len(user_input)>=4:
      self.interactive_set(['s',user_input[1],user_input[3]])
    return False

  def interactive_fit(self, user_input: list[str]) -> bool:
    """Performs fit

    Args:
      user_input: user input 

    Returns:
      true if user should exit interactive session
    """
    self.workspace.saveSnapshot('prefit',','.join(self.param_names))
    if self.fit_range != '':
      fit_result_ptr = self.pdf.fitTo(self.data,ROOT.RooFit.Save(True),
                                      ROOT.RooFit.Range(self.fit_range))
    else:
      fit_result_ptr = self.pdf.fitTo(self.data,ROOT.RooFit.Save(True))
    self.fit_status = fit_result_ptr.status()
    ROOT.free_memory_RooFitResult(fit_result_ptr)
    print('Fit status: '+str(self.fit_status))
    self.plot_callback()
    return False

  def interactive_help(self, user_input: list[str]) -> bool:
    """Prints help for interactive session

    Args:
      user_input: user input 

    Returns:
      true if user should exit interactive session
    """
    print(self.help_info)
    print('Commands include:')
    print('c(onstant) <var> <value>    set a value to constant or not')
    print('f(it)                       attempt a fit')
    print('j(son) <s/l> <fname> [pars] saves or loads parameters to JSON file')
    print('l(ist)                      display values of variables')
    print('n(ext)                      finish fit and proceed to next bin')
    print('p(reset) <s/l> <#> [pars]   save or load parameter preset')
    print('q(uit)                      exit interactive fitting session')
    print('q(uit)!                     exit session without saving result')
    print('ra(ndomize)                 randomize non-constant values')
    print('r(evert)                    revert to prefit parameter values')
    print('s(et) <var> <value>         set variable <var> to <value>')
    print('w(rite) <fname>             write current canvas to a file')
    return False

  def interactive_json(self, user_input: list[str]) -> bool:
    """Reads or writes parameter values to json file

    Args:
      user_input: user input

    Returns:
      true if user should exit interactive session
    """
    if len(user_input)<3:
      print('ERROR: (j)son takes two arguments: json <s/l> <fname>')
      return False
    if not user_input[1] in ['s','save','l','load','S','L']:
      print('ERROR: (j)son second argument must be (s)ave or (l)oad')
      return False
    is_save = (user_input[1] in ['s','save','S'])
    filename = user_input[2] 
    if (not is_save) and (not os.path.exists(filename)):
      print('ERROR: File not found.')
      return False
    if is_save and os.path.exists(filename):
      print('WARNING: file already exists, enter "y" to overwrite')
      user_input_2 = input()
      if user_input_2 != 'y':
        print('Aborting write.')
        return False
    if (is_save):
      param_dict = dict()
      if len(user_input)<4:
        for param_name in self.param_names:
          param_dict[param_name] = self.workspace.var(param_name).getValV()
      else: 
        user_param_names = user_input[3].split(',')
        for param_name in user_param_names:
          if param_name in self.param_names:
            param_dict[param_name] = self.workspace.var(param_name).getValV()
          else:
            print('Warning: parameter {} not known, skipping'.format(
                param_name))
      with open(filename,'w') as output_file:
        output_file.write(json.dumps(param_dict))
    if (not is_save):
      with open(filename,'r') as input_file:
        param_dict = json.loads(input_file.read())
        for param_name in param_dict:
          self.workspace.var(param_name).setVal(param_dict[param_name])
      self.plot_callback()
    return False

  def interactive_list(self, user_input: list[str]) -> bool:
    """Prints variables and current values for interactive session

    Args:
      user_input: user input

    Returns:
      true if user should exit interactive session
    """
    for param_name in self.param_names:
      print(param_name+': '+str(self.workspace.var(param_name).getValV()))
    #workspace.Print('v') 
    return False

  def interactive_next(self, user_input: list[str]) -> bool:
    """Finishes fit and goes to next category

    Args:
      user_input: user input

    Returns:
      true if user should exit interactive session
    """
    if not self.save_output():
      return False
    self.exit_code = 1
    return True

  def interactive_preset(self, user_input: list[str]) -> bool:
    """Saves or loads a parameter preset

    Args:
      user_input: user input

    Returns:
      true if user should exit interactive session
    """
    if len(user_input)<3:
      print('ERROR: preset requires at least 2 arguments')
      return False
    if user_input[1] in ['load','l']:
      self.interactive_json(['j','l','preset_{}.json'.format(user_input[2])])
    elif user_input[1] in ['save','s']:
      if len(user_input)<4:
        print('ERROR: preset s requires at least 3 arguments')
        return False
      param_strs = user_input[3].split(',')
      param_names = ''
      for param_str in param_strs:
        param_name, param_value_str = param_str.split(':')
        if not param_name in self.param_names:
          print('ERROR: no parameter {}'.format(param_name))
          return False
        self.interactive_set(['s',param_name,param_value_str])
        if len(param_names) != 0:
          param_names += ','
        param_names += param_name
      self.interactive_json(['j','s','preset_{}.json'.format(user_input[2]),
                            param_names])
    else:
      print('ERROR: unrecognized argument (must be (s)ave or (l)oad)')

    return False

  def interactive_quit(self, user_input: list[str]) -> bool:
    """Exits interactive fitter

    Args:
      user_input: user input

    Returns:
      true if user should exit interactive session
    """
    if not self.save_output():
      return False
    return True

  def interactive_quitnosave(self, user_input: list[str]) -> bool:
    """Exits interactive fitter without saving

    Args:
      user_input: user input

    Returns:
      true if user should exit interactive session
    """
    return True

  def interactive_randomize(self, user_input: list[str]) -> bool:
    """Randomizes values of nonconstant parameters

    Args:
      user_input: user input

    Returns:
      true if user should exit interactive session
    """
    for param_name in self.param_names:
      if not self.workspace.var(param_name).isConstant():
        param_max = self.workspace.var(param_name).getMax()
        param_min = self.workspace.var(param_name).getMin()
        self.workspace.var(param_name).setVal(random.uniform(param_min, 
                                                             param_max))
      self.plot_callback()
    return False

  def interactive_revert(self, user_input: list[str]) -> bool:
    """Returns parameters to prefit values

    Args:
      user_input: user input

    Returns:
      true if user should exit interactive session
    """
    self.workspace.loadSnapshot('prefit')
    self.plot_callback()
    return False

  def interactive_set(self, user_input: list[str]) -> bool:
    """Sets the value of a parameter

    Args:
      user_input: user input

    Returns:
      true if user should exit interactive session
    """
    if len(user_input)<3:
      print('ERROR: (s)et takes two arguments: set <var> <value>')
      return False
    if not (user_input[1] in self.param_names):
      print('ERROR: unknown parameter '+user_input[1])
      return False
    try:
      float(user_input[2])
      self.workspace.var(user_input[1]).setVal(float(user_input[2]))
      self.plot_callback()
    except ValueError:
      print('ERROR: Unable to cast value, skipping input.')
    return False

  def interactive_write(self, user_input: list[str]) -> bool:
    """Saves canvas with fit plot

    Args:
      user_input: user input

    Returns:
      true if user should exit interactive session
    """
    if len(user_input)<2:
      print('ERROR: (w)rite takes one argument: write <fname>')
      return False
    self.saveplot_callback(user_input[1])
    return False

  def run_interactive(self) -> int:
    """Run interactive fit

    Returns:
      exit status. 0 = exit, 1 = continue to next fit
    """
    function_lookup = {'c': self.interactive_constant,
                       'constant': self.interactive_constant,
                       'f': self.interactive_fit,
                       'fit': self.interactive_fit,
                       'h': self.interactive_help,
                       'help': self.interactive_help,
                       'j': self.interactive_json,
                       'json': self.interactive_json,
                       'l': self.interactive_list,
                       'list': self.interactive_list,
                       'n': self.interactive_next,
                       'next': self.interactive_next,
                       'p': self.interactive_preset,
                       'preset': self.interactive_preset,
                       'q': self.interactive_quit,
                       'quit': self.interactive_quit,
                       'q!': self.interactive_quitnosave,
                       'quit!': self.interactive_quitnosave,
                       'ra': self.interactive_randomize,
                       'randomize': self.interactive_randomize,
                       'r': self.interactive_revert,
                       'revert': self.interactive_revert,
                       's': self.interactive_set,
                       'set': self.interactive_set,
                       'w': self.interactive_write,
                       'write': self.interactive_write}
    exit_loop = False
    fit_status = 0
    while not exit_loop:
      user_input = input('>:')
      user_input = user_input.split()
      if len(user_input)<1:
        continue
      if user_input[0] in function_lookup:
        exit_loop = function_lookup[user_input[0]](user_input)
    return self.exit_code

