"""@package docstring
Package containing utilities for tag-and-probe analyses
"""

import array
import ROOT

#dictionary for converting between array type codes and ROOT codes
array_to_root_type = dict()
array_to_root_type['b'] = 'B'
array_to_root_type['B'] = 'b'
array_to_root_type['h'] = 'S'
array_to_root_type['H'] = 's'
array_to_root_type['i'] = 'I'
array_to_root_type['I'] = 'i'
array_to_root_type['l'] = 'L'
array_to_root_type['L'] = 'l'
array_to_root_type['f'] = 'F'
array_to_root_type['d'] = 'D'

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
  if not var_type in array_to_root_type:
    raise ValueError('Unsupported variable type to write to ROOT file.')
  dummy_tree = ROOT.TTree('tree','')
  var_pointer = array.array(var_type, [var])
  dummy_tree.Branch('var',var_pointer,'var/'+array_to_root_type[var_type])
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
    vert_offset -= 0.8*latex.GetTextSize()


