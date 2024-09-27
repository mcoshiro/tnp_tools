"""@package docstring
Package containing utilities for tag-and-probe analyses
"""

import array
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

CMS_COLORS = [ROOT.TColor.GetColor('#3f90da'), 
              ROOT.TColor.GetColor('#ffa90e'),
              ROOT.TColor.GetColor('#bd1f01'), 
              ROOT.TColor.GetColor('#832db6'), 
              ROOT.TColor.GetColor('#94a4a2'), 
              ROOT.TColor.GetColor('#a96b59'),
              ROOT.TColor.GetColor('#e76300'), 
              ROOT.TColor.GetColor('#b9ac70'),
              ROOT.TColor.GetColor('#717581'), 
              ROOT.TColor.GetColor('#92dadd')]

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
  var_pointer = array.array(var_type, [var])
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

def get_hist_integral_and_error(hist):
  '''
  Returns integral and associated uncertainty of TH1 as a tuple

  hist  TH1 to integrate
  '''
  uncertainty = array.array('d', [0.0])
  integral = hist.IntegralAndError(0,integral.GetXaxis().GetNbins()+1,uncertainty)
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



