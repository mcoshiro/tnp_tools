import os
import shutil
import subprocess

def strip_extension(filename):
  '''
  Function that removes extension from filename if an extension exists

  Param: filename - string of filename
  Returns: string that is filename with extension removed
  '''
  if '.' in filename:
    dot_pos = [i for i in range(len(filename)) if filename[i]=='.']
    return filename[0:dot_pos[-1]]
  else:
    return filename

def merge_pdfs(input_filenames,figures_per_row,output_filename):
  '''Merge PDFs using LaTeX
  input_names      list of strings
  figures_per_row  int, figures per row in output
  output_name      string
  '''
  #calculate figure width
  figure_width = 0.99/figures_per_row
  figure_width_string = '{:.4f}'.format(figure_width)
  #write latex file
  with open('rootpdf_to_png_latexdoc.tex','w') as latex_file:
    latex_file.write('\\documentclass[10pt,oneside]{report}\n')
    latex_file.write('\\usepackage{graphicx,float}\n')
    latex_file.write('\\usepackage[active,tightpage]{preview}\n')
    latex_file.write('\\begin{document}\n')
    latex_file.write('\\begin{preview}\n')
    latex_file.write('\\begin{figure}[H]\n')
    input_index = 1
    for input_filename in input_filenames:
      latex_file.write('\\includegraphics[width='+figure_width_string)
      latex_file.write('\\textwidth]{'+input_filename+'}')
      if (input_index % figures_per_row)==0:
        latex_file.write('\n')
      else:
        latex_file.write('%\n')
      input_index += 1
    latex_file.write('\\end{figure}')
    latex_file.write('\\end{preview}\n')
    latex_file.write('\\end{document}')
  #compile latex document
  subprocess.run(['pdflatex','rootpdf_to_png_latexdoc.tex'])
  #clean up
  os.remove('rootpdf_to_png_latexdoc.tex')
  os.remove('rootpdf_to_png_latexdoc.aux')
  os.remove('rootpdf_to_png_latexdoc.log')
  shutil.copy('rootpdf_to_png_latexdoc.pdf',output_filename)
  os.remove('rootpdf_to_png_latexdoc.pdf')

