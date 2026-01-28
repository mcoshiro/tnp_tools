"""@package docstring
Small script to rename files for redoing fits with different binning after EGM comments
"""

import os
import subprocess

if __name__=='__main__':
  year = '2016APV'
  measurement = 'eliso0p1'
  if os.path.isdir(f'out/hzg_{measurement}_{year}_data_nom'):
    raise RuntimeError('Output folder already exists.')
  fit_cats = ['data_nom','data_altsig','data_altbkg','data_altsigbkg','mc_nom',
              'mc_alt']
  for fit_cat in fit_cats:
    out_dir = f'hzg_{measurement}_{year}_{fit_cat}'
    os.mkdir('out/'+out_dir)
    if fit_cat == 'mc_alt':
      continue
    for pf in ['pass','fail']:
      if measurement = 'elid':
        for ibin in range(32):
          subprocess.run((f'cp out_softlink/{out_dir}/fitinfo_bin{ibin+8}'
              +f'_{pf}.json out/{out_dir}/fitinfo_bin{ibin+4}_{pf}'
              +f'.json').split())
      elif measurement == 'eliso0p1' or measurement == 'eliso0p15':
        for ibin in range(40):
          subprocess.run((f'cp out_softlink/{out_dir}/fitinfo_bin{ibin}'
              +f'_{pf}.json out/{out_dir}/fitinfo_bin{ibin}_{pf}'
              +f'.json').split())
      else:
        raise ValueError('Unsupported measurement type')


