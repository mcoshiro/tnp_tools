"""Hacky script used to copy plots etc. around
"""

import subprocess 

years_run2 = ['2016APV','2016','2017','2018']
years_run3 = ['2022','2022EE','2023','2023BPix']
#measurements = ['elid','eltrig27','eltrig23','eltrig12','eliso0p15','eliso0p1']
measurements = ['elid','eltrig27','eltrig23','eltrig12','eliso0p15','eliso0p1',
                'muid','mutrig24','mutrig17','mutrig8','muiso0p15','muiso0p1']

singleel_names = {'2016APV': 'eltrig27',
                  '2016': 'eltrig27',
                  '2017': 'eltrig32',
                  '2018': 'eltrig32',
                  '2022': 'eltrig30',
                  '2022EE': 'eltrig30',
                  '2023': 'eltrig30',
                  '2023BPix': 'eltrig30',
                  '2023BPixHole': 'eltrig30'}

singlemu_names = {'2016APV': 'mutrig24',
                  '2016': 'mutrig24',
                  '2017': 'mutrig27',
                  '2018': 'mutrig24',
                  '2022': 'mutrig24',
                  '2022EE': 'mutrig24',
                  '2023': 'mutrig24',
                  '2023BPix': 'mutrig24'}

pico_years = {'2016APV': '2016preVFP_UL',
              '2016': '2016postVFP_UL',
              '2017': '2017_UL',
              '2018': '2018_UL',
              '2022': '2022',
              '2022EE': '2022EE',
              '2023': '2023',
              '2023BPix': '2023BPix',
              '2023BPixHole': '2023BPix'}

def fix_measurement_name(name,year):
  if name==measurements[1]:
    return singleel_names[year]
  return name

def make_summary_plots():
  for measurement in measurements:
    for run, years in [('run2',years_run2),('run3',years_run3)]:
      for bins in ['','gap']:
        command = 'python3 /homes/oshiro/analysis/small_phys_utils/scripts/rootpdf_to_png_new.py'
        nrow = 4
        if run=='run3' and bins=='gap':
          nrow = 5
        command += ' -f {} -o {}_{}{}_summary.pdf'.format(nrow,measurement,run,
                                                         bins)
        for year in years:
          meas_name = fix_measurement_name(measurement,year)
          command += (' out/hzg_{0}_{1}/hzg_{0}_{1}_eff_{2}ptbinned.pdf'.format(
              meas_name, year,bins))
        if run=='run3' and bins=='gap':
          command += (' out/hzg_{0}_2023BPixHole/hzg_{0}_2023BPixHole_eff_ptbinned.pdf'.format(
              meas_name, year))
        for year in years:
          meas_name = fix_measurement_name(measurement,year)
          command += (' out/hzg_{0}_{1}/hzg_{0}_{1}_sfpass_{2}ptbinned.pdf'
              .format(meas_name, year,bins))
        if run=='run3' and bins=='gap':
          command += (' out/hzg_{0}_2023BPixHole/hzg_{0}_2023BPixHole_sfpass_ptbinned.pdf'.format(
              meas_name, year))
        subprocess.run(command.split())

def make_web_tarballs():
  all_years = years_run2 + years_run3 + ['2023BPixHole']
  for measurement in measurements:
    for year in all_years:
      meas_name = fix_measurement_name(measurement,year)
      command = ('tar -czvf hzg_{0}_{1}.tgz out/web_hzg_{0}_{1}'.format(
          meas_name,year))
      subprocess.run(command.split())

def make_corr_directory():
  years = years_run2 + years_run3
  subprocess.run('mkdir hzg_corrections'.split())
  for year in years:
    subprocess.run('mkdir hzg_corrections/{}'.format(year).split())
    for measurement in measurements:
      meas_name = measurement
      if meas_name=='muid':
        if not year in ['2023','2023BPix']:
          continue
      elif meas_name=='eltrig27':
        meas_name = singleel_names[year]
      elif meas_name=='mutrig24':
        meas_name = singlemu_names[year]
      corr_type = 'efficiencies'
      if meas_name in ['elid','muid']:
        corr_type = 'scalefactors'
      cmd = ('cp out/hzg_{0}_{1}/hzg_{0}_{1}_{2}.json hzg_corrections/{1}/'
             .format(meas_name,year,corr_type))
      subprocess.run(cmd.split())
  #BPixHole
  year = '2023BPixHole'
  for measurement in measurements:
    if not 'el' in measurement:
      continue
    meas_name = measurement
    if meas_name=='eltrig27':
      meas_name = singleel_names[year]
    corr_type = 'efficiencies'
    if meas_name in ['elid','muid']:
      corr_type = 'scalefactors'
    cmd = ('cp out/hzg_{0}_{1}/hzg_{0}_{1}_{2}.json hzg_corrections/2023BPix/'
           .format(meas_name,year,corr_type))
    subprocess.run(cmd.split())

def copy_an_plots():
  years = years_run2 + years_run3
  for year in years:
    for measurement in measurements:
      meas_name = measurement
      if meas_name=='muid':
        if not year in ['2023','2023BPix']:
          continue
      elif meas_name=='eltrig27':
        meas_name = singleel_names[year]
      elif meas_name=='mutrig24':
        meas_name = singlemu_names[year]
      cmd = (('cp out/hzg_{0}_{1}/hzg_{0}_{1}_eff_ptbinned.pdf '
             +'../AN-22-027/figs/App-Corrections/').format(meas_name,year))
      subprocess.run(cmd.split())

def copy_to_n2p():
  years = years_run2 + years_run3
  for year in years:
    for measurement in measurements:
      meas_name = measurement
      if meas_name=='muid':
        if not year in ['2023','2023BPix']:
          continue
      elif meas_name=='eltrig27':
        meas_name = singleel_names[year]
      elif meas_name=='mutrig24':
        meas_name = singlemu_names[year]
      corr_type = 'efficiencies'
      if meas_name in ['elid','muid']:
        corr_type = 'scalefactors'
      cmd = (('cp out/hzg_{0}_{1}/hzg_{0}_{1}_{2}.json '
             +'../nano2pico/data/zgamma/{3}/')
             .format(meas_name,year,corr_type,pico_years[year]))
      subprocess.run(cmd.split())
  #BPixHole
  year = '2023BPixHole'
  for measurement in measurements:
    if not 'el' in measurement:
      continue
    meas_name = measurement
    if meas_name=='eltrig27':
      meas_name = singleel_names[year]
    corr_type = 'efficiencies'
    if meas_name in ['elid','muid']:
      corr_type = 'scalefactors'
    cmd = (('cp out/hzg_{0}_{1}/hzg_{0}_{1}_{2}.json '
           +'../nano2pico/data/zgamma/{3}/')
           .format(meas_name,year,corr_type,pico_years[year]))
    subprocess.run(cmd.split())

def make_rebin_page():
  years = years_run2 + years_run3
  index_page_text = '<!DOCTYPE html>\n'
  index_page_text += '<html>\n'
  index_page_text += '<head>\n'
  index_page_text += '<link rel="stylesheet" href="../style.css">\n'
  index_page_text += '</head>\n'
  index_page_text += '<body>\n'
  index_page_text += '<h1>Rebinned electron trigger efficiency and scale '
  index_page_text += 'factor summary plots</h1>\n'
  index_page_text += '<p>This page was auto-generated.</p>\n'
  index_page_text += '<p><a href="../index.html">Back to top</a></p>\n'
  for year in years:
    for meas_name in ['eltrig27','eltrig23','eltrig12']:
      if meas_name=='eltrig27':
        meas_name = singleel_names[year]
      meas_desc = 'Electron trigger '+meas_name[-2:]
      index_page_text += f'<h1>{year} {meas_desc}</h1>'
      for plot_type, desc in [
          ('eff_data','Data efficiency as a function of pT and eta.'),
          ('eff_mc','MC efficiency as a function of pT and eta.'),
          ('eff_ptbinned','Data (solid) and MC (dashed) efficiency as a '
                          'function of pT for eta bins (colors).'),
          ('eff_etabinned','Data (solid) and MC (dashed) efficiency as a '
                           'function of eta for pT bins (colors).'),
          ('sfpass','Data/MC efficiency ratio (scale factor) as a function '
                    'of pT and eta.'),
          ('sfpass_unc','Data/MC efficiency ratio (scale factor) uncertainty'
                        'as a function of pT and eta.'),
          ('sfpass_ptbinned','Data/MC efficiency ratio (scale factor) as a '
                             'function of pT for different eta bins '
                             '(colors).'),
          ('sfpass_etabinned','Data/MC efficiency ratio (scale factor) as a '
                              'function of eta for different pT bins '
                              '(colors).'),
          ('sffail','(1-data efficiency)/(1-MC efficiency) as a function '
                    'of pT and eta.'),
          ('sffail_unc','(1-data efficiency)/(1-MC efficiency) uncertainty'
                        'as a function of pT and eta.'),
          ('sffail_ptbinned','(1-data efficiency)/(1-MC efficiency) as a '
                             'function of pT for different eta bins '
                             '(colors).'),
          ('sffail_etabinned','(1-data efficiency)/(1-MC efficiency) as a '
                              'function of eta for different pT bins '
                              '(colors).')]:
        plot_name = f'hzg_{meas_name}_{year}_{plot_type}.png'
        cmd = (f'cp out/hzg_{meas_name}_{year}_rebinned/{plot_name} '
               +'rebinned_page/')
        subprocess.run(cmd.split())
        index_page_text += f'<figure><img src="{plot_name}" width="256"'
        index_page_text += f' height="256"><figcaption>(Above) {desc}'
        index_page_text += '</figcaption></figure>\n'
  index_page_text += '</body>\n'
  index_page_text += '</html>'
  with open('rebinned_page/index.html','w') as html_file:
    html_file.write(index_page_text)

if __name__=='__main__':
  #make_summary_plots()
  #make_web_tarballs()
  #make_corr_directory()
  #copy_an_plots()
  #copy_to_n2p()
  make_rebin_page()
