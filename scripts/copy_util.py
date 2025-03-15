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

if __name__=='__main__':
  #make_summary_plots()
  #make_web_tarballs()
  #make_corr_directory()
  copy_an_plots()
