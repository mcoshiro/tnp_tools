#generate SFs for photons passing MVAID WP80
#  phase space: pt>15, |etasc|<2.5, not in overlap region
#  Fall17v2 (run2) or Winter22 (run3) photon ID WP80

from rms_sf_analyzer import RmsSFAnalyzer
from argparse import ArgumentParser

if __name__=='__main__':

  argument_parser = ArgumentParser(prog='mutrig',
      description='Driver script for muon ID SF measurment.')
  argument_parser.add_argument('-y','--year',choices=['2016APV','2016',
      '2017','2018','2022','2022EE','2023','2023BPix'],default='2016APV')
  args = argument_parser.parse_args()
  
  data_filenames = ['/data2/oshiro/ntuples/2016APV/Run2016B_L1matched.root',
                    '/data2/oshiro/ntuples/2016APV/Run2016C_L1matched.root',
                    '/data2/oshiro/ntuples/2016APV/Run2016D_L1matched.root',
                    '/data2/oshiro/ntuples/2016APV/Run2016E_L1matched.root',
                    '/data2/oshiro/ntuples/2016APV/Run2016F_L1matched.root']
  mc_filenames = ['/data2/oshiro/ntuples/2016APV/DY_NLO_L1matched.root']
  mcalt_filenames = ['/data2/oshiro/ntuples/2016APV/DY_LO_L1matched.root']
  measurement_cut = 'probe_ph_id80'
  measurement_desc = 'Fall17v2 Photon ID WP80'
  preselection = '1'
  year = args.year
  if year=='2018':
    data_filenames = ['/data2/oshiro/ntuples/2018/phwp80skim_data.root']
    mc_filenames = ['/data2/oshiro/ntuples/2018/phwp80skim_dyg.root']
    mcalt_filenames = ['/data2/oshiro/ntuples/2018/phwp80skim_dyg2.root']
  #pt_binning = [15.0,20.0]
  pt_binning = [15.0,20.0,35.0,50.0,80.0] #just for validation with EGM SFs
  eta_binning = [-2.5,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.5]

  analyzer_name = 'hzg_phid_'+year
  analyzer_name = 'hzg_phidvalidate_'+year #just for validation with EGM SFs
  
  phid_analyzer = RmsSFAnalyzer(analyzer_name)
  phid_analyzer.year = year
  phid_analyzer.set_input_files(data_filenames,mc_filenames,mcalt_filenames,
                                'tree')
  phid_analyzer.set_fitting_variable('pair_mass','m_{#mu#mu#gamma} [GeV]')
  phid_analyzer.set_measurement_variable(measurement_cut,measurement_desc)
  phid_analyzer.set_preselection(preselection,preselection,preselection)
  phid_analyzer.add_standard_binning(pt_binning,eta_binning,
                                     'probe_ph_pt','probe_ph_eta')
  phid_analyzer.run_interactive(gamma_add_gauss=True)

