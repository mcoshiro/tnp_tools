#generate SFs for electrons passing H->Zgamma trigger criteria
#  phase space: pt>7, |etasc|<2.5
#  |dz|<1 cm, |dxy|<0.5 cm, fall17v2 WPL (run 2) or HZZ 2018 (run 3)

from rms_sf_analyzer import RmsSFAnalyzer
from argparse import ArgumentParser

if __name__=='__main__':

  argument_parser = ArgumentParser(prog='mutrig',
      description='Driver script for muon ID SF measurment.')
  argument_parser.add_argument('-y','--year',choices=['2016APV','2016',
      '2017','2018','2022','2022EE','2023','2023BPix'],default='2023')
  args = argument_parser.parse_args()
  
  #currently filenames aren't used here since the plots come from spark_tnp
  data_filenames = ['/data2/oshiro/ntuples/2016APV/Run2016B_L1matched.root',
                    '/data2/oshiro/ntuples/2016APV/Run2016C_L1matched.root',
                    '/data2/oshiro/ntuples/2016APV/Run2016D_L1matched.root',
                    '/data2/oshiro/ntuples/2016APV/Run2016E_L1matched.root',
                    '/data2/oshiro/ntuples/2016APV/Run2016F_L1matched.root']
  mc_filenames = ['/data2/oshiro/ntuples/2016APV/DY_NLO_L1matched.root']
  mcalt_filenames = ['/data2/oshiro/ntuples/2016APV/DY_LO_L1matched.root']
  measurement_cut = 'HtoZZMuID'
  measurement_desc = 'H#rightarrow Z#gamma Muon ID'
  preselection = 'muon_preselection'
  analyzer_name = 'muid'
  year = args.year
  pt_binning = [5.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,60.0,120.0,500.0]
  eta_binning = [0.0,0.9,1.2,2.1,2.4]

  analyzer_name = 'hzg_muid_'+year
  
  mutrig_analyzer = RmsSFAnalyzer(analyzer_name)
  mutrig_analyzer.year = year
  mutrig_analyzer.set_input_files(data_filenames,mc_filenames,mcalt_filenames,
                                  'tnpEleTrig/fitter_tree')
  mutrig_analyzer.set_fitting_variable('pair_mass','m_{#mu#Mu} [GeV]')
  mutrig_analyzer.set_measurement_variable(measurement_cut,measurement_desc)
  mutrig_analyzer.set_preselection(preselection,preselection,preselection)
  mutrig_analyzer.add_standard_binning(pt_binning,eta_binning,
                                       'mu_pt','mu_abs_eta')
  mutrig_analyzer.run_interactive()

