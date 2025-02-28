#generate SFs for electrons passing H->Zgamma trigger criteria
#  phase space: pt>7, |etasc|<2.5
#  |dz|<1 cm, |dxy|<0.5 cm, fall17v2 WPL (run 2) or HZZ 2018 (run 3)

from rms_sf_analyzer import RmsSFAnalyzer
from argparse import ArgumentParser

if __name__=='__main__':

  argument_parser = ArgumentParser(prog='muiso',
      description='Driver script for muon Iso SF measurment.')
  argument_parser.add_argument('-y','--year',choices=['2016APV','2016',
      '2017','2018','2022','2022EE','2023','2023BPix'],default='2023')
  argument_parser.add_argument('-c','--cut',choices=['0.1','0.15'],
                               default='0.15')
  args = argument_parser.parse_args()
  
  data_filenames = ['/data2/oshiro/ntuples/'+args.year+'/muisoskim_data0.root',
                    '/data2/oshiro/ntuples/'+args.year+'/muisoskim_data1.root']
  mc_filenames = ['/data2/oshiro/ntuples/'+args.year+'/muisoskim_dynlo0.root',
                  '/data2/oshiro/ntuples/'+args.year+'/muisoskim_dynlo1.root']
  mcalt_filenames = ['/data2/oshiro/ntuples/'+args.year+'/muisoskim_dylo0.root',
                     '/data2/oshiro/ntuples/'+args.year+'/muisoskim_dylo1.root']
  if (args.year=='2023' or args.year=='2023BPix'):
    mcalt_filenames = mc_filenames

  #cut is 0.1 for ttH lep, 0.15 for WH/ZH lep
  measurement_cut = 'probe_mu_miniso<'+args.cut
  measurement_desc = 'I_{mini}<'+args.cut
  preselection = '1'
  year = args.year
  pt_binning = [5.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,60.0,120.0,500.0]
  eta_binning = [0.0,0.9,1.2,2.1,2.4]

  cut_cleanname = args.cut.replace('.','p')
  analyzer_name = 'hzg_muiso'+cut_cleanname+'_'+year
  
  muiso_analyzer = RmsSFAnalyzer(analyzer_name)
  muiso_analyzer.year = year
  muiso_analyzer.set_input_files(data_filenames,mc_filenames,mcalt_filenames,
                                 'tree')
  muiso_analyzer.set_fitting_variable('pair_mass','m_{#mu#mu} [GeV]',
                                      weight_mc='w_lumiyearpu')
  muiso_analyzer.set_measurement_variable(measurement_cut,measurement_desc)
  muiso_analyzer.set_preselection(preselection,preselection,preselection)
  muiso_analyzer.add_standard_binning(pt_binning,eta_binning,
                                      'probe_mu_pt','probe_mu_abseta')
  muiso_analyzer.run_interactive()

