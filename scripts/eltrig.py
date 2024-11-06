#generate SFs for electrons passing H->Zgamma trigger criteria
#  phase space: pt>7, |etasc|<2.5
#  |dz|<1 cm, |dxy|<0.5 cm, fall17v2 WPL (run 2) or HZZ 2018 (run 3)

from rms_sf_analyzer import RmsSFAnalyzer
from argparse import ArgumentParser

if __name__=='__main__':

  argument_parser = ArgumentParser(prog='eltrig',
      description='Driver script for electron trigger SF measurment.')
  argument_parser.add_argument('-y','--year',choices=['2016APV','2016',
      '2017','2018','2022','2022EE','2023','2023BPix'],default='2016APV')
  argument_parser.add_argument('-t','--trig',choices=['singleel','diel23',
      'diel12'],default='singleel')
  args = argument_parser.parse_args()
  
  #default: 2016APV ele27
  data_filenames = ['/data2/oshiro/ntuples/2016APV/Run2016B_L1matched.root',
                    '/data2/oshiro/ntuples/2016APV/Run2016C_L1matched.root',
                    '/data2/oshiro/ntuples/2016APV/Run2016D_L1matched.root',
                    '/data2/oshiro/ntuples/2016APV/Run2016E_L1matched.root',
                    '/data2/oshiro/ntuples/2016APV/Run2016F_L1matched.root']
  mc_filenames = ['/data2/oshiro/ntuples/2016APV/DY_NLO_L1matched.root']
  mcalt_filenames = ['/data2/oshiro/ntuples/2016APV/DY_LO_L1matched.root']
  measurement_cut = 'passHltEle27WPTightGsf'
  measurement_desc = 'HLT_Ele27_WPTight_Gsf'
  preselection = 'tag_Ele_pt>40&&tag_Ele_abseta<2.17&&(tag_Ele_q+el_q)==0&&el_pt>7&&fabs(el_sc_eta)<2.5&&fabs(el_dz)<1.0&&fabs(el_dxy)<0.5&&passingMVA94XwpLooseisoV2'
  analyzer_name = 'eltrig27'
  year = '2016APV'
  pt_binning = [7.0,25.0,26.0,27.0,28.0,29.0,31.0,35.0,50.0,100.0,500.0]
  eta_binning = [-2.5,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.5]
  gappt_binning = [7.0,25.0,35.0,500.0]
  if args.trig=='diel23':
    measurement_cut = 'passHltEle23Ele12CaloIdLTrackIdLIsoVLLeg1L1match'
    measurement_desc = 'Dielectron trigger, 23 GeV leg'
    analyzer_name = 'eltrig23'
    pt_binning = [7.0,21.0,22.0,23.0,24.0,25.0,27.0,35.0,50.0,100.0,500.0]
    eta_binning = [-2.5,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.5]
    gappt_binning = [7.0,23.0,35.0,500.0]
  if args.trig=='diel12':
    measurement_cut = 'passHltEle23Ele12CaloIdLTrackIdLIsoVLLeg2'
    measurement_desc = 'Dielectron trigger, 12 GeV leg'
    analyzer_name = 'eltrig12'
    pt_binning = [7.0,11.0,12.0,13.0,14.0,16.0,20.0,30.0,50.0,100.0,500.0]
    eta_binning = [-2.5,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.5]
    gappt_binning = [7.0,12.0,35.0,500.0]

  analyzer_name = 'hzg_'+analyzer_name+'_'+year
  
  eltrig_analyzer = RmsSFAnalyzer(analyzer_name)
  eltrig_analyzer.year = year
  eltrig_analyzer.set_input_files(data_filenames,mc_filenames,mcalt_filenames,'tnpEleTrig/fitter_tree')
  eltrig_analyzer.set_fitting_variable('pair_mass','m_{ee} [GeV]')
  eltrig_analyzer.set_measurement_variable(measurement_cut,measurement_desc)
  eltrig_analyzer.set_preselection(preselection,preselection,preselection)
  eltrig_analyzer.add_standard_gap_binning(pt_binning,eta_binning,
                                           gappt_binning,'el_pt','el_sc_eta')
  eltrig_analyzer.run_interactive()

