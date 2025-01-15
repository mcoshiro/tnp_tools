#generate SFs for electrons passing H->Zgamma signal criteria:
#  phase space: pt>7, |etasc|<2.5
#  |dz|<1 cm, |dxy|<0.5 cm, fall17v2 WPL (run 2) or HZZ 2018 (run 3)

from rms_sf_analyzer import RmsSFAnalyzer
from argparse import ArgumentParser

if __name__=='__main__':

  argument_parser = ArgumentParser(prog='elid',
      description='Driver script for electron ID SF measurment.')
  argument_parser.add_argument('-y','--year',choices=['2016APV','2016',
      '2017','2018','2022','2022EE','2023','2023BPix','2023BPixHole'],
      default='2016APV')
  args = argument_parser.parse_args()

  #default: 2016APV
  data_filenames = ['/data2/oshiro/ntuples/2016APV/Run2016B_L1matched.root',
                    '/data2/oshiro/ntuples/2016APV/Run2016C_L1matched.root',
                    '/data2/oshiro/ntuples/2016APV/Run2016D_L1matched.root',
                    '/data2/oshiro/ntuples/2016APV/Run2016E_L1matched.root',
                    '/data2/oshiro/ntuples/2016APV/Run2016F_L1matched.root']
  mc_filenames = ['/data2/oshiro/ntuples/2016APV/DY_NLO_L1matched.root']
  mcalt_filenames = ['/data2/oshiro/ntuples/2016APV/DY_LO_L1matched.root']
  measurement_cut = ('passingMVA94XV2wp80')
  measurement_desc = 'Photon ID WP80'
  preselection = 'tag_Ele_pt>40&&tag_Ele_abseta<2.17&&((ph_sc_abseta<1.4442)||(ph_sc_abseta>.1566&&ph_sc_abseta<2.5))'
  year = args.year
  if (year == '2016'):
    file_path = '/eos/cms/store/group/phys_egamma/tnpTuples/rasharma/2021-02-10/UL2016postVFP/merged/'
    data_filenames = [file_path+'Run2016F_L1merged.root',
                      file_path+'Run2016G_L1matched.root',
                      file_path+'Run2016H_L1matched.root']
    mc_filenames = [file_path+'DY_NLO_L1matched.root']
    mcalt_filenames = [file_path+'DY_LO_L1matched.root']
  elif (year == '2017'):
    file_path = '/eos/cms/store/group/phys_egamma/tnpTuples/tomc/2020-05-20/UL2017/merged/'
    data_filenames = [file_path+'Run2017B.root',
                      file_path+'Run2017C.root',
                      file_path+'Run2017D.root',
                      file_path+'Run2017E.root',
                      file_path+'Run2017F.root']
    mc_filenames = [file_path+'DY_NLO.root']
    mcalt_filenames = [file_path+'DY_LO.root']
  elif (year == '2018'):
    file_path = '/eos/cms/store/group/phys_egamma/tnpTuples/tomc/2020-05-20/UL2018/merged/'
    data_filenames = [file_path+'Run2018A.root',
                      file_path+'Run2018B.root',
                      file_path+'Run2018C.root',
                      file_path+'Run2018D.root']
    mc_filenames = [file_path+'DY_NLO.root']
    mcalt_filenames = [file_path+'DY_LO.root']
  elif (year == '2023BPix'):
    preselection += '&&!(el_eta>-1.5&&el_eta<0.0&&el_phi>-1.2&&el_phi<-0.8)'
  elif (year == '2023BPixHole'):
    preselection += '&&(el_eta>-1.5&&el_eta<0.0&&el_phi>-1.2&&el_phi<-0.8)'
  
  elid_analyzer = RmsSFAnalyzer('hzg_phidel_{}'.format(year))
  elid_analyzer.year = year
  elid_analyzer.set_input_files(data_filenames,mc_filenames,mcalt_filenames,
                                'tnpPhoIDs/fitter_tree')
  elid_analyzer.set_fitting_variable('pair_mass','m_{e#gamma} [GeV]')
  elid_analyzer.set_measurement_variable(measurement_cut,measurement_desc)
  elid_analyzer.set_preselection(preselection,preselection,preselection)
  elid_analyzer.add_standard_binning([15.0,20.0,35.0,50.0,80.0],
      [-2.5,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.5],'ph_et',
      'ph_sc_eta')
  elid_analyzer.run_interactive()

