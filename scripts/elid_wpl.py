#generate SFs for electrons passing H->Zgamma signal criteria:
#  phase space: pt>7, |etasc|<2.5
#  |dz|<1 cm, |dxy|<0.5 cm, fall17v2 WPL (run 2) or HZZ 2018 (run 3)

from rms_sf_analyzer import RmsSFAnalyzer
from argparse import ArgumentParser
from bin_utils import add_standard_gap_lohipt_binning

if __name__=='__main__':

  argument_parser = ArgumentParser(prog='elid',
      description='Driver script for electron ID SF measurment.')
  argument_parser.add_argument('-y','--year',choices=['2016APV','2016',
      '2017','2018','2022','2022EE','2023','2023BPix','2023BPixHole'],
      default='2016APV')
  args = argument_parser.parse_args()

  #default: 2016APV
  file_path = ('/eos/cms/store/group/phys_egamma/tnpTuples/rasharma/'
               '2021-02-10/UL2016preVFP/merged/')
  data_filenames = [file_path+'Run2016B_L1matched.root',
                    file_path+'Run2016C_L1matched.root',
                    file_path+'Run2016D_L1matched.root',
                    file_path+'Run2016E_L1matched.root',
                    file_path+'Run2016F_L1matched.root']
  mc_filenames = [file_path+'DY_NLO_L1matched.root']
  mcalt_filenames = [file_path+'DY_LO_L1matched.root']
  measurement_cut = ('fabs(el_dz)<1.0&&fabs(el_dxy)<0.5'
                    +'&&passingMVA94XwpLooseisoV2')
  measurement_desc = 'H#rightarrow Z#gamma electron ID'
  preselection = 'tag_Ele_pt>40&&tag_Ele_abseta<2.17&&(tag_Ele_q+el_q)==0'
  mc_weight = 'totWeight'
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

  elif (year == '2022'):
    file_path = '/eos/cms/store/group/phys_egamma/tnpTuples/jgrassi/2024-03-12/' #Check that the files exist first. As of 22-11-2024, had not merged some of these files
    data_filenames = [file_path+'data_EGamma_2022CD_merged.root']
    mc_filenames = [file_path+'mc_DY_NLO_2022preEE_merged.root']
    mcalt_filenames = [file_path+'mc_DY_LO_2022preEE_merged.root']
    measurement_cut = ('fabs(el_dz)<1.0&&fabs(el_dxy)<0.5'
                    +'&&((fabs(el_eta)<0.8 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9266) || (fabs(el_eta)<0.8 && el_pt>10 && el_hzzMVA>0.3527) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9138) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>10 && el_hzzMVA>0.2601) || (fabs(el_eta)>1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9682) || (fabs(el_eta)>1.479 && el_pt>10 && el_hzzMVA>-0.4963))')
    measurement_desc = 'H#rightarrow Z#gamma electron ID'
    preselection = 'tag_Ele_pt>30&&tag_Ele_abseta<2.17&&(tag_Ele_q+el_q)==0'
    mc_weight = 'getpusf2022(truePU)*weight'
  elif (year == '2022EE'):
    file_path = '/eos/cms/store/group/phys_egamma/tnpTuples/jgrassi/2024-03-12/'
    data_filenames = [file_path+'data_EGamma_2022EFG_merged.root']
    mc_filenames = [file_path+'mc_DY_NLO_2022postEE_merged.root']
    mcalt_filenames = [file_path+'mc_DY_LO_2022postEE_merged.root']
    measurement_cut = ('fabs(el_dz)<1.0&&fabs(el_dxy)<0.5'
                    +'&&((fabs(el_eta)<0.8 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9266) || (fabs(el_eta)<0.8 && el_pt>10 && el_hzzMVA>0.3527) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9138) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>10 && el_hzzMVA>0.2601) || (fabs(el_eta)>1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9682) || (fabs(el_eta)>1.479 && el_pt>10 && el_hzzMVA>-0.4963))')
    measurement_desc = 'H#rightarrow Z#gamma electron ID'
    preselection = 'tag_Ele_pt>30&&tag_Ele_abseta<2.17&&(tag_Ele_q+el_q)==0'
    mc_weight = 'getpusf2022EE(truePU)*weight'
  elif (year == '2023'):
    file_path = "/eos/cms/store/group/phys_egamma/tnpTuples/jgrassi/2024-03-12/"
    data_filenames = [file_path+'data_EGamma_2023C_merged.root']
    mc_filenames = [file_path+'mc_DY_NLO_2023preBPIX_merged.root']
    mcalt_filenames = [file_path+'mc_DY_LO_2023preBPIX_merged.root']
    measurement_cut = ('fabs(el_dz)<1.0&&fabs(el_dxy)<0.5'
                    +'&&((fabs(el_eta)<0.8 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9266) || (fabs(el_eta)<0.8 && el_pt>10 && el_hzzMVA>0.3527) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9138) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>10 && el_hzzMVA>0.2601) || (fabs(el_eta)>1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9682) || (fabs(el_eta)>1.479 && el_pt>10 && el_hzzMVA>-0.4963))')
    measurement_desc = 'H#rightarrow Z#gamma electron ID'
    preselection = 'tag_Ele_pt>30&&tag_Ele_abseta<2.17&&(tag_Ele_q+el_q)==0'
    mc_weight = 'getpusf2023(truePU)*weight'
  elif (year == '2023BPix'):
    file_path = "/eos/cms/store/group/phys_egamma/tnpTuples/jgrassi/2024-03-12/"
    data_filenames = [file_path+'data_EGamma_2023D_merged.root']
    mc_filenames = [file_path+'mc_DY_NLO_2023postBPIX_merged.root']
    mcalt_filenames = [file_path+'mc_DY_LO_2023postBPIX_merged.root']
    measurement_cut = ('fabs(el_dz)<1.0&&fabs(el_dxy)<0.5'
                    +'&&((fabs(el_eta)<0.8 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9266) || (fabs(el_eta)<0.8 && el_pt>10 && el_hzzMVA>0.3527) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9138) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>10 && el_hzzMVA>0.2601) || (fabs(el_eta)>1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9682) || (fabs(el_eta)>1.479 && el_pt>10 && el_hzzMVA>-0.4963))')
    preselection = 'tag_Ele_pt>30&&tag_Ele_abseta<2.17&&(tag_Ele_q+el_q)==0'
    preselection += '&&!(el_eta>-1.5&&el_eta<0.0&&el_phi>-1.2&&el_phi<-0.8)'
    mc_weight = 'getpusf2023BPix(truePU)*weight'
  elif (year == '2023BPixHole'):
    file_path = "/eos/cms/store/group/phys_egamma/tnpTuples/jgrassi/2024-03-12/"
    data_filenames = [file_path+'data_EGamma_2023D_merged.root']
    mc_filenames = [file_path+'mc_DY_NLO_2023postBPIX_merged.root']
    mcalt_filenames = [file_path+'mc_DY_LO_2023postBPIX_merged.root']
    measurement_cut = ('fabs(el_dz)<1.0&&fabs(el_dxy)<0.5'
                    +'&&((fabs(el_eta)<0.8 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9266) || (fabs(el_eta)<0.8 && el_pt>10 && el_hzzMVA>0.3527) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9138) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>10 && el_hzzMVA>0.2601) || (fabs(el_eta)>1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9682) || (fabs(el_eta)>1.479 && el_pt>10 && el_hzzMVA>-0.4963))')
    preselection = 'tag_Ele_pt>30&&tag_Ele_abseta<2.17&&(tag_Ele_q+el_q)==0'
    preselection += '&&(el_eta>-1.5&&el_eta<0.0&&el_phi>-1.2&&el_phi<-0.8)'
    mc_weight = 'getpusf2023BPix(truePU)*weight'

  preselection_mc = preselection + '&&(mcTrue==1)'
  
  elid_analyzer = RmsSFAnalyzer('hzg_elid_{}'.format(year))
  elid_analyzer.year = year
  elid_analyzer.set_input_files(data_filenames,mc_filenames,mcalt_filenames,
                                'tnpEleIDs/fitter_tree')
  elid_analyzer.set_fitting_variable('pair_mass','m_{ee} [GeV]',
                                     weight_mc=mc_weight)
  elid_analyzer.set_measurement_variable(measurement_cut,measurement_desc)
  elid_analyzer.set_preselection(preselection,preselection_mc,preselection)
  if (year != '2023BPixHole'):
    add_standard_gap_lohipt_binning(elid_analyzer, 
        [7.0,15.0,20.0,35.0,50.0,100.0],
        [-2.5,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.5],[7.0,35.0,500.0],'el_pt',
        'el_sc_eta')
  else:
    elid_analyzer.add_standard_binning([7.0,20.0,35.0,50.0,500.0],
        [-1.566,-1.4442,-0.8,0.0],'el_pt','el_sc_eta')
  #elid_analyzer.print_binning()
  elid_analyzer.run_interactive()

