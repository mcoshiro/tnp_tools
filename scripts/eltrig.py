#generate SFs for electrons passing H->Zgamma trigger criteria
#  phase space: pt>7, |etasc|<2.5
#  |dz|<1 cm, |dxy|<0.5 cm, fall17v2 WPL (run 2) or HZZ 2018 (run 3)

from rms_sf_analyzer import RmsSFAnalyzer
from argparse import ArgumentParser

if __name__=='__main__':

  argument_parser = ArgumentParser(prog='eltrig',
      description='Driver script for electron trigger SF measurment.')
  argument_parser.add_argument('-y','--year',choices=['2016APV','2016',
      '2017','2018','2022','2022EE','2023','2023BPix','2023BPixHole'],
      default='2016APV')
  argument_parser.add_argument('-t','--trig',choices=['singleel','diel23',
      'diel12'],default='singleel')
  args = argument_parser.parse_args()
  
  #default: 2016APV ele27
  file_path = '/eos/cms/store/group/phys_egamma/tnpTuples/rasharma/2021-02-10/UL2016preVFP/merged/'
  data_filenames = [file_path+'Run2016B_L1matched.root',
                    file_path+'Run2016C_L1matched.root',
                    file_path+'Run2016D_L1matched.root',
                    file_path+'Run2016E_L1matched.root',
                    file_path+'Run2016F_L1matched.root']
  mc_filenames = [file_path+'DY_NLO_L1matched.root']
  mcalt_filenames = [file_path+'DY_LO_L1matched.root']
  measurement_cut = 'passHltEle27WPTightGsf'
  measurement_desc = 'HLT_Ele27_WPTight_Gsf'
  preselection = ('tag_Ele_pt>40&&tag_Ele_abseta<2.17&&(tag_Ele_q+el_q)==0'
                 +'&&el_pt>7&&fabs(el_sc_eta)<2.5&&fabs(el_dz)<1.0'
                 +'&&fabs(el_dxy)<0.5&&passingMVA94XwpLooseisoV2')
  analyzer_name = 'eltrig27'
  year = args.year
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
    if args.trig=='singleel':
      measurement_cut = 'passHltEle32DoubleEGWPTightGsf&&passEGL1SingleEGOr'
      measurement_desc = 'HLT_Ele32_WPTight (emulated)'
      analyzer_name = 'eltrig32'
      pt_binning = [7.0,31.0,32.0,33.0,34.0,35.0,38.0,45.0,80.0,120.0,500.0]
      eta_binning = [-2.5,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.5]
      gappt_binning = [7.0,32.0,45.0,500.0]
  elif (year == '2018'):
    file_path = '/eos/cms/store/group/phys_egamma/tnpTuples/tomc/2020-05-20/UL2018/merged/'
    data_filenames = [file_path+'Run2018A.root',
                      file_path+'Run2018B.root',
                      file_path+'Run2018C.root',
                      file_path+'Run2018D.root']
    mc_filenames = [file_path+'DY_NLO.root']
    mcalt_filenames = [file_path+'DY_LO.root']
    if args.trig=='singleel':
      measurement_cut = 'passHltEle32WPTightGsf'
      measurement_desc = 'HLT_Ele32_WPTight'
      analyzer_name = 'eltrig32'
      pt_binning = [7.0,31.0,32.0,33.0,34.0,35.0,38.0,45.0,80.0,120.0,500.0]
      eta_binning = [-2.5,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.5]
      gappt_binning = [7.0,32.0,45.0,500.0]
  elif (year == '2022'):
    preselection = ('tag_Ele_pt>40&&tag_Ele_abseta<2.17&&(tag_Ele_q+el_q)==0'
                 +'&&el_pt>7&&fabs(el_sc_eta)<2.5&&fabs(el_dz)<1.0'
                 +'&&fabs(el_dxy)<0.5' + '&&((fabs(el_eta)<0.8 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9266) || (fabs(el_eta)<0.8 && el_pt>10 && el_hzzMVA>0.3527) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9138) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>10 && el_hzzMVA>0.2601) || (fabs(el_eta)>1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9682) || (fabs(el_eta)>1.479 && el_pt>10 && el_hzzMVA>-0.4963))')
    file_path = '/eos/cms/store/group/phys_egamma/tnpTuples/jgrassi/2024-03-12/'
    data_filenames = [file_path+'data_EGamma_2022CD_merged.root']
    mc_filenames = [file_path+'mc_DY_NLO_2022preEE_merged.root']
    mcalt_filenames = [file_path+'mc_DY_LO_2022preEE_merged.root']
    if args.trig=='singleel':
      measurement_cut = 'passHltEle30WPTightGsf'
      measurement_desc = 'HLT_Ele30_WPTight'
      analyzer_name = 'eltrig30'
      pt_binning = [7.0,31.0,32.0,33.0,34.0,35.0,38.0,45.0,80.0,120.0,500.0]
      eta_binning = [-2.5,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.5]
      gappt_binning = [7.0,30.0,45.0,500.0]
  elif (year == '2022EE'):
    preselection = ('tag_Ele_pt>40&&tag_Ele_abseta<2.17&&(tag_Ele_q+el_q)==0'
                 +'&&el_pt>7&&fabs(el_sc_eta)<2.5&&fabs(el_dz)<1.0'
                 +'&&fabs(el_dxy)<0.5' + '&&((fabs(el_eta)<0.8 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9266) || (fabs(el_eta)<0.8 && el_pt>10 && el_hzzMVA>0.3527) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9138) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>10 && el_hzzMVA>0.2601) || (fabs(el_eta)>1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9682) || (fabs(el_eta)>1.479 && el_pt>10 && el_hzzMVA>-0.4963))')
    file_path = '/eos/cms/store/group/phys_egamma/tnpTuples/jgrassi/2024-03-12/'
    data_filenames = [file_path+'data_EGamma_2022EFG_merged.root']
    mc_filenames = [file_path+'mc_DY_NLO_2022postEE_merged.root']
    mcalt_filenames = [file_path+'mc_DY_LO_2022postEE_merged.root']
    if args.trig=='singleel':
      measurement_cut = 'passHltEle30WPTightGsf'
      measurement_desc = 'HLT_Ele30_WPTight'
      analyzer_name = 'eltrig30'
      pt_binning = [7.0,31.0,32.0,33.0,34.0,35.0,38.0,45.0,80.0,120.0,500.0]
      eta_binning = [-2.5,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.5]
      gappt_binning = [7.0,30.0,45.0,500.0]
  elif (year == '2023'):
    preselection = ('tag_Ele_pt>40&&tag_Ele_abseta<2.17&&(tag_Ele_q+el_q)==0'
               +'&&el_pt>7&&fabs(el_sc_eta)<2.5&&fabs(el_dz)<1.0'
               +'&&fabs(el_dxy)<0.5' + '&&((fabs(el_eta)<0.8 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9266) || (fabs(el_eta)<0.8 && el_pt>10 && el_hzzMVA>0.3527) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9138) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>10 && el_hzzMVA>0.2601) || (fabs(el_eta)>1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9682) || (fabs(el_eta)>1.479 && el_pt>10 && el_hzzMVA>-0.4963))')
    file_path = "/eos/cms/store/group/phys_egamma/tnpTuples/jgrassi/2024-03-12/"
    data_filenames = [file_path+'data_EGamma_2023C_merged.root']
    mc_filenames = [file_path+'mc_DY_NLO_2023preBPIX_merged.root']
    mcalt_filenames = [file_path+'mc_DY_LO_2023preBPIX_merged.root']
    if args.trig=='singleel':
      measurement_cut = 'passHltEle30WPTightGsf'
      measurement_desc = 'HLT_Ele30_WPTight'
      analyzer_name = 'eltrig30'
      pt_binning = [7.0,31.0,32.0,33.0,34.0,35.0,38.0,45.0,80.0,120.0,500.0]
      eta_binning = [-2.5,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.5]
      gappt_binning = [7.0,30.0,45.0,500.0]

  elif (year == '2023BPix'):
    preselection = ('tag_Ele_pt>40&&tag_Ele_abseta<2.17&&(tag_Ele_q+el_q)==0'
               +'&&el_pt>7&&fabs(el_sc_eta)<2.5&&fabs(el_dz)<1.0'
               +'&&fabs(el_dxy)<0.5' + '&&((fabs(el_eta)<0.8 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9266) || (fabs(el_eta)<0.8 && el_pt>10 && el_hzzMVA>0.3527) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9138) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>10 && el_hzzMVA>0.2601) || (fabs(el_eta)>1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9682) || (fabs(el_eta)>1.479 && el_pt>10 && el_hzzMVA>-0.4963))')
    preselection += '&&!(el_eta>-1.5&&el_eta<0.0&&el_phi>-1.2&&el_phi<-0.8)'
    file_path = "/eos/cms/store/group/phys_egamma/tnpTuples/jgrassi/2024-03-12/"
    data_filenames = [file_path+'data_EGamma_2023D_merged.root']
    mc_filenames = [file_path+'mc_DY_NLO_2023postBPIX_merged.root']
    mcalt_filenames = [file_path+'mc_DY_LO_2023postBPIX_merged.root']
    if args.trig=='singleel':
      measurement_cut = 'passHltEle30WPTightGsf'
      measurement_desc = 'HLT_Ele30_WPTight'
      analyzer_name = 'eltrig30'
      pt_binning = [7.0,31.0,32.0,33.0,34.0,35.0,38.0,45.0,80.0,120.0,500.0]
      eta_binning = [-2.5,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.5]
      gappt_binning = [7.0,30.0,45.0,500.0]
  elif (year == '2023BPixHole'):
    preselection = ('tag_Ele_pt>40&&tag_Ele_abseta<2.17&&(tag_Ele_q+el_q)==0'
               +'&&el_pt>7&&fabs(el_sc_eta)<2.5&&fabs(el_dz)<1.0'
               +'&&fabs(el_dxy)<0.5' + '&&((fabs(el_eta)<0.8 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9266) || (fabs(el_eta)<0.8 && el_pt>10 && el_hzzMVA>0.3527) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9138) || (fabs(el_eta)>0.8 && el_eta<1.479 && el_pt>10 && el_hzzMVA>0.2601) || (fabs(el_eta)>1.479 && el_pt>5 && el_pt<10 && el_hzzMVA>0.9682) || (fabs(el_eta)>1.479 && el_pt>10 && el_hzzMVA>-0.4963))')
    preselection += '&&(el_eta>-1.5&&el_eta<0.0&&el_phi>-1.2&&el_phi<-0.8)'
    file_path = "/eos/cms/store/group/phys_egamma/tnpTuples/jgrassi/2024-03-12/"
    data_filenames = [file_path+'data_EGamma_2023D_merged.root']
    mc_filenames = [file_path+'mc_DY_NLO_2023postBPIX_merged.root']
    mcalt_filenames = [file_path+'mc_DY_LO_2023postBPIX_merged.root']
    if args.trig=='singleel':
      measurement_cut = 'passHltEle30WPTightGsf'
      measurement_desc = 'HLT_Ele30_WPTight'
      analyzer_name = 'eltrig30'
      pt_binning = [7.0,31.0,32.0,33.0,34.0,35.0,38.0,45.0,80.0,120.0,500.0]
      eta_binning = [-2.5,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.5]
      gappt_binning = [7.0,30.0,45.0,500.0]

  analyzer_name = 'hzg_'+analyzer_name+'_'+year
  preselection_mc = preselection + '&&(mcTrue==1)'
  
  eltrig_analyzer = RmsSFAnalyzer(analyzer_name)
  eltrig_analyzer.year = year
  eltrig_analyzer.set_input_files(data_filenames,mc_filenames,mcalt_filenames,'tnpEleTrig/fitter_tree')
  eltrig_analyzer.set_fitting_variable('pair_mass','m_{ee} [GeV]',
                                       weight_mc='totWeight')
  eltrig_analyzer.set_measurement_variable(measurement_cut,measurement_desc)
  eltrig_analyzer.set_preselection(preselection,preselection_mc,preselection)
  if (year != '2023BPixHole'):
    eltrig_analyzer.add_standard_gap_binning(pt_binning,eta_binning,
                                             gappt_binning,'el_pt','el_sc_eta')
  else:
    eltrig_analyzer.add_standard_binning(pt_binning,[-1.566,-1.4442,-0.8,0.0],
                                         'el_pt','el_sc_eta')
  eltrig_analyzer.run_interactive()


