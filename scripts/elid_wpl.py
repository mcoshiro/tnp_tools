#generate SFs for electrons passing H->Zgamma signal criteria:
#  phase space: pt>7, |etasc|<2.5
#  |dz|<1 cm, |dxy|<0.5 cm, fall17v2 WPL (run 2) or HZZ 2018 (run 3)

from rms_sf_analyzer import RmsSFAnalyzer

data_filenames = ['/data2/oshiro/ntuples/2016APV/Run2016B_L1matched.root',
                  '/data2/oshiro/ntuples/2016APV/Run2016C_L1matched.root',
                  '/data2/oshiro/ntuples/2016APV/Run2016D_L1matched.root',
                  '/data2/oshiro/ntuples/2016APV/Run2016E_L1matched.root',
                  '/data2/oshiro/ntuples/2016APV/Run2016F_L1matched.root']
mc_filenames = ['/data2/oshiro/ntuples/2016APV/DY_NLO_L1matched.root']
mcalt_filenames = ['/data2/oshiro/ntuples/2016APV/DY_LO_L1matched.root']
measurement_cut = 'fabs(el_dz)<1.0&&fabs(el_dxy)<0.5&&passingMVA94XwpLooseisoV2'
preselection = 'tag_Ele_pt>40&&tag_Ele_abseta<2.17&&(tag_Ele_q+el_q)==0'

elid_analyzer = RmsSFAnalyzer('hzg_elid_2016')
elid_analyzer.year = '2016APV'
elid_analyzer.set_input_files(data_filenames,mc_filenames,mcalt_filenames,'tnpEleIDs/fitter_tree')
elid_analyzer.set_fitting_variable('pair_mass','m_{ee} [GeV]')
elid_analyzer.set_measurement_variable(measurement_cut)
elid_analyzer.set_preselection(preselection,preselection,preselection)
#elid_analyzer.add_standard_gap_binning([7.0,15.0,20.0,35.0,50.0,100.0,500.0],
#                                       [-2.5,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.5],
#                                       [7.0,35.0,500.0],
#                                       'el_pt','el_sc_eta')
elid_analyzer.add_standard_gap_binning([20.0,35.0,50.0],
                                       [-2.5,-1.5,0.0,1.5,2.5],
                                       [20.0,40.0,50.0],
                                       'el_pt','el_sc_eta')
elid_analyzer.run_interactive()

