#generate SFs for electrons passing H->Zgamma mini iso criteria
#  phase space: pt>7, |etasc|<2.5
#  |dz|<1 cm, |dxy|<0.5 cm, fall17v2 WPL (run 2) or HZZ 2018 (run 3)
#  and OR of single and dielectron triggers

from rms_sf_analyzer import RmsSFAnalyzer
from bin_utils import add_standard_gap_highpt_binning
from argparse import ArgumentParser

if __name__=='__main__':

  argument_parser = ArgumentParser(prog='muiso',
      description='Driver script for electron Iso SF measurment.')
  argument_parser.add_argument('-y','--year',choices=['2016APV','2016',
      '2017','2018','2022','2022EE','2023','2023BPix','2023BPixHole'],
      default='2016APV')
  argument_parser.add_argument('-c','--cut',choices=['0.1','0.15'],
                               default='0.15')
  args = argument_parser.parse_args()
  
  year = args.year
  if year=='2023BPixHole':
    year = '2023BPix'
  data_filenames = ['/data2/oshiro/ntuples/'+year+'/elisoskim_data0.root',
                    '/data2/oshiro/ntuples/'+year+'/elisoskim_data1.root']
  mc_filenames = ['/data2/oshiro/ntuples/'+year+'/elisoskim_dynlo0.root',
                  '/data2/oshiro/ntuples/'+year+'/elisoskim_dynlo1.root']
  mcalt_filenames = ['/data2/oshiro/ntuples/'+year+'/elisoskim_dylo0.root',
                     '/data2/oshiro/ntuples/'+year+'/elisoskim_dylo1.root']
  if (year=='2023' or year=='2023BPix'):
    mcalt_filenames = mc_filenames

  #cut is 0.1 for ttH lep, 0.15 for WH/ZH lep
  measurement_cut = 'probe_el_miniso<'+args.cut
  measurement_desc = 'I_{mini}<'+args.cut
  preselection = '1'
  pt_binning = [7.0,15.0,20.0,35.0,50.0,100.0]
  eta_binning = [-2.5,-2.0,-1.5,-0.8,0.0,0.8,1.5,2.0,2.5]
  gap_pt_binning = [7.0,35.0,500.0]
  if args.year=='2023BPix':
    preselection = ('!(probe_el_eta>-1.5&&probe_el_eta<0.0'
                    +'&&probe_el_phi>-1.2&&probe_el_phi<-0.8)')
  elif args.year=='2023BPixHole':
    preselection = ('(probe_el_eta>-1.5&&probe_el_eta<0.0'
                    +'&&probe_el_phi>-1.2&&probe_el_phi<-0.8)')

  cut_cleanname = args.cut.replace('.','p')
  analyzer_name = 'hzg_eliso'+cut_cleanname+'_'+args.year
  
  eliso_analyzer = RmsSFAnalyzer(analyzer_name)
  eliso_analyzer.year = args.year
  eliso_analyzer.set_input_files(data_filenames,mc_filenames,mcalt_filenames,
                                 'tree')
  eliso_analyzer.set_fitting_variable('pair_mass','m_{#el#el} [GeV]',
                                      weight_mc='w_lumiyearpu')
  eliso_analyzer.set_measurement_variable(measurement_cut,measurement_desc)
  eliso_analyzer.set_preselection(preselection,preselection,preselection)
  if (args.year != '2023BPixHole'):
    add_standard_gap_highpt_binning(eliso_analyzer, pt_binning,
        eta_binning, gap_pt_binning, 'probe_el_pt', 'probe_el_eta')
  else:
    eliso_analyzer.add_standard_binning([7.0,20.0,35.0,50.0,500.0],
                                        [-1.566,-1.4442,-0.8,0.0],
                                        'probe_el_pt','probe_el_eta')
  eliso_analyzer.run_interactive()

