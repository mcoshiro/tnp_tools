'''Copies plots made by spark_tnp into format tnp_tools can use
'''

from argparse import ArgumentParser
import os
import ROOT

#ROOT.gInterpreter.Declare('''
#
#TH1D* cast_tobject_to_th1d(TObject* var) {
#  return static_cast<TH1D*>(var);
#}
#
#''')

if __name__=='__main__':

  argument_parser = ArgumentParser(prog='muon_plot_converter',
     description='script to convert from muon POG format to tnp_tools format')
  argument_parser.add_argument('-c','--clean',action='store_true')
  argument_parser.add_argument('-p','--parent_directory')
  argument_parser.add_argument('-t','--trigger',choices=['singlemu','dimu17',
      'dimu8'])
  argument_parser.add_argument('-y','--year',choices=['2016APV','2016','2017',
      '2018','2022','2022EE','2023','2023BPix'])
  args = argument_parser.parse_args()

  #set up constants
  pog_year = 'Run2016_UL_HIPM'
  tnp_basename = 'NUM_IsoMu24_or_IsoTkMu24'
  trig_threshold = '24'
  n_pt_bins = 12
  n_eta_bins = 4

  if args.trigger == 'dimu17':
    tnp_basename = 'NUM_Mu17leg'
    trig_threshold = '17'
    n_pt_bins = 11
  if args.trigger == 'dimu8':
    tnp_basename = 'NUM_Mu8leg'
    trig_threshold = '8'
    n_pt_bins = 11

  if args.year == '2016':
    pog_year = 'Run2016_UL'
  elif args.year == '2017':
    pog_year = 'Run2017_UL'
    if args.trigger == 'singlemu':
      trig_threshold = '27'
      tnp_basename = 'NUM_IsoMu27'
      n_pt_bins = 11
  elif args.year=='2018':
    pog_year = 'Run2018_UL'
    if args.trigger == 'singlemu':
      tnp_basename = 'NUM_IsoMu24'

  tnp_basename += '_DEN_HToZGamma_SignalMuons'
  pog_filename = tnp_basename+'_abseta_pt.root'

  pog_data_filename = ('{0}/flat/muon/generalTracks/Z/{1}/{1}/Nominal/{2}'
      .format(args.parent_directory,pog_year,pog_filename))
  pog_nmmc_filename = (('{0}/flat/muon/generalTracks/Z/{1}/DY_madgraph/'
      +'Nominal/{2}').format(args.parent_directory,pog_year,
      pog_filename))
  pog_almc_filename = (('{0}/flat/muon/generalTracks/Z/{1}/DY_MassBinned/'
      +'Nominal/{2}').format(args.parent_directory,pog_year,
      pog_filename))

  out_name = 'hzg_mutrig{}_{}'.format(trig_threshold,args.year)
  tnp_data_filename = 'out/{0}_data_nom/{0}_data_nom.root'.format(out_name)
  tnp_dtas_filename = 'out/{0}_data_altsig/{0}_data_altsig.root'.format(
      out_name)
  tnp_dtab_filename = 'out/{0}_data_altbkg/{0}_data_altbkg.root'.format(
      out_name)
  tnp_dtsb_filename = 'out/{0}_data_altsigbkg/{0}_data_altsigbkg.root'.format(
      out_name)
  tnp_nmmc_filename = 'out/{0}_mc_nom/{0}_mc_nom.root'.format(out_name)
  tnp_almc_filename = 'out/{0}_mc_alt/{0}_mc_alt.root'.format(out_name)

  if args.clean:

    #remove directory structure
    os.system('rm -r out/{}*'.format(out_name))

  else:

    #make directory structure
    os.mkdir('out/{}_data_nom'.format(out_name))
    os.mkdir('out/{}_data_altsig'.format(out_name))
    os.mkdir('out/{}_data_altbkg'.format(out_name))
    os.mkdir('out/{}_data_altsigbkg'.format(out_name))
    os.mkdir('out/{}_mc_nom'.format(out_name))
    os.mkdir('out/{}_mc_alt'.format(out_name))

    #move histograms
    pog_data_file = ROOT.TFile(pog_data_filename,'READ')
    tnp_data_file = ROOT.TFile(tnp_data_filename,'CREATE')
    pog_nmmc_file = ROOT.TFile(pog_nmmc_filename,'READ')
    tnp_nmmc_file = ROOT.TFile(tnp_nmmc_filename,'CREATE')
    pog_almc_file = ROOT.TFile(pog_almc_filename,'READ')
    tnp_almc_file = ROOT.TFile(tnp_almc_filename,'CREATE')
    ibin = 0
    for ipt in range(n_pt_bins):
      for ieta in range(n_eta_bins):
        for pf, PF in [('pass','Pass'),('fail','Fail')]:
          pog_hist_name = '{}_abseta_{}_pt_{}_{}'.format(tnp_basename,ieta+1,
                  ipt+1,PF)
          tnp_hist_name = 'hist_{}_bin{}'.format(pf,ibin)
          data_hist = pog_data_file.Get(pog_hist_name)
          data_hist.SetDirectory(ROOT.nullptr)
          tnp_data_file.WriteObject(data_hist,tnp_hist_name)
          nmmc_hist = pog_nmmc_file.Get(pog_hist_name)
          nmmc_hist.SetDirectory(ROOT.nullptr)
          tnp_nmmc_file.WriteObject(nmmc_hist,tnp_hist_name)
          almc_hist = pog_almc_file.Get(pog_hist_name)
          almc_hist.SetDirectory(ROOT.nullptr)
          tnp_almc_file.WriteObject(almc_hist,tnp_hist_name)
        ibin += 1
    pog_data_file.Close()
    tnp_data_file.Close()
    pog_nmmc_file.Close()
    tnp_nmmc_file.Close()
    pog_almc_file.Close()
    tnp_almc_file.Close()

    #create others by copying files
    os.system('cp {} {}'.format(tnp_data_filename, tnp_dtas_filename))
    os.system('cp {} {}'.format(tnp_data_filename, tnp_dtab_filename))
    os.system('cp {} {}'.format(tnp_data_filename, tnp_dtsb_filename))

