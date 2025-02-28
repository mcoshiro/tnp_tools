#!/usr/bin/env python
'''Script to produce tuples that can be used for tag&probe from pico
'''

from argparse import ArgumentParser
import ROOT

LUMI = {'2016APV' : 19.51,
        '2016' : 16.80,
        '2017' : 41.48,
        '2018' : 59.83,
        '2022' : 8.17,
        '2022EE' : 27.01,
        '2023' : 17.61,
        '2023BPix' : 9.53}

#ROOT JIT'ed C++ definitions
ROOT.gInterpreter.Declare('''
template <class C>
using RVec = ROOT::VecOps::RVec<C>;
using std::vector;
using ROOT::Math::PtEtaPhiMVector;

float dr(float eta1, float eta2, float phi1, float phi2) {
  float dphi = fmod(fabs(phi2-phi1), 2.0*M_PI);
  if (dphi > M_PI) {
    dphi = 2*M_PI-dphi;
  }
  return sqrt((eta2-eta1)*(eta2-eta1)+dphi*dphi);
}

int get_probe_photon_idx(RVec<float> photon_pt, RVec<bool> photon_isScEtaEB,
                         RVec<bool> photon_isScEtaEE, RVec<float> photon_drmin) {
  int idx = -1;
  float max_pt = -999;
  for (unsigned iph = 0; iph < photon_pt.size(); iph++) {
    float pt = photon_pt[iph];
    if (pt > 15.0 && (photon_isScEtaEB[iph] || photon_isScEtaEE[iph])
        && photon_drmin[iph] > 0.3) {
      if (pt > max_pt) {
        idx = static_cast<int>(iph);
        max_pt = pt;
      }
    }
  }
  return idx;
}

int get_probe_photon_idx_mc(RVec<float> photon_pt, RVec<bool> photon_isScEtaEB,
                            RVec<bool> photon_isScEtaEE, RVec<float> photon_drmin,
                            RVec<int> photon_pflavor) {
  int idx = -1;
  float max_pt = -999;
  for (unsigned iph = 0; iph < photon_pt.size(); iph++) {
    float pt = photon_pt[iph];
    if (pt > 15.0 && (photon_isScEtaEB[iph] || photon_isScEtaEE[iph])
        && photon_drmin[iph] > 0.3 && photon_pflavor[iph]==1) {
      if (pt > max_pt) {
        idx = static_cast<int>(iph);
        max_pt = pt;
      }
    }
  }
  return idx;
}

float get_mllg(RVec<float> ll_pt, RVec<float> ll_eta, RVec<float> ll_phi, 
               RVec<float> ll_m, RVec<float> photon_pt, RVec<float> photon_eta,
               RVec<float> photon_phi, int probe_ph_idx) {
  PtEtaPhiMVector ll(ll_pt[0], ll_eta[0], ll_phi[0], ll_m[0]);
  PtEtaPhiMVector ph(photon_pt[probe_ph_idx], photon_eta[probe_ph_idx], 
                     photon_phi[probe_ph_idx], 0.0);
  return (ll+ph).M();
}

''')

def get_year(input_name):
  '''Gets year from input filename
  '''
  for year in ['2016APV','2016','2017','2018','2022EE','2022','2023BPix',
               '2023']:
    if year in input_name:
      return year
  return 'Unknown'

if __name__=='__main__':
  argument_parser = ArgumentParser(prog='nano_to_trigntuple',
      description='Converts NanoAOD to TnP n-tuple')
  argument_parser.add_argument('-i','--input_filename')
  argument_parser.add_argument('-c','--channel',choices=['el','mu'],
                               default='el')
  argument_parser.add_argument('-o','--output_filename')
  args = argument_parser.parse_args()

  year = get_year(args.input_filename)
  is_data = ('data' in args.input_filename)

  for index in ['0','1']:

    output_filename = args.output_filename[:-5]+index+'.root'

    df = ROOT.RDataFrame('tree',args.input_filename)
    if args.channel=='mu':
      df = df.Filter('nmu==2&&nll>=1&&nel==0')
      df = df.Filter('(trig_single_mu&&mu_pt[0]>28)'
                     +'||(trig_double_mu&&mu_pt[0]>20&&mu_pt[1]>10)')
      df = df.Filter('ll_m[0]>50&&ll_m[0]<130')
      if not is_data:
        df = df.Filter('(mu_pflavor[0]==1&&mu_pflavor[1]==1)')
      df = df.Define('probe_mu_pt','mu_pt['+index+']')
      df = df.Define('probe_mu_abseta','fabs(mu_eta['+index+'])')
      df = df.Define('probe_mu_miniso','mu_miniso['+index+']')
      df = df.Define('pair_mass', 'll_m[0]')

      if is_data:
        df = df.Define('w_lumiyearpu','1')
      else:
        df = df.Define('w_year',str(LUMI[year]))
        df = df.Define('w_lumiyearpu','w_lumi*w_year*w_pu')
      df.Snapshot('tree',output_filename,['probe_mu_pt','probe_mu_abseta',
                                          'probe_mu_miniso','pair_mass',
                                          'w_lumiyearpu'])
    elif args.channel=='el':
      df = df.Filter('nel==2&&nll>=1&&nmu==0')
      df = df.Filter('(trig_single_el&&el_pt[0]>35)'
                     +'||(trig_double_el&&el_pt[0]>25&&el_pt[1]>15)')
      df = df.Filter('ll_m[0]>50&&ll_m[0]<130')
      if not is_data:
        df = df.Filter('(el_pflavor[0]==1&&el_pflavor[1]==1)')
      df = df.Define('probe_el_pt','el_pt['+index+']')
      df = df.Define('probe_el_eta','el_eta['+index+']')
      df = df.Define('probe_el_miniso','el_miniso['+index+']')
      df = df.Define('pair_mass', 'll_m[0]')

      if is_data:
        df = df.Define('w_lumiyearpu','1')
      else:
        df = df.Define('w_year',str(LUMI[year]))
        df = df.Define('w_lumiyearpu','w_lumi*w_year*w_pu')
      df.Snapshot('tree',output_filename,['probe_el_pt','probe_el_eta',
                                          'probe_el_miniso','pair_mass',
                                          'w_lumiyearpu'])

