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
  argument_parser.add_argument('-o','--output_filename')
  args = argument_parser.parse_args()

  year = get_year(args.input_filename)
  is_data = ('data' in args.input_filename)

  df = ROOT.RDataFrame('tree',args.input_filename)
  df = df.Filter('nmu==2&nll>=1')
  df = df.Filter('(trig_single_mu&&mu_pt[0]>28)'
                 +'||(trig_double_mu&&mu_pt[0]>20&&mu_pt[1]>10)')
  if not is_data:
    df = df.Define('probe_ph_idx','get_probe_photon_idx_mc(photon_pt,'
                   +'photon_isScEtaEB,photon_isScEtaEE,photon_drmin,'
                   +'photon_pflavor)')
  else:
    df = df.Define('probe_ph_idx','get_probe_photon_idx(photon_pt,'
                   +'photon_isScEtaEB,photon_isScEtaEE,photon_drmin)')
  df = df.Filter('probe_ph_idx != -1')
  df = df.Define('probe_ph_pt','photon_pt[probe_ph_idx]')
  df = df.Define('probe_ph_eta','photon_eta[probe_ph_idx]')
  df = df.Define('probe_ph_abseta','fabs(probe_ph_eta)')
  df = df.Define('probe_ph_idmva','photon_idmva[probe_ph_idx]')
  df = df.Define('probe_ph_id80','photon_id80[probe_ph_idx]')
  df = df.Define('ll_mass','ll_m[0]')
  df = df.Define('pair_mass','get_mllg(ll_pt, ll_eta, ll_phi, ll_m, photon_pt,'
                 +'photon_eta, photon_phi, probe_ph_idx)')
  if is_data:
    df = df.Define('w_lumiyearpu','1')
  else:
    df = df.Define('w_year',str(LUMI[year]))
    df = df.Define('w_lumiyearpu','w_lumi*w_year*w_pu')
  df.Snapshot('tree',args.output_filename,['probe_ph_pt','probe_ph_eta',
                                           'probe_ph_abseta','probe_ph_idmva',
                                           'probe_ph_id80','ll_mass',
                                           'pair_mass','w_lumiyearpu'])

