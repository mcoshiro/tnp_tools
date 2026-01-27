import ROOT
import gc
from root_plot_lib import RplPlot

ROOT.gInterpreter.Declare('''

float get_pair_phi(float el_pt, float el_eta, float el_phi, float tag_Ele_pt,
                   float tag_Ele_eta, float tag_Ele_phi) {
  TLorentzVector tag, probe;
  tag.SetPtEtaPhiM(tag_Ele_pt, tag_Ele_eta, tag_Ele_phi, 0.000511);
  probe.SetPtEtaPhiM(el_pt, el_eta, el_phi, 0.000511);
  return (tag+probe).Phi();
}

''')

if __name__ == '__main__':
  lowpt_barrel_sel = '(tag_Ele_pt>30&&tag_Ele_abseta<2.17&&(tag_Ele_q+el_q)==0&&el_pt>7&&el_pt<15&&el_sc_eta>-0.8&&el_sc_eta<0.8&&!passingMVA94XwpLooseisoV2)'
  #lowpt_barrel_sel = '(tag_Ele_pt>30&&tag_Ele_abseta<2.5&&(tag_Ele_q+el_q)==0&&el_pt>7&&el_pt<15&&el_sc_eta>-0.8&&el_sc_eta<0.8&&!passingMVA94XwpLooseisoV2)'
  data_filename = '/eos/cms/store/group/phys_egamma/tnpTuples/rasharma/2021-02-10/UL2016preVFP/merged/Run2016D_L1matched.root'
  simu_filename = '/eos/cms/store/group/phys_egamma/tnpTuples/rasharma/2021-02-10/UL2016preVFP/merged/DY_NLO_L1matched.root'

  pair_phi_expression = 'get_pair_phi(el_pt, el_eta, el_phi, tag_Ele_pt, tag_Ele_eta, tag_Ele_phi)'
  mt_expression = 'sqrt(2.0*event_met_pfmet*tag_Ele_pt*(1.0-cos(event_met_pfphi-tag_Ele_phi)))'
  ptt_expression = 'sqrt(tag_Ele_pt*tag_Ele_pt+event_met_pfmet*event_met_pfmet+2.0*event_met_pfmet*tag_Ele_pt*cos(event_met_pfphi-tag_Ele_phi))'
  dphi_met_tag_expression = 'fabs(fmod(event_met_pfphi-tag_Ele_phi+3.141593, 2.0*3.141593)-3.141593)'
  dphi_met_probe_expression = 'fabs(fmod(event_met_pfphi-el_phi+3.141593, 2.0*3.141593)-3.141593)'
  dphi_met_pair_expression = 'fabs(fmod(event_met_pfphi-pair_phi+3.141593, 2.0*3.141593)-3.141593)'
  dphi_tag_probe_expression = 'fabs(fmod(tag_Ele_phi-el_phi+3.141593, 2.0*3.141593)-3.141593)'
  dr_tag_probe_expression = 'sqrt((tag_Ele_eta-el_eta)*(tag_Ele_eta-el_eta)+dphi_tag_probe*dphi_tag_probe)'
  df_data = ROOT.RDataFrame('tnpEleIDs/fitter_tree', data_filename)
  df_simu = ROOT.RDataFrame('tnpEleIDs/fitter_tree', simu_filename)
  df_data = df_data.Filter(lowpt_barrel_sel)
  df_simu = df_simu.Filter(lowpt_barrel_sel)
  df_data = df_data.Define('mt',mt_expression)
  df_simu = df_simu.Define('mt',mt_expression)
  df_data = df_data.Define('ptt',ptt_expression)
  df_simu = df_simu.Define('ptt',ptt_expression)
  df_data = df_data.Define('pair_phi',pair_phi_expression)
  df_simu = df_simu.Define('pair_phi',pair_phi_expression)
  df_data = df_data.Define('dphi_met_tag',dphi_met_tag_expression)
  df_simu = df_simu.Define('dphi_met_tag',dphi_met_tag_expression)
  df_data = df_data.Define('dphi_met_probe',dphi_met_probe_expression)
  df_simu = df_simu.Define('dphi_met_probe',dphi_met_probe_expression)
  df_data = df_data.Define('dphi_met_pair',dphi_met_pair_expression)
  df_simu = df_simu.Define('dphi_met_pair',dphi_met_pair_expression)
  df_data = df_data.Define('dphi_tag_probe',dphi_tag_probe_expression)
  df_simu = df_simu.Define('dphi_tag_probe',dphi_tag_probe_expression)
  df_data = df_data.Define('dr_tag_probe',dr_tag_probe_expression)
  df_simu = df_simu.Define('dr_tag_probe',dr_tag_probe_expression)
  df_data = df_data.Filter('event_met_pfmet<60&&tag_Ele_trigMVA>0.995&&ptt>50&&tag_Ele_pt>50')
  df_simu = df_simu.Filter('event_met_pfmet<60&&tag_Ele_trigMVA>0.995&&ptt>50&&tag_Ele_pt>50')
  df_data_nomasscut = df_data
  df_simu_nomasscut = df_simu
  df_data = df_data.Filter('(pair_mass>75&&pair_mass<85)||(pair_mass>95&&pair_mass<105)')
  df_simu = df_simu.Filter('pair_mass>75&&pair_mass<105')
  df_simu_offmass = df_simu.Filter('pair_mass<85||pair_mass>95')
  df_simu_onmass = df_simu.Filter('pair_mass>85&&pair_mass<95')
  #df_data = df_data.Filter('mt<45&&tag_Ele_trigMVA>0.92')
  #df_simu = df_simu.Filter('mt<45&&tag_Ele_trigMVA>0.92')
  plot_params = [('tag_Ele_pt','Tag Electron p_{T} [GeV]', 0.0, 175.0),
                 ('pair_pt', 'Pair p_{T} [GeV]', 0.0, 200.0),
                 ('pair_eta', 'Pair #eta', -5.0, 5.0),
                 ('event_met_pfmet', 'p_{T}^{miss} [GeV]', 0.0, 125.0),
                 ('mt', 'm_{T}(tag) [GeV]', 0.0, 150.0),
                 ('ptt', '|#vec{p}_{T}(tag)+#vec{p}_{T}^{miss}| [GeV]', 
                  0.0, 150.0),
                 ('dphi_met_tag', '#Delta#phi(tag, p_{T}^{miss})', 
                  0.0, 3.1416),
                 ('dphi_met_probe', '#Delta#phi(probe, p_{T}^{miss})', 
                  0.0, 3.1416),
                 ('dphi_met_pair', '#Delta#phi(pair, p_{T}^{miss})', 
                  0.0, 3.1416),
                 ('dr_tag_probe', '#Delta R(tag, probe)', 0.0, 5.0),
                 ('tag_Ele_dxy', 'Tag Electron d_{xy} [cm]', 0.0, 0.025),
                 ('tag_Ele_dz', 'Tag Electron d_{z} [cm]', 0.0, 0.025),
                 ('tag_el_sip', 'Tag SIP', 0.0, 4.0),
                 ('tag_Ele_eta', 'Tag Electron #eta', -2.5, 2.5),
                 ('tag_Ele_trigMVA', 'Tag MVA score', 0.95, 1.0),
                 ('pair_mass', 'm_{ee} [GeV]', 75.0, 105.0)]
  data_hists = []
  simu_hists = []
  for var_name, var_desc, var_lo, var_hi in plot_params:
    if var_name=='pair_mass':
      data_hists.append(df_data_nomasscut.Histo1D((f'data_hist_{var_name}',
          f'data;{var_desc};% Events/bin',40,var_lo,var_hi),var_name))
      simu_hists.append(df_simu_nomasscut.Histo1D((f'simu_hist_{var_name}',
          f'MC;{var_desc};% Events/bin',40,var_lo,var_hi),var_name,
          'totWeight'))
    else:
      data_hists.append(df_data.Histo1D((f'data_hist_{var_name}',
          f'data;{var_desc};% Events/bin',40,var_lo,var_hi),var_name))
      simu_hists.append(df_simu.Histo1D((f'simu_hist_{var_name}',
          f'MC;{var_desc};% Events/bin',40,var_lo,var_hi),var_name,
          'totWeight'))
    #data_hists.append(df_simu_onmass.Histo1D((f'data_hist_{var_name}',
    #    f'Z-peak;{var_desc};% Events/bin',40,var_lo,var_hi),var_name))
    #simu_hists.append(df_simu_offmass.Histo1D((f'simu_hist_{var_name}',
    #    f'Off Z-peak;{var_desc};% Events/bin',40,var_lo,var_hi),var_name,
    #    'totWeight'))
  gc.collect()
  gc.disable()
  plots = []
  for ihist in range(len(data_hists)):
    var_name = plot_params[ihist][0]
    data_hist = data_hists[ihist].GetPtr()
    simu_hist = simu_hists[ihist].GetPtr()
    data_hist.Scale(1.0/data_hist.Integral())
    simu_hist.Scale(1.0/simu_hist.Integral())
    plots.append(RplPlot())
    plots[-1].y_max_lower = 1.75
    plots[-1].y_min_lower = 0.25
    plots[-1].plot_outline(data_hist)
    plots[-1].plot_outline(simu_hist)
    plots[-1].add_ratio(f'data_hist_{var_name}',f'simu_hist_{var_name}')
    plots[-1].y_title = '% Events/bin'
    plots[-1].draw(f'lowpt_datamc_{plot_params[ihist][0]}.pdf')
  gc.enable()
    
