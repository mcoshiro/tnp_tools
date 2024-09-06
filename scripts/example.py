#!/usr/bin/env python3
"""@package docstring
Script demonstrating example usage of TnpAnalyzer
"""

import ROOT

from tnp_analyzer import *
from model_initializers import *

#enable multithreading
ROOT.EnableImplicitMT()

#initialize
my_tnp_analysis = TnpAnalyzer('my_tnp_analysis')

#set input filename(s)
my_tnp_analysis.set_input_files(['/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_lassen_v0/2018/data/merged_zgdata_llg/merged_raw_pico_llg_DoubleMuon_zgdata_llg_nfiles_206.root'],
                                'tree')

#set fitting variable
my_tnp_analysis.set_fitting_variable('llphoton_m[0]','m_{ll#gamma} [GeV]',50,(60,100))

#set measurement variable
my_tnp_analysis.set_measurement_variable('photon_idmva[0]>0.6')

#set preselection
my_tnp_analysis.set_preselection('nmu==2&&nphoton==1&&trig_double_mu','Preselection')

#add binning
my_tnp_analysis.add_nd_binning(
    [TnpAnalyzer.make_nd_bin_dimension('photon_pt[0]','p_{T}^{#gamma} [GeV]',[10.0,20.0,35.0,50.0]),
    TnpAnalyzer.make_nd_bin_dimension('fabs(photon_eta[0])','|#eta_{#gamma}|',[0.0,1.5,2.5])]
    )

#add fit model (shapes)
my_tnp_analysis.add_model('dscb_p_cms',model_initializer_dscb_p_cms)

#run interactive analysis
my_tnp_analysis.run_interactive()
