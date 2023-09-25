
from Utils import *

lumi = 41.42
#f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_2016/"
#f_dir = "/uscms_data/d3/oamram/CMSSW_12_4_0/src/CASE/CASEUtils/H5_maker/Lund_output_files16_june30/"
f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_2017/"
#f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_2018/"

#f_data = h5py.File(f_dir + "SingleMu_merge.h5", "r")
f_ttbar = h5py.File(f_dir + "TT.h5", "r")
f_wjets = h5py.File(f_dir + "QCD_WJets.h5", "r")
f_diboson = h5py.File(f_dir + "diboson.h5", "r")
f_tw = h5py.File(f_dir + "TW.h5", "r")
f_singletop = h5py.File(f_dir + "SingleTop_merge.h5", "r")

jms_corr = 1.0

d_tw = Dataset(f_tw, label = "tW", color = ROOT.kMagenta, jms_corr = jms_corr)
d_wjets = Dataset(f_wjets, label = "W+Jets + QCD", color = ROOT.kGray, jms_corr = jms_corr)
d_diboson = Dataset(f_diboson, label = "Diboson", color = ROOT.kCyan, jms_corr = jms_corr)
d_singletop = Dataset(f_singletop, label = "Single Top", color = ROOT.kMagenta+4, jms_corr = jms_corr)

d_ttbar_w_match = Dataset(f_ttbar, label = "ttbar : W-matched", color = ROOT.kRed, jms_corr =jms_corr)
d_ttbar_t_match = Dataset(f_ttbar, label = "ttbar : t-matched ", color = ROOT.kOrange-3, jms_corr = jms_corr)
d_ttbar_nomatch = Dataset(f_ttbar, label = "ttbar : unmatched", color = ROOT.kGreen+3, jms_corr = jms_corr)

bkgs = [d_ttbar_nomatch, d_ttbar_t_match, d_tw, d_diboson, d_wjets, d_singletop]
#bkgs = [d_ttbar_nomatch, d_ttbar_t_match ]

for d in bkgs:
    print(d.f['event_info'][:100,-1])

keys = sys_weights_map.keys()
keys.remove("nom_weight")
keys.remove("bkg_norm_down")
keys.remove("bkg_norm_up")

for i,sys_name in enumerate(keys):
    print(sys_name)
    for d in bkgs:

        all_sys_weights = d.get_masked('sys_weights')
        weights_nom = d.get_weights()
        sys_idx = sys_weights_map[sys_name]
        weights_sys = weights_nom * all_sys_weights[:, sys_idx]

        print(d.label, np.mean(weights_sys/ weights_nom), np.std(weights_sys/weights_nom))
    print("\n")
