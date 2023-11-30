import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.Utils import *


parser = input_options()
options = parser.parse_args()

print(options)

if(not os.path.exists(options.outdir)): os.system("mkdir %s" % options.outdir)
jet_str = 'CA'

f_sig = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/TagNTrain/data/LundRW/WkkToWRadionToWWW_M3000_Mr400_TuneCP5_13TeV-madgraph-pythia8_TIMBER_Lund.h5", "r")
#f_sig = h5py.File("/uscms_data/d3/oamram/CMSSW_12_4_0/src/CASE/CASEUtils/H5_maker/Wkk_M3500_test.h5", "r")
f_ratio = ROOT.TFile.Open(options.fin)
label = "Radion"
n_prongs = (4,2)
sig_mass = 3000.


j_idx = 0


jetR = 1.0
num_excjets = -1

max_evts = None
#max_evts = 10

d = Dataset(f_sig, label = label, color = ROOT.kRed, dtype = 1)
is_lep = f_sig['event_info'][:,4]
mjj = f_sig['jet_kinematics'][:,0]



j1_m = f_sig['jet_kinematics'][:,5]
j2_m = f_sig['jet_kinematics'][:,9]


j1_pt = f_sig['jet_kinematics'][:,2]
j2_pt = f_sig['jet_kinematics'][:,6]

max_pt = np.maximum(j1_pt, j2_pt)
min_pt = np.minimum(j1_pt, j2_pt)

#cut = (is_lep < 0.5)
cut = (is_lep < 0.5) & (mjj > 0.8*sig_mass) & (mjj < 1.2*sig_mass)
#cut = (mjj > 0.8*sig_mass) & (mjj < 1.2*sig_mass)
print(np.mean(cut))
#jet mass window
cut = cut & (j1_m > 70.) & ( (j2_m > 70.) & (j2_m < 100.))
print(np.mean(cut))
cut = cut & (max_pt > 400.) & (min_pt > 200.)
print(np.mean(cut))

d.apply_cut(cut)

WH_score = f_sig['jet1_extraInfo'][:,8][cut]
print(WH_score.shape)
print(WH_score[:10])

score_cut = WH_score > 0.8

#use gen weights ? 
weights_nom = np.ones_like(WH_score)

weights_rw = copy.deepcopy(weights_nom)[:max_evts]

h_ratio = f_ratio.Get("ratio_nom")
f_ratio.cd('pt_extrap')
rdir = ROOT.gDirectory
#rdir = None

nToys = 100

#Noise used to generated smeared ratio's based on stat unc
rand_noise = np.random.normal(size = (nToys, h_ratio.GetNbinsX(), h_ratio.GetNbinsY(), h_ratio.GetNbinsZ()))
pt_rand_noise = np.random.normal(size = (nToys, h_ratio.GetNbinsY(), h_ratio.GetNbinsZ(), 3))

LP_rw = LundReweighter(jetR = jetR, pt_extrap_dir = rdir, charge_only = options.charge_only)


subjets, splittings, bad_match = d.get_matched_splittings(LP_rw, num_excjets = num_excjets)
print(np.array(subjets[0]).shape, np.array(splittings[0]).shape)
d_LP_weights, d_LP_smeared_weights, d_pt_smeared_weights = d.reweight_LP(LP_rw, h_ratio, num_excjets = num_excjets, 
        max_evts = max_evts,  prefix = "", rand_noise = rand_noise, pt_rand_noise = pt_rand_noise, subjets = subjets, splittings = splittings)

LP_weights = d_LP_weights


#apply weights, keep normalization fixed
old_norm = np.sum(weights_rw)
weights_rw *= d_LP_weights

new_norm = np.sum(weights_rw)

weights_rw *= old_norm / new_norm
LP_smeared_weights = np.array(d_LP_smeared_weights * np.expand_dims(weights_nom, -1) * (old_norm / new_norm))
pt_smeared_weights = np.array(d_pt_smeared_weights * np.expand_dims(weights_nom, -1) * (old_norm / new_norm))


sys_variations = dict()
if(not options.no_sys):
    #sys_list = list(sys_weights_map.keys())
    sys_list = ["sys_tot_up", "sys_tot_down"]
    for sys in sys_list:
        if(sys == 'nom_weight'): continue
        sys_ratio = f_ratio.Get("ratio_" + sys)
        sys_ratio.Print()
        sys_str = sys + "_"

        sys_LP_weights = d.reweight_LP(LP_rw, sys_ratio, num_excjets = num_excjets, prefix = "", 
                max_evts = max_evts, sys_str = sys_str, subjets = subjets, splittings = splittings)
        sys_weights = weights_nom * sys_LP_weights
        rw = np.sum(weights_nom) / np.sum(sys_weights)
        sys_weights *= rw
        sys_variations[sys] = sys_weights


    #vary weights up/down for b-quark subjets by ratio of b-quark to light quark LP
    b_light_ratio = f_ratio.Get("h_bl_ratio")
    bquark_rw = d.reweight_LP(LP_rw, b_light_ratio, num_excjets = num_excjets,  prefix = "", 
            max_evts = max_evts, sys_str = 'bquark', subjets = subjets, splittings = splittings)

    up_bquark_weights = bquark_rw * weights_rw
    down_bquark_weights = (1./ bquark_rw) * weights_rw

    up_bquark_weights *= old_norm / np.sum(up_bquark_weights)
    down_bquark_weights *= old_norm / np.sum(down_bquark_weights)

    sys_variations['bquark_up'] = up_bquark_weights
    sys_variations['bquark_down'] = down_bquark_weights


#compute 'Scalefactor'

eff_nom = np.average(score_cut, weights = weights_nom)
eff_rw = np.average(score_cut, weights = weights_rw)

print("Nom %.3f, RW %.3f" % (eff_nom, eff_rw))


eff_toys = []
pt_eff_toys = []
for i in range(nToys):
    eff = np.average(score_cut, weights = LP_smeared_weights[:,i])
    eff_toys.append(eff)

    eff1 = np.average(score_cut, weights = pt_smeared_weights[:,i])
    pt_eff_toys.append(eff1)

toys_mean = np.mean(eff_toys)
toys_std = np.std(eff_toys)

print("Toys avg %.3f, std dev %.3f" % (toys_mean, toys_std))


pt_toys_mean = np.mean(pt_eff_toys)
pt_toys_std = np.std(pt_eff_toys)

print("Pt variation toys avg %.3f, std dev %.3f" % (pt_toys_mean, pt_toys_std))

#Add systematic differences in quadrature
sys_unc_up = sys_unc_down = 0.
if(not options.no_sys):

    eff_sys_tot_up = np.average(score_cut, weights = sys_variations['sys_tot_up'])
    eff_sys_tot_down = np.average(score_cut, weights = sys_variations['sys_tot_down'])
    SF_sys_unc_up = abs(eff_sys_tot_up - eff_rw)/eff_nom
    SF_sys_unc_down = abs(eff_sys_tot_down - eff_rw)/eff_nom
    SF_sys_unc = (SF_sys_unc_up + SF_sys_unc_down) / 2.0

    eff_bquark_up = np.average(score_cut, weights = sys_variations['bquark_up'])
    eff_bquark_down = np.average(score_cut, weights = sys_variations['bquark_down'])
    SF_bquark_up = abs(eff_bquark_up - eff_rw)/eff_nom
    SF_bquark_down = abs(eff_bquark_down - eff_rw)/eff_nom
    SF_bquark_unc = (SF_bquark_up + SF_bquark_down) /2.0

else: 
    SF_sys_unc = SF_sys_unc_up = SF_sys_unc_down = bquark_unc = 0. 
    


SF = eff_rw / eff_nom
SF_stat_unc = abs(toys_mean - eff_rw)/eff_nom + toys_std /eff_nom
SF_pt_unc = abs(pt_toys_mean - eff_rw)/eff_nom + pt_toys_std /eff_nom

bad_matching_unc = np.mean(bad_match) * SF


print("\n\nSF (WH_score >  0.8 ) is %.2f +/- %.2f  (stat) +/- %.2f (sys) +/- %.2f (pt) +/- %.2f (bquark) +/- %.2f (matching) \n\n"  
        % (SF, SF_stat_unc, SF_sys_unc, SF_pt_unc, SF_bquark_unc, bad_matching_unc))
f_ratio.Close()

#approximate uncertainty on the reweighting for the plots
overall_unc = (SF_stat_unc **2 + SF_sys_unc**2 + SF_pt_unc**2 + SF_bquark_unc**2 + bad_matching_unc**2) **0.5 / SF
print("overall unc %.3f" % overall_unc)


