
from Utils import *
import os



f_sig = h5py.File("../TagNTrain/data/WpToBpT_Wp3000_Bp400_Top170_Zbt_TuneCP5_13TeV-madgraphMLM-pythia8_TIMBER_fix.h5", "r")
num_excjets = 3
jetR = 1.0 


subjet_rw = False
excjet_rw = True

d_sig = Dataset(f_sig, label = "WToBpT")

f_ratio = ROOT.TFile("ttbar_UL_oct14_W_rw_sys_test/ratio.root", "READ")
h_ratio = f_ratio.Get("ratio_nom")
h_ratio.Print()
do_sys_variations = True

loadMax = 200000
numJets = 10000

#pt_min = 200.
pt_min = 0.
pt_max = 99999.

fill_z = False

ratio_range = [0.5, 1.5]

outdir = "Wp_rw_test_oct14/"

tau32_cut = 0.52

if(not os.path.exists(outdir)): os.system("mkdir " + outdir)



jet_kinematics_sig= f_sig['jet_kinematics'][:loadMax]

pt_cut_mask = (jet_kinematics_sig[:,2] > pt_min) & (jet_kinematics_sig[:,2] < pt_max)

print("pt_cut", np.mean(pt_cut_mask))

d_sig.apply_cut(pt_cut_mask)

cut_cands_sig = d_sig.get_masked("jet1_PFCands")[:numJets]
cut_feats_sig = d_sig.get_masked("jet1_extraInfo")[:numJets]
cut_kinematics_sig = d_sig.get_masked("jet_kinematics")[:numJets]

print(cut_cands_sig.shape)






eps = 1e-8


tau21_sig = (cut_feats_sig[:,1] / (cut_feats_sig[:,0] + eps))

tau32_sig = (cut_feats_sig[:,2] / (cut_feats_sig[:,1] + eps))

tau43_sig = (cut_feats_sig[:,3] / (cut_feats_sig[:,2] + eps))

nPF_sig= cut_feats_sig[:,-1]

msd_sig = cut_kinematics_sig[:,5]

n_bins = 12



weights_nom = np.array([1.] * len(tau21_sig))

nToys = 300
#uncs given as a fractional uncertainty
rand_noise = np.random.normal(size = (h_ratio.GetNbinsX(), h_ratio.GetNbinsY(), h_ratio.GetNbinsZ(), nToys))
LP_weights, LP_uncs, LP_smeared_weights = d_sig.reweight_LP(h_ratio, subjet_rw = subjet_rw, fill_z = fill_z, jetR = jetR, num_excjets = num_excjets, 
                     uncs = True, max_evts = numJets, rand_noise = rand_noise)

sys_ratios = []
sys_variations = dict()
if(do_sys_variations):
    for sys in sys_weights_map.keys():
        if(sys == 'nom_weight'): continue
        print(sys)
        sys_ratio = f_ratio.Get("ratio_" + sys)
        sys_ratio.Print()

        #limit to 1 signal for now
        sys_LP_weights, _ = d_sig.reweight_LP(sys_ratio, subjet_rw = subjet_rw, fill_z = fill_z, jetR = jetR, num_excjets = num_excjets, uncs = False, prefix = "2prong", max_evts = numJets)
        sys_weights = weights_nom * sys_LP_weights
        rw = np.sum(weights_nom) / np.sum(sys_weights)
        sys_weights *= rw
        sys_variations[sys] = sys_weights


LP_weights = np.clip(LP_weights, 0., 10.)
LP_weights_stat_unc = LP_uncs

LP_weights_overall_unc = LP_weights_stat_unc**2
for _,sys_weights in sys_variations.items():
    #Up and down variations separate so add half of weight diff in quadrature
    diff_weights = weights_nom - sys_weights
    LP_weights_overall_unc += (diff_weights/2.)**2

LP_weights_overall_unc = np.sqrt(LP_weights_overall_unc)
print(LP_weights_stat_unc[:5])
print(LP_weights_overall_unc[:5])

LP_weights_up = np.clip(LP_weights + LP_weights_overall_unc, 0., 10.)
LP_weights_down = np.clip(LP_weights - LP_weights_overall_unc, 0., 10.)
print(LP_weights[:10])
print(LP_weights_up[:10])
print(LP_weights_down[:10])


weights = [  weights_nom, LP_weights, LP_weights_up, LP_weights_down]
labels = ["Nom", "RW", "RW Up", "RW Down"]
colors = ['gray','black', 'blue', 'red']


make_histogram(LP_weights, "Reweighting factors", 'b', 'Weight', "Lund Plane Reweighting Factors", 20 , h_range = (0., 2.0),
     normalize=False, fname=outdir + "lundPlane_weights.png")

make_histogram(LP_uncs/LP_weights, "Fractional Uncertainties", 'b', 'Weight Statistical Fractional Uncertainty ', "Lund Plane Reweighting Factors Uncertainty", 20,
     normalize=False, fname=outdir + "lundPlane_weights_unc.png", h_range = (0., 1.5))



make_multi_ratio_histogram([tau21_sig]*4, labels, colors, 'Tau21', "Tau21 : Before vs. After Lund Plane Reweighting", n_bins,
     ratio_range = ratio_range, normalize=True, weights = weights, fname=outdir + "tau21_before_vs_after.png")

make_multi_ratio_histogram([tau32_sig]*4, labels, colors, 'tau32', "tau32 : Before vs. After Lund Plane Reweighting", n_bins,
     ratio_range = ratio_range, normalize=True, weights = weights, fname=outdir + "tau32_before_vs_after.png")

make_multi_ratio_histogram([tau43_sig]*4, labels, colors, 'tau43', "tau43 : Before vs. After Lund Plane Reweighting", n_bins,
     ratio_range = ratio_range, normalize=True, weights = weights, fname=outdir + "tau43_before_vs_after.png")

make_multi_ratio_histogram([nPF_sig ]*4, labels, colors, 'nPF', "nPF : Before vs. After Lund Plane Reweighting", n_bins,
     ratio_range = ratio_range, normalize=True, weights = weights, fname=outdir + "nPF_before_vs_after.png")

make_multi_ratio_histogram([msd_sig ]*4, labels, colors, 'mSoftDrop', "mSoftDrop : Before vs. After Lund Plane Reweighting", n_bins,
     ratio_range = ratio_range, normalize=True, weights = weights, fname=outdir + "msd_before_vs_after.png")


cut = tau32_sig < tau32_cut

eff_nom = np.average(cut)
eff_rw = np.average(cut, weights = LP_weights)

print("Nom %.3f, RW %.3f" % (eff_nom, eff_rw))


eff_toys = []
for i in range(nToys):
    eff = np.average(cut, weights = LP_smeared_weights[:,i])
    eff_toys.append(eff)

toys_mean = np.mean(eff_toys)
toys_std = np.std(eff_toys)

print("Toys avg %.3f, std dev %.3f" % (toys_mean, toys_std))

#Add systematic differences in quadrature
sys_unc_up = sys_unc_down = 0.
for sys in sys_variations.keys():
    eff = np.average(cut, weights = sys_variations[sys])
    diff = eff - eff_rw
    if(diff > 0): sys_unc_up += diff**2
    else: sys_unc_down += diff**2
    print("%s %.4f" % (sys,  diff))

sys_unc_up = sys_unc_up**(0.5)
sys_unc_down = sys_unc_down**(0.5)


SF = eff_rw / eff_nom
SF_stat_unc = toys_std / eff_nom
SF_sys_unc_up = sys_unc_up / eff_nom
SF_sys_unc_down = sys_unc_down / eff_nom

print("SF (cut val %.3f ) is %.2f +/- %.2f  (stat) + %.2f - %.2f (sys)"  % (tau32_cut, SF, SF_stat_unc, SF_sys_unc_up, SF_sys_unc_down))
f_ratio.Close()
