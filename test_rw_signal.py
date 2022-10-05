
from Utils import *
import os



f_sig = h5py.File("../TagNTrain/data/WpToBpT_Wp3000_Bp400_Top170_Zbt_TuneCP5_13TeV-madgraphMLM-pythia8_TIMBER_fix.h5", "r")
num_excjets = 3
jetR = 1.0 


subjet_rw = False
excjet_rw = True

d_sig = Dataset(f_sig, label = "WToBpT")

f_ratio = ROOT.TFile("ttbar_UL_oct4_W_rw/ratio.root", "READ")
h_rw = f_ratio.Get("h_ratio")

loadMax = 200000
numJets = 10000

#pt_min = 200.
pt_min = 0.
pt_max = 99999.

fill_z = False

ratio_range = [0.5, 1.5]

outdir = "Wp_rw_test_oct4/"

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





#uncs given as a fractional uncertainty
LP_weights, LP_uncs = d_sig.reweight_LP(h_rw, subjet_rw = subjet_rw, fill_z = fill_z, jetR = jetR, num_excjets = num_excjets, uncs = True, max_evts = numJets)




#for i in range(0, numJets):
#    pf_cands_sig = cut_cands_sig[i]
#
#    if(subjet_rw):
#
#        jet_4vec = convert_4vec(cut_kinematics_sig[i,2:6])
#        boost_vec = fj.PseudoJet(jet_4vec[0], jet_4vec[1], jet_4vec[2], jet_4vec[3])
#        rw = reweight_lund_plane(h_rw, pf_cands_sig,  boost_vec = boost_vec, fill_z =fill_z, dR = jetR, jetR = jetR, maxJets = 4, pt_min = 10.)
#    else:
#
#    LP_weights.append(rw)




LP_weights = np.clip(LP_weights, 0., 10.)

LP_weights_up = LP_weights + LP_uncs*LP_weights
LP_weights_down = LP_weights - LP_uncs*LP_weights

LP_weights_up = np.clip(LP_weights_up, 0., 10.)
LP_weights_down = np.clip(LP_weights_down, 0., 10.)

print(LP_weights[:10])
print(LP_weights_up[:10])
print(LP_weights_down[:10])



weights = [  [1.]*len(tau21_sig), LP_weights, LP_weights_up, LP_weights_down]
labels = ["Nom", "RW", "RW Up", "RW Down"]
colors = ['green','black', 'blue', 'red']


make_histogram(LP_weights, "Reweighting factors", 'b', 'Weight', "Lund Plane Reweighting Factors", 20 , h_range = (0., 2.0),
     normalize=False, fname=outdir + "lundPlane_weights.png")

make_histogram(LP_uncs, "Fractional Uncertainties", 'b', 'Weight Fractional Uncertainty ', "Lund Plane Reweighting Factors Uncertainty", 20,
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
