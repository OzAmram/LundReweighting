
from Utils import *



f_sig = h5py.File("../TagNTrain/data/WpToBpT_Wp3000_Bp400_UL17_TIMBER.h5", "r")

f_ratio = ROOT.TFile("ttbar_rw_ak8_test_v2/ratio.root", "READ")
h_rw = f_ratio.Get("h_ratio")

loadMax = 200000
numJets = 10000

pt_min = 200.
pt_max = 99999.

dr_bin_min = -1.
dr_bin_max = 8.
#y_bin_min = np.log(1./0.5)
#y_bin_max = 20*y_bin_min
#y_label = "ln(1/z)"
y_bin_min = -5
y_bin_max = np.log(pt_max)
y_label = "ln(kt)"
n_bins = 20

jetR = 0.4
fill_z = False

ratio_range = [0.5, 1.5]

outDir = "Wp_rw_test/"

subjet_rw = False


if(subjet_rw):
    n_pt_bins = 5
    pt_bins = array('f', [0., 10., 25., 40., 60., 99999.])
else:
    n_pt_bins = 4
    pt_bins = array('f', [200., 300., 400., 600., 99999.])
    #pt_bins = array('f', [200., 300.])



jet_kinematics_sig= f_sig['jet_kinematics'][:loadMax]
all_cands_sig = f_sig['jet1_PFCands'][:loadMax].astype(np.float64)
jet1_feats_sig = f_sig['jet1_extraInfo'][:loadMax]

pt_cut_mask = (jet_kinematics_sig[:,2] > pt_min) & (jet_kinematics_sig[:,2] < pt_max)

print("pt_cut", np.mean(pt_cut_mask))

cut_cands_sig = all_cands_sig[ pt_cut_mask]
cut_feats_sig = jet1_feats_sig[ pt_cut_mask]
cut_kinematics_sig = jet_kinematics_sig[pt_cut_mask]

print(cut_cands_sig.shape)






eps = 1e-8


tau21_sig = (cut_feats_sig[:,1] / (cut_feats_sig[:,0] + eps))[:numJets ]

tau32_sig = (cut_feats_sig[:,2] / (cut_feats_sig[:,1] + eps))[:numJets ]

tau43_sig = (cut_feats_sig[:,3] / (cut_feats_sig[:,2] + eps))[:numJets]

nPF_sig= cut_feats_sig[:,-1][:numJets]

n_bins = 12



LP_weights = []

for i in range(0, numJets):
    pf_cands_sig = cut_cands_sig[i]

    rw = reweight_lund_plane(h_rw, pf_cands_sig, fill_z = fill_z)
    if(subjet_rw):

        jet_4vec = convert_4vec(cut_kinematics_sig[i,2:6])
        boost_vec = fj.PseudoJet(jet_4vec[0], jet_4vec[1], jet_4vec[2], jet_4vec[3])
        rw = reweight_lund_plane(h_rw, pf_cands_sig,  boost_vec = boost_vec, fill_z =fill_z, dR = jetR, jetR = jetR, maxJets = 4, pt_min = 10.)
    LP_weights.append(rw)




LP_weights = np.clip(LP_weights, 0., 10.)
print(LP_weights[:10])

weights = [  LP_weights, [1.]*len(tau21_sig)]
print(LP_weights.shape, tau21_sig.shape)



make_ratio_histogram([tau21_sig, tau21_sig], ["sig RW", "sig"], ['b', 'r'], 'Tau21', "Tau21 : Before vs. After Lund Plane Reweighting", n_bins,
     ratio_range = ratio_range, normalize=True, weights = weights, save = True, fname=outDir + "tau21_before_vs_after.png")

make_ratio_histogram([tau32_sig, tau32_sig], ["sig RW", "sig"], ['b', 'r'], 'tau32', "tau32 : Before vs. After Lund Plane Reweighting", n_bins,
     ratio_range = ratio_range, normalize=True, weights = weights, save = True, fname=outDir + "tau32_before_vs_after.png")

make_ratio_histogram([tau43_sig, tau43_sig], ["sig RW", "sig"], ['b', 'r'], 'tau43', "tau43 : Before vs. After Lund Plane Reweighting", n_bins,
     ratio_range = ratio_range, normalize=True, weights = weights, save = True, fname=outDir + "tau43_before_vs_after.png")

make_ratio_histogram([nPF_sig, nPF_sig], ["sig RW", "sig"], ['b', 'r'], 'nPF', "nPF : Before vs. After Lund Plane Reweighting", n_bins,
     ratio_range = ratio_range, normalize=True, weights = weights, save = True, fname=outDir + "nPF_before_vs_after.png")
