from Utils import *
import os

parser = input_options()
options = parser.parse_args()

print(options)


def get_dRs(gen_eta_phi, j_4vec):
    dR = np.sqrt(np.square(gen_eta_phi[:,0] - j_4vec[1]) + 
            np.square(ang_dist(gen_eta_phi[:,1], j_4vec[2] )))
    return dR


if(not os.path.exists(options.outdir)): os.system("mkdir %s" % options.outdir)
jet_str = 'CA'

#f_sig = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/CASE_signals/Wkk_M2500_R200_nonUL_test.h5", "r")
#f_sig = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/CASE_signals/Wkk_M3000_R170_2018_UL.h5", "r")
#label = "Radion"
#n_prongs = (4,2)
#sig_mass = 3000.

f_sig = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/CASE_signals/ZpToTp_Zp5000_Tp400_test.h5", "r")
label = "ZpToTpTp"
n_prongs = (5,5)
sig_mass = 5000.

j_idx = 0

sys = ""

jetR = 1.0

max_evts = 100000


d = Dataset(f_sig, label = label, color = ROOT.kRed)
is_lep = f_sig['event_info'][:,4]
mjj = f_sig['jet_kinematics'][:,0]

#cut = (is_lep < 0.5) & (mjj > 0.8*sig_mass) & (mjj < 1.2*sig_mass)
#cut = (mjj > 0.8*sig_mass) & (mjj < 1.2*sig_mass)
#cut = (mjj > 0.8*sig_mass) & (mjj < 1.2*sig_mass)
cut = (is_lep < 0.5)

j1_m = f_sig['jet_kinematics'][:,5]
j2_m = f_sig['jet_kinematics'][:,9]

if(label == 'Radion'):
    cut = cut & (j1_m > 70.) & ( (j2_m > 70.) & (j2_m < 100.))


#cut = (mjj > 0.8*sig_mass) & (mjj < 1.2*sig_mass)
d.apply_cut(cut)

gen_parts = d.get_masked('gen_info')[:max_evts]
pf_cands1 = d.get_masked('jet1_PFCands')[:max_evts].astype(np.float64)
pf_cands2 = d.get_masked('jet2_PFCands')[:max_evts].astype(np.float64)
j1_4vec = d.get_masked('jet_kinematics')[:max_evts][:,2:6].astype(np.float64)
j2_4vec = d.get_masked('jet_kinematics')[:max_evts][:,6:10].astype(np.float64)

n = pf_cands1.shape[0]
n = 10
j_subjets = np.zeros((n, n_prongs[j_idx], 4), dtype = np.float32)

gen_parts_eta_phi = gen_parts[:n,:,1:3]
print(gen_parts_eta_phi.shape)

dists = [[0] * n]
j_closest = [[0] * n]
j_which = [[0] * n]

if(j_idx == 0): 
    cands = pf_cands1
    j_4vec = j1_vec

else : 
    cands = pf_cands2
    j_4vec = j2_vec

nBad = 0
is_bad = []
for i in range(n): 
    dRs = get_dRs(gen_parts_eta_phi[i], j_4vec[i])
    n_prongs_i = np.sum(dRs < 0.8)
    just_outside = (np.sum( (dRs < 1.0) & (dRs > 0.8))) > 0
    n_prongs_i = n_prongs[0]
    if(n_prongs_i > 2):
        j_subjets[i], splittings1 = get_splittings(pf_cands1[i], jetR = jetR, num_excjets = n_prongs_i)
    else:
        j_subjets[i] = [[-999,-999,0,-999]]


    dists[i] = get_subjet_dist(gen_parts_eta_phi[i], j_subjets[i][:,1:3])
    j_closest[i] = np.amin(dists[i], axis = -1)
    j_which[i] = np.argmin(dists[i], axis = -1)
    rs_bad.append(n_prongs_i < 2 or just_outside)



#want closest gen-quark for each subjet 
#(NOT other way around because only doing subjets for one jet, so not all quarks will be matched)


print(j_closest)

print("Avg DeltaR : " + str(np.mean(j_closest, axis = 0)))

j_farthest_per_evt = np.amax(j_closest , axis = -1)


j_correct_sum = np.sum(np.arange(n_prongs[j_idx]))
j_correct_idxs = np.arange(n_prongs[j_idx])
j_idx_sums = np.sum(j_which, axis = -1)
j_unique = np.unique(j_which, axis = 1)

deltaR_cut = 0.2

#find number of subjets matched to same gen particle, discounting ones already failing DeltaR cut
n_repeated = [np.sum(j_closest[i] < deltaR_cut)  - np.unique(j_which[i, j_closest[i] < deltaR_cut]).shape[0] for i in range(j_which.shape[0])]


j_samejet = (j_idx_sums != j_correct_sum)


samejet_frac = np.mean(n_repeated) / j_which.shape[1]

print("\nSubjet matching :")
print("Overall avg %.3f, std %.3f. Frac of subjets over %.1f %.3f, frac of events %.3f" % 
        (np.mean(j_closest), np.std(j_closest), deltaR_cut, np.mean(j_closest), np.mean(j_farthest_per_evt > deltaR_cut)))
print("Same-jet frac %.3f" % (samejet_frac) )
print("Overall bad matching frac %.3f" % ( np.mean(j_closest) + samejet_frac)) 


x = np.clip(j_closest.reshape(-1), 1e-6, 1.0)

make_histogram(x, "", colors = 'blue', xaxis_label = r'$\Delta R$ to closest sub-jet', 
                title = "%s : quark-subjet matching (%i prongs)" % (label, n_prongs[j_idx]), num_bins = 40, normalize = True, fname = options.outdir + label + "_j%i_quark_subjet_matching.png" % j_idx)

j_subjet_pts = j_subjets[:,:,0].reshape(-1)
make_histogram([j_subjet_pts], ["Subjets"], colors = ['blue'], xaxis_label = 'Subjet pt (GeV)', 
                title = "%s : subjet pt (%i prongs)" % (label, n_prongs[j_idx]), num_bins = 40, normalize = True, fname = options.outdir + label + "_j%i_subjet_pt.png" % j_idx)

print("Fraction of subjets with pt > 350 : %.3f" % (np.mean(j_subjet_pts > 350.)))
print("Pt scaling unc would be +/- %.3f" %( np.mean(j_subjet_pts > 350.) * 0.21 * (n_prongs[j_idx]) ))

decay_pdgids = np.abs(gen_parts[:,:,3]).reshape(-1)
make_histogram(decay_pdgids, "", colors = 'blue', xaxis_label = 'Decay Product PdgID', 
                title = "%s : decay products (%i+%i prongs)" % (label, n_prongs[0], n_prongs[1]), num_bins = 16, h_range = (0.5, 16.5), normalize = True, fname = options.outdir + label + "_decay_pdgid.png")

