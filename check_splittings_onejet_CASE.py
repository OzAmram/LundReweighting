from Utils import *
import os

def ang_dist(phi1, phi2):
    phi1 = phi1 % (2. * np.pi)
    phi2 = phi2 % (2. * np.pi)
    dphi = phi1 - phi2
    dphi[dphi < -np.pi] += 2.*np.pi
    dphi[dphi > np.pi] -= 2.*np.pi
    return dphi

def get_dists(q_eta_phis, subjets_eta_phis):
    q_eta_phis = np.expand_dims(q_eta_phis, 2)
    subjets_eta_phis = np.expand_dims(subjets_eta_phis, 1)
    #print(q_eta_phis.shape)
    #print(subjets_eta_phis.shape)
    return np.sqrt(np.square(subjets_eta_phis[:,:,:,0] - q_eta_phis[:,:,:,0]) + 
            np.square(ang_dist(subjets_eta_phis[:,:,:,1], q_eta_phis[:,:,:,1] )))




outdir = "CASE_gen_matching_checks_jan3/"
if(not os.path.exists(outdir)): os.system("mkdir %s" % outdir)
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

subjet_rw = False
excjet_rw = True
sys = ""

jetR = 0.4

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

n = pf_cands1.shape[0]
j_subjets = np.zeros((n, n_prongs[j_idx], 4), dtype = np.float32)

gen_parts_eta_phi = gen_parts[:,:,1:3]

if(j_idx == 0): 
    for i in range(n): j_subjets[i], splittings1 = get_splittings(pf_cands1[i], jetR = jetR, num_excjets = n_prongs[0])
    dists = get_dists(gen_parts_eta_phi, j_subjets[:,:,1:3])
else : 
    for i in range(n): j_subjets[i], splittings2 = get_splittings(pf_cands2[i], jetR = jetR, num_excjets = n_prongs[1])
    dists = get_dists(gen_parts_eta_phi, j_subjets[:,:,1:3])
    #if(i ==0):
        #print(j1_subjets[i])
        #print(splittings1)




print(dists.shape)
#want closest gen-quark for each subjet 
#(NOT other way around because only doing subjets for one jet, so not all quarks will be matched)
dists = np.swapaxes(dists, 1,2)

j_closest = np.amin(dists, axis = -1)
j_which = np.argmin(dists, axis = -1)

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
                title = "%s : quark-subjet matching (%i prongs)" % (label, n_prongs[j_idx]), num_bins = 40, normalize = True, fname = outdir + label + "_j%i_quark_subjet_matching.png" % j_idx)

j_subjet_pts = j_subjets[:,:,0].reshape(-1)
make_histogram([j_subjet_pts], ["Subjets"], colors = ['blue'], xaxis_label = 'Subjet pt (GeV)', 
                title = "%s : subjet pt (%i prongs)" % (label, n_prongs[j_idx]), num_bins = 40, normalize = True, fname = outdir + label + "_j%i_subjet_pt.png" % j_idx)

print("Fraction of subjets with pt > 350 : %.3f" % (np.mean(j_subjet_pts > 350.)))
print("Pt scaling unc would be +/- %.3f" %( np.mean(j_subjet_pts > 350.) * 0.21 * (n_prongs[j_idx]) ))

decay_pdgids = np.abs(gen_parts[:,:,3]).reshape(-1)
make_histogram(decay_pdgids, "", colors = 'blue', xaxis_label = 'Decay Product PdgID', 
                title = "%s : decay products (%i+%i prongs)" % (label, n_prongs[0], n_prongs[1]), num_bins = 16, h_range = (0.5, 16.5), normalize = True, fname = outdir + label + "_decay_pdgid.png")

