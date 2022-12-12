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
    q_eta_phis = np.expand_dims(q_eta_phis, 1)
    return np.sqrt(np.square(subjets_eta_phis[:,:,0] - q_eta_phis[:,:,0]) + 
            np.square(ang_dist(subjets_eta_phis[:,:,1], q_eta_phis[:,:,1] )))




f_ttbar = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_sep29/TT.h5", "r")
outdir = "gen_matching_check_dec7/"
if(not os.path.exists(outdir)): os.system("mkdir %s" % outdir)
jet_str = 'CA'


subjet_rw = False
excjet_rw = True
sys = ""

jetR = 1.0

max_evts = 10000

jms_corr = 0.95

m_cut_top_min = 125.
m_cut_top_max = 225.
pt_cut_top = 500.


m_cut_w_min = 60.
m_cut_w_max = 110.
pt_cut_w = 225.


d_ttbar_w_match = Dataset(f_ttbar, label = "ttbar : W-matched", color = ROOT.kRed, jms_corr =jms_corr)
d_ttbar_top_match = Dataset(f_ttbar, label = "ttbar : t-matched ", color = ROOT.kOrange-3, jms_corr = jms_corr)

ttbar_gen_matching = d_ttbar_w_match.f['gen_parts'][:,0]

w_match_cut = (ttbar_gen_matching  > 0.9) &  (ttbar_gen_matching < 1.1)
t_match_cut = (ttbar_gen_matching  > 1.9) &  (ttbar_gen_matching < 2.1)


d_ttbar_w_match.apply_cut(w_match_cut)
d_ttbar_top_match.apply_cut(t_match_cut)

top_jet_kinematics = d_ttbar_top_match.f['jet_kinematics'][:]
top_msd_cut_mask = (top_jet_kinematics[:,3] * jms_corr > m_cut_top_min) & (top_jet_kinematics[:,3] * jms_corr < m_cut_top_max)
top_pt_cut_mask = top_jet_kinematics[:,0] > pt_cut_top
d_ttbar_top_match.apply_cut(top_msd_cut_mask & top_pt_cut_mask)

w_jet_kinematics = d_ttbar_w_match.f['jet_kinematics'][:]
w_msd_cut_mask = (w_jet_kinematics[:,3] * jms_corr > m_cut_w_min) & (w_jet_kinematics[:,3] * jms_corr < m_cut_w_max)
w_pt_cut_mask = w_jet_kinematics[:,0] > pt_cut_w
d_ttbar_w_match.apply_cut(w_msd_cut_mask & w_pt_cut_mask)





pf_cands_w = d_ttbar_w_match.get_masked('jet1_PFCands')[:max_evts].astype(np.float64)

w_subjets = np.zeros((pf_cands_w.shape[0], 2, 4), dtype = np.float32)

for i,pf_cand_w in enumerate(pf_cands_w):
    w_subjets[i], splittings = get_splittings(pf_cand_w, jetR = jetR, num_excjets = 2)
    #print(w_subjets[i])
    #print(splittings)

gen_parts_w = d_ttbar_w_match.get_masked('gen_parts')[:max_evts]
w_q1_eta_phi = gen_parts_w[:,18:20]
w_q2_eta_phi = gen_parts_w[:,22:24]

w_q1_dists = get_dists(w_q1_eta_phi, w_subjets[:,:,1:3])
w_q2_dists = get_dists(w_q2_eta_phi, w_subjets[:,:,1:3])

w_q1_close = np.amin(w_q1_dists, axis = -1)
w_q2_close = np.amin(w_q2_dists, axis = -1)

w_all_far = (w_q1_close > 0.4) | (w_q2_close > 0.4)

w_q1_which = np.argmin(w_q1_dists, axis = -1)
w_q2_which = np.argmin(w_q2_dists, axis = -1)
w_qs_samejet = (w_q1_which == w_q2_which)

print("\nW-subjet matching :")
print("Q1 dists, avg %.3f, std %.3f. Frac over 0.4 %.3f" % (np.mean(w_q1_close), np.std(w_q1_close), np.mean(w_q1_close > 0.4)))
print("Q2 dists, avg %.3f, std %.3f. Frac over 0.4 %.3f" % (np.mean(w_q2_close), np.std(w_q2_close), np.mean(w_q2_close > 0.4)))
print("Overall frac > 0.4 %.3f" % np.mean(w_all_far))
print("Fraction of quarks matched to same subjet %.4f" % (np.mean(w_qs_samejet)))
print("Overall bad matching frac %.4f" % np.mean(np.mean(w_qs_samejet | w_all_far)) )

w_subjets_pt = w_subjets[:,:,0].reshape(-1)
make_histogram(w_subjets_pt, "W subjets", colors = 'blue', xaxis_label = 'Subjet pt (GeV)', 
                title = "W-jets : subjet pt", num_bins = 40, normalize = True, fname = outdir + "W_subjet_pt.png")


make_histogram([w_q1_close, w_q2_close], ['Q1', 'Q2'], colors = ['red', 'blue'], xaxis_label = r'$\Delta R$ to closest sub-jet', 
                title = 'W-jets : quark-subjet matching %s' % jet_str, num_bins = 40, normalize = True, fname = outdir + "W_quark_%s_matching.png" % jet_str)

jet_kin_top = d_ttbar_top_match.get_masked('jet_kinematics')[:max_evts]
gen_parts_top = d_ttbar_top_match.get_masked('gen_parts')[:max_evts]
pf_cands_top = d_ttbar_top_match.get_masked('jet1_PFCands')[:max_evts].astype(np.float64)

top_subjets = np.zeros((pf_cands_top.shape[0], 3, 4), dtype = np.float32)

for i,pf_cand_top in enumerate(pf_cands_top):
    top_subjets[i], _ = get_splittings(pf_cand_top, jetR = jetR, num_excjets = 3)

top_q1_eta_phi = gen_parts_top[:,18:20]
top_q2_eta_phi = gen_parts_top[:,22:24]
top_b_eta_phi = gen_parts_top[:,26:28]


top_q1_dists = get_dists(top_q1_eta_phi, top_subjets[:,:,1:3])
top_q2_dists = get_dists(top_q2_eta_phi, top_subjets[:,:,1:3])
top_b_dists = get_dists(top_b_eta_phi, top_subjets[:,:,1:3])

top_q1_close = np.amin(top_q1_dists, axis = -1)
top_q2_close = np.amin(top_q2_dists, axis = -1)
top_b_close = np.amin(top_b_dists, axis = -1)

top_q1_which = np.argmin(top_q1_dists, axis = -1)
top_q2_which = np.argmin(top_q2_dists, axis = -1)
top_b_which = np.argmin(top_b_dists, axis = -1)

#if all different should be 0 1 2
top_qs_samejet = (top_q1_which == top_q2_which) | (top_q1_which == top_b_which) |  (top_q2_which == top_b_which)
top_all_far = (top_q1_close > 0.4) | (top_q2_close > 0.4) | (top_b_close > 0.4)

print("\ntop-subjet matching :")
print("Q1 dists, avg %.3f, std %.3f. Frac over 0.4 %.3f" % (np.mean(top_q1_close), np.std(top_q1_close), np.mean(top_q1_close > 0.4)))
print("Q2 dists, avg %.3f, std %.3f. Frac over 0.4 %.3f" % (np.mean(top_q2_close), np.std(top_q2_close), np.mean(top_q2_close > 0.4)))
print("b dists, avg %.3f, std %.3f. Frac over 0.4 %.3f" % (np.mean(top_b_close), np.std(top_b_close), np.mean(top_b_close > 0.4)))
print("Overall frac > 0.4 %.3f" % np.mean(top_all_far))
print("Fraction of quarks matched to same subjet %.4f" % (np.mean(top_qs_samejet)))
print("Overall bad matching frac %.4f" % np.mean(np.mean(top_qs_samejet | top_all_far)) )


make_histogram([top_q1_close, top_q2_close, top_b_close], ['Q1', 'Q2', 'b Quark'], colors = ['red', 'blue', 'green'], xaxis_label = r'$\Delta R$ to closest sub-jet', 
                title = 'top-jets : quark-subjet matching %s' % jet_str, num_bins = 40, normalize = True, fname = outdir + "top_quark_%s_matching.png" % jet_str)

top_subjets_pt = top_subjets[:,:,0].reshape(-1)
print("Frac with pt > 350 : %.3f" % np.mean(top_subjets_pt > 350.))
make_histogram(top_subjets_pt, "top subjets", colors = 'blue', xaxis_label = 'Subjet pt (GeV)', 
                title = "top-jets : subjet pt", num_bins = 40, normalize = True, fname = outdir + "top_subjet_pt.png")
