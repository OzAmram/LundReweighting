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

parser = input_options()
options = parser.parse_args()

print(options)



outdir = options.outdir
f_ttbar = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_may22/TT.h5", "r")
if(not os.path.exists(outdir)): os.system("mkdir %s" % outdir)

pt_max = 1000


dr_bin_min = -1.
dr_bin_max = 8.
#y_bin_min = np.log(1./0.5)
#y_bin_max = 20*y_bin_min
#y_label = "ln(1/z)"
kt_bin_min = -5
kt_bin_max = np.log(pt_max)
z_label = "ln(kt/GeV)"
y_label = "ln(0.8/#Delta)"
n_bins_LP = 20
n_bins = 40

kt_bins = array('f', np.linspace(kt_bin_min, kt_bin_max, num = n_bins_LP+1))
dr_bins = array('f', np.linspace(dr_bin_min, dr_bin_max, num = n_bins_LP+1))


n_pt_bins = 6
pt_bins = array('f', [0., 50., 100., 175., 250., 350., 99999.])

h_qs = ROOT.TH3F("l_LP", "Lund Plane MC", n_pt_bins, pt_bins, n_bins_LP,  dr_bins, n_bins_LP, kt_bins) 
h_bs = ROOT.TH3F("b_LP", "Lund Plane MC", n_pt_bins, pt_bins, n_bins_LP,  dr_bins, n_bins_LP, kt_bins) 
h_ratio = ROOT.TH3F("h_bl_ratio", "Lund Plane MC", n_pt_bins, pt_bins, n_bins_LP,  dr_bins, n_bins_LP, kt_bins) 


h_qs.GetZaxis().SetTitle(z_label)
h_qs.GetYaxis().SetTitle(y_label)

h_bs.GetZaxis().SetTitle(z_label)
h_bs.GetYaxis().SetTitle(y_label)


subjet_rw = False
excjet_rw = True
sys = ""

jetR = 1.0

max_evts = 10000

jms_corr = 1.0

m_cut_top_min = 125.
m_cut_top_max = 225.
#m_cut_top_max = 130.
pt_cut_top = 500.


#for bad matching check
deltaR_cut = 0.2


d_ttbar_top_match = Dataset(f_ttbar, label = "ttbar : t-matched ", color = ROOT.kOrange-3, jms_corr = jms_corr)

ttbar_gen_matching = d_ttbar_top_match.f['gen_parts'][:,0]

t_match_cut = (ttbar_gen_matching  > 1.9) &  (ttbar_gen_matching < 2.1)


d_ttbar_top_match.apply_cut(t_match_cut)

top_jet_kinematics = d_ttbar_top_match.f['jet_kinematics'][:]
top_msd_cut_mask = (top_jet_kinematics[:,3] * jms_corr > m_cut_top_min) & (top_jet_kinematics[:,3] * jms_corr < m_cut_top_max)
top_pt_cut_mask = top_jet_kinematics[:,0] > pt_cut_top
d_ttbar_top_match.apply_cut(top_msd_cut_mask & top_pt_cut_mask)


jet_kin_top = d_ttbar_top_match.get_masked('jet_kinematics')
gen_parts_top = d_ttbar_top_match.get_masked('gen_parts')
pf_cands_top = d_ttbar_top_match.get_masked('jet1_PFCands').astype(np.float64)

top_subjets = np.zeros((pf_cands_top.shape[0], 3, 4), dtype = np.float32)
splittings = [None] * pf_cands_top.shape[0]


rw = LundReweighter(jetR = jetR, charge_only = False)

top_pt = d_ttbar_top_match.get_masked('jet_kinematics')[:,0]

for i,pf_cand_top in enumerate(pf_cands_top):
    top_subjets[i], splittings[i] = rw.get_splittings(pf_cand_top, num_excjets = 3, rescale_subjets = "vec", rescale_val = top_pt[i])

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
top_all_far = (top_q1_close > deltaR_cut) | (top_q2_close > deltaR_cut) | (top_b_close > deltaR_cut)




for i in range(len(top_subjets)):
    if(top_qs_samejet[i] or top_all_far[i]): continue

    for j in range(len(top_subjets[i])):
        if(j == top_q1_which[i] or j == top_q2_which[i]): 
            #fill light quark LP
            rw.fill_lund_plane(h_qs, subjets = np.array(top_subjets[i]).reshape(-1), splittings = splittings[i], subjet_idx = j)
        else: 
            #b quark LP
            rw.fill_lund_plane(h_bs, subjets = np.array(top_subjets[i]).reshape(-1), splittings = splittings[i], subjet_idx = j)

h_qs.Print()
h_bs.Print()


f_out = ROOT.TFile.Open(outdir + "ratio.root", "UPDATE")
w_qs = f_out.Get("mc_nom")

for i in range(1,h_qs.GetNbinsX()+1):
    h_qs_clone1 = h_qs.Clone("h_bkg_clone%i" %i)
    h_qs_W_clone1 = w_qs.Clone("h_w_clone%i" %i)
    h_bs_clone1 = h_bs.Clone("h_mc_clone%i"% i)

    h_qs_clone1.GetXaxis().SetRange(i,i)
    h_qs_W_clone1.GetXaxis().SetRange(i,i)
    h_bs_clone1.GetXaxis().SetRange(i,i)


    h_qs_proj = h_qs_clone1.Project3D("zy")
    h_qs_W_proj = h_qs_W_clone1.Project3D("zy")
    h_bs_proj = h_bs_clone1.Project3D("zy")


    h_qs_proj.SetTitle("MC Light Quarks Pt %i-%i (N = %.0f)" % (pt_bins[i-1], pt_bins[i], h_qs_proj.Integral()))
    h_qs_W_proj.SetTitle("MC Light Quarks (from W) Pt %i-%i (N = %.0f)" % (pt_bins[i-1], pt_bins[i], h_qs_proj.Integral()))
    h_bs_proj.SetTitle("MC b Quarks Pt %i-%i (N = %.0f)" % (pt_bins[i-1], pt_bins[i], h_bs_proj.Integral()))


    print(h_qs_proj.Integral(), h_bs_proj.Integral())

    h_qs_proj.Scale(1./h_qs_proj.Integral())
    h_qs_W_proj.Scale(1./h_qs_W_proj.Integral())
    h_bs_proj.Scale(1./h_bs_proj.Integral())

    h_ratio_proj = h_bs_proj.Clone("h_ratio%i" % i)
    h_ratio_proj.Divide(h_qs_proj)
    h_ratio_proj.SetTitle("Ratio b/light Pt %i-%i " % (pt_bins[i-1], pt_bins[i]) )
    cleanup_ratio(h_ratio_proj, h_min =0., h_max = 2.0)


    h_ratio2_proj = h_bs_proj.Clone("h_ratio%i_v2" % i)
    h_ratio2_proj.Divide(h_qs_W_proj)
    h_ratio2_proj.SetTitle("Ratio b/light (from W) Pt %i-%i " % (pt_bins[i-1], pt_bins[i]) )
    cleanup_ratio(h_ratio2_proj, h_min =0., h_max = 2.0)


    #copy_proj(i, h_ratio_proj, h_ratio)
    #use light quarks from W
    copy_proj(i, h_ratio2_proj, h_ratio)


    c1 = ROOT.TCanvas("c1", "", 1000, 1000)
    h_qs_proj.Draw("colz")
    c1.SetRightMargin(0.2)
    c1.Print(outdir + "lundPlane_light_quarks_bin%i.png" % i)

    c1 = ROOT.TCanvas("c1", "", 1000, 1000)
    h_qs_W_proj.Draw("colz")
    c1.SetRightMargin(0.2)
    c1.Print(outdir + "lundPlane_light_quarks_fromW_bin%i.png" % i)



    c2 = ROOT.TCanvas("c2", "", 1000, 1000)
    h_bs_proj.Draw("colz")
    c2.SetRightMargin(0.2)
    c2.Print(outdir + "lundPlane_b_quarks_bin%i.png" % i)


    c3 = ROOT.TCanvas("c3", "", 1000, 1000)
    h_ratio2_proj.Draw("colz")
    c3.SetRightMargin(0.2)
    c3.Print(outdir + "lundPlane_b_light_ratio_bin%i.png" % i)



    diffs = []
    w_diffs = []
    for j in range(h_qs_proj.GetNbinsX()+1):
        for k in range(h_qs_proj.GetNbinsY()+1):
            c1 = h_qs_proj.GetBinContent(j,k)
            cw = h_qs_W_proj.GetBinContent(j,k)
            c2 = h_bs_proj.GetBinContent(j,k)

            ew = h_qs_W_proj.GetBinError(j,k)
            e1 = h_qs_proj.GetBinError(j,k)
            e2 = h_bs_proj.GetBinError(j,k)

            if(  c1 > 1e-6 and c2 > 1e-6 and  (e1/c1) < 0.2 and (e2/c2) < 0.2):
                frac_diff = abs(c1 - c2)/(c1 + c2)
                frac_diff = np.clip(frac_diff, 0.001, 100)
                diffs.append(frac_diff)

            if(  c1 > 1e-6 and cw > 1e-6 and  (e1/c1) < 0.1 and (ew/cw) < 0.1):
                frac_diff_w = abs(c1 - cw)/(c1 + cw)
                w_diffs.append(frac_diff_w)

    print("Pt bin %i, avg diff %.2f" % (i, np.mean(diffs)))
    print("Pt bin %i, avg W diff %.4f" % (i, np.mean(w_diffs)))

h_ratio.Write()
f_out.Print()
f_out.Close()

