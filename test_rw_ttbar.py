from Utils import *
import os


def cleanup_hist(h):
    if(type(h) == ROOT.TH3F):
        for i in range(h.GetNbinsX()+1):
            for j in range(h.GetNbinsY()+1):
                for k in range(h.GetNbinsZ()+1):
                    c = h.GetBinContent(i,j,k)
                    if( c < 0): 
                        h.SetBinContent(i,j,k, 0.)
                        h.SetBinError(i,j,k, 0.)
    elif(type(h) == ROOT.TH2F):
        for i in range(h.GetNbinsX()+1):
            for j in range(h.GetNbinsY()+1):
                c = h.GetBinContent(i,j)
                if( c < 0): 
                    h.SetBinContent(i,j, 0.)
                    h.SetBinError(i,j, 0.)

def get_unc_hist(h):
    h_unc = h.Clone(h.GetName() + "_unc")
    for i in range(1, h.GetNbinsX() + 1):
        for j in range(1, h.GetNbinsY() + 1):
            err = h.GetBinError(i,j)
            cont = h.GetBinContent(i,j)

            if(cont > 0):
                h_unc.SetBinContent(i,j, err/cont)
                h_unc.SetBinError(i,j, 0.)
            else:
                h_unc.SetBinContent(i,j, 0.)
    return h_unc


def copy_proj(bin_i, h_ratio_proj, h_ratio):
    for j in range(h_ratio.GetNbinsY()):
        for k in range(h_ratio.GetNbinsZ()):
            h_ratio.SetBinContent(bin_i,j,k, h_ratio_proj.GetBinContent(j,k))
            h_ratio.SetBinError(bin_i,j,k, h_ratio_proj.GetBinError(j,k))
    return

#NON UL (2018B only)
lumi = 6.90
f_data = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/CASEUtils/H5_maker/ttbar_output_files_v2/SingleMu_2018C.h5", "r")
f_ttbar = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/CASEUtils/H5_maker/ttbar_output_files_v2/TTToSemiLep_2018.h5", "r")

#UL
#lumi = 59.74
#f_data = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files/SingleMu_2018_total.h5", "r")
#f_ttbar = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files/TTbar_semilep_2018.h5", "r")

f_bkg = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/CASEUtils/H5_maker/ttbar_output_files_v2/QCD_WJets_merged.h5", "r")


subjet_rw = False
excjet_rw = True
outdir = "ttbar_test_UL_sep2/"


if(not os.path.exists(outdir)): os.system("mkdir " + outdir)


pt_max = 1000

norm = False

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
n_bins = 12

kt_bins = array('f', np.linspace(kt_bin_min, kt_bin_max, num = n_bins_LP+1))

dr_bins = array('f', np.linspace(dr_bin_min, dr_bin_max, num = n_bins_LP+1))



fill_z = False
jetR = -1.0

num_excjets = -1

if(subjet_rw):
    jetR = 0.4
    n_pt_bins = 5
    pt_bins = array('f', [0., 10., 25., 40., 60., 99999.])
elif(excjet_rw):
    jetR = 1.0 #idk if this matters ? 
    n_pt_bins = 5
    num_excjets = 2
    pt_bins = array('f', [0., 50., 100., 200., 300., 450.])
else:
    n_pt_bins = 4
    pt_bins = array('f', [200., 300., 400., 600., 99999.])
    #pt_bins = array('f', [200., 300.])


ratio_range = [0.5, 1.5]
h_mc = ROOT.TH3F("lp_mc", "Lund Plane MC", n_pt_bins, pt_bins, n_bins_LP,  dr_bins, n_bins_LP, kt_bins) 
h_bkg = ROOT.TH3F("lp_bkg", "Lund Plane Bkg", n_pt_bins, pt_bins, n_bins_LP,  dr_bins, n_bins_LP, kt_bins) 
h_data = ROOT.TH3F("lp_data", "Lund Plane Data", n_pt_bins, pt_bins, n_bins_LP, dr_bins, n_bins_LP, kt_bins) 




jet_kinematics_data= f_data['jet_kinematics'][()]
pf_cands_data = f_data['jet1_PFCands'][()].astype(np.float64)
jet1_feats_data = f_data['jet1_extraInfo'][()]
event_info_data = f_data['event_info'][()]


jet_kinematics_ttbar= f_ttbar['jet_kinematics'][()]
pf_cands_ttbar = f_ttbar['jet1_PFCands'][()].astype(np.float64)
jet1_feats_ttbar = f_ttbar['jet1_extraInfo'][()]
event_info_ttbar = f_ttbar['event_info'][()]
weights_ttbar_raw = f_ttbar['norm_weights'][()].reshape(-1)


jet_kinematics_bkg= f_bkg['jet_kinematics'][()]
pf_cands_bkg = f_bkg['jet1_PFCands'][()].astype(np.float64)
jet1_feats_bkg = f_bkg['jet1_extraInfo'][()]
event_info_bkg = f_bkg['event_info'][()]
label_bkg = f_bkg['truth_label'][()]
weights_bkg_raw = f_bkg['norm_weights'][()].reshape(-1)

#soft drop mass cut at 10 gev
msd_cut_ttbar = jet_kinematics_ttbar[:,3] > 10.
msd_cut_data = jet_kinematics_data[:,3] > 10.
msd_cut_bkg = jet_kinematics_bkg[:,3] > 10.
#msd_cut_ttbar = jet_kinematics_ttbar[:,3] > 0.
#msd_cut_data = jet_kinematics_data[:,3] > 0.

#pt_cut_ttbar = (jet_kinematics_ttbar[:,0] > 200.) & (jet_kinematics_ttbar[:,0]  <  300.)
#pt_cut_data = (jet_kinematics_data[:,0] > 200.) & (jet_kinematics_data[:,0]  <  300.)


pt_cut_ttbar = (jet_kinematics_ttbar[:,0] > 200.)
pt_cut_data = (jet_kinematics_data[:,0] > 200.)
pt_cut_bkg = (jet_kinematics_bkg[:,0] > 200.)

cut_data = msd_cut_data & pt_cut_data
cut_ttbar = msd_cut_ttbar & pt_cut_ttbar
cut_bkg = msd_cut_bkg & pt_cut_bkg




jet_kinematics_data = jet_kinematics_data[cut_data]
pf_cands_data = pf_cands_data[cut_data]
jet1_feats_data = jet1_feats_data[cut_data]
event_info_data = event_info_data[cut_data]


jet_kinematics_ttbar = jet_kinematics_ttbar[cut_ttbar]
pf_cands_ttbar = pf_cands_ttbar[cut_ttbar]
jet1_feats_ttbar = jet1_feats_ttbar[cut_ttbar]
event_info_ttbar = event_info_ttbar[cut_ttbar]
weights_ttbar_raw = weights_ttbar_raw[cut_ttbar] * lumi
print(pf_cands_ttbar.shape, jet_kinematics_ttbar.shape)


jet_kinematics_bkg = jet_kinematics_bkg[cut_bkg]
pf_cands_bkg = pf_cands_bkg[cut_bkg]
jet1_feats_bkg = jet1_feats_bkg[cut_bkg]
event_info_bkg = event_info_bkg[cut_bkg]
label_bkg = label_bkg[cut_bkg].reshape(-1)
weights_bkg_raw = weights_bkg_raw[cut_bkg] * lumi

QCD_mask = (label_bkg == 0)
WJets_mask = (label_bkg == -1)


num_ttbar = pf_cands_ttbar.shape[0]
num_data = pf_cands_data.shape[0]
num_bkg = pf_cands_data.shape[0]

print("Num data %i. Num ttbar MC %i " % (num_data, num_ttbar))

print(weights_bkg_raw.shape)
print(label_bkg.shape)


#print(weights_ttbar_raw[:5])
#print(np.mean(weights_ttbar_raw))
print("%i data, %.0f ttbar %.0f wjets %.0f qcd" % ( num_data, np.sum(weights_ttbar_raw), np.sum(weights_bkg_raw[WJets_mask]), np.sum(weights_bkg_raw[QCD_mask])))
normalization = num_data  / (np.sum(weights_ttbar_raw) + np.sum(weights_bkg_raw))
print("normalization", normalization)

weights_ttbar = weights_ttbar_raw * normalization
weights_bkg = weights_bkg_raw * normalization

weights_nom = [weights_ttbar, [1.]*num_data]

print(np.sum(weights_nom[0]), np.sum(weights_nom[1]))



eps = 1e-8


tau21_ttbar = (jet1_feats_ttbar[:,1] / (jet1_feats_ttbar[:,0] + eps))
tau32_ttbar = (jet1_feats_ttbar[:,2] / (jet1_feats_ttbar[:,1] + eps))
tau43_ttbar = (jet1_feats_ttbar[:,3] / (jet1_feats_ttbar[:,2] + eps))
nPF_ttbar= jet1_feats_ttbar[:,6]


tau21_bkg = (jet1_feats_bkg[:,1] / (jet1_feats_bkg[:,0] + eps))
tau32_bkg = (jet1_feats_bkg[:,2] / (jet1_feats_bkg[:,1] + eps))
tau43_bkg = (jet1_feats_bkg[:,3] / (jet1_feats_bkg[:,2] + eps))
nPF_bkg= jet1_feats_bkg[:,6]



tau21_data = (jet1_feats_data[:,1] / (jet1_feats_data[:,0] + eps))
tau32_data = (jet1_feats_data[:,2] / (jet1_feats_data[:,1] + eps))
tau43_data = (jet1_feats_data[:,3] / (jet1_feats_data[:,2] + eps))
nPF_data= jet1_feats_data[:,6]

root_colors = [ROOT.kBlue, ROOT.kRed]

weights_nom = [weights_bkg, weights_ttbar]

make_stack_ratio_histogram(data = tau21_data, entries = [tau21_bkg, tau21_ttbar], weights = weights_nom, labels = ["QCD + WJets", "ttbar"], 
    colors = root_colors, axis_label = 'Tau21',  title = 'Tau21 : No Reweighting', num_bins = n_bins, normalize = False, ratio_range = (0.5, 1.5), fname = outdir + 'tau21_ratio_before.png' )

make_stack_ratio_histogram(data = tau32_data, entries = [tau32_bkg, tau32_ttbar], weights = weights_nom, labels = ["QCD + WJets", "ttbar"], 
    colors = root_colors, axis_label = 'Tau32',  title = 'Tau32 : No Reweighting', num_bins = n_bins, normalize = False, ratio_range = (0.5, 1.5), fname = outdir + 'tau32_ratio_before.png' )


make_stack_ratio_histogram(data = tau43_data, entries = [tau43_bkg, tau43_ttbar], weights = weights_nom, labels = ["QCD + WJets", "ttbar"], 
    colors = root_colors, axis_label = 'Tau43',  title = 'Tau43 : No Reweighting', num_bins = n_bins, normalize = False, ratio_range = (0.5, 1.5), fname = outdir + 'tau43_ratio_before.png' )

make_stack_ratio_histogram(data = nPF_data, entries = [nPF_bkg, nPF_ttbar], weights = weights_nom, labels = ["QCD + WJets", "ttbar"], 
    colors = root_colors, axis_label = 'nPF',  title = 'nPF : No Reweighting', num_bins = n_bins, normalize = False, ratio_range = (0.5, 1.5), fname = outdir + 'nPF_ratio_before.png' )


make_stack_ratio_histogram(data = jet_kinematics_data[:,3], entries = [jet_kinematics_bkg[:,3], jet_kinematics_ttbar[:,3]], 
        weights = weights_nom , labels = ["QCD + WJets", "ttbar"], colors = root_colors, axis_label = 'AK8 Softdrop Mass',  title = 'Softdrop Mass: No Reweighting', num_bins = n_bins, 
    normalize = False, ratio_range = (0.5, 1.5), h_range = (0., 300.), fname = outdir + 'jet_mass_ratio_before.png' )

make_stack_ratio_histogram(data = jet_kinematics_data[:,0], entries = [jet_kinematics_bkg[:,0], jet_kinematics_ttbar[:,0]], 
        weights = weights_nom , labels = ["QCD + WJets", "ttbar"], colors = root_colors, axis_label = 'AK8 Jet pT',  title = 'AK8 Jet pT : No Reweighting', num_bins = n_bins, 
    normalize = False, ratio_range = (0.5, 1.5), h_range = (200., 1000.), fname = outdir + 'jet_pt_ratio_before.png' )


make_stack_ratio_histogram(data = event_info_data[:,1], entries = [event_info_bkg[:,1], event_info_ttbar[:,1]], 
        weights = weights_nom , labels = ["QCD + WJets", "ttbar"], colors = root_colors, axis_label = 'MET',  title = 'MET', num_bins = n_bins, 
    normalize = False, ratio_range = (0.5, 1.5), fname = outdir + 'MET_ratio_before.png' )



#make_ratio_histogram([tau21_ttbar[:num_rw], tau21_data], ["ttbar MC", "Data"], ['b', 'r'], 'Tau21', "Tau21 : No Reweighting", n_bins,
#     ratio_range = ratio_range, normalize=norm,  weights = weights_nom, save = True, fname=outdir + "tau21_ratio_before.png")
#
#make_ratio_histogram([tau32_ttbar[:num_rw], tau32_data], ["ttbar MC", "Data"], ['b', 'r'], 'tau32', "tau32 : No Reweighting", n_bins,
#     ratio_range = ratio_range, normalize=norm,  weights = weights_nom, save = True, fname=outdir + "tau32_ratio_before.png")
#
#make_ratio_histogram([tau43_ttbar[:num_rw], tau43_data], ["ttbar MC", "Data"], ['b', 'r'], 'tau43', "tau43 : No Reweighting", n_bins,
#     ratio_range = ratio_range, normalize=norm,  weights = weights_nom, save = True, fname=outdir + "tau43_ratio_before.png")
#
#make_ratio_histogram([nPF_ttbar[:num_rw], nPF_data], ["ttbar MC", "Data"], ['b', 'r'], 'nPF', "nPF : No Reweighting", n_bins,
#     ratio_range = ratio_range, normalize=norm,  weights = weights_nom, save = True, fname=outdir + "nPF_ratio_before.png")
#
#make_ratio_histogram([jet_kinematics_ttbar[:num_rw,0], jet_kinematics_data[:,0]], ["ttbar MC", "Data"], ['b', 'r'], 'AK8 Jet pT', "AK8 Jet PT : No Reweighting", n_bins,
#     ratio_range = ratio_range, normalize=norm,  weights = weights_nom, save = True, fname=outdir + "jet_pt_ratio_before.png")
#
#make_ratio_histogram([jet_kinematics_ttbar[:num_rw, 3], jet_kinematics_data[:, 3]], ["ttbar MC", "Data"], ['b', 'r'], 'AK8 Softdrop Mass', "Softdrop Mass : No Reweighting", n_bins,
#     h_range = (0., 300.), ratio_range = ratio_range, normalize=norm,  weights = weights_nom, save = True, fname=outdir + "jet_mass_ratio_before.png")
#
#make_ratio_histogram([event_info_ttbar[:num_rw,1], event_info_data[:,1]], ["ttbar MC", "Data"], ['b', 'r'], 'MET', "MET", n_bins,
#     ratio_range = ratio_range, normalize=norm,  weights = weights_nom, save = True, fname=outdir + "MET_ratio_before.png")
#

LP_weights = []

for i,pf_cand in enumerate(pf_cands_ttbar):
    weight =weights_ttbar[i]
    if(subjet_rw):
        pt_eta_phi_m_vec = jet_kinematics_ttbar[i]
        jet_4vec = convert_4vec(pt_eta_phi_m_vec)
        boost_vec = fj.PseudoJet(jet_4vec[0], jet_4vec[1], jet_4vec[2], jet_4vec[3])
        fill_lund_plane(h_mc, pf_cand,  boost_vec = boost_vec, fill_z =fill_z, dR = jetR, jetR = jetR, weight = weight)
    else: fill_lund_plane(h_mc, pf_cand, fill_z = fill_z, jetR = jetR, num_excjets = num_excjets, weight = weight)


for i,pf_cand in enumerate(pf_cands_bkg):
    weight =weights_bkg[i]
    if(subjet_rw):
        pt_eta_phi_m_vec = jet_kinematics_bkg[i]
        jet_4vec = convert_4vec(pt_eta_phi_m_vec)
        boost_vec = fj.PseudoJet(jet_4vec[0], jet_4vec[1], jet_4vec[2], jet_4vec[3])
        fill_lund_plane(h_bkg, pf_cand,  boost_vec = boost_vec, fill_z =fill_z, dR = jetR, jetR = jetR, weight = weight)
    else: fill_lund_plane(h_bkg, pf_cand, fill_z = fill_z, jetR = jetR, num_excjets = num_excjets, weight = weight)


for i,pf_cand in enumerate(pf_cands_data):
    weight = 1.
    if(subjet_rw):
        pt_eta_phi_m_vec = jet_kinematics_data[i]
        jet_4vec = convert_4vec(pt_eta_phi_m_vec)
        boost_vec = fj.PseudoJet(jet_4vec[0], jet_4vec[1], jet_4vec[2], jet_4vec[3])
        fill_lund_plane(h_data, pf_cand,  boost_vec = boost_vec, fill_z =fill_z, dR = jetR, jetR = jetR, weight = weight)
    else: fill_lund_plane(h_data, pf_cand, fill_z = fill_z, jetR = jetR, num_excjets = num_excjets, weight = weight)



h_mc.GetZaxis().SetTitle(z_label)
h_mc.GetYaxis().SetTitle(y_label)
h_bkg.GetZaxis().SetTitle(z_label)
h_bkg.GetYaxis().SetTitle(y_label)
h_data.GetZaxis().SetTitle(z_label)
h_data.GetYaxis().SetTitle(y_label)


h_ratio = h_data.Clone("h_ratio")
h_ratio.SetTitle("(Data - Bkg ) / TTbar MC")

h_data.Print()

h_data_sub = h_data.Clone("h_data_sub")
h_data_sub.Add(h_bkg, -1.)
h_data_sub.Print()

cleanup_hist(h_data_sub)
h_data_sub.Print()

default = ROOT.TStyle("Default","Default Style");
default.cd()
ROOT.gROOT.SetStyle('Default')
ROOT.gStyle.SetOptStat(0) # To display the mean and RMS:   SetOptStat("mr")

for i in range(1, h_data.GetNbinsX() + 1):
    h_bkg_clone = h_bkg.Clone("h_bkg_clone%i" %i)
    h_mc_clone = h_mc.Clone("h_mc_clone%i"% i)
    #h_data_clone = h_data.Clone("h_data_clone%i" %i)
    h_data_clone = h_data_sub.Clone("h_data_clone%i" %i)




    h_mc_clone.GetXaxis().SetRange(i,i)
    h_bkg_clone.GetXaxis().SetRange(i,i)
    h_data_clone.GetXaxis().SetRange(i,i)

    cleanup_hist(h_mc_clone)
    cleanup_hist(h_bkg_clone)
    cleanup_hist(h_data_clone)

    h_mc_proj = h_mc_clone.Project3D("zy")
    h_bkg_proj = h_bkg_clone.Project3D("zy")
    h_data_proj = h_data_clone.Project3D("zy")

    h_mc_proj.Print()
    h_bkg_proj.Print()
    h_data_proj.Print()

    


    h_bkg_proj.Scale(1./h_bkg_proj.Integral())

    data_norm = h_data_proj.Integral()
    h_data_proj.Scale(1./data_norm)
    h_mc_proj.Scale(1./h_mc_proj.Integral())
    h_ratio_proj = h_data_proj.Clone("h_ratio_proj%i" %i)
    h_ratio_proj.Divide(h_mc_proj)

    copy_proj(i, h_ratio_proj, h_ratio)


    h_bkg_proj.SetTitle("Bkg MC pT %.0f - %.0f" % (pt_bins[i-1], pt_bins[i]))
    h_mc_proj.SetTitle("TTbar MC pT %.0f - %.0f" % (pt_bins[i-1], pt_bins[i]))
    h_data_proj.SetTitle("Data - Bkg pT %.0f - %.0f (N = %.0f)" % (pt_bins[i-1], pt_bins[i], data_norm))
    h_ratio_proj.SetTitle("Ratio pT %.0f - %.0f (N = %.0f)" % (pt_bins[i-1], pt_bins[i], data_norm))

    c_mc = ROOT.TCanvas("c", "", 1000, 800)
    h_mc_proj.Draw("colz")
    c_mc.Print(outdir + "lundPlane_bin%i_MC.png" % i)


    c_bkg = ROOT.TCanvas("c", "", 1000, 800)
    h_bkg_proj.Draw("colz")
    c_bkg.Print(outdir + "lundPlane_bin%i_bkg.png" % i)

    c_data = ROOT.TCanvas("c", "", 1000, 800)
    h_data_proj.Draw("colz")
    c_data.Print(outdir + "lundPlane_bin%i_data.png" %i )



    c_ratio = ROOT.TCanvas("c", "", 1000, 800)
    h_ratio_proj.Draw("colz")
    c_ratio.Print(outdir + "lundPlane_bin%i_ratio.png" % i)

    h_ratio_unc = get_unc_hist(h_ratio_proj)
    c_ratio_unc = ROOT.TCanvas("c_unc", "", 800, 800)
    h_ratio_unc.SetTitle("Ratio pT %.0f - %.0f (N = %.0f) Relative Unc." % (pt_bins[i-1], pt_bins[i], data_norm))
    h_ratio_unc.Draw("colz")
    c_ratio_unc.Print(outdir + "lundPlane_bin%i_ratio_unc.png" % i)
    h_ratio_unc.Reset()


f_out = ROOT.TFile.Open(outdir + "ratio.root", "RECREATE")
h_ratio.Write()
f_out.Close()


LP_weights = []
LP_uncs = []
for i,pf_cand in enumerate(pf_cands_ttbar):
    if(subjet_rw):
        pt_eta_phi_m_vec = jet_kinematics_ttbar[i]
        jet_4vec = convert_4vec(pt_eta_phi_m_vec)
        boost_vec = fj.PseudoJet(jet_4vec[0], jet_4vec[1], jet_4vec[2], jet_4vec[3])
        rw, unc = reweight_lund_plane(h_ratio, pf_cand,  boost_vec = boost_vec, fill_z =fill_z, dR = jetR, jetR = jetR, uncs = True)
    else: 
        rw, unc = reweight_lund_plane(h_ratio, pf_cand, fill_z = fill_z, jetR = jetR, num_excjets = num_excjets, uncs = True)

    LP_weights.append(rw)
    if(rw >= 1e-6):
        LP_uncs.append(unc/rw)
    else:
        LP_uncs.append(0.)


LP_weights = np.clip(LP_weights, 0., 10.)
LP_uncs = np.clip(LP_uncs, 0., 1.5)
mean = np.mean(LP_weights)
LP_weights /= mean
print(LP_weights[:10])
print(LP_uncs[:10])

weights_rw = [  weights_bkg, LP_weights * weights_ttbar]


make_histogram(LP_weights, "Reweighting factors", 'b', 'Weight', "Lund Plane Reweighting Factors", n_bins,
     normalize=norm, fname=outdir + "lundPlane_weights.png")

make_histogram(LP_uncs, "Fractional Uncertainties", 'b', 'Weight Fractional Uncertainty ', "Lund Plane Reweighting Factors Uncertainty", 20,
     normalize=norm, fname=outdir + "lundPlane_weights_unc.png")

make_stack_ratio_histogram(data = tau21_data, entries = [tau21_bkg, tau21_ttbar], weights = weights_rw, labels = ["QCD + WJets", "ttbar"], 
    colors = root_colors, axis_label = 'Tau21',  title = 'Tau21 :  LP Reweighted', num_bins = n_bins, normalize = False, ratio_range = (0.5, 1.5), fname = outdir + 'tau21_ratio_after.png' )

make_stack_ratio_histogram(data = tau32_data, entries = [tau32_bkg, tau32_ttbar], weights = weights_rw, labels = ["QCD + WJets", "ttbar"], 
    colors = root_colors, axis_label = 'Tau32',  title = 'Tau32 :  LP Reweighted', num_bins = n_bins, normalize = False, ratio_range = (0.5, 1.5), fname = outdir + 'tau32_ratio_after.png' )


make_stack_ratio_histogram(data = tau43_data, entries = [tau43_bkg, tau43_ttbar], weights = weights_rw, labels = ["QCD + WJets", "ttbar"], 
    colors = root_colors, axis_label = 'Tau43',  title = 'Tau43 :  LP Reweighted', num_bins = n_bins, normalize = False, ratio_range = (0.5, 1.5), fname = outdir + 'tau43_ratio_after.png' )

make_stack_ratio_histogram(data = nPF_data, entries = [nPF_bkg, nPF_ttbar], weights = weights_rw, labels = ["QCD + WJets", "ttbar"], 
    colors = root_colors, axis_label = 'nPF',  title = 'nPF :  LP Reweighted', num_bins = n_bins, normalize = False, ratio_range = (0.5, 1.5), fname = outdir + 'nPF_ratio_after.png' )


make_stack_ratio_histogram(data = jet_kinematics_data[:,3], entries = [jet_kinematics_bkg[:,3], jet_kinematics_ttbar[:,3]], 
        weights = weights_rw , labels = ["QCD + WJets", "ttbar"], colors = root_colors, axis_label = 'AK8 Softdrop Mass',  title = 'Softdrop Mass:  LP Reweighted', num_bins = n_bins, 
    normalize = False, ratio_range = (0.5, 1.5), h_range = (0., 300.), fname = outdir + 'jet_mass_ratio_after.png' )


#make_ratio_histogram([tau21_ttbar[num_rw:], tau21_data], ["ttbar MC", "Data"], ['b', 'r'], 'Tau21', "Tau21 : After Reweighting", n_bins,
#     ratio_range = ratio_range, normalize=norm,  weights = weights_rw, save = True, fname=outdir + "tau21_ratio_after.png")
#
#make_ratio_histogram([tau32_ttbar[num_rw:], tau32_data], ["ttbar MC", "Data"], ['b', 'r'], 'tau32', "tau32 : After Reweighting", n_bins,
#     ratio_range = ratio_range, normalize=norm,  weights = weights_rw, save = True, fname=outdir + "tau32_ratio_after.png")
#
#make_ratio_histogram([tau43_ttbar[num_rw:], tau43_data], ["ttbar MC", "Data"], ['b', 'r'], 'tau43', "tau43 : After Reweighting", n_bins,
#     ratio_range = ratio_range, normalize=norm,  weights = weights_rw, save = True, fname=outdir + "tau43_ratio_after.png")
#
#make_ratio_histogram([nPF_ttbar[num_rw:], nPF_data], ["ttbar MC", "Data"], ['b', 'r'], 'nPF', "nPF : After Reweighting", n_bins,
#     ratio_range = ratio_range, normalize=norm,  weights = weights_rw, save = True, fname=outdir + "nPF_ratio_after.png")
#
#make_ratio_histogram([jet_kinematics_ttbar[num_rw:, 3], jet_kinematics_data[:, 3]], ["ttbar MC", "Data"], ['b', 'r'], 'AK8 Softdrop Mass', "Softdrop Mass : After Reweighting", n_bins,
#        h_range = (0., 300.),    ratio_range = ratio_range, normalize=norm,  weights = weights_rw, save = True, fname=outdir + "jet_mass_ratio_after.png")


