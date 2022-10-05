from Utils import *
import os

def cleanup_ratio(h, h_min = 0., h_max = 2.):
    for i in range(1, h.GetNbinsX() + 1):
        for j in range(1, h.GetNbinsY() + 1):
            cont = h.GetBinContent(i,j)
            cont = max(h_min, min(cont, h_max))
            h.SetBinContent(i,j,cont)
    #h.GetZAxis().SetRangeUser(h_min, h_max);



#UL
lumi = 59.74
f_data = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_sep29/SingleMu_2018_merge.h5", "r")
f_ttbar = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_sep29/TT.h5", "r")
f_wjets = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_sep29/QCD_WJets.h5", "r")
#f_wjets = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_sep29/WJets_merge.h5", "r")
f_diboson = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_sep29/diboson.h5", "r")
f_tw = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_sep29/TW.h5", "r")
f_singletop = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_sep29/SingleTop_merge.h5", "r")


f_ratio = "ttbar_UL_oct4_W_rw/ratio.root"
#f_ratio = "ttbar_UL_oct3_W_rw_PS_FSR_down_squared/ratio.root"



subjet_rw = False
excjet_rw = True
sys = ""
outdir = "ttbar_UL_oct4_top/"

norm = True

jms_corr = 0.95

m_cut_min = 125.
m_cut_max = 225.
pt_cut = 500.

if(not os.path.exists(outdir)): os.system("mkdir " + outdir)

d_data = Dataset(f_data, is_data = True)

d_tw = Dataset(f_tw, label = "tW", color = ROOT.kMagenta, jms_corr = jms_corr)
d_wjets = Dataset(f_wjets, label = "W+Jets + QCD", color = ROOT.kGray, jms_corr = jms_corr)
d_diboson = Dataset(f_diboson, label = "Diboson", color = ROOT.kCyan, jms_corr = jms_corr)
d_singletop = Dataset(f_singletop, label = "Single Top", color = ROOT.kMagenta+4, jms_corr = jms_corr)


d_ttbar_w_match = Dataset(f_ttbar, label = "ttbar : W-matched", color = ROOT.kRed, jms_corr =jms_corr)
d_ttbar_t_match = Dataset(f_ttbar, label = "ttbar : t-matched ", color = ROOT.kOrange-3, jms_corr = jms_corr)
d_ttbar_nomatch = Dataset(f_ttbar, label = "ttbar : unmatched", color = ROOT.kGreen+3, jms_corr = jms_corr)

ttbar_gen_matching = d_ttbar_w_match.f['gen_parts'][:,0]

#0 is unmatched, 1 is W matched, 2 is top matched
nomatch_cut = ttbar_gen_matching < 0.1
w_match_cut = (ttbar_gen_matching  > 0.9) &  (ttbar_gen_matching < 1.1)
t_match_cut = (ttbar_gen_matching  > 1.9) &  (ttbar_gen_matching < 2.1)

d_ttbar_w_match.apply_cut(w_match_cut)
d_ttbar_t_match.apply_cut(t_match_cut)
d_ttbar_nomatch.apply_cut(nomatch_cut)


sigs = [d_ttbar_t_match]
bkgs = [d_ttbar_nomatch, d_ttbar_w_match, d_tw,  d_diboson, d_wjets, d_singletop]



if(len(sys) != 0):
    for d in sigs: 
        d.apply_sys(sys)
        d.sys_power = 2.0



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
n_bins = 20

kt_bins = array('f', np.linspace(kt_bin_min, kt_bin_max, num = n_bins_LP+1))

dr_bins = array('f', np.linspace(dr_bin_min, dr_bin_max, num = n_bins_LP+1))


fill_z = False
jetR = 1.0

num_excjets = 3


ratio_range = [0.5, 1.5]


jet_kinematics_data= f_data['jet_kinematics'][()]
msd_cut_data = (jet_kinematics_data[:,3] > m_cut_min) & (jet_kinematics_data[:,3] < m_cut_max)
pt_cut_data = jet_kinematics_data[:,0] > pt_cut
d_data.apply_cut(msd_cut_data & pt_cut_data)
d_data.compute_obs()

for d in (bkgs + sigs):

    d.norm_factor = lumi

    jet_kinematics = d.f['jet_kinematics'][:]
    msd_cut_mask = (jet_kinematics[:,3] * jms_corr > m_cut_min) & (jet_kinematics[:,3] * jms_corr < m_cut_max)
    pt_cut_mask = jet_kinematics[:,0] > pt_cut
    d.apply_cut(msd_cut_mask & pt_cut_mask)
    d.compute_obs()



#print("Num data %i. Num ttbar MC %i " % (d_data.n(), d_ttbar.n()))


num_data = np.sum(d_data.get_weights())
num_ttbar_nomatch = np.sum(d_ttbar_nomatch.get_weights())
num_ttbar_w_match = np.sum(d_ttbar_w_match.get_weights())
num_ttbar_t_match = np.sum(d_ttbar_t_match.get_weights())
num_ttbar_tot = num_ttbar_nomatch + num_ttbar_w_match + num_ttbar_t_match
num_tw = np.sum(d_tw.get_weights())

tot_bkg = 0.
for d in (d_diboson, d_wjets, d_singletop):
    tot_bkg += np.sum(d.get_weights())
print("%i data, %.0f ttbar (%.0f unmatched, %.0f W matched, %.0f t matched), %.0f tW %.0f bkg" % ( num_data, num_ttbar_tot,num_ttbar_nomatch, 
                                                                                          num_ttbar_w_match, num_ttbar_t_match, num_tw, tot_bkg))
normalization = num_data  / (num_ttbar_tot + num_tw + tot_bkg)
print("normalization", normalization)

if(norm):
    for d in (bkgs + sigs):
        d.norm_factor *= normalization



obs = ["tau21", "tau32", "tau43", "nPF", "mSoftDrop", "pt"]

colors = []
weights_nom = []
labels = []
for d in (bkgs + sigs):
    colors.append(d.color)
    weights_nom.append(d.get_weights())
    labels.append(d.label)

for l in obs:
    a = []
    for d in (bkgs + sigs):
        a.append(getattr(d, l))

    if(l == 'mSoftDrop'): 
        h_range = (m_cut_min, m_cut_max)
        n_bins_ = n_bins
    elif(l == 'nPF'): 
        h_range = (0.5,120.5)
        n_bins_ = 40
    elif(l == 'pt'): 
        h_range = (pt_cut, 1200.)
        n_bins_ = n_bins
    else: 
        n_bins_ = n_bins
        h_range = None
    make_multi_sum_ratio_histogram(data = getattr(d_data, l), entries = a, weights = weights_nom, labels = labels, h_range = h_range, drawSys = False, stack = False,
            colors = colors, axis_label = l,  title = l + " : No Reweighting", num_bins = n_bins, normalize = False, ratio_range = (0.5, 1.5), fname = outdir + l + '_ratio_before.png' )



weights_rw = copy.deepcopy(weights_nom)

f = ROOT.TFile.Open(f_ratio)
h_ratio = f.Get("h_ratio")

LP_weights = []
LP_uncs = []
for i,d in enumerate(sigs):
    print("Reweighting ", d.f, i)
    d_LP_weights, d_LP_uncs = d.reweight_LP(h_ratio, subjet_rw = subjet_rw, fill_z = fill_z, jetR = jetR, num_excjets = num_excjets, uncs = True, prefix = "3prong")
    LP_weights.append(d_LP_weights)
    LP_uncs.append(d_LP_uncs)

    #apply weights, keep normalization fixed
    old_norm = np.sum(weights_rw[len(bkgs) + i])
    weights_rw[len(bkgs) + i] *= d_LP_weights
    new_norm = np.sum(weights_rw[len(bkgs) + i])
    
    weights_rw[len(bkgs) + i] *= old_norm / new_norm


f.Close()


#print(LP_weights[0][:10])
#print(LP_uncs[0][:10])

make_histogram(LP_weights[0], "Reweighting factors", 'b', 'Weight', "Lund Plane Reweighting Factors", 20 , h_range = (0., 2.0),
     normalize=False, fname=outdir + "lundPlane_weights.png")

make_histogram(LP_uncs[0], "Fractional Uncertainties", 'b', 'Weight Fractional Uncertainty ', "Lund Plane Reweighting Factors Uncertainty", 20,
     normalize=False, fname=outdir + "lundPlane_weights_unc.png", h_range = (0., 1.5))

for l in obs:
    a = []
    for d in (bkgs + sigs):
        a.append(getattr(d, l))
    if(l == 'mSoftDrop'): 
        h_range = (m_cut_min, m_cut_max)
        n_bins_ = n_bins
    elif(l == 'nPF'): 
        h_range = (0.5,120.5)
        n_bins_ = 40
    elif(l == 'pt'): 
        h_range = (pt_cut, 1200.)
        n_bins_ = n_bins
    else: 
        n_bins_ = n_bins
        h_range = None
    make_multi_sum_ratio_histogram(data = getattr(d_data, l), entries = a, weights = weights_rw, labels = labels, h_range = h_range, drawSys = False, stack = False,
            colors = colors, axis_label = l,  title = l + " : LP Reweighting", num_bins = n_bins, normalize = False, ratio_range = (0.5, 1.5), fname = outdir + l + '_ratio_after.png' )

