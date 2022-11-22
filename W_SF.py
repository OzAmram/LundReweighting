from Utils import *
import os




#NON UL (2018B only)
#lumi = 6.90
#f_data = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/CASEUtils/H5_maker/ttbar_output_files_v2/SingleMu_2018C.h5", "r")
#f_ttbar = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/CASEUtils/H5_maker/ttbar_output_files_v2/TTToSemiLep_2018.h5", "r")
#f_bkg = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/CASEUtils/H5_maker/ttbar_output_files_v2/QCD_WJets_merged.h5", "r")
#
#UL
lumi = 59.74
#f_data = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files/SingleMu_2018_merge.h5", "r")
#f_ttbar = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files/TTToSemiLep_sep21.h5", "r")
#f_wjets = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files/WJets_merge.h5", "r")
#f_diboson = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files/diboson.h5", "r")
#f_tw = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files/TW.h5", "r")
#f_singletop = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files/SingleTop_merge.h5", "r")

f_data = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_sep29/SingleMu_2018_merge.h5", "r")
f_ttbar = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_sep29/TT.h5", "r")
f_wjets = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_sep29/QCD_WJets.h5", "r")
#f_wjets = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_sep29/WJets.h5", "r")
f_diboson = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_sep29/diboson.h5", "r")
f_tw = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_sep29/TW.h5", "r")
f_singletop = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_sep29/SingleTop_merge.h5", "r")




subjet_rw = False
excjet_rw = True
outdir = "ttbar_UL_W_SF_nov7_kt_sys/"
f_ratio = ROOT.TFile.Open("ttbar_UL_nov7_W_rw_kt_sys/ratio.root")
sys = ""
prefix = "2prong_kt"



#tau21_cut = 0.2700 # med
tau21_cut = 0.3452 # loose
DeepAK8_cut = 0.479

do_sys_variations = True

norm = True

jms_corr = 0.95

m_cut_min = 60.
m_cut_max = 110.
pt_cut = 225.

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


sigs = [d_ttbar_w_match]
bkgs = [d_ttbar_nomatch, d_ttbar_t_match, d_tw, d_diboson, d_wjets, d_singletop]


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
n_bins = 40

kt_bins = array('f', np.linspace(kt_bin_min, kt_bin_max, num = n_bins_LP+1))

dr_bins = array('f', np.linspace(dr_bin_min, dr_bin_max, num = n_bins_LP+1))


fill_z = False
jetR = 1.0

num_excjets = -1

if(subjet_rw):
    jetR = 0.4
    n_pt_bins = 5
    pt_bins = array('f', [0., 10., 25., 40., 60., 99999.])
elif(excjet_rw):
    jetR = 1.0 #idk if this matters ? 
    n_pt_bins = 6
    num_excjets = 2
    pt_bins = array('f', [0., 50., 100., 200., 300., 450., 99999.])
else:
    n_pt_bins = 4
    pt_bins = array('f', [200., 300., 400., 600., 99999.])
    #pt_bins = array('f', [200., 300.])




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



obs = ["tau21", "tau32", "tau43", "nPF", "mSoftDrop", "pt", "DeepAK8_W_MD"]

colors = []
weights_nom = []
uncs_nom = []
labels = []
for d in (bkgs + sigs):
    colors.append(d.color)
    weights_nom.append(d.get_weights())
    uncs_nom.append(np.zeros_like(weights_nom[-1]))
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
        h_range = (pt_cut, 800.)
        n_bins_ = n_bins
    else: 
        n_bins_ = n_bins
        h_range = None
    make_multi_sum_ratio_histogram(data = getattr(d_data, l), entries = a, weights = weights_nom, labels = labels, h_range = h_range, drawSys = False, stack = False,
            colors = colors, axis_label = l,  title = l + " : No Reweighting", num_bins = n_bins_, normalize = False, ratio_range = (0.5, 1.5), fname = outdir + l + '_ratio_before.png' )

weights_rw = copy.deepcopy(weights_nom)
uncs_rw = copy.deepcopy(uncs_nom)


h_ratio = f_ratio.Get("ratio_nom")

#Noise used to generated smeared ratio's based on stat unc
nToys = 100
rand_noise = np.random.normal(size = (h_ratio.GetNbinsX(), h_ratio.GetNbinsY(), h_ratio.GetNbinsZ(), nToys))

h_ratio = f_ratio.Get("ratio_nom")

d = sigs[0]

print("Reweighting ", d.f)
sig_idx = len(bkgs)
d_LP_weights, d_LP_uncs, d_LP_smeared_weights = d.reweight_LP(h_ratio, subjet_rw = subjet_rw, fill_z = fill_z, jetR = jetR, num_excjets = num_excjets, uncs = True, prefix = prefix, 
        rand_noise = rand_noise)
LP_weights = d_LP_weights
LP_uncs = d_LP_uncs/d_LP_weights

#cap at 100% uncertainty
d_LP_uncs = np.minimum(d_LP_uncs, d_LP_weights)

#apply weights, keep normalization fixed
old_norm = np.sum(weights_rw[sig_idx])
weights_rw[sig_idx] *= d_LP_weights
new_norm = np.sum(weights_rw[sig_idx])
weights_rw[sig_idx] *= old_norm / new_norm
uncs_rw [sig_idx] = weights_nom[sig_idx] * d_LP_uncs * (old_norm / new_norm)
LP_smeared_weights = np.array(d_LP_smeared_weights * np.expand_dims(weights_nom[sig_idx], -1) * (old_norm / new_norm))


sys_ratios = []
sys_variations = dict()
if(do_sys_variations):
    sys_list = list(sys_weights_map.keys())
    for sys in sys_list:
        if(sys == 'nom_weight'): continue
        sys_ratio = f_ratio.Get("ratio_" + sys)
        sys_ratio.Print()

        #limit to 1 signal for now...
        d = sigs[0]
        sys_LP_weights, _ = d.reweight_LP(sys_ratio, subjet_rw = subjet_rw, fill_z = fill_z, jetR = jetR, num_excjets = num_excjets, uncs = False, prefix = prefix)
        sys_weights = weights_nom[sig_idx] * sys_LP_weights
        rw = np.sum(weights_nom[sig_idx]) / np.sum(sys_weights)
        sys_weights *= rw
        sys_variations[sys] = sys_weights


make_histogram(LP_weights[0], "Reweighting factors", 'b', 'Weight', "Lund Plane Reweighting Factors", 20 , h_range = (0., 2.0),
     normalize=False, fname=outdir + "lundPlane_weights.png")

make_histogram(LP_uncs[0], "Fractional Uncertainties", 'b', 'Weight Fractional Uncertainty ', "Lund Plane Reweighting Factors Uncertainty", 20,
     normalize=False, fname=outdir + "lundPlane_weights_unc.png", h_range = (0., 1.5))


#compute 'Scalefactor'
cut_tau21 = d_ttbar_w_match.tau21 < tau21_cut
cut_DeepAK8 = d_ttbar_w_match.DeepAK8_W_MD > DeepAK8_cut

cuts = [cut_tau21, cut_DeepAK8]
cut_names = ["Tau21", "DeepAK8 W MD"]
cut_vals = [tau21_cut, DeepAK8_cut]

for idx,cut in enumerate(cuts):

    print(cut.shape, weights_nom[sig_idx].shape, LP_smeared_weights[:,idx].shape)
    eff_nom = np.average(cut, weights = weights_nom[sig_idx])
    eff_rw = np.average(cut, weights = weights_rw[sig_idx])
    print("Running %s" % cut_names[idx])

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
    if(do_sys_variations):
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

    print(i, len(cut_names), len(cut_vals))
    print("\n\nSF %s (cut val %.3f ) is %.2f +/- %.2f  (stat) + %.2f - %.2f (sys) \n\n"  % (cut_names[idx], cut_vals[idx], SF, SF_stat_unc, SF_sys_unc_up, SF_sys_unc_down))

    #approximate uncertainty on the reweighting for the plots
    overall_unc = (SF_stat_unc **2 + (0.5 * SF_sys_unc_up + 0.5 * SF_sys_unc_down)**2) **0.5 / SF
    print("overall unc %.3f" % overall_unc)

    uncs_rw[len(bkgs)] = overall_unc * weights_rw[len(bkgs)]


f_ratio.Close()


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
        h_range = (pt_cut, 800.)
        n_bins_ = n_bins
    else: 
        n_bins_ = n_bins
        h_range = None
    make_multi_sum_ratio_histogram(data = getattr(d_data, l), entries = a, weights = weights_rw, labels = labels, h_range = h_range, drawSys = False, stack = False,
            colors = colors, axis_label = l,  title = l + " : LP Reweighting", num_bins = n_bins_, normalize = False, ratio_range = (0.5, 1.5), fname = outdir + l + '_ratio_after.png' )



