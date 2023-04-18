from Utils import *
import os

parser = input_options()
options = parser.parse_args()

print(options)

#UL
lumi = 59.74
f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_jan31/"

f_data = h5py.File(f_dir + "SingleMu_2018_merge.h5", "r")
f_ttbar = h5py.File(f_dir + "TT.h5", "r")
f_wjets = h5py.File(f_dir + "QCD_WJets.h5", "r")
f_diboson = h5py.File(f_dir + "diboson.h5", "r")
f_tw = h5py.File(f_dir + "TW.h5", "r")
f_singletop = h5py.File(f_dir + "SingleTop_merge.h5", "r")



#f_ratio = ROOT.TFile.Open("ttbar_UL_jan20_W_rw_kt_sys/ratio.root")
f_ratio = ROOT.TFile.Open(options.fin)


#for SF computation
tau32_cut = 0.52



outdir = options.outdir
if(not os.path.exists(outdir)): os.system("mkdir %s" % outdir)
do_sys_variations = True

max_evts = None

norm = True

jms_corr = 0.95

m_cut_min = 125.
#m_cut_max = 130.
m_cut_max = 225.
pt_cut = 500.

if(not os.path.exists(outdir)): os.system("mkdir " + outdir)

d_data = Dataset(f_data, is_data = True)

d_tw = Dataset(f_tw, label = "tW", color = ROOT.kMagenta, jms_corr = jms_corr)
d_wjets = Dataset(f_wjets, label = "W+Jets + QCD", color = ROOT.kGray, jms_corr = jms_corr)
d_diboson = Dataset(f_diboson, label = "Diboson", color = ROOT.kCyan, jms_corr = jms_corr)
d_singletop = Dataset(f_singletop, label = "Single Top", color = ROOT.kMagenta+4, jms_corr = jms_corr)


d_ttbar_w_match = Dataset(f_ttbar, label = "ttbar : W-matched", color = ROOT.kRed, jms_corr =jms_corr, dtype = 2)
d_ttbar_t_match = Dataset(f_ttbar, label = "ttbar : t-matched ", color = ROOT.kOrange-3, jms_corr = jms_corr, dtype = 3)
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

#jetR = 0.4

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
        h_range = (pt_cut, 1200.)
        n_bins_ = n_bins
    else: 
        n_bins_ = n_bins
        h_range = None
    make_multi_sum_ratio_histogram(data = getattr(d_data, l), entries = a, weights = weights_nom, labels = labels, h_range = h_range, drawSys = False, stack = False, draw_chi2 = True,
            colors = colors, axis_label = l,  title = l + " : No Reweighting", num_bins = n_bins, normalize = False, ratio_range = (0.4, 1.6), fname = outdir + l + '_ratio_before.png' )



weights_rw = copy.deepcopy(weights_nom)
uncs_rw = copy.deepcopy(uncs_nom)

h_ratio = f_ratio.Get("ratio_nom")
f_ratio.cd('pt_extrap')
rdir = ROOT.gDirectory

nToys = 100

#Noise used to generated smeared ratio's based on stat unc
rand_noise = np.random.normal(size = (nToys, h_ratio.GetNbinsX(), h_ratio.GetNbinsY(), h_ratio.GetNbinsZ()))
pt_rand_noise = np.random.normal(size = (nToys, h_ratio.GetNbinsY(), h_ratio.GetNbinsZ(), 3))

LP_rw = LundReweighter(jetR = jetR, pt_extrap_dir = rdir, charge_only = options.charge_only)


d_sig = sigs[0]

sig_idx = len(bkgs)
print("Reweighting ", d.f )
subjets, splittings, bad_match = d_sig.get_matched_splittings(LP_rw, num_excjets = num_excjets)
d_LP_weights, _, d_LP_smeared_weights, d_pt_smeared_weights = d_sig.reweight_LP(LP_rw, h_ratio, num_excjets = num_excjets, uncs = False, prefix = "", 
        rand_noise = rand_noise, pt_rand_noise = pt_rand_noise, subjets = subjets, splittings = splittings)

LP_weights = d_LP_weights

#apply weights, keep normalization fixed
old_norm = np.sum(weights_rw[sig_idx])
weights_rw[sig_idx] *= d_LP_weights

new_norm = np.sum(weights_rw[sig_idx])


weights_rw[sig_idx] *= old_norm / new_norm
LP_smeared_weights = np.array(d_LP_smeared_weights * np.expand_dims(weights_nom[sig_idx], -1) * (old_norm / new_norm))
pt_smeared_weights = np.array(d_pt_smeared_weights * np.expand_dims(weights_nom[sig_idx], -1) * (old_norm / new_norm))


sys_variations = dict()
if(do_sys_variations):
    #sys_list = list(sys_weights_map.keys())
    sys_list = ['sys_tot_up', 'sys_tot_down']
    for sys in sys_list:
        if(sys == 'nom_weight'): continue
        sys_ratio = f_ratio.Get("ratio_" + sys)
        sys_ratio.Print()
        sys_str = sys + "_"

        sys_LP_weights, _ = d_sig.reweight_LP(LP_rw, sys_ratio, num_excjets = num_excjets, uncs = False, prefix = "", 
                sys_str = sys_str, subjets = subjets, splittings = splittings)

        sys_weights = weights_nom[sig_idx] * sys_LP_weights
        rw = np.sum(weights_nom[sig_idx]) / np.sum(sys_weights)
        sys_weights *= rw
        sys_variations[sys] = sys_weights

    b_light_ratio = f_ratio.Get("h_bl_ratio")
    bquark_rw, _ = d.reweight_LP(LP_rw, b_light_ratio, num_excjets = num_excjets, uncs = False, prefix = "", 
            max_evts = max_evts, sys_str = 'bquark', subjets = subjets, splittings = splittings)

    up_bquark_weights = bquark_rw * weights_rw[sig_idx]
    down_bquark_weights = (1./ bquark_rw) * weights_rw[sig_idx]


    up_bquark_weights *= old_norm / np.sum(up_bquark_weights)
    down_bquark_weights *= old_norm / np.sum(down_bquark_weights)

    sys_variations['bquark_up'] = up_bquark_weights
    sys_variations['bquark_down'] = down_bquark_weights




make_histogram(LP_weights, "Reweighting factors", 'b', 'Weight', "Lund Plane Reweighting Factors", 20 , h_range = (0., 2.0),
     normalize=False, fname=outdir + "lundPlane_weights.png")


#compute 'Scalefactor'
cut = d_ttbar_t_match.tau32 < tau32_cut

eff_nom = np.average(cut, weights = weights_nom[sig_idx])
eff_rw = np.average(cut, weights = weights_rw[sig_idx])

print("Nom %.3f, RW %.3f" % (eff_nom, eff_rw))


eff_toys = []
pt_eff_toys = []
for i in range(nToys):
    eff = np.average(cut, weights = LP_smeared_weights[:,i])
    eff_toys.append(eff)

    eff1 = np.average(cut, weights = pt_smeared_weights[:,i])
    pt_eff_toys.append(eff1)

toys_mean = np.mean(eff_toys)
toys_std = np.std(eff_toys)

print("Toys avg %.3f, std dev %.3f" % (toys_mean, toys_std))

pt_toys_mean = np.mean(pt_eff_toys)
pt_toys_std = np.std(pt_eff_toys)

print("Pt variation toys avg %.3f, std dev %.3f" % (pt_toys_mean, pt_toys_std))

#Add systematic differences in quadrature
sys_unc_up = sys_unc_down = 0.
if(do_sys_variations):

    eff_sys_tot_up = np.average(cut, weights = sys_variations['sys_tot_up'])
    eff_sys_tot_down = np.average(cut, weights = sys_variations['sys_tot_down'])
    SF_sys_unc_up = abs(eff_sys_tot_up - eff_rw)/eff_nom
    SF_sys_unc_down = abs(eff_sys_tot_down - eff_rw)/eff_nom
    SF_sys_unc = (SF_sys_unc_up + SF_sys_unc_down) / 2.0

    eff_bquark_up = np.average(cut, weights = sys_variations['bquark_up'])
    eff_bquark_down = np.average(cut, weights = sys_variations['bquark_down'])
    SF_bquark_up = abs(eff_bquark_up - eff_rw)/eff_nom
    SF_bquark_down = abs(eff_bquark_down - eff_rw)/eff_nom
    SF_bquark_unc = (SF_bquark_up + SF_bquark_down) /2.0


SF = eff_rw / eff_nom
SF_stat_unc = abs(toys_mean - eff_rw)/eff_nom + toys_std /eff_nom
SF_pt_unc = abs(pt_toys_mean - eff_rw)/eff_nom + pt_toys_std /eff_nom

bad_matching_unc = np.mean(bad_match) * SF

print("\n\nSF (cut val %.2f ) is %.2f +/- %.2f  (stat) +/- %.2f (sys) +/- %.2f (pt) +/- %.2f (bquark) +/- %.2f (matching) \n\n"  
        % (tau32_cut, SF, SF_stat_unc, SF_sys_unc, SF_pt_unc, SF_bquark_unc, bad_matching_unc))

f_ratio.Close()

#approximate uncertainty on the reweighting for the plots
overall_unc = (SF_stat_unc **2 + SF_sys_unc**2 + SF_pt_unc**2 + SF_bquark_unc**2 + bad_matching_unc**2) **0.5 / SF
print("overall unc %.3f" % overall_unc)

uncs_rw[len(bkgs)] = overall_unc * weights_rw[len(bkgs)]



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
    make_multi_sum_ratio_histogram(data = getattr(d_data, l), entries = a, weights = weights_rw, labels = labels, uncs = uncs_rw, h_range = h_range, drawSys = False, stack = False, draw_chi2 = True,
            colors = colors, axis_label = l,  title = l + " : LP Reweighting", num_bins = n_bins, normalize = False, ratio_range = (0.4, 1.6), fname = outdir + l + '_ratio_after.png' )



