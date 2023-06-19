from Utils import *
import os

parser = input_options()
options = parser.parse_args()

print(options)

ROOT.TGaxis.SetMaxDigits(3);



#NON UL (2018B only)
#lumi = 6.90
#f_data = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/CASEUtils/H5_maker/ttbar_output_files_v2/SingleMu_2018C.h5", "r")
#f_ttbar = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/CASEUtils/H5_maker/ttbar_output_files_v2/TTToSemiLep_2018.h5", "r")
#f_bkg = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/CASEUtils/H5_maker/ttbar_output_files_v2/QCD_WJets_merged.h5", "r")
#
#UL
lumi = 59.74
year = 2018
f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_may22/"

f_data = h5py.File(f_dir + "SingleMu_2018_merge.h5", "r")
f_ttbar = h5py.File(f_dir + "TT.h5", "r")
f_wjets = h5py.File(f_dir + "QCD_WJets.h5", "r")
f_diboson = h5py.File(f_dir + "diboson.h5", "r")
f_tw = h5py.File(f_dir + "TW.h5", "r")
f_singletop = h5py.File(f_dir + "SingleTop_merge.h5", "r")




#f_ratio = ROOT.TFile.Open("ttbar_UL_feb14_W_rw/ratio.root")
#f_ratio = ROOT.TFile.Open("ttbar_UL_feb14_W_rw_chg_only/ratio.root")
f_ratio = ROOT.TFile.Open(options.fin)
outdir = options.outdir
prefix = ""
drawSys = False


ratio_range = [0.2, 1.8]

#tau21_cut = 0.2700 # med
tau21_cut_val = 0.3452 # loose
DeepAK8_MD_cut_val = 0.479
DeepAK8_cut_val = 0.762

do_sys_variations = not (options.no_sys)

norm = True

#jms_corr = 0.95

jms_corr = 1.0

#m_cut_min = 60.
#m_cut_max = 110.
m_cut_min = 50.
m_cut_max = 250.
#m_cut_max = 65.
pt_cut = 225.

if(not os.path.exists(outdir)): os.system("mkdir " + outdir)

d_data = Dataset(f_data, is_data = True)

d_tw = Dataset(f_tw, label = "tW", color = ROOT.kYellow-7, jms_corr = jms_corr, dtype = -1)
d_wjets = Dataset(f_wjets, label = "W+Jets + QCD", color = ROOT.kOrange-3, jms_corr = jms_corr)
d_diboson = Dataset(f_diboson, label = "Diboson", color = ROOT.kCyan, jms_corr = jms_corr)
d_singletop = Dataset(f_singletop, label = "Single Top", color = ROOT.kMagenta-1, jms_corr = jms_corr)


d_ttbar_w_match = Dataset(f_ttbar, label = "t#bar{t} : W-matched", color = ROOT.kRed-7, jms_corr =jms_corr, dtype = 2)
d_ttbar_t_match = Dataset(f_ttbar, label = "t#bar{t} : t-matched ", color = ROOT.kBlue-7, jms_corr = jms_corr, dtype = 3)
d_ttbar_nomatch = Dataset(f_ttbar, label = "t#bar{t} : unmatched", color = ROOT.kGreen-6, jms_corr = jms_corr)


ttbar_gen_matching = d_ttbar_w_match.f['gen_parts'][:,0]

#0 is unmatched, 1 is W matched, 2 is top matched
nomatch_cut = ttbar_gen_matching < 0.1
w_match_cut = (ttbar_gen_matching  > 0.9) &  (ttbar_gen_matching < 1.1)
t_match_cut = (ttbar_gen_matching  > 1.9) &  (ttbar_gen_matching < 2.1)

d_ttbar_w_match.apply_cut(w_match_cut)
d_ttbar_t_match.apply_cut(t_match_cut)
d_ttbar_nomatch.apply_cut(nomatch_cut)


sigs = [d_ttbar_w_match]
bkgs = [d_diboson, d_singletop, d_wjets,  d_ttbar_t_match, d_ttbar_nomatch, d_tw]
tw_idx = len(bkgs)-1





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


jetR = 1.0

num_excjets = -1

n_pt_bins = 6
num_excjets = 2
pt_bins = array('f', [0., 50., 100., 200., 300., 450., 99999.])




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



obs = ["tau21", "tau32", "tau43", "nPF", "mSoftDrop", "pt", "DeepAK8_W_MD", "DeepAK8_W"]

obs_attrs = {
        'mSoftDrop' : (m_cut_min, m_cut_max, n_bins, "m_{SD}"),
        'tau21' : (0.05, 0.8, 20, "#tau_{21}"),
        'tau32' : (0.4, 1.0, 20, "#tau_{32}"),
        'tau43' : (0.6, 1.0, 20, "#tau_{43}"),
        'nPF' : (20.5, 80.5, 30, "Num. PF Cands."),
        'pt' : (pt_cut, 1200., 20, "p_{T}"),
        'DeepAK8_W' : (0., 1., 20, "DeepAK8 (W vs. QCD)"),
        'DeepAK8_W_MD' : (0., 1., 20, "DeepAK8-MD (W vs. QCD)"),
        }


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
    a_data = getattr(d_data, l)

    low,high, nbins_, label = obs_attrs.get(l, (None, None, 20, l))

    make_multi_sum_ratio_histogram(data = a_data, entries = a, weights = weights_nom, labels = labels, h_range = (low, high), drawSys = False, stack = True, draw_chi2 = True,
            year = year, colors = colors, axis_label = label,  title = l + " : No Reweighting", 
            num_bins = nbins_, normalize = False, ratio_range = ratio_range, fname = outdir + l + '_ratio_before.png' )

weights_rw = copy.deepcopy(weights_nom)

#10% normalization unc on bkgs
uncs = [0.1] * len(bkgs + sigs)


h_ratio = f_ratio.Get("ratio_nom")

if('pt_extrap' in f_ratio.GetListOfKeys()):
    rdir = f_ratio.GetDirectory("pt_extrap")
    rdir.cd()

else: 
    print("NO Pt extrapolation")
    rdir = None

#Noise used to generated smeared ratio's based on stat unc
nToys = 100
rand_noise = np.random.normal(size = (nToys, h_ratio.GetNbinsX(), h_ratio.GetNbinsY(), h_ratio.GetNbinsZ()))
pt_rand_noise = np.random.normal(size = (nToys, h_ratio.GetNbinsY(), h_ratio.GetNbinsZ(), 3))

h_ratio = f_ratio.Get("ratio_nom")

LP_rw = LundReweighter(jetR = jetR, pt_extrap_dir = rdir, charge_only = options.charge_only)

d_tw_LP_weights = d_tw.reweight_LP(LP_rw, h_ratio, num_excjets = num_excjets,  prefix = "",)

d_sig = sigs[0]
print("Reweighting ", d.f)
sig_idx = len(bkgs)

subjets, splittings, bad_match, deltaRs = d_sig.get_matched_splittings(LP_rw, num_excjets = num_excjets, return_dRs = True)

subjet_responses = []
jet_responses = []

gen_parts_raw = d_sig.get_masked('gen_parts')[:]
top = gen_parts_raw[:,1:5]
antitop = gen_parts_raw[:,5:9]
W = gen_parts_raw[:,9:13]
antiW = gen_parts_raw[:,13:17]
q1 = gen_parts_raw[:,17:20]
q2 = gen_parts_raw[:,21:24]
b = gen_parts_raw[:,25:28]
gen_parts = np.stack([q1, q2], axis = 1)

j_4vec = d_sig.get_masked('jet_kinematics')[:,:4].astype(np.float64)

for i,sjs in enumerate(subjets):
    if(bad_match[i]): continue

    if(deltaR(W[i], j_4vec[i]) < deltaR(antiW[i], j_4vec[i])):
        jet_responses.append(j_4vec[i][0] / W[i][0])
    else:
        jet_responses.append(j_4vec[i][0] / antiW[i][0])

    subjet_responses.append(d_sig.get_pt_response(gen_parts[i], subjets[i]))

make_histogram(np.array(subjet_responses).reshape(-1), "W subjets", 'b', 'Subjet pt / gen pt', "Subjet pt response ", 20 , h_range = (0.5, 1.5),
     normalize=True, fname=outdir + "subjet_response.png", mean_std = True)


make_histogram(np.array(jet_responses).reshape(-1), "W jets", 'b', 'Jet pt / gen pt', "Jet pt response ", 20 , h_range = (0.5, 1.5),
     normalize=True, fname=outdir + "jet_response.png", mean_std = True)


d_LP_weights, d_LP_smeared_weights, d_pt_smeared_weights = d_sig.reweight_LP(LP_rw, h_ratio, num_excjets = num_excjets, prefix = "", 
        rand_noise = rand_noise, pt_rand_noise = pt_rand_noise, subjets = subjets, splittings = splittings)


#apply weights, keep normalization fixed
old_norm = np.sum(weights_rw[sig_idx])
weights_rw[sig_idx] *= d_LP_weights
new_norm = np.sum(weights_rw[sig_idx])
weights_rw[sig_idx] *= old_norm / new_norm


old_norm = np.sum(weights_rw[tw_idx])
weights_rw[tw_idx] *= d_tw_LP_weights
new_norm = np.sum(weights_rw[tw_idx])
weights_rw[tw_idx] *= old_norm / new_norm


LP_smeared_weights = np.array(d_LP_smeared_weights * np.expand_dims(weights_nom[sig_idx], -1) * (old_norm / new_norm))
pt_smeared_weights = np.array(d_pt_smeared_weights * np.expand_dims(weights_nom[sig_idx], -1) * (old_norm / new_norm))


sys_ratios = []
sys_variations = dict()
if(do_sys_variations):
    #sys_list = list(sys_weights_map.keys())
    sys_list = ['sys_tot_up', 'sys_tot_down']
    for sys in sys_list:
        if(sys == 'nom_weight'): continue
        sys_ratio = f_ratio.Get("ratio_" + sys)
        sys_ratio.Print()
        sys_str = sys + "_"

        sys_LP_weights = d_sig.reweight_LP(LP_rw, sys_ratio, num_excjets = num_excjets, prefix = "", 
                sys_str = sys_str, subjets = subjets, splittings = splittings)
        sys_weights = weights_nom[sig_idx] * sys_LP_weights
        rw = np.sum(weights_nom[sig_idx]) / np.sum(sys_weights)
        sys_weights *= rw
        sys_variations[sys] = sys_weights


clip_weights = np.clip(d_LP_weights, 0., 5.)
make_histogram(clip_weights, "Reweighting factors", 'b', 'Weight', "Lund Plane Reweighting Factors", 20 , h_range = (0., 5.0),
     normalize=False, fname=outdir + "lundPlane_weights.png")

#Save subjet pts and deltaR
subjet_pts =  []
deltaRs = np.reshape(deltaRs, -1)

for i,sjs in enumerate(subjets):
    for sj in sjs:
        subjet_pts.append(sj[0])
    
num_bins = 40
pt_bins = array('d', np.linspace(0., 800., num_bins + 1))
dR_bins = array('d', np.linspace(0., 0.8, num_bins + 1))

h_subjet_pts = make_root_hist(data = subjet_pts, name = 'h_W_subjetpt', num_bins = num_bins, bins = pt_bins)
h_dRs = make_root_hist(data = deltaRs, name = 'h_W_dRs', num_bins = num_bins, bins = dR_bins)
f_ptout = ROOT.TFile.Open(outdir + "subjet_pt_dR.root", "RECREATE")
h_subjet_pts.Write()
h_dRs.Write()
f_ptout.Close()


#compute 'Scalefactor'
cut_tau21 = d_ttbar_w_match.tau21 < tau21_cut_val
cut_DeepAK8_MD = d_ttbar_w_match.DeepAK8_W_MD > DeepAK8_MD_cut_val
cut_DeepAK8 = d_ttbar_w_match.DeepAK8_W > DeepAK8_cut_val

cuts = [cut_tau21, cut_DeepAK8_MD, cut_DeepAK8]
cut_names = ["Tau21", "DeepAK8 W MD", "DeepAK8 W"]
cut_vals = [tau21_cut_val, DeepAK8_MD_cut_val, DeepAK8_cut_val]

for idx,cut in enumerate(cuts):

    eff_nom = np.average(cut, weights = weights_nom[sig_idx])
    eff_rw = np.average(cut, weights = weights_rw[sig_idx])
    print("Running %s" % cut_names[idx])

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
    SF_sys_unc = 0.
    if(do_sys_variations):

        eff_sys_tot_up = np.average(cut, weights = sys_variations['sys_tot_up'])
        eff_sys_tot_down = np.average(cut, weights = sys_variations['sys_tot_down'])
        SF_sys_unc_up = (eff_sys_tot_up - eff_rw)/eff_nom
        SF_sys_unc_down = (eff_sys_tot_down - eff_rw)/eff_nom
        SF_sys_unc = (SF_sys_unc_up + SF_sys_unc_down) / 2.0

        #for sys in sys_variations.keys():
        #    eff = np.average(cut, weights = sys_variations[sys])
        #    diff = eff - eff_rw
        #    if(diff > 0): sys_unc_up += diff**2
        #    else: sys_unc_down += diff**2
        #    print("%s %.4f" % (sys,  diff))

        #sys_unc_up = sys_unc_up**(0.5)
        #sys_unc_down = sys_unc_down**(0.5)


    SF = eff_rw / eff_nom
    SF_stat_unc = abs(toys_mean - eff_rw)/eff_nom + toys_std /eff_nom
    SF_pt_unc = abs(pt_toys_mean - eff_rw)/eff_nom + pt_toys_std /eff_nom


    bad_matching_unc = np.mean(bad_match) * SF

    print("\n\nSF %s (cut val %.3f ) is %.2f +/- %.2f  (stat) +/- %.2f (sys) +/- %.2f (pt) +/- %.2f (matching) \n\n"  
            % (cut_names[idx], cut_vals[idx], SF, SF_stat_unc, SF_sys_unc, SF_pt_unc, bad_matching_unc))

    #approximate uncertainty on the reweighting for the plots
    overall_unc = (SF_stat_unc **2 + SF_sys_unc**2 + SF_pt_unc**2 + bad_matching_unc**2) **0.5 / SF
    print("overall unc %.3f" % overall_unc)

    uncs[len(bkgs)] = overall_unc

    #apply same unc to tW
    uncs[len(bkgs)-1] = overall_unc


f_ratio.Close()


for l in obs:
    a = []
    for d in (bkgs + sigs):
        a.append(getattr(d, l))
    a_data = getattr(d_data, l)

    low,high, nbins_, label = obs_attrs.get(l, (None, None, 20, l))
    make_multi_sum_ratio_histogram(data = a_data, entries = a, weights = weights_rw, labels = labels, h_range = (low, high), stack = True, uncs = uncs, drawSys = drawSys, draw_chi2 = True,
            year = year, colors = colors, axis_label = label,  title = l + " : LP Reweighting", num_bins = nbins_, 
            normalize = False, ratio_range = ratio_range, fname = outdir + l + '_ratio_after.png' )



