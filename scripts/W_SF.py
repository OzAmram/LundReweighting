import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.Utils import *

parser = input_options()
options = parser.parse_args()

print(options)

ROOT.TGaxis.SetMaxDigits(3);



#UL
if(options.year == 2018):
    lumi = 59.74
    year = 2018
    f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_2018/"
    f_ratio_name = "data/ratio_2018.root"

elif(options.year == 2017):
    lumi = 41.42
    year = 2017
    f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_2017/"
    f_ratio_name = "data/ratio_2017.root"

elif(options.year == 2016):
    year = 2016
    lumi = 16.8 + 19.5
    f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_2016/"
    f_ratio_name = "data/ratio_2016.root"
else:
    exit(1)


f_data = h5py.File(f_dir + "SingleMu_merge.h5", "r")
f_ttbar = h5py.File(f_dir + "TT.h5", "r")
f_wjets = h5py.File(f_dir + "QCD_WJets.h5", "r")
f_diboson = h5py.File(f_dir + "diboson.h5", "r")
f_tw = h5py.File(f_dir + "TW.h5", "r")
f_singletop = h5py.File(f_dir + "SingleTop_merge.h5", "r")




if(options.fin != ""): f_ratio_name = options.fin
f_ratio = ROOT.TFile.Open(f_ratio_name)
outdir = options.outdir
prefix = ""
drawSys = False
if(not os.path.exists(outdir)): os.system("mkdir " + outdir)


ratio_range = [0.2, 1.8]

#tau21_cut = 0.2700 # med
tau21_cut_val = 0.3452 # loose
DeepAK8_MD_cut_val = 0.704
DeepAK8_cut_val = 0.918
ParticleNet_cut_val = 0.94 # 1\% mistag

do_sys_variations = not (options.no_sys)

norm = True

#jms_corr = 0.95

jms_corr = 1.0

#m_cut_min = 80.
#m_cut_max = 81.
m_cut_min = 70.
m_cut_max = 110.
pt_cut = 225.


d_data = Dataset(f_data, is_data = True)

d_tw = Dataset(f_tw, label = "tW", color = ROOT.kYellow-7, jms_corr = jms_corr, dtype = -1)
d_wjets = Dataset(f_wjets, label = "W+Jets + QCD", color = ROOT.kOrange-3, jms_corr = jms_corr)
d_diboson = Dataset(f_diboson, label = "Diboson", color = ROOT.kCyan, jms_corr = jms_corr)
d_singletop = Dataset(f_singletop, label = "Single Top", color = ROOT.kMagenta-1, jms_corr = jms_corr)


d_ttbar_w_match = Dataset(f_ttbar, label = "t#bar{t} : W-matched", color = ROOT.kRed-7, jms_corr =jms_corr, dtype = 2)
d_ttbar_t_match = Dataset(f_ttbar, label = "t#bar{t} : t-matched ", color = ROOT.kBlue-7, jms_corr = jms_corr, dtype = 3)
d_ttbar_nomatch = Dataset(f_ttbar, label = "t#bar{t} : unmatched", color = ROOT.kGreen-6, jms_corr = jms_corr)

d_tw_w_match = Dataset(f_tw, label = "tW : W-matched", color = ROOT.kMagenta, jms_corr = jms_corr)
d_tw_nomatch = Dataset(f_tw, label = "tW : unmatched", color = ROOT.kMagenta+3, jms_corr = jms_corr)


ttbar_gen_matching = d_ttbar_w_match.f['gen_parts'][:,0]
tW_gen_matching = d_tw_w_match.f['gen_parts'][:,0]

#0 is unmatched, 1 is W matched, 2 is top matched
nomatch_cut = ttbar_gen_matching < 0.1
w_match_cut = (ttbar_gen_matching  > 0.9) &  (ttbar_gen_matching < 1.1)
t_match_cut = (ttbar_gen_matching  > 1.9) &  (ttbar_gen_matching < 2.1)

d_ttbar_w_match.apply_cut(w_match_cut)
d_ttbar_t_match.apply_cut(t_match_cut)
d_ttbar_nomatch.apply_cut(nomatch_cut)

tW_w_match_cut = (tW_gen_matching  > 0.9) &  (tW_gen_matching < 1.1)
d_tw_w_match.apply_cut(tW_w_match_cut)
d_tw_nomatch.apply_cut(~tW_w_match_cut)


sigs = [d_ttbar_w_match, d_tw_w_match]
bkgs = [d_ttbar_nomatch, d_ttbar_t_match, d_tw_nomatch, d_diboson, d_wjets, d_singletop]

tw_idx = len(bkgs)-1


pt_max = 1000

num_excjets = -1




jet_kinematics_data= d_data.get_masked('jet_kinematics')
msd_cut_data = (jet_kinematics_data[:,3] > m_cut_min) & (jet_kinematics_data[:,3] < m_cut_max)
pt_cut_data = jet_kinematics_data[:,0] > pt_cut
d_data.compute_kinematics()
mu_b_dr_cut = d_data.dR_mu_bjet > dR_mu_bjet_cut
d_data.apply_cut(msd_cut_data & pt_cut_data & mu_b_dr_cut)
d_data.compute_obs()

for d in (bkgs + sigs):

    d.norm_factor = lumi
    d.compute_kinematics()

    jet_kinematics = d.get_masked('jet_kinematics')
    msd_cut_mask = (jet_kinematics[:,3] * jms_corr > m_cut_min) & (jet_kinematics[:,3] * jms_corr < m_cut_max)
    pt_cut_mask = jet_kinematics[:,0] > pt_cut
    mu_b_dr_cut = d.dR_mu_bjet > dR_mu_bjet_cut
    d.apply_cut(msd_cut_mask & pt_cut_mask & mu_b_dr_cut)
    d.compute_obs()



#print("Num data %i. Num ttbar MC %i " % (d_data.n(), d_ttbar.n()))


num_data = np.sum(d_data.get_weights())
num_ttbar_nomatch = np.sum(d_ttbar_nomatch.get_weights())
num_ttbar_w_match = np.sum(d_ttbar_w_match.get_weights())
num_ttbar_t_match = np.sum(d_ttbar_t_match.get_weights())
num_ttbar_tot = num_ttbar_nomatch + num_ttbar_w_match + num_ttbar_t_match
num_tw_nomatch = np.sum(d_tw_nomatch.get_weights())
num_tw_w_match = np.sum(d_tw_w_match.get_weights())
num_tw = num_tw_nomatch + num_tw_w_match

tot_bkg = 0.
for d in (d_diboson, d_wjets, d_singletop):
    tot_bkg += np.sum(d.get_weights())
print("%i data, %.0f ttbar (%.0f unmatched, %.0f W matched, %.0f t matched), %.0f tW (%0.f unmatched %.0f W matched) %.0f bkg" % ( num_data, 
            num_ttbar_tot,num_ttbar_nomatch, num_ttbar_w_match, num_ttbar_t_match, 
            num_tw, num_tw_nomatch, num_tw_w_match, tot_bkg))
normalization = num_data  / (num_ttbar_tot + num_tw + tot_bkg)
print("normalization", normalization)

if(norm):
    for d in (bkgs + sigs):
        d.norm_factor *= normalization



obs = ["tau21", "tau32", "tau43", "nPF", "mSoftDrop", "pt", "DeepAK8_W_MD", "DeepAK8_W", "ParticleNet_W"]

obs_attrs = {
        'mSoftDrop' : (70, 110, 20, "m_{SD} [GeV] ", "Events / 2 GeV") if m_cut_max < 200 else (50, 230, 45, "m_{SD} [GeV]", "Events / 4 GeV"),
        'tau21' : (0.05, 0.8, 25, "#tau_{21}", "Events / 0.03"),
        'tau32' : (0.4, 0.95, 25, "#tau_{32}", "Events / 0.022)"),
        'tau43' : (0.6, 0.96, 18, "#tau_{43}", "Events / 0.02"),
        'nPF' : (20.5, 80.5, 30, "Num. PF Cands.", "Events / 2" ),
        'pt' : (225, 825., 20, "p_{T}", "Events / 30 GeV"),
        'DeepAK8_W' : (0., 1., 20, "DeepAK8 (W vs. QCD)", "Events / 0.05"),
        'DeepAK8_W_MD' : (0., 1., 20, "DeepAK8-MD (W vs. QCD)", "Events / 0.05"),
        'ParticleNet_W' : (0., 1., 20, "ParticleNet (W vs. QCD)", "Events / 0.05"),
        }


colors = []
weights_nom = []
labels = []
sig_idx = -2
for d in (bkgs + sigs):
    colors.append(d.color)
    weights_nom.append(d.get_weights())
    labels.append(d.label)

for l in obs:
    a = []
    for d in (bkgs + sigs):
        a.append(getattr(d, l))
    a_data = getattr(d_data, l)

    low,high, nbins_, label, ylabel = obs_attrs.get(l, (None, None, 20, l, ""))

    make_multi_sum_ratio_histogram(data = a_data, entries = a, weights = weights_nom, labels = labels, h_range = (low, high), drawSys = False, stack = True, draw_chi2 = True,
            year = year, colors = colors, axis_label = label,  title = l + " : No Reweighting",
            num_bins = nbins_, normalize = False, ratio_range = ratio_range, fname = outdir + l + '_ratio_before.png' )

weights_rw = copy.deepcopy(weights_nom)

#10% normalization unc on bkgs
uncs = [0.1] * len(bkgs + sigs)


LP_rw = LundReweighter(f_ratio = f_ratio, charge_only = options.charge_only)

tw_LP_weights = d_tw_w_match.reweight_all(LP_rw)

for key in tw_LP_weights.keys():
    if('nom' in key or 'up' in key or 'down' in key):
        if(isinstance(tw_LP_weights[key], np.ndarray)) : tw_LP_weights[key] *= d_tw_w_match.get_weights()
weights_rw[-1] = tw_LP_weights['nom']


LP_weights = d_ttbar_w_match.reweight_all(LP_rw)

make_histogram(LP_weights['nom'], "Reweighting factors", 'b', 'Weight', "Lund Plane Reweighting Factors", 20 , h_range = (0., 5.0),
     normalize=False, fname=outdir + "lundPlane_weights.png")

for key in LP_weights.keys():
    if('nom' in key or 'up' in key or 'down' in key):
        if(isinstance(LP_weights[key], np.ndarray)) : LP_weights[key] *= d_ttbar_w_match.get_weights()
weights_rw[-2] = LP_weights['nom']

for i in range(len(weights_nom)):
    print(i, np.sum(weights_nom[i]), np.sum(weights_rw[i]))

subjets = LP_weights['subjet_pts']
weights_nom2 = np.repeat(weights_nom[sig_idx], 2)
weights_rw2 = np.repeat(LP_weights['nom'], 2)
make_histogram([subjets, subjets], labels = ["Original", "Reweighted"], colors = ['black', 'blue'], 
        xaxis_label = 'Subjet pt (GeV)',  weights = [weights_nom2, weights_rw2], fname = options.outdir + "subjet_pts.png")

#Fraction of prongs that are not well matched to subjets (want this to be low)
print("Bad match frac %.2f" % np.mean(LP_weights['bad_match']))
#Fraction of prongs that are still not well matched after reclustering with varied number of prongs
print("Reclustered bad match frac %.2f" % np.mean(LP_weights['reclust_still_bad_match']))

#compute 'Scalefactor'
cut_tau21 = d_ttbar_w_match.tau21 < tau21_cut_val
cut_DeepAK8_MD = d_ttbar_w_match.DeepAK8_W_MD > DeepAK8_MD_cut_val
cut_DeepAK8 = d_ttbar_w_match.DeepAK8_W > DeepAK8_cut_val
cut_ParticleNet = d_ttbar_w_match.ParticleNet_W > ParticleNet_cut_val

cuts = [cut_tau21, cut_DeepAK8_MD, cut_DeepAK8, cut_ParticleNet]
cut_names = ["Tau21", "DeepAK8 W MD", "DeepAK8 W", "ParticleNet"]
cut_vals = [tau21_cut_val, DeepAK8_MD_cut_val, DeepAK8_cut_val, ParticleNet_cut_val]

f_SFs = open(options.outdir + "SFs.txt", "w")

for idx,cut in enumerate(cuts):

    eff_nom = np.average(cut, weights = weights_nom[sig_idx])
    eff_rw = np.average(cut, weights = LP_weights['nom'])
    print("Running %s" % cut_names[idx])

    print("Nom %.3f, RW %.3f" % (eff_nom, eff_rw))
    sf_nom = eff_rw/eff_nom

    nToys = LP_weights['stat_vars'].shape[1]
    eff_toys = []
    pt_eff_toys = []
    for i in range(nToys):
        eff = np.average(cut, weights = LP_weights['stat_vars'][:,i])
        eff_toys.append(eff)

        eff1 = np.average(cut, weights = LP_weights['pt_vars'][:,i])
        pt_eff_toys.append(eff1)

    #Compute stat and pt uncertainty based on variation in the toys
    toys_mean = np.mean(eff_toys)
    toys_std = np.std(eff_toys)
    pt_toys_mean = np.mean(pt_eff_toys)
    pt_toys_std = np.std(pt_eff_toys)

    eff_stat_unc = (abs(toys_mean - eff_rw)  + toys_std) /eff_nom
    eff_pt_unc = (abs(pt_toys_mean - eff_rw) + pt_toys_std) /eff_nom

    print("Stat variation toys eff. avg %.3f, std dev %.3f" % (toys_mean, toys_std))
    print("Pt variation toys eff. avg %.3f, std dev %.3f" % (pt_toys_mean, pt_toys_std))

    #Add systematic differences in quadrature
    sys_keys = ['sys', 'bquark', 'prongs', 'unclust']
    sys_uncs = dict()

    for sys in sys_keys: sys_uncs[sys] = [0.,0.]
    if(do_sys_variations):

        #Compute difference in efficiency due to weight variations as uncertainty
        def get_uncs(score_cut, weights_up, weights_down, eff_baseline):
            eff_up =  np.average(score_cut, weights = weights_up)
            eff_down =  np.average(score_cut, weights = weights_down)

            unc_up = eff_up - eff_baseline
            unc_down = eff_down - eff_baseline 
            return unc_up, unc_down

        for sys in sys_keys:
            unc_up, unc_down = get_uncs(cut, LP_weights[sys + '_up'], LP_weights[sys + '_down'], eff_rw)
            print(sys, unc_up, unc_down)
            sys_uncs[sys] = [unc_up/eff_nom, unc_down/eff_nom]

    cut_str = "Cut %s, %.2f \n" % (cut_names[idx], cut_vals[idx])



    SF_str = "SF is %.2f +/- %.2f (stat) +/- %.2f (pt)" % (sf_nom, eff_stat_unc, eff_pt_unc )
    tot_unc_up = tot_unc_down = eff_stat_unc**2 + eff_pt_unc**2

    for sys in sys_keys:
        SF_str += " %.2f/%.2f (%s)" % (sys_uncs[sys][0], sys_uncs[sys][1], sys)
        up_var = max(sys_uncs[sys][0], sys_uncs[sys][1])
        down_var = min(sys_uncs[sys][0], sys_uncs[sys][1])
        tot_unc_up += up_var**2
        tot_unc_down += down_var**2

    tot_unc_up = tot_unc_up**0.5
    tot_unc_down = tot_unc_down**0.5

    SF_str += "\n Overall %.2f +%.2f/-%.2f \n\n"  % (sf_nom, tot_unc_up, tot_unc_down)

    print(cut_str)
    print(SF_str)
    f_SFs.write(cut_str)
    f_SFs.write(SF_str)

    tot_unc = (abs(tot_unc_up) + (tot_unc_down))/2.0
    uncs[sig_idx] = tot_unc
    #apply same unc to tW
    uncs[-1] = tot_unc


f_SFs.close()
f_ratio.Close()


for l in obs:
    a = []
    for d in (bkgs + sigs):
        a.append(getattr(d, l))
    a_data = getattr(d_data, l)

    low,high, nbins_, label, ylabel = obs_attrs.get(l, (None, None, 20, l, ""))
    make_multi_sum_ratio_histogram(data = a_data, entries = a, weights = weights_rw, labels = labels, h_range = (low, high), stack = True, uncs = uncs, drawSys = drawSys, draw_chi2 = True,
            year = year, colors = colors, axis_label = label,  title = l + " : LP Reweighting", num_bins = nbins_, 
            normalize = False, ratio_range = ratio_range, fname = outdir + l + '_ratio_after.png' )



