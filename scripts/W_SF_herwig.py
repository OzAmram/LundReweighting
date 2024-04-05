import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.Utils import *


parser = input_options()
parser.add_argument("--topSF", default=False, action='store_true',  help="Top SF")
options = parser.parse_args()

print(options)

#UL
lumi = 59.74
f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_gen/"

no_bkg = True

f_powheg = h5py.File(f_dir + "TT_powheg.h5", "r")
f_herwig = h5py.File(f_dir + "TT_herwig.h5", "r")


f_ratio_name = ""
if(options.fin != ""): f_ratio_name = options.fin
f_ratio = ROOT.TFile.Open(f_ratio_name)
h_ratio = f_ratio.Get("ratio_nom")

outdir = options.outdir

do_sys_variations = not (options.no_sys)

tau21_thresholds = [0.2, 0.27, 0.35]
tau32_thresholds = [0.3, 0.4, 0.5]


num_excjets = 2


if(not os.path.exists(outdir)): os.system("mkdir " + outdir)

d_powheg_w_match = Dataset(f_powheg, label = "powheg : W-matched", color = ROOT.kRed, is_data = True)
d_powheg_t_match = Dataset(f_powheg, label = "powheg : t-matched ", color = ROOT.kOrange-3, is_data = True)

d_herwig_w_match = Dataset(f_herwig, label = "herwig : W-matched", color = ROOT.kRed, is_data = True)
d_herwig_t_match = Dataset(f_herwig, label = "herwig : t-matched ", color = ROOT.kOrange-3, is_data = True)


powheg_gen_matching = d_powheg_w_match.f['gen_parts'][:,0]
#0 is unmatched, 1 is W matched, 2 is top matched
powheg_w_match_cut = (powheg_gen_matching  > 0.9) &  (powheg_gen_matching < 1.1)
powheg_t_match_cut = (powheg_gen_matching  > 1.9) &  (powheg_gen_matching < 2.1)

d_powheg_w_match.apply_cut(powheg_w_match_cut)
d_powheg_t_match.apply_cut(powheg_t_match_cut)


herwig_gen_matching = d_herwig_w_match.f['gen_parts'][:,0]
herwig_w_match_cut = (herwig_gen_matching  > 0.9) &  (herwig_gen_matching < 1.1)
herwig_t_match_cut = (herwig_gen_matching  > 1.9) &  (herwig_gen_matching < 2.1)

d_herwig_w_match.apply_cut(herwig_w_match_cut)
d_herwig_t_match.apply_cut(herwig_t_match_cut)

if(options.topSF):
    d_powheg, d_herwig = d_powheg_t_match, d_herwig_t_match
    num_excjets = 3
    thresholds = tau32_thresholds
    title="Top-matched"
    obs = 'tau32'
    pt_cut = 500.
else:
    d_powheg, d_herwig = d_powheg_w_match, d_herwig_w_match
    num_excjets = 2
    thresholds = tau21_thresholds
    title="W-matched"
    obs = 'tau21'
    pt_cut = 225.


if('pt_extrap' in f_ratio.GetListOfKeys()):
    rdir = f_ratio.GetDirectory("pt_extrap")
    rdir.cd()
else: 
    print("NO Pt extrapolation")
    rdir = None

for d in [d_powheg, d_herwig]:
    jet_kinematics = d.f['jet_kinematics'][:]
    pt_cut_mask = jet_kinematics[:,0] > pt_cut
    d.apply_cut(pt_cut_mask)
    d.compute_obs()
    weights = d.get_weights()
    d.nom_weights = weights

weights_nom = d_powheg.nom_weights
print("%i powheg, %i herwig evts" % (len(d_powheg.nom_weights), len(d_herwig.nom_weights)))


#Noise used to generated smeared ratio's based on stat unc
nToys = 100
rand_noise = np.random.normal(size = (nToys, h_ratio.GetNbinsX(), h_ratio.GetNbinsY(), h_ratio.GetNbinsZ()))
pt_rand_noise = np.random.normal(size = (nToys, h_ratio.GetNbinsY(), h_ratio.GetNbinsZ(), 3))

LP_rw = LundReweighter(jetR = jetR, pt_extrap_dir = rdir, charge_only = options.charge_only)

reclust_nom, reclust_prongs_up, reclust_prongs_down = d_powheg.get_all_matched_splittings(LP_rw, num_excjets = num_excjets)

LP_weights, stat_smeared_weights, pt_smeared_weights = d_powheg.reweight_all(LP_rw, h_ratio, reclust_objs = reclust_nom, num_excjets = num_excjets, 
        rand_noise = rand_noise, pt_rand_noise = pt_rand_noise)

LP_weights_prong_up, LP_weights_prong_down, LP_weights_match_up, LP_weights_match_down = d_powheg.reweight_all_up_down_prongs(
        LP_rw, h_ratio = h_ratio, reclust_prongs_up = reclust_prongs_up, reclust_prongs_down = reclust_prongs_down, nom_weights = LP_weights)

sys_ratios = []
sys_variations = dict()
LP_weights_sys_up = LP_weights_sys_down = b_weights_up = b_weights_down = np.ones_like(weights_nom)
if(do_sys_variations):
    #sys_list = list(sys_weights_map.keys())
    sys_list = ['sys_tot_up', 'sys_tot_down']
    for sys in sys_list:

        if(f_ratio.GetListOfKeys().Contains("ratio_" + sys)):
            sys_ratio = f_ratio.Get("ratio_" + sys)
            sys_str = sys + "_"

            sys_LP_weights = d_powheg.reweight_all(LP_rw, sys_ratio, num_excjets = num_excjets, prefix = "", 
                    sys_str = sys_str, reclust_objs = reclust_nom)
            sys_weights = weights_nom * sys_LP_weights
            rw = np.sum(weights_nom) / np.sum(sys_weights)
            sys_weights *= rw
        else:
            print("Missing sys %s!" % sys)
            sys_weights = LP_weights

        if('up' in sys): LP_weights_sys_up = sys_weights
        else: LP_weights_sys_down = sys_weights

    if(f_ratio.GetListOfKeys().Contains("h_bl_ratio")):
        b_light_ratio = f_ratio.Get("h_bl_ratio")
        bquark_rw = d_powheg.reweight_all(LP_rw, b_light_ratio, num_excjets = num_excjets, prefix = "", sys_str = 'bquark', reclust_objs = reclust_nom)
    else: 
        bquark_rw = np.ones_like(LP_weights)

    b_weights_down = (1./ bquark_rw) * LP_weights
    b_weights_up = bquark_rw * LP_weights


all_badmatches = [ro.badmatch for ro in reclust_nom]
print("Bad match frac is %.2f" % np.mean(all_badmatches))



#The nominal Lund Plane correction event weights
LP_weights = LP_rw.normalize_weights(LP_weights) * weights_nom 

#Toy variations for stat and pt uncertainties
stat_smeared_weights = LP_rw.normalize_weights(stat_smeared_weights) * weights_nom.reshape(-1, 1)
pt_smeared_weights = LP_rw.normalize_weights(pt_smeared_weights) * weights_nom.reshape(-1,1)

#Systematic up/down variations
LP_weights_sys_up = LP_rw.normalize_weights(LP_weights_sys_up) * weights_nom
LP_weights_sys_down = LP_rw.normalize_weights(LP_weights_sys_down) * weights_nom

LP_weights_match_up = LP_rw.normalize_weights(LP_weights_match_up) * weights_nom
LP_weights_match_down = LP_rw.normalize_weights(LP_weights_match_down) * weights_nom

LP_weights_prong_up = LP_rw.normalize_weights(LP_weights_prong_up) * weights_nom
LP_weights_prong_down = LP_rw.normalize_weights(LP_weights_prong_down) * weights_nom

diffs = LP_weights_match_up - LP_weights
print("Diffs")
print(np.mean(diffs > 0.1))
print(diffs[diffs > 0.1][:10])


#b quark systematic up/down variations
b_weights_up = LP_rw.normalize_weights(b_weights_up) * weights_nom
b_weights_down = LP_rw.normalize_weights(b_weights_down) * weights_nom


make_histogram(LP_weights, "Reweighting factors", 'b', 'Weight', "Lund Plane Reweighting Factors", 20 , h_range = (0., 5.0),
     normalize=False, fname=outdir + "lundPlane_weights.png")

#Save subjet pts and deltaR
subjet_pts =  []

for i,ro in enumerate(reclust_nom):
    sjs = ro.subjet
    for sj in sjs:
        subjet_pts.append(sj[0])
    
print("Min subjet pt %.2f " % np.amin(subjet_pts))
num_bins = 40
pt_bins = array('d', np.linspace(0., 800., num_bins + 1))

h_subjet_pts = make_root_hist(data = subjet_pts, name = 'h_W_subjetpt', num_bins = num_bins, bins = pt_bins)

#compute 'Scalefactor'
powheg_cuts = [getattr(d_powheg, obs) < thresh for thresh in thresholds]
herwig_cuts = [getattr(d_herwig, obs) < thresh for thresh in thresholds]

print("%i Powheg evts, %i Herwig" % (len(getattr(d_powheg, obs)), len(getattr(d_herwig, obs))))

f_effs = open(options.outdir + "Effs.txt", "w")

for idx in range(len(powheg_cuts)):

    pow_cut = powheg_cuts[idx]
    her_cut = herwig_cuts[idx]

    eff_nom = np.average(pow_cut, weights = weights_nom)
    eff_rw = np.average(pow_cut, weights = LP_weights)

    eff_herwig = np.average(her_cut)

    print("Powheg: Nom %.3f, RW %.3f" % (eff_nom, eff_rw))
    print("Herwig: %.3f " % (eff_herwig))


    eff_toys = []
    pt_eff_toys = []
    for i in range(nToys):
        eff = np.average(pow_cut, weights = stat_smeared_weights[:,i])
        eff_toys.append(eff)

        eff1 = np.average(pow_cut, weights = pt_smeared_weights[:,i])
        pt_eff_toys.append(eff1)

    #Compute stat and pt uncertainty based on variation in the toys
    toys_mean = np.mean(eff_toys)
    toys_std = np.std(eff_toys)
    pt_toys_mean = np.mean(pt_eff_toys)
    pt_toys_std = np.std(pt_eff_toys)

    eff_stat_unc = (abs(toys_mean - eff_rw)  + toys_std) 
    eff_pt_unc = (abs(pt_toys_mean - eff_rw) + pt_toys_std)

    print("Stat variation toys eff. avg %.3f, std dev %.3f" % (toys_mean, toys_std))
    print("Pt variation toys eff. avg %.3f, std dev %.3f" % (pt_toys_mean, pt_toys_std))

    #Add systematic differences in quadrature
    sys_LPunc_up = sys_LPunc_down = sys_prongs_up = sys_prongs_down = sys_matching_up = sys_matching_down = b_unc_up = b_unc_down = 0.
    if(do_sys_variations):

        #Compute difference in efficiency due to weight variations as uncertainty
        def get_uncs(score_cut, weights_up, weights_down, eff_baseline):
            eff_up =  np.average(score_cut, weights = weights_up)
            eff_down =  np.average(score_cut, weights = weights_down)
            print(eff_up, eff_down, eff_baseline)

            unc_up = eff_up - eff_baseline
            unc_down = eff_down - eff_baseline 
            return unc_up, unc_down


        #Compute efficiency of systematic variations
        sys_LPunc_up, sys_LPunc_down = get_uncs(pow_cut, LP_weights_sys_up, LP_weights_sys_down, eff_rw)
        print("prong")
        sys_prong_up, sys_prong_down = get_uncs(pow_cut, LP_weights_prong_up, LP_weights_prong_down, eff_rw)
        print("match")
        sys_matching_up, sys_matching_down = get_uncs(pow_cut, LP_weights_match_up, LP_weights_match_down, eff_rw)

        b_unc_up, b_unc_down = get_uncs(pow_cut, b_weights_up, b_weights_down, eff_rw)


    f_effs.write("Thresh %.2f" % thresholds[idx])

    eff_str = ("Calibrated efficiency  is %.2f +/- %.2f  (stat) +/- %.2f (pt) %.2f/%.2f (sys) %.2f/%.2f (bquark) %.2f/%.2f (prongs)  %.2f/%.2f (matching) \n"  % 
        (eff_rw, eff_stat_unc, eff_pt_unc, sys_LPunc_up, sys_LPunc_down, b_unc_up, b_unc_down, sys_prongs_up, sys_prongs_down, sys_matching_up, sys_matching_down ))

    tot_unc_up = (eff_stat_unc**2 + eff_pt_unc**2 + max(sys_LPunc_up, sys_LPunc_down)**2 + max(b_unc_up, b_unc_down)**2 + 
            max(sys_prongs_up, sys_prongs_down)**2 + max(sys_matching_up, sys_matching_down)**2)**0.5

    tot_unc_down = (eff_stat_unc**2 + eff_pt_unc**2 + min(sys_LPunc_up, sys_LPunc_down)**2 + min(b_unc_up, b_unc_down)**2 + 
            min(sys_prongs_up, sys_prongs_down)**2 + min(sys_matching_up, sys_matching_down)**2)**0.5

    eff_str += "Herwig %.2f, Reweighted Powheg %.2f +%.2f/-%.2f \n"  % (eff_herwig, eff_rw, tot_unc_up, tot_unc_down)


    print(eff_str)
    f_effs.write(eff_str)
    #approximate uncertainty on the reweighting for the plots


f_effs.close()
f_ratio.Close()

obs_attrs = {
        'mSoftDrop' : (50, 230, 45, "m_{SD} [GeV]", "Events / 4 GeV"),
        'tau21' : (0.05, 0.8, 25, "#tau_{21}", "Events  "),
        'tau32' : (0.1, 0.95, 25, "#tau_{32}", "Events "),
        'tau43' : (0.6, 0.96, 18, "#tau_{43}", "Events "),
        'nPF' : (0.5, 100.5, 50, "Num. PF Cands.", "Events " ),
        'pt' : (pt_cut, 825., 20, "p_{T}", "Events "),
        }
labels = ['herwig', 'pythia', 'pythia, reweighted']
colors = ['b', 'r', 'green']



for l in obs_attrs.keys():
    a = []

    low,high, nbins_, label, ylabel = obs_attrs.get(l, (None, None, 20, l, ""))

    obs = [getattr(d_herwig, l), getattr(d_powheg, l), getattr(d_powheg, l)]
    weights = [d_herwig.nom_weights, d_powheg.nom_weights, LP_weights]
    sys_weights = None

    make_multi_ratio_histogram(obs, weights = weights, sys_weights = sys_weights, labels = labels, colors = colors, axis_label = label, num_bins = nbins_, h_range = (low, high),
            normalize = True, ratio_range = (0.5, 1.5), title = title, fname = outdir +title + '_' + l + "_cmp.png" )

