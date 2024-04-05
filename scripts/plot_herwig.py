import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.Utils import *


parser = input_options()
options = parser.parse_args()

print(options)

#UL
lumi = 59.74
f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_gen/"

no_bkg = True

f_powheg = h5py.File(f_dir + "TT_powheg.h5", "r")
f_herwig = h5py.File(f_dir + "TT_herwig.h5", "r")



outdir = options.outdir
sys = ""
#CA_prefix = "2prong"
CA_prefix = ""
charge_only = False


do_sys_variations = False
do_plot = True

norm = True

jms_corr = 1.0

m_cut_min = 60.
m_cut_max = 110.
pt_cut = 225.

num_excjets = 2


if(not os.path.exists(outdir)): os.system("mkdir " + outdir)

d_powheg_w_match = Dataset(f_powheg, label = "powheg : W-matched", color = ROOT.kRed, is_data = True)
d_powheg_t_match = Dataset(f_powheg, label = "powheg : t-matched ", color = ROOT.kOrange-3, is_data = True)
d_powheg_nomatch = Dataset(f_powheg, label = "powheg : unmatched", color = ROOT.kGreen+3, is_data = True)

d_herwig_w_match = Dataset(f_herwig, label = "herwig : W-matched", color = ROOT.kRed, is_data = True)
d_herwig_t_match = Dataset(f_herwig, label = "herwig : t-matched ", color = ROOT.kOrange-3, is_data = True)
d_herwig_nomatch = Dataset(f_herwig, label = "herwig : unmatched", color = ROOT.kGreen+3, is_data = True)


powheg_gen_matching = d_powheg_w_match.f['gen_parts'][:,0]
#0 is unmatched, 1 is W matched, 2 is top matched
powheg_nomatch_cut = powheg_gen_matching < 0.1
powheg_w_match_cut = (powheg_gen_matching  > 0.9) &  (powheg_gen_matching < 1.1)
powheg_t_match_cut = (powheg_gen_matching  > 1.9) &  (powheg_gen_matching < 2.1)

d_powheg_w_match.apply_cut(powheg_w_match_cut)
d_powheg_t_match.apply_cut(powheg_t_match_cut)
d_powheg_nomatch.apply_cut(powheg_nomatch_cut)


herwig_gen_matching = d_herwig_w_match.f['gen_parts'][:,0]
herwig_nomatch_cut = herwig_gen_matching < 0.1
herwig_w_match_cut = (herwig_gen_matching  > 0.9) &  (herwig_gen_matching < 1.1)
herwig_t_match_cut = (herwig_gen_matching  > 1.9) &  (herwig_gen_matching < 2.1)

d_herwig_w_match.apply_cut(herwig_w_match_cut)
d_herwig_t_match.apply_cut(herwig_t_match_cut)
d_herwig_nomatch.apply_cut(herwig_nomatch_cut)


d_all = [d_powheg_w_match, d_powheg_t_match, d_powheg_nomatch,
         d_herwig_w_match, d_herwig_t_match, d_herwig_nomatch]

sys_keys = ['PS_ISR_up', 'PS_ISR_down', 'PS_FSR_up', 'PS_FSR_down']

for d in (d_all):
    jet_kinematics = d.f['jet_kinematics'][:]
    #msd_cut_mask = (jet_kinematics[:,3] * jms_corr > m_cut_min) & (jet_kinematics[:,3] * jms_corr < m_cut_max)
    pt_cut_mask = jet_kinematics[:,0] > pt_cut
    d.apply_cut(pt_cut_mask)
    d.compute_obs()
    weights = d.get_weights()
    d.nom_weights = weights
    if('powheg' in d.label):
        all_sys_weights = d.get_masked('sys_weights')
        d_sys_weights = []
        for sys in sys_keys:
            sys_idx = sys_weights_map[sys]
            sys_weight = all_sys_weights[:,sys_idx] * weights 
            sys_weight /= (np.mean(sys_weight) / np.mean(weights)) #correct norm
            d_sys_weights.append(sys_weight)

        d.sys_weights = d_sys_weights


obs_attrs = {
        'mSoftDrop' : (50, 230, 45, "jet mass [GeV]", "Events / 4 GeV"),
        'tau21' : (0.05, 0.8, 25, "#tau_{21}", "Events  "),
        'tau32' : (0.1, 0.95, 25, "#tau_{32}", "Events "),
        'tau43' : (0.6, 0.96, 18, "#tau_{43}", "Events "),
        'nPF' : (0.5, 100.5, 50, "Num. PF Cands.", "Events " ),
        'pt' : (pt_cut, 825., 20, "p_{T}", "Events "),
        }

pairs = [[d_powheg_w_match, d_herwig_w_match], 
        [d_powheg_t_match, d_herwig_t_match], 
        [d_powheg_nomatch, d_herwig_nomatch]]

titles = ["W-Matched", "Top-Matched", "Unmatched"]
colors = ['b', 'r']
labels = ['powheg', 'herwig']

for l in obs_attrs.keys():
    a = []

    low,high, nbins_, label, ylabel = obs_attrs.get(l, (None, None, 20, l, ""))
    for i,pair in enumerate(pairs):
        obs = [getattr(pair[0], l), getattr(pair[1], l)]
        weights = [pair[0].nom_weights, pair[1].nom_weights]
        sys_weights = [pair[0].sys_weights, []]


        make_multi_ratio_histogram(obs, weights = weights, sys_weights = sys_weights, labels = labels, colors = colors, axis_label = l, num_bins = nbins_, h_range = (low,high),
                normalize = True, ratio_range = (0.5, 1.5), title = titles[i], fname = outdir +titles[i] +"_"  + l + "_cmp.png" )

        #make_multi_sum_ratio_histogram(data = getattr(d_herwig, l), entries = a, weights = weights_rw, labels = labels, h_range = h_range, drawSys = False, stack = False, data_label = "herwig",
                #colors = colors, axis_label = l,  title = l + " : LP Reweighting", num_bins = n_bins_, normalize = True, ratio_range = (0.5, 1.5), fname = outdir + l + '_ratio_after.png' )


m_cut_min = 60.
m_cut_max = 110.
#Check non-match frac
for d in [d_powheg_w_match, d_powheg_nomatch, d_herwig_w_match, d_herwig_nomatch]:
    jet_kinematics = d.f['jet_kinematics'][:]
    msd_cut_mask = (jet_kinematics[:,3]  > m_cut_min) & (jet_kinematics[:,3]  < m_cut_max)
    d.apply_cut(msd_cut_mask)

powheg_tot = np.sum(d_powheg_w_match.get_weights()) + np.sum(d_powheg_nomatch.get_weights())
powheg_nomatch_frac = np.sum(d_powheg_nomatch.get_weights())/ powheg_tot

herwig_tot = np.sum(d_herwig_w_match.get_weights()) + np.sum(d_herwig_nomatch.get_weights())
herwig_nomatch_frac = np.sum(d_herwig_nomatch.get_weights())/ herwig_tot
print("In W region: Powheg no-match frac is %.3f, Herwig is %.3f" % (powheg_nomatch_frac, herwig_nomatch_frac))

