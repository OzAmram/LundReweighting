import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.Utils import *

def make_herwig_ratio_histogram(entries = None, labels = None, colors = None, axis_label = None, title = None, num_bins=10, normalize = False, h_range = None, first_like_data = True,
        weights = None, fname="", ratio_range = -1, errors = False, logy = False, max_rw = 5, sys_weights = None):
    h_type= 'step'
    alpha = 1.
    fontsize = 22
    label_size = 18
    lw = 3
    fig = plt.figure(figsize = fig_size)
    gs = gridspec.GridSpec(2,1, height_ratios = [3,1])
    ax0 =  plt.subplot(gs[0])

    if(h_range is None):
        low = np.amin(entries[0])
        high = np.amax(entries[0])
    else:
        low,high = h_range

    data,bins = np.histogram(entries[0], bins = num_bins, range=(low,high), weights = weights[0])
    data_uncs = np.sqrt(data)

    if(normalize):
        norm = np.sum(data) * ((high - low)/num_bins) #counts * bin width
        data /= norm
        data_uncs /=norm

    ns, bins, patches  = ax0.hist(entries[1:], bins=num_bins, range=(low,high), color=colors[1:], alpha=alpha,label=labels[1:len(entries)], 
            density = normalize, weights = weights[1:], histtype=h_type, linewidth = lw)

    bincenters = 0.5*(bins[1:]+bins[:-1]) 
    ax0.errorbar(bincenters, data, yerr=data_uncs, fmt='ko', markerfacecolor = colors[0], ecolor = colors[0], markeredgecolor = colors[0], label = labels[0])

    plt.xlim(h_range)

    if(logy): plt.yscale("log")
    plt.title(title, fontsize=fontsize)

    bin_size = bins[1] - bins[0]
    bincenters = 0.5*(bins[1:]+bins[:-1]) 

    ax1 = plt.subplot(gs[1])
    ax1.errorbar(bincenters, np.ones_like(bincenters), yerr = data_uncs/data, markerfacecolor = colors[0], ecolor = colors[0], markeredgecolor = colors[0], fmt='ko')

    ratios = []
    for i in range(len(ns)):
        ratio =  np.clip(data, 1e-8, None)/ np.clip(ns[i], 1e-8, None)
        ratio = np.append(ratio, ratio[-1])
        ratios.append(ratio)

        plt.step(bins, ratio, color = colors[i+1], linewidth = lw, where = 'post')

    leg = ax0.legend(loc='best', fontsize = 14)

    #draw sys unc band
    if(sys_weights is not None):
        for j in range(len(sys_weights)): # separate set for each observable
            sys_w = sys_weights[j]
            if(len(sys_w) == 0): continue
            nom = ns[j]

            uncs_up = np.zeros(nom.shape)
            uncs_down = np.zeros(nom.shape)
            for i in range(len(sys_w)):
                weights_up = sys_w[i][0]
                weights_down = sys_w[i][1]

                ns_sys_up, _ = np.histogram(entries[j+1], bins = bins, weights = weights_up, density =normalize )
                ns_sys_down, _ = np.histogram(entries[j+1], bins = bins, weights = weights_down, density =normalize )

                ns_up, ns_down = np.maximum(ns_sys_up, ns_sys_down), np.minimum(ns_sys_up, ns_sys_down)

                uncs_up += (ns_up - nom)**2
                uncs_down += (ns_down - nom)**2

            uncs_up = np.sqrt(uncs_up)
            uncs_down = np.sqrt(uncs_down)

            vals_up = nom + uncs_up
            vals_down = nom - uncs_down

            #compute ratio unc
            ratio_up =  np.clip(data, 1e-8, None)/ np.clip(vals_up, 1e-8, None) 
            ratio_down =  np.clip(data, 1e-8, None)/ np.clip(vals_down, 1e-8, None) 
            
            #Pad with dummy val
            vals_up = np.append(vals_up, vals_up[-1])
            vals_down = np.append(vals_down, vals_down[-1])

            ratio_up = np.append(ratio_up, ratio_up[-1])
            ratio_down = np.append(ratio_down, ratio_down[-1])

            ax0.fill_between(bins, vals_down, vals_up, color = colors[j+1], alpha = 0.5, step = 'post')
            ax1.fill_between(bins, ratio_down, ratio_up, color = colors[j+1], alpha = 0.5, step = 'post')


    ax1.set_ylabel("Herwig/Pythia", fontsize= label_size)
    ax1.set_xlabel(axis_label, fontsize = label_size)

    plt.xlim([low, high])

    if(type(ratio_range) == list or type(ratio_range) == tuple):
        plt.ylim(ratio_range[0], ratio_range[1])
    else:
        if(ratio_range > 0):
            plt.ylim([1-ratio_range, 1+ratio_range])

    plt.grid(axis='y')


    if(fname != ""): 
        plt.savefig(fname)
        print("saving fig %s" %fname)

    return 



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

outdir = options.outdir

do_sys_variations = not (options.no_sys)

tau21_thresholds = [0.18, 0.24, 0.3]
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



for d in [d_powheg, d_herwig]:
    jet_kinematics = d.f['jet_kinematics'][:]
    pt_cut_mask = jet_kinematics[:,0] > pt_cut
    d.apply_cut(pt_cut_mask)
    d.compute_obs()
    d.nom_weights = d.get_weights()

nom_weights = d_powheg.nom_weights
print("%i powheg, %i herwig evts" % (len(d_powheg.nom_weights), len(d_herwig.nom_weights)))

LP_rw = LundReweighter(jetR = jetR, f_ratio = f_ratio, charge_only = options.charge_only)

LP_weights = d_powheg.reweight_all(LP_rw, do_sys_weights = do_sys_variations)

for key in LP_weights.keys():
    if('nom' in key or 'up' in key or 'down' in key):
        LP_weights[key] *= nom_weights


print("Bad match frac %.2f" % np.mean(LP_weights['bad_match']))
print("Reclustered bad match frac %.2f" % np.mean(LP_weights['reclust_still_bad_match']))


make_histogram(LP_weights['nom'], "Reweighting factors", 'b', 'Weight', "Lund Plane Reweighting Factors", 20 , h_range = (0., 5.0),
     normalize=False, fname=outdir + "lundPlane_weights.png")
quantiles = np.arange(0.,1.0, 0.1)
threshs = np.quantile(LP_weights['nom'], quantiles)
print("Quantiles are ", threshs)



#Save subjet pts and deltaR
subjet_pts =  LP_weights['subjet_pts']

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

    eff_nom = np.average(pow_cut, weights = nom_weights)
    eff_rw = np.average(pow_cut, weights = LP_weights['nom'])

    eff_herwig = np.average(her_cut)

    print("Powheg: Nom %.3f, RW %.3f" % (eff_nom, eff_rw))
    print("Herwig: %.3f " % (eff_herwig))


    nToys = LP_weights['stat_vars'].shape[1]
    eff_toys = []
    pt_eff_toys = []
    for i in range(nToys):
        eff = np.average(pow_cut, weights = LP_weights['stat_vars'][:,i])
        eff_toys.append(eff)

        eff1 = np.average(pow_cut, weights = LP_weights['pt_vars'][:,i])
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
    sys_keys = ['sys', 'bquark', 'prongs', 'matching', 'unclust']
    sys_uncs = dict()

    diffs_up = np.abs(LP_weights['nom'] - LP_weights['prongs_up'])
    diffs_down = np.abs(LP_weights['nom'] - LP_weights['prongs_down'])
    filt = diffs_down > 0.1
    print(diffs_up[filt][:5], diffs_down[filt][:5])
    print("Avg up down", np.mean(diffs_up), np.mean(diffs_down))

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
            unc_up, unc_down = get_uncs(pow_cut, LP_weights[sys + '_up'], LP_weights[sys + '_down'], eff_rw)
            print(sys, unc_up, unc_down)
            sys_uncs[sys] = [unc_up, unc_down]


    f_effs.write("Thresh %.2f" % thresholds[idx])

    eff_str = "Calibrated efficiency  is %.2f +/- %.2f (stat) +/- %.2f (pt)" % (eff_rw, eff_stat_unc, eff_pt_unc )
    tot_unc_up = tot_unc_down = eff_stat_unc**2 + eff_pt_unc**2

    for sys in sys_keys:
        eff_str += " %.2f/%.2f (%s)" % (sys_uncs[sys][0], sys_uncs[sys][1], sys)
        up_var = max(sys_uncs[sys][0], sys_uncs[sys][1])
        down_var = min(sys_uncs[sys][0], sys_uncs[sys][1])
        tot_unc_up += up_var**2
        tot_unc_down += down_var**2



    tot_unc_up = tot_unc_up**0.5
    tot_unc_down = tot_unc_down**0.5

    eff_str += "\n Herwig %.2f, Reweighted Powheg %.2f +%.2f/-%.2f \n"  % (eff_herwig, eff_rw, tot_unc_up, tot_unc_down)


    print(eff_str)
    f_effs.write(eff_str)
    #approximate uncertainty on the reweighting for the plots


f_effs.close()
f_ratio.Close()

#Plotting
tau21_start = 0.2 if options.topSF else 0.05
tau32_start = 0.15 if options.topSF else 0.3
tau43_start = 0.4 if options.topSF else 0.6
obs_attrs = {
        'mSoftDrop' : (50, 230, 45, "m_{SD} [GeV]", "Events / 4 GeV"),
        'tau21' : (tau21_start, 0.8, 12, r"$\tau_{21}$", "Events  "),
        'tau32' : (tau32_start, 0.9, 15, r"$\tau_{32}$", "Events "),
        'tau43' : (tau43_start, 0.96, 12, r"$\tau_{43}$", "Events "),
        'nPF' : (0.5, 100.5, 50, "Num. PF Cands.", "Events " ),
        'pt' : (pt_cut, 825., 20, r"$p_{T}$", "Events "),
        }
labels = ['herwig', 'pythia', 'pythia, reweighted']
colors = [c_red, c_lightblue, c_purple]



hist_weights = [d_herwig.nom_weights, d_powheg.nom_weights, LP_weights['nom']]
sys_list = ['sys', 'bquark', 'prongs', 'matching', 'unclust']
hist_sys_weights = [ [], [[LP_weights[sys+"_up"], LP_weights[sys+"_down"]] for sys in sys_list]]


for l in obs_attrs.keys():
    print(l)
    a = []

    low,high, nbins_, label, ylabel = obs_attrs.get(l, (None, None, 20, l, ""))
    obs = [getattr(d_herwig, l), getattr(d_powheg, l), getattr(d_powheg, l)]

    make_herwig_ratio_histogram(obs, weights = hist_weights, sys_weights = hist_sys_weights, first_like_data = True, 
            labels = labels, colors = colors, axis_label = label, num_bins = nbins_, h_range = (low, high),
            normalize = True, ratio_range = (0.5, 1.5), title = title, fname = outdir +title + '_' + l + "_cmp.png" )

