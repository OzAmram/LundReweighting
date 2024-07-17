import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.Utils import *


parser = input_options()
parser.add_argument("--LPorder", default=1, type=int,  help="LP max order")
parser.add_argument("--reco", default=False, action='store_true',  help="Reco level")
parser.add_argument("--YtoHH", default=False, action='store_true',  help="YtoHH signal")
options = parser.parse_args()

tau43_thresholds = [0.45, 0.55, 0.65, 0.75]
DeepH4q_thresholds = [0.3, 0.4, 0.6, 0.7]
tau54_thresholds = [0.55, 0.65, 0.75]

pt_cut = 1000

print(options)

if(options.YtoHH):
    sig_name = "YtoHH"
    title = r"H$\to t\bar{t}$ (6 pronged)"
    if(options.reco):
        obs = 'DeepAK8_H4q'
        thresholds = DeepH4q_thresholds
    else:
        obs = 'tau54'
        thresholds = tau54_thresholds
else:
    sig_name = "Wkk"
    title = r"R $\to WW$ (4 pronged)"
    if(options.reco):
        obs = 'DeepAK8_H4q'
        thresholds = DeepH4q_thresholds
    else:
        obs = 'tau43'
        thresholds = tau43_thresholds


#UL
lumi = 59.74
if(not options.reco):
    f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_herwig/"
    f_pythia = h5py.File(f_dir + "%s_pythia.h5" %sig_name, "r")
    f_herwig = h5py.File(f_dir + "%s_herwig.h5" %sig_name, "r")
elif(options.YtoHH):
    f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/"
    f_herwig = h5py.File(f_dir + "Lund_output_files_herwig/YtoHH_herwig_reco.h5", "r")
    f_pythia = h5py.File(f_dir + "Lund_output_files_herwig/YtoHH_pythia_reco.h5", "r")
else:
    f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/"
    f_herwig = h5py.File(f_dir + "Lund_output_files_herwig/Wkk_herwig_reco.h5", "r")
    f_pythia = h5py.File(f_dir + "Lund_output_files_herwig/Wkk_pythia_reco.h5", "r")

    #fname = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/TagNTrain/data/LundRW/WkkToWRadionToWWW_M3000_Mr400_TuneCP5_13TeV-madgraph-pythia8_TIMBER_Lund.h5"
    #f_pythia = h5py.File(fname, "r")

print(f_herwig, f_pythia)





f_ratio_name = ""
if(options.fin != ""): f_ratio_name = options.fin
f_ratio = ROOT.TFile.Open(f_ratio_name)

outdir = options.outdir

do_sys_variations = not (options.no_sys)




max_evts = None




if(not os.path.exists(outdir)): os.system("mkdir " + outdir)

d_pythia = Dataset(f_pythia, label = "pythia", color = ROOT.kRed,  dtype=1, gen=not options.reco)
d_herwig = Dataset(f_herwig, label = "herwig", color = ROOT.kRed, dtype=1, gen=not options.reco)





for d in [d_pythia, d_herwig]:
    print(list(d.f.keys()))
    jet_kinematics = d.f['jet_kinematics'][:]
    gen_info = d.f['gen_info'][:]
    if(not options.reco): gen_pdg_id = np.abs(gen_info[:,2:,3])
    else: gen_pdg_id = np.abs(gen_info[:,:,3])
    all_had = np.all(gen_pdg_id < 6, axis=-1)
    d.apply_cut(all_had)
    if( not options.reco): pt_cut_mask = jet_kinematics[:,0] > pt_cut
    else: pt_cut_mask = jet_kinematics[:,2] > pt_cut
    d.apply_cut(pt_cut_mask)
    d.compute_obs()
    d.nom_weights = d.get_weights()[:max_evts]

nom_weights = d_pythia.nom_weights
print("%i pythia, %i herwig evts" % (len(d_pythia.nom_weights), len(d_herwig.nom_weights)))

LP_rw = LundReweighter(jetR = jetR, f_ratio = f_ratio, charge_only = options.charge_only, LP_order = options.LPorder )

LP_weights = d_pythia.reweight_all(LP_rw, do_sys_weights = do_sys_variations, max_evts = max_evts)


for key in LP_weights.keys():
    if('nom' in key or 'up' in key or 'down' in key or 'vars' in key):
        if(not isinstance(LP_weights[key], np.ndarray)): continue
        if('vars' in key): LP_weights[key] *= nom_weights.reshape(-1,1)
        else: LP_weights[key] *= nom_weights


print("Bad match frac %.2f" % np.mean(LP_weights['bad_match']))
print("Reclustered bad match frac %.2f" % np.mean(LP_weights['reclust_still_bad_match']))


make_histogram(LP_weights['nom'], "Reweighting factors", 'b', 'Weight', "Lund Plane Reweighting Factors", 20 , h_range = (0., 5.0),
     normalize=False, fname=outdir + "lundPlane_weights.png")

print(LP_weights['n_prongs'].shape, nom_weights.shape, LP_weights['nom'].shape)
make_histogram([LP_weights['n_prongs'], LP_weights['n_prongs']], labels = ['original', 'reweighted'], xaxis_label = "nProngs", colors = ['b', 'r'], num_bins=6 , h_range = (0.5, 6.5),
     normalize=True, weights = [nom_weights, LP_weights['nom']], fname=outdir + "nProngs.png")

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
pythia_cuts = [getattr(d_pythia, obs)[:max_evts] < thresh for thresh in thresholds] if 'Deep' not in obs else  [getattr(d_pythia, obs)[:max_evts] > thresh for thresh in thresholds]
herwig_cuts = [getattr(d_herwig, obs)[:max_evts] < thresh for thresh in thresholds] if 'Deep' not in obs else  [getattr(d_herwig, obs)[:max_evts] > thresh for thresh in thresholds]

f_effs = open(options.outdir + "Effs.txt", "w")

for idx in range(len(pythia_cuts)):

    pow_cut = pythia_cuts[idx]
    her_cut = herwig_cuts[idx]

    eff_nom = np.average(pow_cut, weights = nom_weights)
    eff_rw = np.average(pow_cut, weights = LP_weights['nom'])

    eff_herwig = np.average(her_cut)

    print("pythia: Nom %.3f, RW %.3f" % (eff_nom, eff_rw))
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
    sys_keys = ['sys', 'bquark', 'prongs', 'unclust', 'distortion']
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

    eff_str += "\n Herwig %.2f, Reweighted pythia %.2f +%.2f/-%.2f \n"  % (eff_herwig, eff_rw, tot_unc_up, tot_unc_down)


    print(eff_str)
    f_effs.write(eff_str)
    #approximate uncertainty on the reweighting for the plots


f_effs.close()
f_ratio.Close()

#Plotting
tau21_start = 0.05
tau32_start = 0.1
tau43_start = 0.25
tau54_start = 0.35
obs_attrs = {
        'mSoftDrop' : (200, 500, 30, r"$m_\mathrm{SD}$ [GeV]", "Events / 4 GeV", 'upper left'),
        'tau21' : (tau21_start, 0.8, 12, r"$\tau_{21}$", "Events  ", 'upper left'),
        'tau32' : (tau32_start, 0.9, 15, r"$\tau_{32}$", "Events ", 'upper left'),
        'tau43' : (tau43_start, 0.96, 12, r"$\tau_{43}$", "Events ", 'upper left'),
        'nPF' : (0.5, 100.5, 50, "Num. PF Cands.", "Events ", 'upper left' ),
        'pt' : (300., 825., 20, r"$p_{T}$", "Events ", 'upper left'),
        }

if(not options.reco):
    obs_attrs['tau54'] = (tau54_start, 1.05, 12, r"$\tau_{54}$", "Events "),
    obs_attrs['tau65'] =  (tau54_start, 1.05, 12, r"$\tau_{65}$", "Events "),

else:
    obs_attrs['ParticleNet_W'] = (0., 1., 15, r"ParticleNet W Tag Score", "Events ", 'upper left')
    obs_attrs['ParticleNet_H4q'] = (0., 1., 15, r"ParticleNet H4q Tag Score", "Events ", 'upper left')
    obs_attrs['DeepAK8_W'] = (0., 1., 15, r"DeepAK8 W Tag Score", "Events ", 'upper left')
    obs_attrs['DeepAK8_W_MD'] = (0., 1., 15, r"DeepAK8 W MD Tag Score", "Events ", 'upper left')
    obs_attrs['DeepAK8_H4q'] = (0., 1., 15, r"DeepAK8 H4q Tag Score", "Events ", 'upper left')
labels = ['herwig', 'pythia', 'pythia, reweighted']
colors = [c_red, c_lightblue, c_purple]



hist_weights = [d_herwig.nom_weights, d_pythia.nom_weights, LP_weights['nom']]
sys_list = ['sys', 'bquark', 'prongs', 'unclust', 'distortion']
hist_sys_weights = [ [], [[LP_weights[sys+"_up"], LP_weights[sys+"_down"]] for sys in sys_list]]
hist_stat_weights = [ [], [LP_weights['stat_vars'], LP_weights['pt_vars']]]

#hist_sys_weights = None
#hist_stat_weights = None


for l in obs_attrs.keys():
    print(l)
    a = []

    low,high, nbins_, label, ylabel, leg_loc = obs_attrs.get(l, (None, None, 20, l, "", "best"))
    obs = [getattr(d_herwig, l)[:max_evts], getattr(d_pythia, l)[:max_evts], getattr(d_pythia, l)[:max_evts]]

    make_herwig_ratio_histogram(obs, weights = hist_weights, sys_weights = hist_sys_weights, stat_weights = hist_stat_weights, first_like_data = True, 
            labels = labels, colors = colors, axis_label = label, num_bins = nbins_, h_range = (low, high), leg_loc = leg_loc,
            normalize = True, ratio_range = (0.5, 1.5), title = title, fname = outdir + l + "_cmp.png" )

