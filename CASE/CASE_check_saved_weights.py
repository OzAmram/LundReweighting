from Utils import *
import os



#Setup 
parser = input_options()
options = parser.parse_args()

outdir = options.outdir
if(not os.path.exists(outdir)): os.system("mkdir %s" % outdir)
jet_str = 'CA'

#fname = "test_signal_CASE.h5"
fname = "/eos/uscms/store/user/oamram/case/sig_files/LundRW/XToYYprimeTo4Q_MX3000_MY170_MYprime170_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER_Lund.h5"
label = "XtoYY"

f_sig = h5py.File(fname, "r")
d = Dataset(f_sig, label = label, color = ROOT.kRed, dtype = 1)

print("input file", fname)
d.compute_obs()

obs = ["tau21", "tau32", "tau43", "nPF", "mSoftDrop", "pt"]
n_bins = 20



#Cut we will check the efficiency of
tag_obs = 'tau21'
score_thresh = 0.34


#Only consider limited sample of events for simplicity
#max_evts = 2000
max_evts = None


score = getattr(d, tag_obs)[:max_evts]

print("%i events" % score.shape[0])
score_cut = score < score_thresh

weights_nom = np.ones_like(score)
weights_rw = copy.deepcopy(weights_nom)


#get saved lund weight info
weights_rw = d.f['lund_weights'][:max_evts]
print("Mean RW", np.mean(weights_rw))
print("First weights", weights_rw[:10])

#stat and pt extrapolation toy weights
weights_stat = d.f['lund_weights_stat_var'][:max_evts]
weights_pt = d.f['lund_weights_pt_var'][:max_evts]

#multiply sys variations  by nominal weight
weights_sys = d.f['lund_weights_sys_var'][:max_evts] * np.expand_dims(weights_rw, -1)

#matching unc (single number)
bad_match_frac = d.f['lund_weights_matching_unc'][0]


#compute nominal 'Scalefactor' for our selection
eff_nom = np.average(score_cut, weights = weights_nom)
eff_rw = np.average(score_cut, weights = weights_rw)

print("Nom %.3f, RW %.3f" % (eff_nom, eff_rw))
SF = eff_rw / eff_nom

#Plot distribution of Lund weights
make_histogram(weights_rw, "Reweighting factors", 'b', 'Weight', "Lund Plane Reweighting Factors", 20 , h_range = (0., 2.0),
     normalize=False, fname=outdir + "lundPlane_weights.png")

#Compute selection efficiency of pt and stat variations
eff_toys = []
pt_eff_toys = []
nToys = 100
for i in range(nToys):
    eff = np.average(score_cut, weights = weights_stat[:,i])
    eff_toys.append(eff)

    eff1 = np.average(score_cut, weights = weights_pt[:,i])
    pt_eff_toys.append(eff1)

toys_mean = np.mean(eff_toys)
toys_std = np.std(eff_toys)

print("Stat toys avg %.3f, std dev %.3f" % (toys_mean, toys_std))

pt_toys_mean = np.mean(pt_eff_toys)
pt_toys_std = np.std(pt_eff_toys)

print("Pt variation toys avg %.3f, std dev %.3f" % (pt_toys_mean, pt_toys_std))

SF_stat_unc = abs(toys_mean - eff_rw)/eff_nom + toys_std /eff_nom
SF_pt_unc = abs(pt_toys_mean - eff_rw)/eff_nom + pt_toys_std /eff_nom


#Compute systematics for two up/down variations
sys_unc_up = sys_unc_down = 0.

#Regular sys unc. on Lund Plane
eff_sys_tot_up = np.average(score_cut, weights = weights_sys[:,0])
eff_sys_tot_down = np.average(score_cut, weights = weights_sys[:,1])
SF_sys_unc_up = abs(eff_sys_tot_up - eff_rw)/eff_nom
SF_sys_unc_down = abs(eff_sys_tot_down - eff_rw)/eff_nom
SF_sys_unc = (SF_sys_unc_up + SF_sys_unc_down) / 2.0

#b-quark specific variations
eff_bquark_up = np.average(score_cut, weights = weights_sys[:,2])
eff_bquark_down = np.average(score_cut, weights = weights_sys[:,3])
SF_bquark_up = abs(eff_bquark_up - eff_rw)/eff_nom
SF_bquark_down = abs(eff_bquark_down - eff_rw)/eff_nom
SF_bquark_unc = (SF_bquark_up + SF_bquark_down) /2.0

#fraction of evts with bad match, take as fractional unc on SF
bad_matching_unc =  bad_match_frac * SF

print("\n\nSF (%s val %.2f ) is %.2f +/- %.2f  (stat) +/- %.2f (sys) +/- %.2f (pt) +/- %.2f (bquark) +/- %.2f (matching) \n\n"  
        % (tag_obs, score_thresh, SF, SF_stat_unc, SF_sys_unc, SF_pt_unc, SF_bquark_unc, bad_matching_unc))


overall_unc = (SF_stat_unc **2 + SF_sys_unc**2 + SF_pt_unc**2 + SF_bquark_unc**2 + bad_matching_unc**2) **0.5 
print("overall unc +/- %.3f" % overall_unc)



#Make plots of different observables before vs after reweighting
#also include up/down shape variations from the sys variation 
weights = [  weights_nom, weights_rw, weights_sys[:,0], weights_sys[:,1]]
labels = ["Nom", "RW", "RW Sys. Up", "RW Sys. Down"]
colors = ['gray','black', 'blue', 'red']
ratio_range = [0.5, 1.5]


for l in obs:

    if(l == 'nPF'): 
        h_range = (0.5,120.5)
        n_bins_ = 40
    else: 
        n_bins_ = n_bins
        h_range = None

    x = getattr(d,l)[:max_evts]

    make_multi_ratio_histogram([x]*4, labels, colors, l, "%s : Before vs. After Lund Plane Reweighting" % l, n_bins,
         ratio_range = ratio_range, normalize=True, weights = weights, fname=outdir + "%s_before_vs_after.png" %l )

