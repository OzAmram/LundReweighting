import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.Utils import *



#Setup 
parser = input_options()
options = parser.parse_args()

outdir = options.outdir
if(not os.path.exists(outdir)): os.system("mkdir %s" % outdir)

#fname = "test_signal_CASE.h5"
#fname = "test_signal_CASE.h5"
#fname = "/eos/uscms/store/user/oamram/case/sig_files/LundRW/XToYYprimeTo4Q_MX3000_MY170_MYprime170_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER_Lund.h5"
fname = options.fin
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

sf_stat_unc = abs(toys_mean - eff_rw)/eff_nom + toys_std /eff_nom
sf_pt_unc = abs(pt_toys_mean - eff_rw)/eff_nom + pt_toys_std /eff_nom

sys_keys = ['sys', 'bquark', 'prongs', 'unclust' , 'distortion']
sys_uncs = dict()

for sys in sys_keys: sys_uncs[sys] = [0.,0.]

#Compute difference in efficiency due to weight variations as uncertainty
def get_uncs(score_cut, weights_up, weights_down, eff_baseline):
    eff_up =  np.average(score_cut, weights = weights_up)
    eff_down =  np.average(score_cut, weights = weights_down)

    unc_up = eff_up - eff_baseline
    unc_down = eff_down - eff_baseline 
    return unc_up, unc_down

for j,sys in enumerate(sys_keys):
    unc_up, unc_down = get_uncs(score_cut, weights_sys[:,j], weights_sys[:,j+1], eff_rw)
    sys_uncs[sys] = [unc_up/eff_nom, unc_down/eff_nom]


sf_nom = eff_rw / eff_nom

SF_str = "SF is %.2f +/- %.2f (stat) +/- %.2f (pt)" % (sf_nom, sf_stat_unc, sf_pt_unc )
tot_unc_up = tot_unc_down = sf_stat_unc**2 + sf_pt_unc**2

for sys in sys_keys:
    SF_str += " %.2f/%.2f (%s)" % (sys_uncs[sys][0], sys_uncs[sys][1], sys)
    up_var = max(sys_uncs[sys][0], sys_uncs[sys][1])
    down_var = min(sys_uncs[sys][0], sys_uncs[sys][1])
    tot_unc_up += up_var**2
    tot_unc_down += down_var**2

tot_unc_up = tot_unc_up**0.5
tot_unc_down = tot_unc_down**0.5

SF_str += "\n Overall %.2f +%.2f/-%.2f \n\n"  % (sf_nom, tot_unc_up, tot_unc_down)
print(SF_str)


#Make plots of different observables before vs after reweighting
#also include up/down shape variations from the sys variation 
weights = [  weights_nom, weights_rw, weights_sys[:,0], weights_sys[:,1],  weights_sys[:,4], weights_sys[:,5], weights_sys[:,6], weights_sys[:,7], weights_sys[:,8], weights_sys[:,9],]
labels = ["Nom", "RW", "RW Sys. Up", "RW Sys. Down", "RW prongs up", "RW prongs down", "RW unclust up", "RW unclust down", "RW distort up", "RW distort down"]
colors = ['gray','black', 'blue', 'blue', 'red', 'red', 'green', 'green', 'purple', 'purple']
ratio_range = [0.5, 1.5]


for l in obs:

    if(l == 'nPF'): 
        h_range = (0.5,120.5)
        n_bins_ = 40
    else: 
        n_bins_ = n_bins
        h_range = None

    x = getattr(d,l)[:max_evts]

    make_multi_ratio_histogram([x]*len(weights), labels, colors, l, "%s : Before vs. After Lund Plane Reweighting" % l, n_bins,
         ratio_range = ratio_range, normalize=True, weights = weights, fname=outdir + "%s_before_vs_after.png" %l )

