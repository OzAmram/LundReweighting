from Utils import *
import os



parser = input_options()
options = parser.parse_args()

print(options)


outdir = options.outdir
if(not os.path.exists(outdir)): os.system("mkdir %s" % outdir)
jet_str = 'CA'

fname = "test_signal_CASE.h5"


#label = "YtoHH"
label = "XtoYY"
#label = "ZpToTpTp"

#tag_obs = 'tau43'
#score_thresh = 0.65
tag_obs = 'tau21'
score_thresh = 0.34


f_ratio = ROOT.TFile.Open(options.fin)
f_sig = h5py.File(fname, "r")
j_idx = 0

jetR = 1.0
num_excjets = -1

#max_evts = 5000
max_evts = 2000
#max_evts = None

d = Dataset(f_sig, label = label, color = ROOT.kRed, dtype = 1)
is_lep = f_sig['event_info'][:,4]
mjj = f_sig['jet_kinematics'][:,0]

print("input file", fname)

j1_m = f_sig['jet_kinematics'][:,5]
j2_m = f_sig['jet_kinematics'][:,9]


j1_pt = f_sig['jet_kinematics'][:,2]
j2_pt = f_sig['jet_kinematics'][:,6]

max_pt = np.maximum(j1_pt, j2_pt)
min_pt = np.minimum(j1_pt, j2_pt)

cut = (is_lep < 0.5)
print(np.mean(cut))

d.apply_cut(cut)
d.compute_obs()

obs = ["tau21", "tau32", "tau43", "nPF", "mSoftDrop", "pt"]
n_bins = 20

score = getattr(d, tag_obs)[:max_evts]

print("%i events" % score.shape[0])




score_cut = score < score_thresh

weights_nom = np.ones_like(score)


weights_rw = copy.deepcopy(weights_nom)


weights_rw = d.f['lund_weights'][:max_evts]
weights_stat = d.f['lund_weights_stat_var'][:max_evts]
weights_pt = d.f['lund_weights_pt_var'][:max_evts]
weights_sys = d.f['lund_weights_sys_var'][:max_evts] * np.expand_dims(weights_rw, -1)


#compute 'Scalefactor'

eff_nom = np.average(score_cut, weights = weights_nom)
eff_rw = np.average(score_cut, weights = weights_rw)

print("Nom %.3f, RW %.3f" % (eff_nom, eff_rw))


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

print("Toys avg %.3f, std dev %.3f" % (toys_mean, toys_std))


pt_toys_mean = np.mean(pt_eff_toys)
pt_toys_std = np.std(pt_eff_toys)

print("Pt variation toys avg %.3f, std dev %.3f" % (pt_toys_mean, pt_toys_std))


make_histogram(weights_rw, "Reweighting factors", 'b', 'Weight', "Lund Plane Reweighting Factors", 20 , h_range = (0., 2.0),
     normalize=False, fname=outdir + "lundPlane_weights.png")

#Add systematic differences in quadrature
sys_unc_up = sys_unc_down = 0.
if(not options.no_sys):

    eff_sys_tot_up = np.average(score_cut, weights = weights_sys[:,0])
    eff_sys_tot_down = np.average(score_cut, weights = weights_sys[:,1])
    SF_sys_unc_up = abs(eff_sys_tot_up - eff_rw)/eff_nom
    SF_sys_unc_down = abs(eff_sys_tot_down - eff_rw)/eff_nom
    SF_sys_unc = (SF_sys_unc_up + SF_sys_unc_down) / 2.0

    eff_bquark_up = np.average(score_cut, weights = weights_sys[:,2])
    eff_bquark_down = np.average(score_cut, weights = weights_sys[:,3])
    SF_bquark_up = abs(eff_bquark_up - eff_rw)/eff_nom
    SF_bquark_down = abs(eff_bquark_down - eff_rw)/eff_nom
    SF_bquark_unc = (SF_bquark_up + SF_bquark_down) /2.0

else: 
    SF_sys_unc = SF_sys_unc_up = SF_sys_unc_down = bquark_unc = 0. 



SF = eff_rw / eff_nom
SF_stat_unc = abs(toys_mean - eff_rw)/eff_nom + toys_std /eff_nom
SF_pt_unc = abs(pt_toys_mean - eff_rw)/eff_nom + pt_toys_std /eff_nom

#fraction of evts with bad match, take as fractional unc on SF
bad_matching_unc = d.f['lund_weights_matching_unc'][0] * SF

print("\n\nSF (%s val %.2f ) is %.2f +/- %.2f  (stat) +/- %.2f (sys) +/- %.2f (pt) +/- %.2f (bquark) +/- %.2f (matching) \n\n"  
        % (tag_obs, score_thresh, SF, SF_stat_unc, SF_sys_unc, SF_pt_unc, SF_bquark_unc, bad_matching_unc))

#approximate uncertainty on the reweighting for the plots
overall_unc = (SF_stat_unc **2 + SF_sys_unc**2 + SF_pt_unc**2 + SF_bquark_unc**2 + bad_matching_unc**2) **0.5 
print("overall unc %.3f" % overall_unc)



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

