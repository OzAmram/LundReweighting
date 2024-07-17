import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.Utils import *
import os



parser = input_options()
options = parser.parse_args()

print(options)


outdir = options.outdir
if(not os.path.exists(outdir)): os.system("mkdir %s" % outdir)

#fname = "/eos/uscms/store/user/oamram/case/sig_files/WkkToWRadionToWWW_M3000_Mr400_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5"
#fname = "/eos/uscms/store/user/oamram/case/sig_files/YtoHH_Htott_Y3000_H400_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5"
#fname = "/eos/uscms/store/user/oamram/case/sig_files/ZpToTpTp_Zp5000_Tp400_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5"
fname = "/eos/uscms/store/user/oamram/case/sig_files/XToYYprimeTo4Q_MX3000_MY170_MYprime170_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5"
#fname = "test_signal_CASE.h5"



#label = "YtoHH"
#h_label = label
#label = "Wkk"
#h_label = label
#label = "ZpToTpTp"
#h_label = 'Zp'
label = "XtoYY"
h_label = 'XYY'

#tag_obs = 'tau43'
#score_thresh = 0.65
tag_obs = 'tau21'
score_thresh = 0.34


f_ratio = ROOT.TFile.Open(options.fin)
f_sig = h5py.File(fname, "r")
j_idx = 0

jetR = 1.0
num_excjets = -1

max_evts = 20000
#max_evts = 1000
#max_evts = None

d = Dataset(f_sig, label = label, color = ROOT.kRed, dtype = 1)
is_lep = f_sig['event_info'][:,4]
mjj = f_sig['jet_kinematics'][:,0]

print("input file", fname)

cut = (is_lep < 0.5)
print(np.mean(cut))

d.apply_cut(cut)
d.compute_obs()

obs = ["tau21", "tau32", "tau43", "nPF", "mSoftDrop", "pt"]
n_bins = 20

score = getattr(d, tag_obs)[:max_evts]

print("%i events" % score.shape[0])

score_cut = score < score_thresh

nom_weights = np.ones_like(score)

LP_rw = LundReweighter(f_ratio = f_ratio, charge_only = options.charge_only)

#Compute reweighting factors and all systematic variations
LP_weights = d.reweight_all(LP_rw, max_evts = max_evts)


#multiply Lund plane weights with nominal event weights
for key in LP_weights.keys():
    if('nom' in key or 'up' in key or 'down' in key):
        if(isinstance(LP_weights[key], np.ndarray)) : LP_weights[key] *= nom_weights


#Fraction of prongs that are not well matched to subjets (want this to be low)
print("Bad match frac %.2f" % np.mean(LP_weights['bad_match']))
#Fraction of prongs that are still not well matched after reclustering with varied number of prongs
print("Reclustered bad match frac %.2f" % np.mean(LP_weights['reclust_still_bad_match']))


make_histogram(LP_weights['nom'], "Reweighting factors", 'b', 'Weight', "Lund Plane Reweighting Factors", 40 , h_range = (0., LP_rw.max_rw + 0.1),
     normalize=False, fname=outdir + "lundPlane_weights.png")


frac_low_edge = np.mean( LP_weights['nom'] < np.amin(LP_weights['nom']) * 1.01)
frac_high_edge = np.mean( LP_weights['nom'] > np.amax(LP_weights['nom']) * 0.99)

print(np.amin(LP_weights['nom']), np.amax(LP_weights['nom']))

print("Frac of weights low edge %.3f, high edge %.3f" % (frac_low_edge, frac_high_edge))


#Efficiency of the cut in nominal MC
eff_nom = np.average(score_cut, weights = nom_weights)

#Efficiency of the cut after the Lund Plane reweighting
eff_rw = np.average(score_cut, weights = LP_weights['nom'])

#Nominal 'scale factor'
SF = eff_rw / eff_nom

print("Nominal efficiency %.3f, Corrected efficiency %.3f, SF (corrected / nom) %.3f" % (eff_nom, eff_rw, SF))


######  Compute uncertainties on the efficiency from the various weight variations ##############

#statistical and pt extrapolation uncertainties derived from 100 variations of the weights 
#take std dev to determine unc

nToys = LP_weights['stat_vars'].shape[1]
eff_toys = []
pt_eff_toys = []
for i in range(nToys):
    eff = np.average(score_cut, weights = LP_weights['stat_vars'][:,i])
    eff_toys.append(eff)

    eff1 = np.average(score_cut, weights = LP_weights['pt_vars'][:,i])
    pt_eff_toys.append(eff1)

#Compute stat and pt uncertainty based on variation in the toys
toys_mean = np.mean(eff_toys)
toys_std = np.std(eff_toys)
pt_toys_mean = np.mean(pt_eff_toys)
pt_toys_std = np.std(pt_eff_toys)

#if mean of toys is biased, also include it as an unc (should be zero)
stat_unc = (abs(toys_mean - eff_rw)  + toys_std) / eff_nom
pt_unc = (abs(pt_toys_mean - eff_rw) + pt_toys_std) / eff_nom

print("Stat variation toys eff. avg %.3f, std dev %.3f" % (toys_mean, toys_std))
print("Pt variation toys eff. avg %.3f, std dev %.3f" % (pt_toys_mean, pt_toys_std))

#Other systematics come from up/down variations of the weights
sys_keys = ['sys', 'bquark', 'prongs', 'unclust', 'distortion']
sys_uncs = dict()

for sys in sys_keys: sys_uncs[sys] = [0.,0.]

#Compute difference in efficiency due to weight variations as uncertainty
def get_uncs(cut, weights_up, weights_down, eff_baseline):
    eff_up =  np.average(cut, weights = weights_up)
    eff_down =  np.average(cut, weights = weights_down)

    unc_up = eff_up - eff_baseline
    unc_down = eff_down - eff_baseline 
    return unc_up, unc_down

for sys in sys_keys:
    unc_up, unc_down = get_uncs(score_cut, LP_weights[sys + '_up'], LP_weights[sys + '_down'], eff_rw)
    sys_uncs[sys] = [unc_up/eff_nom, unc_down/eff_nom]


eff_str = "Obs is %s, cut is %.2f \n" %  (tag_obs, score_thresh)
#Print uncertainty breakdown
eff_str += "Original %.2f, Calibrated %.2f \n"  % (eff_nom, eff_rw)
eff_str += "SF  is %.2f +/- %.2f (stat) +/- %.2f (pt)" % (SF, stat_unc, pt_unc )
tot_unc_up = tot_unc_down = stat_unc**2 + pt_unc**2

for sys in sys_keys:
    eff_str += " %.2f/%.2f (%s)" % (sys_uncs[sys][0], sys_uncs[sys][1], sys)
    up_var = max(sys_uncs[sys][0], sys_uncs[sys][1])
    down_var = min(sys_uncs[sys][0], sys_uncs[sys][1])
    tot_unc_up += up_var**2
    tot_unc_down += down_var**2



tot_unc_up = tot_unc_up**0.5
tot_unc_down = tot_unc_down**0.5

#Print final calibrated efficiency and total uncertaintiy
eff_str += "\n SF is %.2f +%.2f/-%.2f \n"  % (SF, tot_unc_up, tot_unc_down)

print(eff_str)
f_effs = open(options.outdir + "Effs.txt", "w")
f_effs.write(eff_str)
f_effs.close()

f_ratio.Close()



#Save subjet pts and deltaR
#deltaRs = np.reshape(deltaRs, -1)


num_bins = 40
pt_bins = array('d', np.linspace(0., 800., num_bins + 1))
dR_bins = array('d', np.linspace(0., 0.8, num_bins + 1))

h_subjet_pts = make_root_hist(data = LP_weights['subjet_pts'], name = 'h_%s_subjetpt' %h_label, num_bins = num_bins, bins = pt_bins)
f_ptout = ROOT.TFile.Open(outdir + "subjet_pt_dR.root", "RECREATE")
h_subjet_pts.Write()
#h_dRs = make_root_hist(data = deltaRs, name = 'h_%s_dRs' %h_label, num_bins = num_bins, bins = dR_bins)
#h_dRs.Write()
f_ptout.Close()

make_histogram([LP_weights['subjet_pts']], ["Subjets"], colors = ['blue'], xaxis_label = 'Subjet pt (GeV)', 
                title = "%s : subjet pt " % (label), num_bins = 40, normalize = True, fname = options.outdir + label + "_subjet_pt.png")

print(np.amin(LP_weights['subjet_pts']))
print("Fraction of subjets with pt > 350 : %.3f" % (np.mean( np.array(LP_weights['subjet_pts']) > 350.)))



#weights = [  weights_nom, weights_rw, sys_variations['sys_tot_up'], sys_variations['sys_tot_down']]
#labels = ["Nom", "RW", "RW Sys. Up", "RW Sys. Down"]
#colors = ['gray','black', 'blue', 'red']

weights = [ nom_weights, LP_weights['nom']]
labels = ["Nom", "RW"]
colors = ['gray','black']
ratio_range = [0.5, 1.5]


for l in obs:

    if(l == 'nPF'): 
        h_range = (0.5,120.5)
        n_bins_ = 40
    else: 
        n_bins_ = n_bins
        h_range = None

    x = getattr(d,l)[:max_evts]

    print(x.shape, weights[0].shape)

    make_multi_ratio_histogram([x]*2, labels, colors, l, "%s : Before vs. After Lund Plane Reweighting" % l, n_bins,
         ratio_range = ratio_range, normalize=True, weights = weights, fname=outdir + "%s_before_vs_after.png" %l )

