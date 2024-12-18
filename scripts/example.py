import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.Utils import *
""" An example how to use the Lund Plane reweighting  code. 
Computes reweighting factors for an example signal. 
Uses them to compute the efficiency and uncertainty of a given substructure cut"""

parser = input_options()
options = parser.parse_args()

######################## Setup 

#Input file 
fname = "data/example_signal.h5"
f_ratio_name = 'data/ratio_2018.root'

f_sig = h5py.File(fname, "r")
f_ratio = ROOT.TFile.Open(f_ratio_name)

#Class to help read input dataset 
#NOTE using this class or data format is not required for using the method
d = Dataset(f_sig, dtype = 1)
d.compute_obs()

#The cut we will compute a SF for 'tau21 < 0.34'
tag_obs = 'tau21'
score_thresh = 0.34



#Main class for reweighting utilities
LP_rw = LundReweighter(f_ratio = f_ratio)

#use small number of events for the example, but generally advisable to use larger sample
max_evts = 5000
score = getattr(d, tag_obs)[:max_evts]
score_cut = score < score_thresh

################### Compute reweighting factors

#PF candidates in the AK8 jet for each event (list of (px,py,pz,E))
pf_cands = d.get_masked("jet1_PFCands").astype(np.float64)[:max_evts]
print("Running on %i events" % len(pf_cands))

#Generator level quarks from hard process 
gen_parts = d.get_masked('gen_info')[:max_evts]
gen_parts_eta_phi = gen_parts[:,:,1:3]
gen_parts_pdg_ids = gen_parts[:,:,3]

#4 vector of AK8 jets we are calibrating
ak8_jets = d.get_masked('jet_kinematics')[:max_evts][:,2:6].astype(np.float64)

#Nominal event weights of the MC, (assume every event is weight '1' for this example)
nom_weights = np.ones(len(pf_cands))

#Compute reweighting factors and all systematic variations
LP_weights = LP_rw.get_all_weights(pf_cands, gen_parts_eta_phi, ak8_jets, gen_parts_pdg_ids = gen_parts_pdg_ids)


#multiply Lund plane weights with nominal event weights
for key in LP_weights.keys():
    if('nom' in key or 'up' in key or 'down' in key):
        if(isinstance(LP_weights[key], np.ndarray)) : LP_weights[key] *= nom_weights


#Fraction of prongs that are not well matched to subjets (want this to be low)
print("Bad match frac %.2f" % np.mean(LP_weights['bad_match']))
#Fraction of prongs that are still not well matched after reclustering with varied number of prongs
print("Reclustered bad match frac %.2f" % np.mean(LP_weights['reclust_still_bad_match']))


###### Use weights to compute efficiency of a cut

#Efficiency of the cut in nominal MC
eff_nom = np.average(score_cut, weights = nom_weights)

#Efficiency of the cut after the Lund Plane reweighting
eff_rw = np.average(score_cut, weights = LP_weights['nom'])

#Nominal 'scale factor'
SF = eff_rw / eff_nom

print("Nominal efficiency %.3f, Corrected efficiency %.3f, SF (corrected / nom) %.3f" % (eff_nom, eff_rw, SF))

#NOTE, because there is kinematic dependence to the correction, it is better to use corrected efficiency computed 
#separately for each MC sample rather than a single 'SF'


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
eff_stat_unc = (abs(toys_mean - eff_rw)  + toys_std) 
eff_pt_unc = (abs(pt_toys_mean - eff_rw) + pt_toys_std)

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
    sys_uncs[sys] = [unc_up, unc_down]


#Print uncertainty breakdown
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

#Print final calibrated efficiency and total uncertaintiy
eff_str += "\n Original %.2f, Calibrated %.2f +%.2f/-%.2f \n"  % (eff_nom, eff_rw, tot_unc_up, tot_unc_down)

print(eff_str)
f_ratio.Close()
