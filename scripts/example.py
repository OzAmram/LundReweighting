import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.Utils import *
""" An example how to use the Lund Plane reweighting  code """

parser = input_options()
options = parser.parse_args()

######################## Setup 

#Input file 
fname = "data/example_signal.h5"
#File containing data/MC Lund Plane ratio
f_ratio_name = 'data/ratio_2018.root'

f_sig = h5py.File(fname, "r")
f_ratio = ROOT.TFile.Open(f_ratio_name)

#Class to help read input dataset 
d = Dataset(f_sig, dtype = 1)
d.compute_obs()

#The cut we will compute a SF for 'tau21 < 0.34'
tag_obs = 'tau21'
score_thresh = 0.34


#nominal data/MC Lund plane ratio (3d histogram)
h_ratio = f_ratio.Get("ratio_nom")
#systematic variations
h_ratio_sys_up = f_ratio.Get("ratio_sys_tot_up")
h_ratio_sys_down = f_ratio.Get("ratio_sys_tot_down")
#MC ratio of b to light quarks
b_light_ratio = f_ratio.Get("h_bl_ratio")


#directory of pt extrapolation fits
f_ratio.cd('pt_extrap')
rdir = ROOT.gDirectory

#Main class for reweighting utilities
LP_rw = LundReweighter(pt_extrap_dir = rdir)

max_evts = 1000
score = getattr(d, tag_obs)[:max_evts]
score_cut = score < score_thresh


#Number of toys for statistical and pt extrapolation uncertainties
nToys = 100
#Noise vectors used to to generated generate the toys
#NOTE the same vector has to be used for the whole sample/signal file for the toys to be consistent 
rand_noise = np.random.normal(size = (nToys, h_ratio.GetNbinsX(), h_ratio.GetNbinsY(), h_ratio.GetNbinsZ()))
pt_rand_noise = np.random.normal(size = (nToys, h_ratio.GetNbinsY(), h_ratio.GetNbinsZ(), 3))


################### Compute reweighting factors

#PF candidates in the AK8 jet
pf_cands = d.get_masked("jet1_PFCands").astype(np.float64)[:max_evts]

#Generator level quarks from hard process
gen_parts = d.get_masked('gen_info')[:max_evts]
gen_parts_eta_phi = gen_parts[:,:,1:3]
gen_parts_pdg_ids = gen_parts[:,:,3]

B_PDG_ID = 5

ak8_jets = d.get_masked('jet_kinematics')[:max_evts][:,2:6].astype(np.float64)

#Nominal event weights of the MC, assume every event is weight '1' for this example
weights_nom = np.ones(max_evts)

LP_weights = []
LP_weights_sys_up = []
LP_weights_sys_down = []
stat_smeared_weights = []
pt_smeared_weights = []
b_weights_up = []
b_weights_down = []
bad_matches = []


for i,cands in enumerate(pf_cands):

    #Get the subjets, splittings and checking matching based on PF candidates in the jet and gen-level quarks
    subjets, splittings, bad_match, deltaRs = LP_rw.get_splittings_and_matching(cands, gen_parts_eta_phi[i], ak8_jets[i])

    #Gets the nominal LP reweighting factor for this event and statistical + pt extrapolation toys
    LP_weight, stat_smeared_weight, pt_smeared_weight = LP_rw.reweight_lund_plane(h_rw = h_ratio, subjets = subjets, splittings = splittings,
            rand_noise = rand_noise, pt_rand_noise = pt_rand_noise, )
    #Now get systematic variations
    LP_weight_sys_up,_,_ = LP_rw.reweight_lund_plane(h_rw = h_ratio_sys_up, subjets = subjets, splittings = splittings)
    LP_weight_sys_down,_,_ = LP_rw.reweight_lund_plane(h_rw = h_ratio_sys_down, subjets = subjets, splittings = splittings)


    #compute special systematic for subjets matched to b quarks
    #not needed if signal does not specifically produce b quark subjets
    gen_bs = [j for j in range(len(gen_parts_pdg_ids[i])) if abs(gen_parts_pdg_ids[i,j]) == B_PDG_ID]

    if(len(gen_bs) == 0): b_rw = 1.0
    else:
        dists = get_subjet_dist(gen_parts_eta_phi[i,gen_bs,1:3], np.array(subjets)[:,1:3])

        b_matches = []
        #which subjet is each quark closest to
        j_closest = np.amin(dists, axis = 0)
        j_which = np.argmin(dists, axis = 0)
        b_matches = np.unique(j_which[j_closest < deltaR_cut])

        #reweight only subjets matched to b quarks
        if(len(b_matches) > 0):
            b_subjet = [subjet[j] for j in range(len(subjet)) if j in b_matches]
            b_split  = [split[j]  for j in range(len(split)) if split[j][0] in b_matches]

            b_rw, _,_   = LP_rw.reweight_lund_plane(b_light_ratio, subjets = b_subjet, splittings = b_split, sys_str = 'bquark')

        else: b_rw = 1.0

    b_weights_up.append(LP_weight * b_rw)
    b_weights_down.append(LP_weight / b_rw)


    LP_weights.append(LP_weight)
    stat_smeared_weights.append(stat_smeared_weight)
    pt_smeared_weights.append(pt_smeared_weight)

    LP_weights_sys_up.append(LP_weight_sys_up)
    LP_weights_sys_down.append(LP_weight_sys_down)
    bad_matches.append(bad_match)



############### Normalize weights to preserve normalization of the MC sample

#The nominal Lund Plane correction event weights
LP_weights = LP_rw.normalize_weights(LP_weights) * weights_nom 

#Toy variations for stat and pt uncertainties
stat_smeared_weights = LP_rw.normalize_weights(stat_smeared_weights) * weights_nom.reshape(max_evts, 1)
pt_smeared_weights = LP_rw.normalize_weights(pt_smeared_weights) * weights_nom.reshape(max_evts,1)

#Systematic up/down variations
LP_weights_sys_up = LP_rw.normalize_weights(LP_weights_sys_up) * weights_nom
LP_weights_sys_down = LP_rw.normalize_weights(LP_weights_sys_down) * weights_nom


#b quark systematic up/down variations
b_weights_up = LP_rw.normalize_weights(b_weights_up) * weights_nom
b_weights_down = LP_rw.normalize_weights(b_weights_down) * weights_nom



############### Compute efficiences and uncertainties


#Efficiency of the cut in nominal MC
eff_nom = np.average(score_cut, weights = weights_nom)

#Efficiency of the cut after the Lund Plane reweighting
eff_rw = np.average(score_cut, weights = LP_weights)

#Nominal 'scale factor'
SF = eff_rw / eff_nom

print("Nominal efficiency %.3f, Corrected efficiency %.3f, SF (corrected / nom) %.3f" % (eff_nom, eff_rw, SF))

#NOTE, better to use corrected efficiency computed separately for each sample rather than a single 'SF'


#Compute efficiency for each of the stat/pt toys
eff_toys = []
pt_eff_toys = []
for i in range(nToys):
    eff = np.average(score_cut, weights = stat_smeared_weights[:,i])
    eff_toys.append(eff)

    eff1 = np.average(score_cut, weights = pt_smeared_weights[:,i])
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


#Compute difference in efficiency due to weight variations as uncertainty
def get_uncs(score_cut, weights_up, weights_down, eff_baseline):
    eff_up =  np.average(score_cut, weights = weights_up)
    eff_down =  np.average(score_cut, weights = weights_down)

    unc_up = eff_up - eff_baseline
    unc_down = eff_down - eff_baseline 
    return unc_up, unc_down


#Compute efficiency of systematic variations
sys_unc_up, sys_unc_down = get_uncs(score_cut, LP_weights_sys_up, LP_weights_sys_down, eff_rw)
b_unc_up, b_unc_down = get_uncs(score_cut, b_weights_up, b_weights_down, eff_rw)


#matching uncertainty, taken as a fractional uncertainty on efficiency
bad_match_frac = np.mean(bad_matches)
bad_match_unc = bad_match_frac * eff_rw


############ Results
print("\n\nCalibrated efficiency  is %.2f +/- %.2f  (stat) +/- %.2f (pt) %.2f/%.2f (sys) %.2f/%.2f (bquark) +/- %.2f (matching)  \n\n"  % 
        (eff_rw, eff_stat_unc, eff_pt_unc, sys_unc_up, sys_unc_down, b_unc_up, b_unc_down, bad_match_unc))
f_ratio.Close()
