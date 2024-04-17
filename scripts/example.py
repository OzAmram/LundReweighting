import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.Utils import *
""" An example how to use the Lund Plane reweighting  code """

parser = input_options()
options = parser.parse_args()

######################## Setup 

#Input file 
#fname = "data/example_signal.h5"
fname = "../TagNTrain/data/LundRW/YtoHH_Htott_Y5000_H400_TuneCP5_13TeV-madgraph-pythia8_TIMBER_Lund.h5"
#File containing data/MC Lund Plane ratio
f_ratio_name = 'data/ratio_2018.root'

f_sig = h5py.File(fname, "r")
f_ratio = ROOT.TFile.Open(f_ratio_name)

#Class to help read input dataset 
d = Dataset(f_sig, dtype = 1)
d.compute_obs()

#The cut we will compute a SF for 'tau21 < 0.34'
#tag_obs = 'tau21'
#score_thresh = 0.34
tag_obs = 'tau43'
score_thresh = 0.65



#Main class for reweighting utilities
LP_rw = LundReweighter(f_ratio)

max_evts = 10000
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
LP_weights_prong_up = []
LP_weights_prong_down = []
LP_weights_matching_up = []
LP_weights_matching_down = []
stat_smeared_weights = []
pt_smeared_weights = []
b_weights_up = []
b_weights_down = []

partial_merge = 0
bad_matches = 0


for i,cands in enumerate(pf_cands):

    #Get the subjets and splittings  based on PF candidates in the jet and gen-level quarks
    #Also compute variations if needed by changing the number of subjets in the reclustering
    reclust_nom, reclust_prongs_up, reclust_prongs_down = LP_rw.get_splittings_and_matching(cands, gen_parts_eta_phi[i], ak8_jets[i])

    if(reclust_nom.badmatch): bad_matches +=1
    if((reclust_prongs_up is not None and reclust_prongs_up.from_prongs_up) or 
       (reclust_prongs_down is not None and reclust_prongs_down.from_prongs_down)): partial_merge +=1

    #Gets the nominal LP reweighting factor for this event and statistical + pt extrapolation toys
    LP_weight, stat_smeared_weight, pt_smeared_weight = LP_rw.reweight_lund_plane(h_rw = h_ratio, reclust_obj = reclust_nom, rand_noise = rand_noise, pt_rand_noise = pt_rand_noise, )

    LP_weight_prong_up, LP_weight_prong_down, LP_weight_match_up, LP_weight_match_down = LP_rw.get_up_down_prongs_weights(h_rw = h_ratio, 
            reclust_prongs_up = reclust_prongs_up, reclust_prongs_down = reclust_prongs_down, nom_weights = LP_weight)


    #Now get systematic variations due to systemtatic uncertainties on LP
    LP_weight_sys_up,_,_ = LP_rw.reweight_lund_plane(h_rw = h_ratio_sys_up, reclust_obj = reclust_nom)
    LP_weight_sys_down,_,_ = LP_rw.reweight_lund_plane(h_rw = h_ratio_sys_down, reclust_obj = reclust_nom)


    #compute special systematic for subjets matched to b quarks
    #not needed if signal does not specifically produce b quark subjets
    gen_bs = [j for j in range(len(gen_parts_pdg_ids[i])) if abs(gen_parts_pdg_ids[i,j]) == B_PDG_ID]

    if(len(gen_bs) == 0): b_rw = 1.0
    else:
        dists = get_subjet_dist(gen_parts_eta_phi[i,gen_bs,:], np.array(reclust_nom.subjet)[:,1:3])

        deltaR_cut = 0.2
        b_matches = []
        #which subjet is each quark closest to
        j_closest = np.amin(dists, axis = 0)
        j_which = np.argmin(dists, axis = 0)
        b_matches = np.unique(j_which[j_closest < deltaR_cut])

        #reweight only subjets matched to b quarks
        if(len(b_matches) > 0):
            b_subjet = [reclust_nom.subjet[j] for j in range(len(reclust_nom.subjet)) if j in b_matches]
            b_split  = [reclust_nom.split[j]  for j in range(len(reclust_nom.split)) if reclust_nom.split[j][0] in b_matches]

            b_rw, _,_   = LP_rw.reweight_lund_plane(b_light_ratio, subjets = b_subjet, splittings = b_split, sys_str = 'bquark')

        else: b_rw = 1.0

    b_weights_up.append(LP_weight * b_rw)
    b_weights_down.append(LP_weight / b_rw)

    LP_weights.append(LP_weight)

    LP_weights_prong_up.append(LP_weight_prong_up)
    LP_weights_prong_down.append(LP_weight_prong_down)

    LP_weights_match_up.append(LP_weight_match_up)
    LP_weights_match_down.append(LP_weight_match_down)


    stat_smeared_weights.append(stat_smeared_weight)
    pt_smeared_weights.append(pt_smeared_weight)

    LP_weights_sys_up.append(LP_weight_sys_up)
    LP_weights_sys_down.append(LP_weight_sys_down)





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



#TODO
#print(partial_merge, bad_matches)
#partial_merge_frac = float(partial_merge) / len(LP_weights)
partial_merge_frac = 99999999
bad_match_frac = np.mean(reclust

print("Fraction of jets with a partially merged prong is %.3f and fraction with bad matches %.3f" % (partial_merge_frac, bad_match_frac))
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
sys_LPunc_up, sys_LPunc_down = get_uncs(score_cut, LP_weights_sys_up, LP_weights_sys_down, eff_rw)
sys_merge_up, sys_merge_down = get_uncs(score_cut, LP_weights_merge_up, LP_weights_merge_down, eff_rw)
sys_matching_up, sys_matching_down = get_uncs(score_cut, LP_weights_match_up, LP_weights_match_down, eff_rw)

b_unc_up, b_unc_down = get_uncs(score_cut, b_weights_up, b_weights_down, eff_rw)



############ Results
print("\n\nCalibrated efficiency  is %.2f +/- %.2f  (stat) +/- %.2f (pt) %.2f/%.2f (sys) %.2f/%.2f (bquark) %.2f/%.2f (merging)  %.2f/%.2f (matching) \n\n"  % 
        (eff_rw, eff_stat_unc, eff_pt_unc, sys_LPunc_up, sys_LPunc_down, b_unc_up, b_unc_down, sys_prongs_up, sys_prongs_down, sys_matching_up, sys_matching_down ))
f_ratio.Close()
