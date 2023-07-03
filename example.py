from Utils import *
import os
""" An example how to use the Lund Plane reweighting  code """

parser = input_options()
options = parser.parse_args()

######################## Setup 

#Input file 
fname = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/TagNTrain/data/YtoHH_Htott_Y3000_H400_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5"
f_sig = h5py.File(fname, "r")
#Class to help read input dataset 
d = Dataset(f_sig, dtype = 1)
d.compute_obs()

#The cut we will compute a SF for 'tau43 < 0.65'
tag_obs = 'tau43'
score_thresh = 0.65

#File containing data/MC Lund Plane ratio
f_ratio = ROOT.TFile.Open(options.fin)

#nominal data/MC Lund plane ratio (3d histogram)
h_ratio = f_ratio.Get("ratio_nom")
#systematic variations
h_ratio_sys_up = f_ratio.Get("ratio_sys_tot_up")
h_ratio_sys_down = f_ratio.Get("ratio_sys_tot_down")

#directory of pt extrapolation fits
f_ratio.cd('pt_extrap')
rdir = ROOT.gDirectory

#Main class for reweighting utilities
LP_rw = LundReweighter(pt_extrap_dir = rdir)

#How many subjets to recluster into, 
#This should be obtained for each event based on a gen-level match. But just use a fixed quantity for now
num_excjets = 6

max_evts = 2000
score = getattr(d, tag_obs)[:max_evts]
score_cut = score < score_thresh


#Number of toys for statistical and pt extrapolation uncertainties
nToys = 100
#Noise vectors used to to generated generate the toys
#NOTE the same vector has to be used for the whole sample/signal file for the toys to be consistent 
rand_noise = np.random.normal(size = (nToys, h_ratio.GetNbinsX(), h_ratio.GetNbinsY(), h_ratio.GetNbinsZ()))
pt_rand_noise = np.random.normal(size = (nToys, h_ratio.GetNbinsY(), h_ratio.GetNbinsZ(), 3))


################### Compute reweighting factors
pf_cands = d.get_masked("jet1_PFCands").astype(np.float64)[:max_evts]

#Nominal event weights of the MC, assume every event is weight '1' for this example
weights_nom = np.ones(max_evts)

LP_weights = []
LP_weights_sys_up = []
LP_weights_sys_down = []
stat_smeared_weights = []
pt_smeared_weights = []

for cands in pf_cands:
    #Get the subjets and splittings for this jet from the PF candidates
    subjets, splittings = LP_rw.get_splittings(cands, num_excjets = num_excjets)

    #First call gets the nominal LP reweighting factor and statistical + pt extrapolation toys
    LP_weight, stat_smeared_weight, pt_smeared_weight = LP_rw.reweight_lund_plane(h_rw = h_ratio, subjets = subjets, splittings = splittings,
            rand_noise = rand_noise, pt_rand_noise = pt_rand_noise, )
    #Now get systematic variations
    LP_weight_sys_up,_,_ = LP_rw.reweight_lund_plane(h_rw = h_ratio_sys_up, subjets = subjets, splittings = splittings)
    LP_weight_sys_down,_,_ = LP_rw.reweight_lund_plane(h_rw = h_ratio_sys_down, subjets = subjets, splittings = splittings)


    LP_weights.append(LP_weight)
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



############### Compute efficiences and an example 'scale factor'


#Efficiency of the cut in nominal MC
eff_nom = np.average(score_cut, weights = weights_nom)

#Efficiency of the cut after the Lund Plane reweighting
eff_rw = np.average(score_cut, weights = LP_weights)

#Nominal 'scale factor'
SF = eff_rw / eff_nom

print("Nom %.3f, RW %.3f, SF %.3f" % (eff_nom, eff_rw, SF))


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

SF_stat_unc = (abs(toys_mean - eff_rw)  + toys_std) /eff_nom
SF_pt_unc = (abs(pt_toys_mean - eff_rw) + pt_toys_std) /eff_nom

print("Stat variation toys avg %.3f, std dev %.3f" % (toys_mean, toys_std))
print("Pt variation toys avg %.3f, std dev %.3f" % (pt_toys_mean, pt_toys_std))

#Compute efficiency of systematic variations
eff_sys_up =  np.average(score_cut, weights = LP_weights_sys_up)
eff_sys_down =  np.average(score_cut, weights = LP_weights_sys_up)

sys_unc_up = abs(eff_rw - eff_sys_up)
sys_unc_down = abs(eff_rw - eff_sys_down)


############ Results
print("\n\nSF  is %.2f +/- %.2f  (stat) +/- %.2f (pt) +%.2f/-%.2f (sys)  \n\n"  % (SF, SF_stat_unc, SF_pt_unc, sys_unc_up, sys_unc_down, ))
f_ratio.Close()
#NOTE This example does not include the subjet matching uncertainty which is generally dominant, see the AN / documentation for how to compute it
