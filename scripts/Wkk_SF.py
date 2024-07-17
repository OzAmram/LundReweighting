import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.Utils import *


parser = input_options()
options = parser.parse_args()

print(options)

if(not os.path.exists(options.outdir)): os.system("mkdir %s" % options.outdir)
jet_str = 'CA'

#f_sig = h5py.File("/eos/uscms/store/user/oamram/case/sig_files/WkkToWRadionToWWW_M3000_Mr170_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5", "r")
f_sig = h5py.File("/uscms_data/d3/oamram/CMSSW_12_4_0/src/CASE/CASEUtils/H5_maker/Wkk_M3000_R170.h5", "r")
f_ratio = ROOT.TFile.Open(options.fin)
label = "Radion"
n_prongs = (4,2)
sig_mass = 3000.


j_idx = 0


jetR = 1.0
num_excjets = -1

max_evts = 10000
#max_evts = 10

d = Dataset(f_sig, label = label, color = ROOT.kRed, dtype = 1)
is_lep = f_sig['event_info'][:,4]
mjj = f_sig['jet_kinematics'][:,0]



j1_m = f_sig['jet_kinematics'][:,5]
j2_m = f_sig['jet_kinematics'][:,9]


j1_pt = f_sig['jet_kinematics'][:,2]
j2_pt = f_sig['jet_kinematics'][:,6]

max_pt = np.maximum(j1_pt, j2_pt)
min_pt = np.minimum(j1_pt, j2_pt)

#cut = (is_lep < 0.5)
cut = (is_lep < 0.5) & (mjj > 0.8*sig_mass) & (mjj < 1.2*sig_mass)
#cut = (mjj > 0.8*sig_mass) & (mjj < 1.2*sig_mass)
print(np.mean(cut))
#jet mass window
cut = cut & (j1_m > 70.) & ( (j2_m > 70.) & (j2_m < 100.))
print(np.mean(cut))
cut = cut & (max_pt > 400.) & (min_pt > 200.)
print(np.mean(cut))

d.apply_cut(cut)

print(f_sig['jet1_extraInfo'].shape)
WH_score = f_sig['jet1_extraInfo'][:,8][cut]
print(WH_score.shape)
print(WH_score[:10])

score_cut = WH_score > 0.8

#use gen weights ? 
weights_nom = np.ones_like(WH_score)

LP_rw = LundReweighter(f_ratio = f_ratio, charge_only = options.charge_only)

LP_weights = d.reweight_all(LP_rw)

weights_rw = LP_weights['nom']

print("Bad match frac %.2f" % np.mean(LP_weights['bad_match']))
#Fraction of prongs that are still not well matched after reclustering with varied number of prongs
print("Reclustered bad match frac %.2f" % np.mean(LP_weights['reclust_still_bad_match']))


#compute 'Scalefactor'

eff_nom = np.average(score_cut, weights = weights_nom)
eff_rw = np.average(score_cut, weights = weights_rw)

sf_nom = eff_rw/eff_nom

print("Nom %.3f, RW %.3f" % (eff_nom, eff_rw))


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

eff_stat_unc = (abs(toys_mean - eff_rw)  + toys_std) /eff_nom
eff_pt_unc = (abs(pt_toys_mean - eff_rw) + pt_toys_std) /eff_nom

print("Stat variation toys eff. avg %.3f, std dev %.3f" % (toys_mean, toys_std))
print("Pt variation toys eff. avg %.3f, std dev %.3f" % (pt_toys_mean, pt_toys_std))

#Add systematic differences in quadrature
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

for sys in sys_keys:
    unc_up, unc_down = get_uncs(score_cut, LP_weights[sys + '_up'], LP_weights[sys + '_down'], eff_rw)
    print(sys, unc_up, unc_down)
    sys_uncs[sys] = [unc_up/eff_nom, unc_down/eff_nom]



SF_str = "SF (WH_Score > 0.8) is %.2f +/- %.2f (stat) +/- %.2f (pt)" % (sf_nom, eff_stat_unc, eff_pt_unc )
tot_unc_up = tot_unc_down = eff_stat_unc**2 + eff_pt_unc**2

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
f_SFs = open(options.outdir + "SFs.txt", "w")
f_SFs.write(SF_str)
f_SFs.close()

f_ratio.Close()
