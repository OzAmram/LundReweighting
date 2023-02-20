from Utils import *
import os

def ang_dist(phi1, phi2):
    phi1 = phi1 % (2. * np.pi)
    phi2 = phi2 % (2. * np.pi)
    dphi = phi1 - phi2
    dphi[dphi < -np.pi] += 2.*np.pi
    dphi[dphi > np.pi] -= 2.*np.pi
    return dphi

def get_dists(q_eta_phis, subjets_eta_phis):
    q_eta_phis = np.expand_dims(q_eta_phis, 2)
    subjets_eta_phis = np.expand_dims(subjets_eta_phis, 1)
    #print(q_eta_phis.shape)
    #print(subjets_eta_phis.shape)
    return np.sqrt(np.square(subjets_eta_phis[:,:,:,0] - q_eta_phis[:,:,:,0]) + 
            np.square(ang_dist(subjets_eta_phis[:,:,:,1], q_eta_phis[:,:,:,1] )))




outdir = "Wkk_SF_jan3/"
if(not os.path.exists(outdir)): os.system("mkdir %s" % outdir)
jet_str = 'CA'

f_sig = h5py.File("/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/CASE_signals/Wkk_M3000_R170_2018_UL.h5", "r")
#f_sig = h5py.File("/uscms_data/d3/oamram/CMSSW_12_4_0/src/CASE/CASEUtils/H5_maker/Wkk_M3500_test.h5", "r")
f_ratio = ROOT.TFile.Open("ttbar_UL_nov7_W_rw_kt_sys/ratio_test.root")
label = "Radion"
n_prongs = (4,2)
sig_mass = 3000.


j_idx = 0
do_sys_variations = True

subjet_rw = False
excjet_rw = True

jetR = 0.4
num_excjets = 4
fill_z = False

max_evts = 100000

d = Dataset(f_sig, label = label, color = ROOT.kRed)
is_lep = f_sig['event_info'][:,4]
mjj = f_sig['jet_kinematics'][:,0]

print(mjj[:10])
print(is_lep[:10])



j1_m = f_sig['jet_kinematics'][:,5]
j2_m = f_sig['jet_kinematics'][:,9]


j1_pt = f_sig['jet_kinematics'][:,2]
j2_pt = f_sig['jet_kinematics'][:,6]

max_pt = np.maximum(j1_pt, j2_pt)
min_pt = np.minimum(j1_pt, j2_pt)

cut = (is_lep < 0.5)
#cut = (is_lep < 0.5) & (mjj > 0.8*sig_mass) & (mjj < 1.2*sig_mass)
#cut = (mjj > 0.8*sig_mass) & (mjj < 1.2*sig_mass)
print(np.mean(cut))
#jet mass window
cut = cut & (j1_m > 70.) & ( (j2_m > 70.) & (j2_m < 100.))
print(np.mean(cut))
cut = cut & (max_pt > 400.) & (min_pt > 200.)
print(np.mean(cut))

d.apply_cut(cut)

WH_score = f_sig['jet1_extraInfo'][:,8][cut]
print(WH_score.shape)
print(WH_score[:10])

WH_cut = WH_score > 0.8

#use gen weights ? 
weights_nom = np.ones_like(WH_score)


weights_rw = copy.deepcopy(weights_nom)

h_ratio = f_ratio.Get("ratio_nom")

nToys = 100

#Noise used to generated smeared ratio's based on stat unc
rand_noise = np.random.normal(size = (h_ratio.GetNbinsX(), h_ratio.GetNbinsY(), h_ratio.GetNbinsZ(), nToys))

d_LP_weights, d_LP_uncs, d_LP_smeared_weights = d.reweight_LP(h_ratio, subjet_rw = subjet_rw, fill_z = fill_z, jetR = jetR, num_excjets = num_excjets, uncs = True, prefix = "4prong", 
        rand_noise = rand_noise)

LP_weights = d_LP_weights


#apply weights, keep normalization fixed
old_norm = np.sum(weights_rw)
weights_rw *= d_LP_weights

new_norm = np.sum(weights_rw)

weights_rw *= old_norm / new_norm
LP_smeared_weights = np.array(d_LP_smeared_weights * np.expand_dims(weights_nom, -1) * (old_norm / new_norm))

print("smeared", LP_smeared_weights.shape)


sys_variations = dict()
if(do_sys_variations):
    #sys_list = list(sys_weights_map.keys())
    sys_list = []
    sys_list.append("sys_tot_up")
    sys_list.append("sys_tot_down")
    for sys in sys_list:
        if(sys == 'nom_weight'): continue
        sys_ratio = f_ratio.Get("ratio_" + sys)
        sys_ratio.Print()

        sys_LP_weights, _ = d.reweight_LP(sys_ratio, subjet_rw = subjet_rw, fill_z = fill_z, jetR = jetR, num_excjets = num_excjets, uncs = False, prefix = "3prong", max_evts = max_evts)
        sys_weights = weights_nom * sys_LP_weights
        rw = np.sum(weights_nom) / np.sum(sys_weights)
        sys_weights *= rw
        sys_variations[sys] = sys_weights




#make_histogram(LP_weights, "Reweighting factors", 'b', 'Weight', "Lund Plane Reweighting Factors", 20 , h_range = (0., 2.0),
     #normalize=False, fname=outdir + "lundPlane_weights.png")

#compute 'Scalefactor'

eff_nom = np.average(WH_cut, weights = weights_nom)
eff_rw = np.average(WH_cut, weights = weights_rw)

print("Nom %.3f, RW %.3f" % (eff_nom, eff_rw))


eff_toys = []
for i in range(nToys):
    eff = np.average(WH_cut, weights = LP_smeared_weights[:,i])
    eff_toys.append(eff)

toys_mean = np.mean(eff_toys)
toys_std = np.std(eff_toys)

print("Toys avg %.3f, std dev %.3f" % (toys_mean, toys_std))

#Add systematic differences in quadrature
sys_unc_up = sys_unc_down = 0.
if(do_sys_variations):

    eff_sys_tot_up = np.average(WH_cut, weights = sys_variations['sys_tot_up'])
    eff_sys_tot_down = np.average(WH_cut, weights = sys_variations['sys_tot_down'])
    SF_sys_unc_up = (eff_sys_tot_up - eff_rw)/eff_nom
    SF_sys_unc_down = (eff_sys_tot_down - eff_rw)/eff_nom
    SF_sys_unc = max(SF_sys_unc_up, SF_sys_unc_down)


    #for sys in sys_variations.keys():
    #    eff = np.average(WH_cut, weights = sys_variations[sys])
    #    diff = eff - eff_rw
    #    if(diff > 0): sys_unc_up += diff**2
    #    else: sys_unc_down += diff**2
    #    print("%s %.4f" % (sys,  diff))

    #sys_unc_up = sys_unc_up**(0.5)
    #sys_unc_down = sys_unc_down**(0.5)
    


SF = eff_rw / eff_nom
SF_stat_unc = abs(toys_mean - eff_rw)/eff_nom + toys_std /eff_nom

print("\n\nSF (WH_cut val 0.8 ) is %.2f +/- %.2f  (stat) +/- %.2f (sys) \n\n"  % (SF, SF_stat_unc, SF_sys_unc))
f_ratio.Close()

#approximate uncertainty on the reweighting for the plots
overall_unc = (SF_stat_unc **2 + (0.5 * SF_sys_unc_up + 0.5 * SF_sys_unc_down)**2) **0.5 / SF
print("overall unc %.3f" % overall_unc)

