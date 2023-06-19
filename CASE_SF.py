from Utils import *
import os



parser = input_options()
options = parser.parse_args()

print(options)


outdir = options.outdir
if(not os.path.exists(outdir)): os.system("mkdir %s" % outdir)
jet_str = 'CA'

fname = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/TagNTrain/data/YtoHH_Htott_Y3000_H400_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5"
#fname = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/TagNTrain/data/ZpToTpTp_Zp5000_Tp400_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5"
#fname = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/TagNTrain/data/XToYYprimeTo4Q_MX3000_MY170_MYprime170_narrow_TuneCP5_13TeV-madgraph-pythia8_TIMBER.h5"
#fname = "test_signal_CASE.h5"



label = "YtoHH"
h_label = label
#label = "ZpToTpTp"
#h_label = 'Zp'
#label = "XtoYY"
#h_label = 'XYY'

tag_obs = 'tau43'
score_thresh = 0.65
#tag_obs = 'tau21'
#score_thresh = 0.34


f_ratio = ROOT.TFile.Open(options.fin)
f_sig = h5py.File(fname, "r")
j_idx = 0

jetR = 1.0
num_excjets = -1

max_evts = 20000
#max_evts = 100
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

#for l in obs:
#
#    if(l == 'nPF'): 
#        h_range = (0.5,120.5)
#        n_bins_ = 40
#    else: 
#        n_bins_ = n_bins
#        h_range = None
#
#    make_histogram(getattr(d,l), l, 'b', l, l, n_bins_, h_range = h_range, normalize=True, fname=outdir + l + "_before.png")
#
#    #make_histogram(data = getattr(d, l), labels = labels, h_range = h_range, drawSys = False, stack = False,
#            #colors = colors, axis_label = l,  title = l + " : No Reweighting", num_bins = n_bins, normalize = False, ratio_range = (0.5, 1.5), fname = outdir + l + '_ratio_before.png' )

score = getattr(d, tag_obs)[:max_evts]

print("%i events" % score.shape[0])




score_cut = score < score_thresh

weights_nom = np.ones_like(score)


weights_rw = copy.deepcopy(weights_nom)

h_ratio = f_ratio.Get("ratio_nom")
f_ratio.cd('pt_extrap')
rdir = ROOT.gDirectory
#rdir = None

nToys = 100

#Noise used to generated smeared ratio's based on stat unc
rand_noise = np.random.normal(size = (nToys, h_ratio.GetNbinsX(), h_ratio.GetNbinsY(), h_ratio.GetNbinsZ()))
pt_rand_noise = np.random.normal(size = (nToys, h_ratio.GetNbinsY(), h_ratio.GetNbinsZ(), 3))

LP_rw = LundReweighter(jetR = jetR, pt_extrap_dir = rdir, charge_only = options.charge_only)


subjets, splittings, bad_match, deltaRs = d.get_matched_splittings(LP_rw, num_excjets = num_excjets, max_evts = max_evts, return_dRs = True)


d_LP_weights, d_LP_smeared_weights, d_pt_smeared_weights = d.reweight_LP(LP_rw, h_ratio, num_excjets = num_excjets, 
        max_evts = max_evts, prefix = "", rand_noise = rand_noise, pt_rand_noise = pt_rand_noise, subjets = subjets, splittings = splittings)






LP_weights = d_LP_weights
print(LP_weights[:10])


#apply weights, keep normalization fixed
old_norm = np.sum(weights_rw)
weights_rw *= d_LP_weights

new_norm = np.sum(weights_rw)

weights_rw *= old_norm / new_norm
LP_smeared_weights = np.array(d_LP_smeared_weights * np.expand_dims(weights_nom, -1) * (old_norm / new_norm))
pt_smeared_weights = np.array(d_pt_smeared_weights * np.expand_dims(weights_nom, -1) * (old_norm / new_norm))


print("MEAN weight", np.mean(weights_rw))


make_histogram(weights_rw, "Reweighting factors", 'b', 'Weight', "Lund Plane Reweighting Factors", 40 , h_range = (0., LP_rw.max_rw + 0.1),
     normalize=False, fname=outdir + "lundPlane_weights.png")


frac_low_edge = np.mean( weights_rw < np.amin(weights_rw) * 1.01)
frac_high_edge = np.mean( weights_rw > np.amax(weights_rw) * 0.99)

print(np.amin(weights_rw), np.amax(weights_rw))
print(weights_rw[:10])

print("Frac of weights low edge %.3f, high edge %.3f" % (frac_low_edge, frac_high_edge))


sys_variations = dict()
if(not options.no_sys):
    #sys_list = list(sys_weights_map.keys())
    sys_list = ["sys_tot_up", "sys_tot_down"]
    for sys in sys_list:
        sys_ratio = f_ratio.Get("ratio_" + sys)
        sys_ratio.Print()
        sys_str = sys + "_"


        sys_LP_weights = d.reweight_LP(LP_rw, sys_ratio, num_excjets = num_excjets, prefix = "", 
                max_evts = max_evts, sys_str = sys_str, subjets = subjets, splittings = splittings)
        sys_weights = weights_nom * sys_LP_weights
        rw = np.sum(weights_nom) / np.sum(sys_weights)
        sys_weights *= rw
        sys_variations[sys] = sys_weights

    #vary weights up/down for b-quark subjets by ratio of b-quark to light quark LP
    b_light_ratio = f_ratio.Get("h_bl_ratio")
    bquark_rw = d.reweight_LP(LP_rw, b_light_ratio, num_excjets = num_excjets, prefix = "", 
            max_evts = max_evts, sys_str = 'bquark', subjets = subjets, splittings = splittings)

    up_bquark_weights = bquark_rw * weights_rw
    down_bquark_weights = (1./ bquark_rw) * weights_rw

    up_bquark_weights *= old_norm / np.sum(up_bquark_weights)
    down_bquark_weights *= old_norm / np.sum(down_bquark_weights)

    sys_variations['bquark_up'] = up_bquark_weights
    sys_variations['bquark_down'] = down_bquark_weights


#Save subjet pts and deltaR
subjet_pts =  []
deltaRs = np.reshape(deltaRs, -1)

for i,sjs in enumerate(subjets):
    for sj in sjs:
        subjet_pts.append(sj[0])

#n_splits = []
#for splits in splittings:
#    count_jet_idx = 0
#    n_split = 0
#    for jet_idx,delta,kt in splits:
#        if(jet_idx == count_jet_idx):
#            n_split += 1
#        else:
#            n_splits.append(n_split)
#            count_jet_idx = jet_idx
#            n_split = 1
#    n_splits.append(n_split)
#
#print(n_splits[:5])
        

#subjet_pts = np.array(subjet_pts).reshape(-1)
#n_splits = np.array(n_splits).reshape(-1)
#print(subjet_pts.shape, n_splits.shape)
#make_scatter_plot(subjet_pts, n_splits, "b", ["subjet pt", "num splittings"], fname = options.outdir + "subjet_pt_vs_nsplit.png")
    
num_bins = 40
pt_bins = array('d', np.linspace(0., 800., num_bins + 1))
dR_bins = array('d', np.linspace(0., 0.8, num_bins + 1))

h_subjet_pts = make_root_hist(data = subjet_pts, name = 'h_%s_subjetpt' %h_label, num_bins = num_bins, bins = pt_bins)
h_dRs = make_root_hist(data = deltaRs, name = 'h_%s_dRs' %h_label, num_bins = num_bins, bins = dR_bins)
f_ptout = ROOT.TFile.Open(outdir + "subjet_pt_dR.root", "RECREATE")
h_subjet_pts.Write()
h_dRs.Write()
f_ptout.Close()

make_histogram([subjet_pts], ["Subjets"], colors = ['blue'], xaxis_label = 'Subjet pt (GeV)', 
                title = "%s : subjet pt " % (label), num_bins = 40, normalize = True, fname = options.outdir + label + "_subjet_pt.png")

subjet_pts = np.array(subjet_pts).reshape(-1)
print(np.amin(subjet_pts))
print("Fraction of subjets with pt > 350 : %.3f" % (np.mean( subjet_pts > 350.)))




#compute 'Scalefactor'

eff_nom = np.average(score_cut, weights = weights_nom)
eff_rw = np.average(score_cut, weights = weights_rw)

print("Nom %.3f, RW %.3f" % (eff_nom, eff_rw))


eff_toys = []
pt_eff_toys = []
for i in range(nToys):
    eff = np.average(score_cut, weights = LP_smeared_weights[:,i])
    eff_toys.append(eff)

    eff1 = np.average(score_cut, weights = pt_smeared_weights[:,i])
    pt_eff_toys.append(eff1)

toys_mean = np.mean(eff_toys)
toys_std = np.std(eff_toys)

print("Toys avg %.3f, std dev %.3f" % (toys_mean, toys_std))


pt_toys_mean = np.mean(pt_eff_toys)
pt_toys_std = np.std(pt_eff_toys)

print("Pt variation toys avg %.3f, std dev %.3f" % (pt_toys_mean, pt_toys_std))

#Add systematic differences in quadrature
sys_unc_up = sys_unc_down = 0.
if(not options.no_sys):

    eff_sys_tot_up = np.average(score_cut, weights = sys_variations['sys_tot_up'])
    eff_sys_tot_down = np.average(score_cut, weights = sys_variations['sys_tot_down'])
    SF_sys_unc_up = abs(eff_sys_tot_up - eff_rw)/eff_nom
    SF_sys_unc_down = abs(eff_sys_tot_down - eff_rw)/eff_nom
    SF_sys_unc = (SF_sys_unc_up + SF_sys_unc_down) / 2.0

    eff_bquark_up = np.average(score_cut, weights = sys_variations['bquark_up'])
    eff_bquark_down = np.average(score_cut, weights = sys_variations['bquark_down'])
    SF_bquark_up = abs(eff_bquark_up - eff_rw)/eff_nom
    SF_bquark_down = abs(eff_bquark_down - eff_rw)/eff_nom
    SF_bquark_unc = (SF_bquark_up + SF_bquark_down) /2.0

else: 
    SF_sys_unc = SF_sys_unc_up = SF_sys_unc_down = bquark_unc = 0. 


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
SF_pt_unc = abs(pt_toys_mean - eff_rw)/eff_nom + pt_toys_std /eff_nom

print(bad_match[:10])
#fraction of evts with bad match, take as fractional unc on SF
bad_matching_unc = np.mean(bad_match) * SF

print("\n\nSF (%s val %.2f ) is %.2f +/- %.2f  (stat) +/- %.2f (sys) +/- %.2f (pt) +/- %.2f (bquark) +/- %.2f (matching) \n\n"  
        % (tag_obs, score_thresh, SF, SF_stat_unc, SF_sys_unc, SF_pt_unc, SF_bquark_unc, bad_matching_unc))
f_ratio.Close()

#approximate uncertainty on the reweighting for the plots
overall_unc = (SF_stat_unc **2 + SF_sys_unc**2 + SF_pt_unc**2 + SF_bquark_unc**2 + bad_matching_unc**2) **0.5 
print("overall unc %.3f" % overall_unc)



weights = [  weights_nom, weights_rw, sys_variations['sys_tot_up'], sys_variations['sys_tot_down']]
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

    print(x.shape, weights[0].shape)

    make_multi_ratio_histogram([x]*4, labels, colors, l, "%s : Before vs. After Lund Plane Reweighting" % l, n_bins,
         ratio_range = ratio_range, normalize=True, weights = weights, fname=outdir + "%s_before_vs_after.png" %l )

