import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from Utils import *


parser = input_options()
options = parser.parse_args()

print(options)

#UL
#lumi = 59.74
#f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_2018/"

lumi = 41.42
f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_2017/"

#lumi = 16.8 + 19.5
#f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_2016/"

f_data = h5py.File(f_dir + "SingleMu_merge.h5", "r")
f_ttbar = h5py.File(f_dir + "TT.h5", "r")
f_wjets = h5py.File(f_dir + "QCD_WJets.h5", "r")
f_diboson = h5py.File(f_dir + "diboson.h5", "r")
f_tw = h5py.File(f_dir + "TW.h5", "r")
f_singletop = h5py.File(f_dir + "SingleTop_merge.h5", "r")



outdir = options.outdir
sys = ""
charge_only = False


do_sys_variations = not options.no_sys
do_plot = False

norm = True

jms_corr = 1.0

m_cut_min = 60.
m_cut_max = 110.
#m_cut_min = 80.
#m_cut_max = 81.
pt_cut = 225.

num_excjets = 2


if(not os.path.exists(outdir)): os.system("mkdir " + outdir)

d_data = Dataset(f_data, is_data = True)

d_tw = Dataset(f_tw, label = "tW", color = ROOT.kMagenta, jms_corr = jms_corr)
d_wjets = Dataset(f_wjets, label = "W+Jets + QCD", color = ROOT.kGray, jms_corr = jms_corr)
d_diboson = Dataset(f_diboson, label = "Diboson", color = ROOT.kCyan, jms_corr = jms_corr)
d_singletop = Dataset(f_singletop, label = "Single Top", color = ROOT.kMagenta+4, jms_corr = jms_corr)


d_ttbar_w_match = Dataset(f_ttbar, label = "ttbar : W-matched", color = ROOT.kRed, jms_corr =jms_corr)
d_ttbar_t_match = Dataset(f_ttbar, label = "ttbar : t-matched ", color = ROOT.kOrange-3, jms_corr = jms_corr)
d_ttbar_nomatch = Dataset(f_ttbar, label = "ttbar : unmatched", color = ROOT.kGreen+3, jms_corr = jms_corr)

#TODO 
d_wjets.norm_unc = 0.06
d_diboson.norm_unc = 0.06
d_singletop.norm_unc = 0.06

d_ttbar_nomatch.norm_unc = 0.06

ttbar_gen_matching = d_ttbar_w_match.f['gen_parts'][:,0]

#0 is unmatched, 1 is W matched, 2 is top matched
nomatch_cut = ttbar_gen_matching < 0.1
w_match_cut = (ttbar_gen_matching  > 0.9) &  (ttbar_gen_matching < 1.1)
t_match_cut = (ttbar_gen_matching  > 1.9) &  (ttbar_gen_matching < 2.1)

d_ttbar_w_match.apply_cut(w_match_cut)
d_ttbar_t_match.apply_cut(t_match_cut)
d_ttbar_nomatch.apply_cut(nomatch_cut)


sigs = [d_ttbar_w_match]
bkgs = [d_ttbar_nomatch, d_ttbar_t_match, d_tw, d_diboson, d_wjets, d_singletop]


ratio_range = [0.5, 1.5]
h_mc = ROOT.TH3F("mc_nom", "Lund Plane MC", n_pt_bins, pt_bins, n_bins_LP,  dr_bins, n_bins_LP, kt_bins) 
h_bkg = ROOT.TH3F("bkg_nom", "Lund Plane Bkg", n_pt_bins, pt_bins, n_bins_LP,  dr_bins, n_bins_LP, kt_bins) 
h_data = ROOT.TH3F("data", "Lund Plane Data", n_pt_bins, pt_bins, n_bins_LP, dr_bins, n_bins_LP, kt_bins) 

h_data_subjets = ROOT.TH1F("data_subjet_pts", "data subjet pts", n_pt_bins, pt_bins)
h_bkg_subjets = ROOT.TH1F("bkg_subjet_pts", "bkg subjet pts", n_pt_bins, pt_bins)
h_mc_subjets = ROOT.TH1F("mc_subjet_pts", "mc subjet pts", n_pt_bins, pt_bins)


h_mc.GetZaxis().SetTitle(z_label)
h_mc.GetYaxis().SetTitle(y_label)
h_bkg.GetZaxis().SetTitle(z_label)
h_bkg.GetYaxis().SetTitle(y_label)
h_data.GetZaxis().SetTitle(z_label)
h_data.GetYaxis().SetTitle(y_label)


bkg_sys_variations = dict()
sig_sys_variations = dict()
if(do_sys_variations):
    keys = sys_weights_map.keys()
    keys.remove("nom_weight")
    for sys in keys: 
        bkg_sys_variations[sys] = h_bkg.Clone(h_bkg.GetName().replace("nom",sys))
        if(sys in sig_sys): sig_sys_variations[sys] = h_mc.Clone(h_mc.GetName().replace("nom",sys))



jet_kinematics_data= f_data['jet_kinematics'][()]
msd_cut_data = (jet_kinematics_data[:,3] > m_cut_min) & (jet_kinematics_data[:,3] < m_cut_max)
pt_cut_data = jet_kinematics_data[:,0] > pt_cut
d_data.apply_cut(msd_cut_data & pt_cut_data)
d_data.compute_obs()

for d in (bkgs + sigs):

    d.norm_factor = lumi

    jet_kinematics = d.f['jet_kinematics'][:]
    msd_cut_mask = (jet_kinematics[:,3] * jms_corr > m_cut_min) & (jet_kinematics[:,3] * jms_corr < m_cut_max)
    pt_cut_mask = jet_kinematics[:,0] > pt_cut
    d.apply_cut(msd_cut_mask & pt_cut_mask)
    d.compute_obs()



#print("Num data %i. Num ttbar MC %i " % (d_data.n(), d_ttbar.n()))


num_data = np.sum(d_data.get_weights())
num_ttbar_nomatch = np.sum(d_ttbar_nomatch.get_weights())
num_ttbar_w_match = np.sum(d_ttbar_w_match.get_weights())
num_ttbar_t_match = np.sum(d_ttbar_t_match.get_weights())
num_ttbar_tot = num_ttbar_nomatch + num_ttbar_w_match + num_ttbar_t_match
num_tw = np.sum(d_tw.get_weights())

tot_bkg = 0.
for d in (d_diboson, d_wjets, d_singletop):
    tot_bkg += np.sum(d.get_weights())
print("%i data, %.0f ttbar (%.0f unmatched, %.0f W matched, %.0f t matched), %.0f tW %.0f bkg" % ( num_data, num_ttbar_tot,num_ttbar_nomatch, 
                                                                                          num_ttbar_w_match, num_ttbar_t_match, num_tw, tot_bkg))
normalization = num_data  / (num_ttbar_tot + num_tw + tot_bkg)
print("normalization", normalization)

if(norm):
    for d in (bkgs + sigs):
        d.norm_factor *= normalization



obs = ["tau21", "tau32", "tau43", "nPF", "mSoftDrop", "pt", "DeepAK8_W_MD"]

colors = []
weights_nom = []
labels = []
for d in (bkgs + sigs):
    colors.append(d.color)
    weights_nom.append(d.get_weights())
    labels.append(d.label)

if(do_plot):
    for l in obs:
        a = []
        for d in (bkgs + sigs):
            a.append(getattr(d, l))

        if(l == 'mSoftDrop'): 
            h_range = (m_cut_min, m_cut_max)
            n_bins_ = n_bins
        elif(l == 'nPF'): 
            h_range = (0.5,120.5)
            n_bins_ = 40
        elif(l == 'pt'): 
            h_range = (pt_cut, 800.)
            n_bins_ = n_bins
        else: 
            n_bins_ = n_bins
            h_range = None
        make_multi_sum_ratio_histogram(data = getattr(d_data, l), entries = a, weights = weights_nom, labels = labels, h_range = h_range, drawSys = False, stack = False,
                colors = colors, axis_label = l,  title = l + " : No Reweighting", num_bins = n_bins_, normalize = False, ratio_range = (0.5, 1.5), fname = outdir + l + '_ratio_before.png' )




LP_rw = LundReweighter(jetR = jetR, charge_only = options.charge_only)

d_data.subjets = d_data.fill_LP(LP_rw, h_data,  num_excjets = num_excjets,  rescale_subjets = "vec" )

for d in sigs:
    d.subjets = d.fill_LP(LP_rw, h_mc,  num_excjets = num_excjets, sys_variations = sig_sys_variations,  rescale_subjets = "vec" )

for d in bkgs:
    d.subjets = d.fill_LP(LP_rw, h_bkg, num_excjets = num_excjets, sys_variations = bkg_sys_variations, rescale_subjets = "vec")


for d in ([d_data] + sigs + bkgs): 
    d.subjet_pt = []
    if(d is d_data): h_subjets = h_data_subjets
    elif(d in sigs): h_subjets = h_mc_subjets
    elif(d in bkgs): h_subjets = h_bkg_subjets

    weights = d.get_weights()

    for idx,sjs in enumerate(d.subjets):
        sj_pts = []
        for sj in sjs: 
            sj_pts.append(sj[0])
            h_subjets.Fill(sj[0], weights[idx])
        d.subjet_pt.append(sj_pts)

obs.append("subjet_pt")

default = ROOT.TStyle("Default","Default Style");
default.cd()
ROOT.gROOT.SetStyle('Default')
ROOT.gStyle.SetOptStat(0) # To display the mean and RMS:   SetOptStat("mr")

f_out = ROOT.TFile.Open(outdir + "ratio.root", "RECREATE")
nom_ratio = LP_rw.make_LP_ratio(h_data, h_bkg,  h_mc,  h_data_subjets, h_bkg_subjets, h_mc_subjets, pt_bins = pt_bins, outdir = outdir, save_plots = True)
nom_ratio.SetName("ratio_nom")
nom_ratio.Write()
h_mc.Write()
h_bkg.Write()
h_data.Write()

if(do_sys_variations):
    keys = sys_weights_map.keys()
    keys.remove("nom_weight")
    for i,sys_name in enumerate(keys):
        print(sys_name)
        h_bkg_subjets_sys = h_bkg_subjets.Clone(sys_name + "_bkg_ptnorm")
        h_bkg_subjets_sys.Reset()

        sys_idx = sys_weights_map[sys_name]

        #compute pt distributions for this sys
        for d in bkgs:
            if('norm' in sys): #normalization uncs split by process
                process = sys.split("_")[0]
                if(process in d.label and 'up' in sys): weights_sys = nom_weights * (1. + d.norm_unc)
                elif(process in d.label and 'down' in sys): weights_sys = nom_weights  * (1. - d.norm_unc)
                else: weights_sys = nom_weights
            else:  
                all_sys_weights = d.get_masked('sys_weights')
                weights_sys = d.get_weights() * all_sys_weights[:, sys_idx]
            fill_hist(h_bkg_subjets_sys, d.subjet_pt, weights_sys)

        if(sys_name in sig_sys): #some systematics effect signal
            h_mc_subjets_sys = h_mc_subjets.Clone(sys_name + "_mc_ptnorm")
            h_mc_subjets_sys.Reset()
            for d in sigs:
                all_sys_weights = d.get_masked('sys_weights')
                weights_sys = d.get_weights() * all_sys_weights[:, sys_idx]
                fill_hist(h_mc_subjets_sys, d.subjet_pt, weights_sys)

            h_mc_sys = sig_sys_variations[sys_name]

        else:
            h_mc_subjets_sys = h_mc_subjets.Clone("nom_clone")
            h_mc_sys = h_mc



        h_bkg_sys = bkg_sys_variations[sys_name]

        #Some systematics only for bkgs not signal (want denom of ratio to be consistent)

        sys_ratio = LP_rw.make_LP_ratio(h_data, h_bkg_sys, h_mc_sys, h_data_subjets, h_bkg_subjets_sys, h_mc_subjets_sys, pt_bins = pt_bins)
        sys_ratio.SetName("ratio_" + sys_name)
        #sys_ratio.Print("range") 
        sys_ratio.Write()





if(do_plot):
    weights_rw = copy.deepcopy(weights_nom)

    LP_weights = []
    for i,d in enumerate(sigs):
        d_LP_weights  = d.reweight_LP(LP_rw, nom_ratio,  num_excjets = num_excjets )
        LP_weights.append(d_LP_weights)

        weights_rw[len(bkgs) + i] *= d_LP_weights


    make_histogram(LP_weights[0], "Reweighting factors", 'b', 'Weight', "Lund Plane Reweighting Factors", 20 , h_range = (0., 2.0),
         normalize=False, fname=outdir + "lundPlane_weights.png")

    for l in obs:
        a = []
        for i,d in enumerate(bkgs + sigs):
            a.append(getattr(d, l))

        if(l == 'mSoftDrop'): 
            h_range = (m_cut_min, m_cut_max)
            n_bins_ = n_bins
        elif(l == 'nPF'): 
            h_range = (0.5,120.5)
            n_bins_ = 40
        elif(l == 'pt'): 
            h_range = (pt_cut, 800.)
            n_bins_ = n_bins
        elif(l == 'subjet_pt'): 
            h_range = (0., 800.)
            n_bins_ = n_bins
        else: 
            n_bins_ = n_bins
            h_range = None

        make_multi_sum_ratio_histogram(data = getattr(d_data, l), entries = a, weights = weights_rw, labels = labels, h_range = h_range, drawSys = False, stack = False,
                colors = colors, axis_label = l,  title = l + " : LP Reweighting", num_bins = n_bins_, normalize = False, ratio_range = (0.5, 1.5), fname = outdir + l + '_ratio_after.png' )


f_out.Close()
