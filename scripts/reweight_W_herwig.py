import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.Utils import *


parser = input_options()
parser.add_argument("--noBkg", default=False, action='store_true',  help="Ideal case, no unmerged bkg")
options = parser.parse_args()

print(options)

#UL
lumi = 59.74
f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_gen/"


f_powheg = h5py.File(f_dir + "TT_powheg.h5", "r")
f_herwig = h5py.File(f_dir + "TT_herwig.h5", "r")



outdir = options.outdir
sys = ""
#CA_prefix = "2prong"
CA_prefix = ""
charge_only = False


do_sys_variations = True
do_plot = True

norm = True


m_cut_min = 60.
m_cut_max = 110.
#m_cut_min = 80.
#m_cut_max = 81.
pt_cut = 225.

num_excjets = 2


if(not os.path.exists(outdir)): os.system("mkdir " + outdir)

d_powheg_w_match = Dataset(f_powheg, label = "powheg : W-matched", color = ROOT.kRed, is_data = True)
d_powheg_t_match = Dataset(f_powheg, label = "powheg : t-matched ", color = ROOT.kOrange-3, is_data = True)
d_powheg_nomatch = Dataset(f_powheg, label = "powheg : unmatched", color = ROOT.kGreen+3, is_data = True)

d_herwig = Dataset(f_herwig, label = "herwig", color = ROOT.kRed, is_data = True)


d_powheg_nomatch.norm_unc = 0.06

powheg_gen_matching = d_powheg_w_match.f['gen_parts'][:,0]

#0 is unmatched, 1 is W matched, 2 is top matched
nomatch_cut = powheg_gen_matching < 0.1
w_match_cut = (powheg_gen_matching  > 0.9) &  (powheg_gen_matching < 1.1)
t_match_cut = (powheg_gen_matching  > 1.9) &  (powheg_gen_matching < 2.1)

d_powheg_w_match.apply_cut(w_match_cut)
d_powheg_t_match.apply_cut(t_match_cut)
d_powheg_nomatch.apply_cut(nomatch_cut)


sigs = [d_powheg_w_match]
bkgs = [d_powheg_nomatch, d_powheg_t_match] 
if(options.noBkg):
    bkgs = []
    herwig_gen_matching = d_herwig.f['gen_parts'][:,0]
    w_match_cut_herwig = (herwig_gen_matching  > 0.9) &  (herwig_gen_matching < 1.1)
    d_herwig.apply_cut(w_match_cut_herwig)


ratio_range = [0.5, 1.5]
h_mc = ROOT.TH3F("mc_nom", "Lund Plane MC", n_pt_bins, pt_bins, n_bins_LP,  dr_bins, n_bins_LP, kt_bins) 
h_bkg = ROOT.TH3F("bkg_nom", "Lund Plane Bkg", n_pt_bins, pt_bins, n_bins_LP,  dr_bins, n_bins_LP, kt_bins) 
h_herwig = ROOT.TH3F("herwig", "Lund Plane herwig", n_pt_bins, pt_bins, n_bins_LP, dr_bins, n_bins_LP, kt_bins) 

h_herwig_subjets = ROOT.TH1F("herwig_subjet_pts", "herwig subjet pts", n_pt_bins, pt_bins)
h_bkg_subjets = ROOT.TH1F("bkg_subjet_pts", "bkg subjet pts", n_pt_bins, pt_bins)
h_mc_subjets = ROOT.TH1F("mc_subjet_pts", "mc subjet pts", n_pt_bins, pt_bins)


h_mc.GetZaxis().SetTitle(z_label)
h_mc.GetYaxis().SetTitle(y_label)
h_bkg.GetZaxis().SetTitle(z_label)
h_bkg.GetYaxis().SetTitle(y_label)
h_herwig.GetZaxis().SetTitle(z_label)
h_herwig.GetYaxis().SetTitle(y_label)


bkg_sys_variations = dict()
sig_sys_variations = dict()

#only Parton shower uncs apply
sys_keys = ['PS_ISR_up', 'PS_ISR_down', 'PS_FSR_up', 'PS_FSR_down', 'unmatched_norm_up', 'unmatched_norm_down']
if(do_sys_variations):
    for sys in sys_keys: 
        bkg_sys_variations[sys] = h_bkg.Clone(h_bkg.GetName().replace("nom",sys))



for d in (bkgs + sigs + [d_herwig]):

    jet_kinematics = d.f['jet_kinematics'][:]
    msd_cut_mask = (jet_kinematics[:,3]  > m_cut_min) & (jet_kinematics[:,3]  < m_cut_max)
    pt_cut_mask = jet_kinematics[:,0] > pt_cut
    d.apply_cut(msd_cut_mask & pt_cut_mask)
    d.compute_obs()




num_herwig = np.sum(d_herwig.get_weights())

num_powheg_nomatch = np.sum(d_powheg_nomatch.get_weights())
num_powheg_w_match = np.sum(d_powheg_w_match.get_weights())
num_powheg_t_match = np.sum(d_powheg_t_match.get_weights())
num_powheg_tot = num_powheg_nomatch + num_powheg_w_match + num_powheg_t_match

print("%i herwig, %.0f powheg (%.0f unmatched, %.0f W matched, %.0f t matched)" % ( num_herwig, num_powheg_tot,num_powheg_nomatch, 
                                                                                          num_powheg_w_match, num_powheg_t_match))

norm = num_herwig / (num_powheg_w_match if options.noBkg else num_powheg_tot)

#normalize two samples to match
for d in (bkgs + sigs):
    d.norm_factor *= norm


obs = ["tau21", "tau32", "tau43", "nPF", "mSoftDrop", "pt" ]

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
        o_data = getattr(d_herwig, l)
        for d in (bkgs + sigs):
            a.append(getattr(d, l))
        if(l == 'mSoftDrop'): 
            h_range = (m_cut_min, m_cut_max)
            n_bins_ = n_bins
            l = 'mass'
        elif(l == 'nPF'): 
            h_range = (0.5,120.5)
            n_bins_ = 40
        elif(l == 'pt'): 
            h_range = (pt_cut, 800.)
            n_bins_ = n_bins
        else: 
            n_bins_ = n_bins
            h_range = None
        make_multi_sum_ratio_histogram(data = o_data, entries = a, weights = weights_nom, labels = labels, h_range = h_range, drawSys = False, stack = False, data_label = "herwig",
                colors = colors, axis_label = l,  title = l + " : No Reweighting", num_bins = n_bins_, normalize = True, ratio_range = (0.5, 1.5), fname = outdir + l + '_ratio_before.png' )




LP_rw = LundReweighter(jetR = jetR, charge_only = options.charge_only)

d_herwig.subjets = d_herwig.fill_LP(LP_rw, h_herwig,  num_excjets = num_excjets, prefix = CA_prefix, rescale_subjets = "vec" )

for d in sigs:
    d.subjets = d.fill_LP(LP_rw, h_mc,  num_excjets = num_excjets, sys_variations = sig_sys_variations, prefix = CA_prefix,  rescale_subjets = "vec" )

for d in bkgs:
    d.subjets = d.fill_LP(LP_rw, h_bkg, num_excjets = num_excjets, sys_variations = bkg_sys_variations, prefix = CA_prefix,  rescale_subjets = "vec")


for d in ([d_herwig] + sigs + bkgs): 
    d.subjet_pt = []
    if(d is d_herwig): h_subjets = h_herwig_subjets
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
nom_ratio = LP_rw.make_LP_ratio(h_herwig, h_bkg,  h_mc,  h_herwig_subjets, h_bkg_subjets, h_mc_subjets, pt_bins = pt_bins, outdir = outdir, save_plots = True)
nom_ratio.SetName("ratio_nom")
nom_ratio.Write()
h_mc.Write()
h_bkg.Write()
h_herwig.Write()

if(do_sys_variations and not options.noBkg):
    for i,sys_name in enumerate(sys_keys):
        print(sys_name)
        h_bkg_subjets_sys = h_bkg_subjets.Clone(sys_name + "_bkg_ptnorm")
        h_bkg_subjets_sys.Reset()
        h_mc_subjets_sys = h_mc_subjets.Clone(sys_name + "_mc_ptnorm")

        sys_idx = sys_weights_map[sys_name]

        #compute pt distributions for this sys
        for d in bkgs:
            all_sys_weights = d.get_masked('sys_weights')

            if('norm' in sys): #normalization uncs split by process
                process = sys.split("_")[0]
                if(process in d.label and 'up' in sys): weights_sys = d.get_weights() * (1. + d.norm_unc)
                elif(process in d.label and 'down' in sys): weights_sys = d.get_weights()  * (1. - d.norm_unc)
                else: weights_sys = d.get_weights()
            else:
                weights_sys = d.get_weights() * all_sys_weights[:, sys_idx]

            fill_hist(h_bkg_subjets_sys, d.subjet_pt, weights_sys)



        h_bkg_sys = bkg_sys_variations[sys_name]
        h_mc_sys = h_mc

        sys_ratio = LP_rw.make_LP_ratio(h_herwig, h_bkg_sys, h_mc_sys, h_herwig_subjets, h_bkg_subjets_sys, h_mc_subjets_sys, pt_bins = pt_bins)
        sys_ratio.SetName("ratio_" + sys_name)
        #sys_ratio.Print("range") 
        sys_ratio.Write()





if(do_plot):
    weights_rw = copy.deepcopy(weights_nom)

    LP_weights = []
    for i,d in enumerate(sigs):
        d_LP_weights  = d.reweight_LP(LP_rw, nom_ratio,  num_excjets = num_excjets, prefix = CA_prefix)
        LP_weights.append(d_LP_weights)

        weights_rw[len(bkgs) + i] *= d_LP_weights


    make_histogram(LP_weights[0], "Reweighting factors", 'b', 'Weight', "Lund Plane Reweighting Factors", 20 , h_range = (0., 2.0),
         normalize=False, fname=outdir + "lundPlane_weights.png")

    for l in obs:
        a = []
        for i,d in enumerate(bkgs + sigs):
            a.append(getattr(d, l))

        o_data = getattr(d_herwig, l)
        if(l == 'mSoftDrop'): 
            h_range = (m_cut_min, m_cut_max)
            n_bins_ = n_bins
            l = 'mass'
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

        make_multi_sum_ratio_histogram(data = o_data, entries = a, weights = weights_rw, labels = labels, h_range = h_range, drawSys = False, stack = False, data_label = "herwig",
                colors = colors, axis_label = l,  title = l + " : LP Reweighting", num_bins = n_bins_, normalize = True, ratio_range = (0.5, 1.5), fname = outdir + l + '_ratio_after.png' )


f_out.Close()
