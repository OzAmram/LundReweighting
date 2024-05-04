import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.Utils import *

parser = input_options()
options = parser.parse_args()

print(options)

if(not os.path.exists(options.outdir)): os.system("mkdir " + options.outdir)

ROOT.TGaxis.SetMaxDigits(3);



#UL
if(options.year == 2018):
    lumi = 59.74
    year = 2018
    f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_2018/"

elif(options.year == 2017):
    lumi = 41.42
    year = 2017
    f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_2017/"

elif(options.year == 2016):
    year = 2016
    lumi = 16.8 + 19.5
    f_dir = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_2016/"
else:
    exit(1)


f_data = h5py.File(f_dir + "SingleMu_merge.h5", "r")
f_ttbar = h5py.File(f_dir + "TT.h5", "r")
f_wjets = h5py.File(f_dir + "QCD_WJets.h5", "r")
f_diboson = h5py.File(f_dir + "diboson.h5", "r")
f_tw = h5py.File(f_dir + "TW.h5", "r")
f_singletop = h5py.File(f_dir + "SingleTop_merge.h5", "r")

d_data = Dataset(f_data, is_data = True)

d_tw = Dataset(f_tw, label = "tW", color = ROOT.kYellow-7,  dtype = -1)
d_wjets = Dataset(f_wjets, label = "W+Jets + QCD", color = ROOT.kOrange-3) 
d_diboson = Dataset(f_diboson, label = "Diboson", color = ROOT.kCyan )
d_singletop = Dataset(f_singletop, label = "Single Top", color = ROOT.kMagenta-1 )



d_ttbar_w_match = Dataset(f_ttbar, label = "t#bar{t} : W-matched", color = ROOT.kRed-7,  dtype = 2)
d_ttbar_t_match = Dataset(f_ttbar, label = "t#bar{t} : t-matched ", color = ROOT.kBlue-7, dtype = 3)
d_ttbar_nomatch = Dataset(f_ttbar, label = "t#bar{t} : unmatched", color = ROOT.kGreen-6, )


ttbar_gen_matching = d_ttbar_w_match.f['gen_parts'][:,0]

#0 is unmatched, 1 is W matched, 2 is top matched
nomatch_cut = ttbar_gen_matching < 0.1
w_match_cut = (ttbar_gen_matching  > 0.9) &  (ttbar_gen_matching < 1.1)
t_match_cut = (ttbar_gen_matching  > 1.9) &  (ttbar_gen_matching < 2.1)

d_ttbar_w_match.apply_cut(w_match_cut)
d_ttbar_t_match.apply_cut(t_match_cut)
d_ttbar_nomatch.apply_cut(nomatch_cut)

samples = [d_diboson, d_singletop, d_wjets, d_ttbar_t_match, d_ttbar_nomatch, d_ttbar_w_match, d_tw]

m_cut_min = 50.
m_cut_max = 250.
#m_cut_max = 65.
pt_cut = 225.

for d in (samples + [d_data]):

    if(d is not d_data): d.norm_factor = lumi

    jet_kinematics = d.f['jet_kinematics'][:]
    msd_cut_mask = (jet_kinematics[:,3] > m_cut_min) & (jet_kinematics[:,3] < m_cut_max)
    pt_cut_mask = jet_kinematics[:,0] > pt_cut
    d.apply_cut(msd_cut_mask & pt_cut_mask)
    d.compute_kinematics()
    d.apply_cut(d.dR_mu_bjet > 0.1)
    d.compute_kinematics()
    d.compute_obs()


obs = ["mSoftDrop", "mu_pt", "mu_eta", "mu_phi", "met_pt", "w_cand_pt", "bjet_pt", "bjet_eta", "bjet_phi", 
       "ak8_pt", "ak8_eta", "ak8_phi", "dphi_mu_bjet", "dR_mu_bjet", "dphi_mu_ak8", "dR_mu_ak8"]

pi = 3.14159
obs_attrs = {
        'mu_pt' : (60, 500., 20, "Muon p_{T} [GeV]", "Events / bin"),
        'mu_eta' : (-2.4, 2.4, 20, "Muon #eta", "Events / bin"),
        'mu_phi' : (-pi, pi, 20, "Muon #phi", "Events / bin"),
        'met_pt' : (50, 500., 20, "MET [GeV]", "Events / bin"),
        'w_cand_pt' : (80, 500., 20, "Leptonic W Cand. p_{T} [GeV]", "Events / bin"),
        'bjet_pt' : (50, 500., 20, "AK4 jet p_{T} [GeV]", "Events / bin"),
        'bjet_eta' : (-2.4, 2.4, 20, "AK4 jet  #eta", "Events / bin"),
        'bjet_phi' : (-pi, pi, 20, "AK4 jet #phi", "Events / bin"),
        'ak8_eta' : (-2.4, 2.4, 20, "AK8 jet  #eta", "Events / bin"),
        'ak8_phi' : (-pi, pi, 20, "AK8 jet #phi", "Events / bin"),
        'mSoftDrop' : (60, 110, 25, "m_{SD} [GeV] ", "Events / 2 GeV") if m_cut_max < 200 else (50, 230, 45, "m_{SD} [GeV]", "Events / 4 GeV"),
        'ak8_pt' : (225, 825., 20, "AK8 Jet p_{T} [GeV]", "Events / 30 GeV"),
        'dphi_mu_ak8' : (0, pi, 20, "Muon - AK8 jet |#Delta#phi|", "Events / bin"),
        'dphi_mu_bjet' : (0, pi, 20, "Muon - AK4 jet |#Delta#phi|", "Events / bin"),
        'dR_mu_ak8' : (0, 4.0, 40, "Muon - AK8 jet |#DeltaR|", "Events / bin"),
        'dR_mu_bjet' : (0, 4.0, 40, "Muon - AK4 jet |#DeltaR|", "Events / bin"),
        }

colors = []
weights_nom = []
labels = []
for d in (samples):
    colors.append(d.color)
    weights_nom.append(d.get_weights())
    labels.append(d.label)


ratio_range = [0.2, 1.8]
for l in obs:
    a = []
    for d in (samples):
        a.append(getattr(d, l))
    a_data = getattr(d_data, l)

    low,high, nbins_, label, ylabel = obs_attrs.get(l, (None, None, 20, l, ""))

    make_multi_sum_ratio_histogram(data = a_data, entries = a, weights = weights_nom, labels = labels, h_range = (low, high), drawSys = False, stack = True, draw_chi2 = False,
            year = year, colors = colors, axis_label = label,  title = l,
            num_bins = nbins_, normalize = False, ratio_range = ratio_range, fname = options.outdir + l + '.png' )
