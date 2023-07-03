from Utils import *
import os

parser = input_options()
options = parser.parse_args()

mode = 'full'

ratio_range = (0.2, 1.8)
if(not os.path.exists(options.outdir)):
    os.system("mkdir " + options.outdir)


incl_after = True
if(mode == 'top'):
    dirs = ["top_SF_2016_july1/", "top_SF_2017_july1/", "top_SF_june16/"]
    labels = [
            "Single Top",
            "W+Jets + QCD",
            "tW",
            "t#bar{t} : W-matched",
            "t#bar{t} : unmatched",
            "t#bar{t} : t-matched",
        ]
    colors = [
        ROOT.kMagenta-1,
        ROOT.kOrange-3,
        ROOT.kYellow-7,
        ROOT.kRed -7,
        ROOT.kGreen-6,
        ROOT.kBlue -7, 
        ]


    obs_attrs = {
            'mSoftDrop' : (120, 210, 30, "m_{SD}"),
            'tau21' : (0.05, 0.8, 20, "#tau_{21}"),
            'tau32' : (0.2, 1.0, 20, "#tau_{32}"),
            'tau43' : (0.6, 1.0, 20, "#tau_{43}"),
            'nPF' : (20.5, 120.5, 25, "Num. PF Cands."),
            'pt' : (500, 1200., 20, "p_{T}"),
            }
    draw_chi2 = True


elif(mode == "full"):
    dirs = ["W_fullrange_2016/", "W_fullrange_2017/", "W_fullrange_2018/"]
    labels = [
            "Diboson",
            "Single Top",
            "W+Jets + QCD",
            "t#bar{t} : t-matched",
            "t#bar{t} : unmatched",
            "tW",
            "t#bar{t} : W-matched",
        ]

    colors = [
        ROOT.kCyan,
        ROOT.kMagenta-1,
        ROOT.kOrange-3,
        ROOT.kBlue -7, 
        ROOT.kGreen-6,
        ROOT.kYellow-7,
        ROOT.kRed -7,
        ]

    obs_attrs = {
            'mSoftDrop' : (50, 210, 30, "m_{SD}"),
            }
    draw_chi2 = False
    incl_after = False

else:
    dirs = ["W_SF_2016_july1/", "W_SF_2017_july1/", "W_SF_june16/"]
    labels = [
            "Diboson",
            "Single Top",
            "W+Jets + QCD",
            "t#bar{t} : t-matched",
            "t#bar{t} : unmatched",
            "tW",
            "t#bar{t} : W-matched",
        ]

    colors = [
        ROOT.kCyan,
        ROOT.kMagenta-1,
        ROOT.kOrange-3,
        ROOT.kBlue -7, 
        ROOT.kGreen-6,
        ROOT.kYellow-7,
        ROOT.kRed -7,
        ]

    obs_attrs = {
            'mSoftDrop' : (60, 110, 30, "m_{SD}"),
            'tau21' : (0.05, 0.8, 20, "#tau_{21}"),
            'tau32' : (0.4, 1.0, 20, "#tau_{32}"),
            'tau43' : (0.6, 1.0, 20, "#tau_{43}"),
            'nPF' : (20.5, 80.5, 30, "Num. PF Cands."),
            'pt' : (225, 1200., 20, "p_{T}"),
            'DeepAK8_W' : (0., 1., 20, "DeepAK8 (W vs. QCD)"),
            'DeepAK8_W_MD' : (0., 1., 20, "DeepAK8-MD (W vs. QCD)"),
            }
    draw_chi2 = True



for obs in obs_attrs.keys():
    plt_name = obs + "_ratio_"
    
    h_tot_after = h_data = None

    hists = {}
    for dir_in in dirs:
        if(incl_after):
            f_after = ROOT.TFile.Open(dir_in + plt_name + "after.root")
            h_tot_after_file = f_after.Get("h_tot")
            if(h_tot_after is None):
                h_tot_after = h_tot_after_file
                h_tot_after.SetDirectory(0)
                h_tot_after.SetName("Total MC: LP RW")
                h_tot_after.SetFillStyle(0)
            else:
                h_tot_after.Add(h_tot_after_file)
            f_after.Close()

        f_before = ROOT.TFile.Open(dir_in + plt_name + "before.root")

        h_data_file = f_before.Get("data")
        if(h_data is None):
            h_data = h_data_file
            h_data.SetDirectory(0)
        else:
            h_data.Add(h_data_file)

        for l in labels:
            h = f_before.Get("h"+l)
            if(isinstance(h, ROOT.TH1)): 
                if(l in hists.keys()):
                    hists[l].Add(h)
                else:
                    hists[l] = h
                    h.SetDirectory(0)
        f_before.Close()





    h_tot_before = hists.items()[0][1].Clone("h_tot_before")
    h_tot_before.Reset()
    for k,h in hists.items(): 
        h_tot_before.Add(h)
    fname = options.outdir + plt_name + "comb.png"

    xstart, xstop, nbins, label = obs_attrs[obs]

    hist_list = [hists[l] for l in labels]
    outlines = [h_tot_after] if incl_after else []

    makeCan("temp", fname, [h_data], bkglist = [hist_list], totlist = [h_tot_before], signals = outlines, colors = colors, bkgNames = labels, titles = [""], logy = False, xtitle = label,
        drawSys = False, ratio_range = ratio_range, stack = True, draw_chi2 = draw_chi2, prelim = True, year = -1)

