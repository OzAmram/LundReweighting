from Utils import *
import os

parser = input_options()
options = parser.parse_args()

ext = ".png"


ratio_range = (0.2, 1.8)
if(not os.path.exists(options.outdir)):
    os.system("mkdir " + options.outdir)


incl_after = True
if(options.mode == 'top'):
    dirs = ["top_SF_2016_july14/", "top_SF_2017_july14/", "top_SF_2018_july14/"]
    labels = [
            #"Single Top",
            "W+Jets + QCD",
            "tW",
            "t#bar{t} : W-matched",
            "t#bar{t} : unmatched",
            "t#bar{t} : t-matched",
        ]
    colors = [
        #ROOT.kMagenta-1,
        ROOT.kOrange-3,
        ROOT.kYellow-7,
        ROOT.kRed -7,
        ROOT.kGreen-6,
        ROOT.kBlue -7, 
        ]


    obs_attrs = {
            'mSoftDrop' : (125, 225, 25, "m_{SD} [GeV]", "Events / 4 GeV"),
            'tau21' : (0.05, 0.8, 15, "#tau_{21}", "Events / 0.05" ),
            'tau32' : (0.2, 0.95, 15, "#tau_{32}", "Events / 0.05"),
            'tau43' : (0.6, 0.96, 18, "#tau_{43}", "Events / 0.02"),
            'nPF' : (20.5, 120.5, 25, "Num. PF Cands.", "Events / 4"),
            'pt' : (500, 1200., 20, "p_{T}", ""),
            }
    draw_chi2 = True


elif(options.mode == "full"):
    dirs = ["W_fullrange_2016/", "W_fullrange_2017/", "W_fullrange_2018/"]
    labels = [
            #"Diboson",
            #"Single Top",
            "W+Jets + QCD",
            "t#bar{t} : t-matched",
            "t#bar{t} : unmatched",
            "tW",
            "t#bar{t} : W-matched",
        ]

    colors = [
        #ROOT.kCyan,
        #ROOT.kMagenta-1,
        ROOT.kOrange-3,
        ROOT.kBlue -7, 
        ROOT.kGreen-6,
        ROOT.kYellow-7,
        ROOT.kRed -7,
        ]

    obs_attrs = {
            'mSoftDrop' : (50, 230, 45, "m_{SD} [GeV]", "Events / 4 GeV"),
            }
    draw_chi2 = False
    incl_after = False

else:
    dirs = ["W_SF_2016_july14/", "W_SF_2017_july14/", "W_SF_2018_july14/"]
    labels = [
            #"Diboson",
            #"Single Top",
            "W+Jets + QCD",
            "t#bar{t} : t-matched",
            "t#bar{t} : unmatched",
            "tW",
            "t#bar{t} : W-matched",
        ]

    colors = [
        #ROOT.kCyan,
        #ROOT.kMagenta-1,
        ROOT.kOrange-3,
        ROOT.kBlue -7, 
        ROOT.kGreen-6,
        ROOT.kYellow-7,
        ROOT.kRed -7,
        ]

    obs_attrs = {
            'mSoftDrop' : (60, 110, 25, "m_{SD} [GeV] ", "Events / 2 GeV"),
            'tau21' : (0.05, 0.8, 25, "#tau_{21}", "Events / 0.03"),
            'tau32' : (0.4, 0.95, 25, "#tau_{32}", "Events / 0.022"),
            'tau43' : (0.6, 0.96, 18, "#tau_{43}", "Events / 0.02"),
            'nPF' : (20.5, 80.5, 30, "Num. PF Cands.", "Events / 2" ),
            'pt' : (225, 825., 20, "p_{T}", "Events / 30 GeV"),
            'DeepAK8_W' : (0., 1., 20, "DeepAK8 (W vs. QCD)", "Events / 0.05"),
            'DeepAK8_W_MD' : (0., 1., 20, "DeepAK8-MD (W vs. QCD)", "Events / 0.05"),
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
                h_tot_after.SetName("Total Reweighted Sim.")
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
    fname = options.outdir + plt_name + "comb" + ext

    xstart, xstop, nbins, label, ylabel = obs_attrs[obs]


    hist_list = [hists[l] for l in labels]
    #for h in hist_list:
    #    h.GetXaxis().SetRangeUser(xstart, xstop)
    #h_tot_after.GetXaxis().SetRangeUser(xstart, xstop)
    #h_tot_before.GetXaxis().SetRangeUser(xstart, xstop)
    #h_data.GetXaxis().SetRangeUser(xstart, xstop)
    outlines = [h_tot_after] if incl_after else []
    labels_new = [l if "Single" not in l else "Single t" for l in labels]
    


    makeCan("temp", fname, [h_data], bkglist = [hist_list], totlist = [h_tot_before], signals = outlines, colors = colors, bkgNames = labels_new, titles = [""], logy = False, xtitle = label,
        ytitle = ylabel, drawSys = False, ratio_range = ratio_range, stack = True, draw_chi2 = draw_chi2, prelim = True, year = -1)

