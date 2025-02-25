import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.Utils import *
from utils.PlotUtils import *

out_dir = "plots/subjet_pt_dR_feb13/"

dsets = [
        ("plots/dRs/W_subjet_pt_dR.root", "W"),
        ("plots/dRs/top_subjet_pt_dR.root", "top"),
        ("plots/dRs/XYY_subjet_pt_dR.root", "XYY"),
        ("plots/dRs/Zp_subjet_pt_dR.root", "Zp"),
        ("plots/dRs/YtoHH_subjet_pt_dR.root", "YtoHH"),
        ]

colors = {
        "W": CMS_red,
        "top": CMS_lightblue,
        "XYY": CMS_purple,
        "Zp": CMS_brown,
        "YtoHH": CMS_orange,
        }

leg_labels = {
        "W":"t#bar{t} (W-matched)" ,
        "top":"t#bar{t} (t-matched)",
        #"XYY":"X #rightarrow YY, Y #rightarrow q#bar{q} ",
        #"Zp":"Z' #rightarrow T'T', T' #rightarrow tZ #rightarrow 6q ",
        #"YtoHH":"Y #rightarrow HH, H #rightarrow t#bar{t} #rightarrow 6q ",
        "XYY":"Y #rightarrow q#bar{q} ",
        "Zp":"T' #rightarrow tZ #rightarrow 5q ",
        "YtoHH":"H #rightarrow t#bar{t} #rightarrow 6q ",
        }

ext = ".pdf"


if(not os.path.exists(out_dir)):
    os.system("mkdir " + out_dir)

d_pts = []
d_dRs = []

dr_max = 0.
pt_max = 0.

for fname, label in dsets:
    print(fname)
    f = ROOT.TFile.Open(fname)
    h_subjet_pt = f.Get("h_%s_subjetpt" % label)
    h_dR = f.Get("h_%s_dRs" % label)
    h_dR.GetXaxis().SetRangeUser(0., 0.4)
    h_subjet_pt.SetDirectory(0)
    h_dR.SetDirectory(0)


    h_subjet_pt.Scale(1./h_subjet_pt.Integral())
    h_dR.Scale(1./h_dR.Integral())

    dr_max = max(h_dR.GetMaximum(), dr_max)
    pt_max = max(h_subjet_pt.GetMaximum(), pt_max)

    f.Close()

    d_pts.append((label, h_subjet_pt))
    d_dRs.append((label, h_dR))

x_start = 0.6
x_size = 0.3
y_size = 0.2 
y_end = 0.9


setTDRStyle()
c1 = ROOT.TCanvas("c1", "c1", 1400, 1200)
leg1 = ROOT.TLegend(x_start,y_end - y_size,x_start + x_size,y_end)
first = True


for label, h in d_pts:
    h.SetLineColor(colors[label])
    h.SetLineWidth(6)
    leg1.AddEntry(h, leg_labels[label], "l")
    if(first):
        h.SetMaximum(1.5 * pt_max)
        mLS = 0.07
        mTS = 0.05
        h.GetYaxis().SetTitleOffset(1.5)
        h.GetXaxis().SetTitleOffset(1.2)
        h.GetYaxis().SetTitle('Arbitrary Units')
        h.GetXaxis().SetTitle('Subjet p_{T}')
        #h.GetYaxis().SetLabelSize(mLS)
        #h.GetXaxis().SetLabelSize(mLS)
        h.GetYaxis().SetTitleSize(mTS)
        h.GetXaxis().SetTitleSize(mTS)
        h.GetYaxis().SetNdivisions(505)
        h.Draw("hist")
        first = False
    else:
        h.Draw("hist same")

leg1.SetBorderSize(0)
leg1.Draw()
extraText = "Simulation"
CMS_loc = 11
period = 0
CMS_lumi(c1, period, CMS_loc, writeExtraText=True, extraText=extraText)
c1.Print(out_dir + "pts" + ext)

x_start = 0.35
x_size = 0.45
y_size = 0.2 
y_end = 0.83

setTDRStyle()
c2 = ROOT.TCanvas("c2", "c2", 1400, 1200)
leg2 = ROOT.TLegend(x_start,y_end - y_size,x_start + x_size,y_end)
leg2.SetNColumns(2)
leg2.SetColumnSeparation(0.3)
leg2.SetTextSize(0.045)
leg2.SetTextFont(42)
first = True


for label, h in d_dRs:
    h.SetLineColor(colors[label])
    h.SetLineWidth(4)
    leg2.AddEntry(h, leg_labels[label], "l")
    if(first):
        h.SetMaximum(5.0 * dr_max)
        mLS = 0.07
        mTS = 0.065
        h.GetYaxis().SetTitleOffset(1.2)
        h.GetXaxis().SetTitleOffset(1.0)
        h.GetYaxis().SetTitle('Arbitrary units')
        h.GetXaxis().SetTitle('#DeltaR (gen quark, subjet)')
        #h.GetYaxis().SetLabelSize(mLS)
        #h.GetXaxis().SetLabelSize(mLS)
        h.GetYaxis().SetTitleSize(mTS)
        h.GetXaxis().SetTitleSize(mTS)
        h.GetYaxis().SetNdivisions(505)
        h.Draw("hist")
        first = False
    else:
        h.Draw("hist same")

ROOT.gPad.RedrawAxis();
c2.SetLogy()

c2.SetRightMargin(0.05)
leg2.SetBorderSize(0)
leg2.Draw()
CMS_lumi(c2, period, CMS_loc, writeExtraText=True, extraText=extraText, extraTextRight=True)
c2.Print(out_dir + "dRs" + ext)
