from PlotUtils import *
from Utils import *
import os
import sys

out_dir = "subjet_pt_dR_july14/"

dsets = [
        ("W_SF_may5/subjet_pt_dR.root", "W"),
        ("top_SF_may5/subjet_pt_dR.root", "top"),
        ("XYY_SF_may10/subjet_pt_dR.root", "XYY"),
        ("Zp_SF_may10/subjet_pt_dR.root", "Zp"),
        ("YToHH_SF_may10/subjet_pt_dR.root", "YtoHH"),
        ]

colors = {
        "W": ROOT.kRed-7,
        "top": ROOT.kBlue-7,
        "XYY": ROOT.kMagenta-1,
        "Zp": ROOT.kGreen-6,
        "YtoHH": ROOT.kOrange-3,
        }

leg_labels = {
        "W":"t#bar{t} : W-matched" ,
        "top":"t#bar{t} : t-matched",
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

tdrstyle.setTDRStyle()
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
CMS_lumi.writeExtraText = True
CMS_lumi.extraText = "Simulation Preliminary"
CMS_loc = 11
period = 0
CMS_lumi.CMS_lumi(c1, period, CMS_loc)
c1.Print(out_dir + "pts" + ext)


tdrstyle.setTDRStyle()
c2 = ROOT.TCanvas("c2", "c2", 1400, 1200)
leg2 = ROOT.TLegend(x_start,y_end - y_size,x_start + x_size,y_end)
first = True


for label, h in d_dRs:
    h.SetLineColor(colors[label])
    h.SetLineWidth(6)
    leg2.AddEntry(h, leg_labels[label], "l")
    if(first):
        h.SetMaximum(3.5 * dr_max)
        mLS = 0.07
        mTS = 0.05
        h.GetYaxis().SetTitleOffset(1.5)
        h.GetXaxis().SetTitleOffset(1.2)
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

c2.SetLogy()

c2.SetRightMargin(0.05)
leg2.SetBorderSize(0)
leg2.Draw()
CMS_lumi.CMS_lumi(c2, period, CMS_loc)
c2.Print(out_dir + "dRs" + ext)
