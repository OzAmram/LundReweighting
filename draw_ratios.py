from PlotUtils import *
from Utils import *
import os
import sys

def style(h):

    mLS = 0.07
    mTS = 0.2

    x_label = "ln(0.8/#Delta)"
    y_label = "ln(kt/GeV)"

    #h.GetXaxis().SetTitleOffset(0.)
    #h.GetYaxis().SetTitleOffset(0.)
    h.GetXaxis().SetTitleSize(0.05)
    h.GetYaxis().SetTitleSize(0.05)
    h.GetZaxis().SetTitleSize(0.05)
    h.GetXaxis().SetTitle(x_label)
    h.GetYaxis().SetTitle(y_label)
    h.GetXaxis().CenterTitle(True)
    h.GetYaxis().CenterTitle(True)
    h.GetZaxis().CenterTitle(True)
    h.GetXaxis().SetRangeUser(0,8)
    h.GetYaxis().SetRangeUser(-4, 4)
    #_3dh.GetYaxis().SetLabelSize(mLS)
    #_3dh.GetXaxis().SetLabelSize(mLS)



def add_uncs(h1, h2):
    for i in range(1, h1.GetNbinsX() + 1):
        for j in range(1, h1.GetNbinsY() + 1):
            cont1 = h1.GetBinContent(i,j)
            cont2 = h2.GetBinContent(i,j)
            new_cont = (cont1**2 + cont2**2)**0.5
            h1.SetBinContent(i,j, new_cont)

def get_sys_unc_hist(h, up, down):
    h_unc = h.Clone(h.GetName() + "_sys_unc")
    for i in range(1, h.GetNbinsX() + 1):
        for j in range(1, h.GetNbinsY() + 1):
            cont = h.GetBinContent(i,j)
            sys_up = up.GetBinContent(i,j)
            sys_down = down.GetBinContent(i,j)
            err = 0.5 * (abs(cont - sys_up) + abs(cont - sys_down) )

            eps = 1e-4
            if(cont  < eps):
                h_unc.SetBinContent(i,j,0)
            else:
                h_unc.SetBinContent(i,j, err/cont)
                h_unc.SetBinError(i,j, 0.)
    return h_unc

parser = input_options()
options = parser.parse_args()

out_dir = options.outdir
if(not os.path.exists(out_dir)): os.system("mkdir " + out_dir)
tdrstyle.setTDRStyle()

fnames = ["W_RW_june9/ratio.root", "W_RW_2017_june30/ratio.root", "W_RW_2016_july1/ratio.root"]
#fnames = ["W_RW_june9/ratio.root", "W_RW_june9/ratio.root", "W_RW_june9/ratio.root"]
weights = [59.74, 41.4, 35.9]
#fname = "W_RW_2017_june30/ratio.root"

h_3d = h_sys_up = h_sys_down = None
for i,fname in enumerate(fnames):
    f = ROOT.TFile.Open(fname)
    h_3d_ = f.Get("ratio_nom")
    h_sys_up_ = f.Get("ratio_sys_tot_up")
    h_sys_down_ = f.Get("ratio_sys_tot_down")

    val = h_3d_.GetBinContent(3,5,11)
    err = h_3d_.GetBinError(3,5,11)
    val_sys = h_sys_up_.GetBinContent(3,5,11)
    print(val, err, val_sys)

    if(h_3d is None):
        h_3d = h_3d_
        h_sys_up = h_sys_up_
        h_sys_down = h_sys_down_

        h_3d.SetDirectory(0)
        h_sys_up.SetDirectory(0)
        h_sys_down.SetDirectory(0)


        h_3d.Scale(weights[i])
        h_sys_up.Scale(weights[i])
        h_sys_down.Scale(weights[i])
    else:
        h_3d.Add(h_3d_, weights[i])
        h_sys_up.Add(h_sys_up_, weights[i])
        h_sys_down.Add(h_sys_down_, weights[i])
    f.Close()

weight_sum = np.sum(weights)
h_3d.Scale(1./weight_sum)
h_sys_up.Scale(1./weight_sum)
h_sys_down.Scale(1./weight_sum)

val = h_3d.GetBinContent(3,5,11)
err = h_3d.GetBinError(3,5,11)
val_sys = h_sys_up.GetBinContent(3,5,11)
print(val, err, val_sys)


    


#set colz colors
ROOT.gStyle.SetPalette(ROOT.kViridis)


pt_bin_labels = [ "0-50", "50-100", "100-175", "175-250", "250-350", "350+"]

pt_bins = array('f', [0., 50., 100., 175., 250., 350., 99999.])


latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextAngle(0)
latex.SetTextColor(ROOT.kBlack)    
latex.SetTextSize(0.04)

latex.SetTextFont(42)
latex.SetTextAlign(21) 


for i in range(1,h_3d.GetNbinsX()+1):
    h_3d_clone = h_3d.Clone("h%i"%i)
    h_3d_clone.GetXaxis().SetRange(i,i)
    h_ratio_proj = h_3d_clone.Project3D("zy")
    h_ratio_unc = get_unc_hist(h_ratio_proj)

    h_sys_up_clone = h_sys_up.Clone("s_up%i"%i)
    h_sys_up_clone.GetXaxis().SetRange(i,i)
    h_sys_up_proj = h_sys_up_clone.Project3D("zy")

    h_sys_down_clone = h_sys_down.Clone("s_down%i"%i)
    h_sys_down_clone.GetXaxis().SetRange(i,i)
    h_sys_down_proj = h_sys_down_clone.Project3D("zy")

    h_ratio_sys_unc = get_sys_unc_hist(h_ratio_proj, h_sys_up_proj, h_sys_down_proj)
    add_uncs(h_ratio_unc, h_ratio_sys_unc)
    

    cleanup_ratio(h_ratio_proj, h_min =0., h_max = 2.0)
    cleanup_ratio(h_ratio_unc, h_min = 0., h_max = 1.0)


    style(h_ratio_proj)
    style(h_ratio_unc)

    CMS_lumi.writeExtraText = True
    CMS_loc = 0
    period = -1


    posX = 0.58
    posY = 0.88
    shiftY = 0.04

    ROOT.gStyle.SetHatchesLineWidth(1)
    ROOT.gStyle.SetHatchesSpacing(0.7)

    h_ratio_unc.SetLineColor(ROOT.kBlack)
    h_ratio_unc.SetFillColor(ROOT.kGray+1)
    h_ratio_unc.SetFillStyle(3354)
    h_ratio_unc.SetMarkerStyle(20)
    h_ratio_unc.SetMarkerSize(0.01)

    h_ratio_sys_unc.SetLineColor(ROOT.kWhite)
    h_ratio_sys_unc.SetFillColor(ROOT.kMagenta)
    h_ratio_sys_unc.SetFillStyle(3345)
    h_ratio_sys_unc.SetMarkerStyle(20)
    h_ratio_sys_unc.SetMarkerSize(0.01)


    h_ratio_proj.SetFillColor(ROOT.kTeal-7)
    h_ratio_proj.SetLineColor(ROOT.kTeal-7)
    c_ratio = ROOT.TCanvas("c_unc", "", 1000, 800)
    h_ratio_proj.Draw("colz")
    h_ratio_unc.Draw("BOX same")
    #h_ratio_sys_unc.Draw("BOX same")

    ratio_plot_label = "#bf{Data/MC Ratio}" 
    unc_plot_label = "#bf{Data/MC Ratio Frac. Unc.}"
    subj_label = "#bf{Subjet p_{T} %s GeV}" % (pt_bin_labels[i-1])

    #latex.DrawLatex(posX, posY, ratio_plot_label)
    #h_ratio_proj.GetZaxis().SetTitle(ratio_plot_label)

    latex.DrawLatex(posX, posY, subj_label)

    leg = ROOT.TLegend(0.47, 0.75, 0.75, 0.85)
    leg.SetTextSize(0.03)
    leg.AddEntry(h_ratio_proj, "Data/MC Ratio", "f")
    leg.AddEntry(h_ratio_unc, "Uncertainty", "f")
    leg.SetBorderSize(0)
    leg.Draw()

    c_ratio.SetRightMargin(0.2)
    CMS_lumi.CMS_lumi(c_ratio, period, CMS_loc)
    c_ratio.Print(out_dir + "lundPlane_bin%i_ratio.png" % i)


    #c_ratio_unc = ROOT.TCanvas("c_unc", "", 800, 800)
    #h_ratio_unc.Draw("colz")
    #latex.DrawLatex(posX, posY, unc_plot_label)
    #latex.DrawLatex(posX, posY-shiftY, subj_label)
    #c_ratio_unc.SetRightMargin(0.2)
    #CMS_lumi.CMS_lumi(c_ratio_unc, period, CMS_loc)
    #c_ratio_unc.Print(out_dir + "lundPlane_bin%i_ratio_unc.png" % i)



