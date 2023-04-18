from Utils import *
import os
import numpy as np
from array import array

def create_func(order):
    expr = ""
    for i in range(order + 1):
        expr += "[%i]  * x ^ (%i) " % (i,i)
        if(i < order): expr += "+ "
    print(expr)
    f = ROOT.TF1("f%i"%order, expr)
    #f.SetParameter(0, 1.0)
    return f



parser = input_options()
options = parser.parse_args()

print(options)
if(not os.path.exists(options.outdir)): os.system("mkdir %s" % options.outdir)


f_ratio = ROOT.TFile.Open(options.fin, "UPDATE")
h_nom = f_ratio.Get("ratio_nom")
h_up = f_ratio.Get("ratio_sys_tot_up")
h_down = f_ratio.Get("ratio_sys_tot_down")
f_ratio.mkdir("pt_extrap")
f_ratio.cd("pt_extrap")
f_ratio.Write()





diffs = []
x_bin1 = 2
x_bin2 = 5
base_order = 1

for j in range(1,h_nom.GetNbinsY()+1):
    for k in range(1,h_nom.GetNbinsZ()+1):
        c1 = h_nom.GetBinContent(x_bin1,j,k)
        c2 = h_nom.GetBinContent(x_bin2,j,k)
        e1 = h_nom.GetBinError(x_bin1,j,k)
        e2 = h_nom.GetBinError(x_bin2,j,k)
        if(  c1 > 1e-6 and c2 > 1e-6 and  abs(e1/c1) < 0.5 and (e2/c2) < 0.5):
            frac_diff = 2. * abs(c1 - c2)/(c1 + c2)
            diffs.append(frac_diff)

print("AVERAGE DIFF: %.3f" % np.mean(diffs))
ROOT.gStyle.SetOptFit(1110)

for h in [h_nom, h_up, h_down]:

    func_name_base = "func_"
    if(h is h_up): func_name_base += "sys_tot_up_"
    if(h is h_down): func_name_base += "sys_tot_down_"

    for j in range(1,h.GetNbinsY()+1):
        for k in range(1,h.GetNbinsZ()+1):
    #for j in [11]:
        #for k in [7]:
            x = array('d')
            ex = array('d')
            y = array('d')
            ey = array('d')
            nbin = 0
            for i in range(1,h.GetNbinsX()+1):
                eps = 1e-6

                c1 = h.GetBinContent(i,j,k)
                e1 = h.GetBinError(i,j,k)

                #unmeasured bins get ratio of 1 with 100% unc
                if( c1 < eps):
                    c1 = 1.0
                    e1 = 1.0
                else:
                    nbin += 1

                center = h.GetXaxis().GetBinCenter(i)
                width = h.GetXaxis().GetBinWidth(i)
                if(center > 300):
                    #center overflow bin on range of actual data
                    center = 400
                    width = 50


                x.append(center)
                ex.append(width)
                y.append(c1)
                ey.append(e1)

            if(len(x) > 1):

                g = ROOT.TGraphErrors(len(x), x, y, ex, ey)
                order = base_order
                chi2_prev = 999999

                #keep adding params based on F-test
                thresh = 0.01
                while(True):

                    func = create_func(order)
                    fit_res = g.Fit(func, "0 S +")
                    chi2_new = fit_res.Chi2()

                    F_num = max((chi2_prev - chi2_new), 0)
                    eps = 1e-8
                    F_denom = chi2_new/(nbin - order + eps)
                    F = F_num / (F_denom + eps)
                    F_prob = 1. - ROOT.TMath.FDistI(F, 1, nbin - order)

                    print(order, chi2_prev, chi2_new, F_prob)

                    if( (nbin - order) <= 1 or chi2_new / (nbin - order) < 1.1):
                        break
                    elif(order == base_order or (F_prob < thresh)):
                        order +=1
                        chi2_prev = chi2_new
                    else:
                        order -=1
                        break

                func = g.GetFunction("f%i" % order)
                f_ratio.cd("pt_extrap")
                func.SetName(func_name_base + ("%i_%i" % (j,k)) )
                func.Write()

                if(h is h_nom):
                    fout  = options.outdir + "fit_%i_%i.png" % (j, k)
                else:
                    fout  = options.outdir + "fit_%s_%i_%i.png" % (func_name_base, j, k)

                c = ROOT.TCanvas("c", "", 800, 800)
                g.SetTitle("Fit Bin %i,%i" %(j, k))
                g.GetYaxis().SetTitle("Correction Factor")
                g.GetYaxis().CenterTitle()
                g.GetXaxis().SetTitle("Subjet pT")
                g.GetXaxis().CenterTitle()
                g.Draw("AP")
                func.Draw("same")
                c.Print(fout)

f_ratio.Write()
f_ratio.Close()
