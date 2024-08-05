import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.Utils import *
import numpy as np
from array import array

def transform(val):
    return 1./ val

def create_func(order, label = ""):
    expr = ""
    for i in range(order + 1):
        expr += "[%i]  * x ^ (%i) " % (i,i)
        if(i < order): expr += "+ "
    print(expr)
    f = ROOT.TF1("%s_%i"% (label, order), expr)
    #f.SetParameter(0, 1.0)
    return f



parser = input_options()
parser.add_argument("--herwig", default=False, action='store_true',  help="Pythia herwig ratio (less sys)")
options = parser.parse_args()

non_const_plots = False

print(options)
if(not os.path.exists(options.outdir)): os.system("mkdir %s" % options.outdir)
else: os.system("rm %s/*" % options.outdir)


f_ratio = ROOT.TFile.Open(options.fin, "UPDATE")
h_nom = f_ratio.Get("ratio_nom")
h_up = f_ratio.Get("ratio_sys_tot_up")
h_down = f_ratio.Get("ratio_sys_tot_down")
f_ratio.mkdir("pt_extrap")
f_ratio.cd("pt_extrap")
f_ratio.Write()





max_order = 1
write_out = True
diffs = []
x_bin1 = 2
x_bin2 = 5
base_order = 0

ROOT.gStyle.SetOptFit(1110)
order_dict = {0:0, 1:0, 2:0, 3:0, 4:0}
lin_orders = []

ratio_max = 20.0

for h in [h_nom, h_up, h_down]:
#for h in [h_nom]:

    func_name_base = "func_"
    if(not isinstance(h, ROOT.TH3)): 
        print("Skipping")
        continue
    if(h is h_up): func_name_base += "sys_tot_up_"
    if(h is h_down): func_name_base += "sys_tot_down_"

    for j in range(1,h.GetNbinsY()+1):
        for k in range(1,h.GetNbinsZ()+1):
    #for j in [10]:
        #for k in [6]:
            x = array('d')
            ex_up = array('d')
            ex_down = array('d')
            y = array('d')
            ey = array('d')
            nbin = 0
            for i in range(1,h.GetNbinsX()+1):
                eps = 1e-4

                c1 = h.GetBinContent(i,j,k)
                e1 = h.GetBinError(i,j,k)

                c1 = min(c1, ratio_max)
                e1 = min(e1, ratio_max)

                #unmeasured bins get ratio of 1 with 100% unc
                if( c1 < eps):
                    c1 = 1.0
                    e1 = 1.0
                else:
                    nbin += 1

                center = h.GetXaxis().GetBinCenter(i)
                width = h.GetXaxis().GetBinWidth(i)
                if(center >= 290):
                    #center overflow bin on range of actual data
                    center = 350
                    width = 50


                x_val = transform(center)
                xmin = transform(max(center - width/2., 20.))
                xmax = transform(center + width/2.)
                x_err_up = abs(x_val - xmax)
                x_err_down = abs(x_val - xmin)

                if(xmin > xmax): x_err_up, x_err_down = x_err_down, x_err_up

                x.append(x_val)
                #ex_up.append(x_err_up)
                #ex_down.append(x_err_down)
                ex_up.append(0.0)
                ex_down.append(0.0)
                y.append(c1)
                ey.append(e1)

            if(len(x) > 1):
                print(x,y,ex_up, ey)

                g = ROOT.TGraphAsymmErrors(len(x), x, y, ex_down, ex_up, ey, ey)
                order = base_order
                chi2_prev = 999999

                #keep adding params based on F-test
                thresh = 0.05
                fit_label = "%s_bin%i_%i" %  (h.GetName(), j,k)

                while(True):
                    g_clone = g.Clone(g.GetName() + "clone")

                    func = create_func(order, fit_label)
                    fit_res = g_clone.Fit(func, "0 S E EX0 +")
                    chi2_new = fit_res.Chi2()

                    #corrs = fit_res.GetCorrelationMatrix()
                    #corrs.Print()

                    F_num = max((chi2_prev - chi2_new), 0)
                    eps = 1e-8
                    F_denom = chi2_new/(nbin - order + eps)
                    F = F_num / (F_denom + eps)
                    F_prob = 1. - ROOT.TMath.FDistI(F, 1, nbin - order)

                    chi2_prob = ROOT.Math.chisquared_cdf_c(chi2_new, nbin - order)
                    print(nbin, order, chi2_prev, chi2_new, F_prob, chi2_prob)

                    if(order == base_order or (F_prob < thresh)):
                        #This order is preferred
                        #if( (nbin - order) <= 1 or chi2_new / (nbin - order) < 1.1):
                        if( (nbin - order) <= 1 or chi2_prob > 0.3 or order >= max_order):
                        #if( (nbin - order) <= 1):
                            #stop now
                            break

                        else:
                            #try higher order
                            order +=1
                            chi2_prev = chi2_new
                    else:
                        #this order not preferred, go back one
                        order -=1
                        break


                if(h is h_nom and nbin > 2):
                    order_dict [order] += 1
                if(order == 1):
                    lin_orders.append((j,k))


                if(write_out):
                    func = create_func(order, fit_label)
                    fit_res = g.Fit(func, "0 S +")
                    cov = fit_res.GetCovarianceMatrix()
                    covar_val = ROOT.TMatrixDRow(cov, 0)(1)
                    covar = ROOT.TParameter("float")(func_name_base +"covar_%i_%i" % (j,k), covar_val)

                    print('covariance :', covar_val)

                    f_ratio.cd("pt_extrap")
                    func.SetName(func_name_base + ("%i_%i" % (j,k)) )
                    covar.Write()
                    func.Write()

                    if( (not non_const_plots) or order >=1):

                        if(h is h_nom):
                            fout  = options.outdir + "fit_%i_%i.png" % (j, k)
                        else:
                            continue
                            fout  = options.outdir + "fit_%s_%i_%i.png" % (func_name_base, j, k)

                        c = ROOT.TCanvas("c", "", 800, 800)
                        g.SetTitle("Fit Bin %i,%i" %(j, k))
                        g.GetYaxis().SetTitle("Correction Factor")
                        g.GetYaxis().CenterTitle()
                        g.GetXaxis().SetTitle("1 / Subjet pT")
                        g.GetXaxis().CenterTitle()
                        g.GetXaxis().SetRangeUser(0., 0.1)
                        g.Draw("AP")
                        func.Draw("same")
                        c.Print(fout)
                
    if(h is h_nom):
        h_clone1 = h_nom.Clone("clone" )
        h_clone1.GetXaxis().SetRange(i,i)
        h_slopes = h_clone1.Project3D("zy")
        h_slopes.Reset()


        for (j,k) in lin_orders:

            f_str = "pt_extrap/func_%i_%i" % (j,k)
            func = f_ratio.Get(f_str)
            slope = func.GetParameter(1)
            h_slopes.SetBinContent(k,j, slope)

        c = ROOT.TCanvas("c", "", 800, 800)
        h_slopes.Draw("colz")
        h_slopes.SetTitle("Slopes of pt Extrap Fits")
        c.Print(options.outdir + "slope_summary.png")

print("Summary of functional orders:")
print(order_dict)
print("Linear order keys")
print(lin_orders)
f_ratio.Write()
f_ratio.Close()




