import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.Utils import *
parser = input_options()
options = parser.parse_args()


def get_neigh_avg(h, i,j,k, max_j, min_k):
    #compute average value of neighbors
    avg = 0.
    w_sum = 0.
    pairs = [(j-1, k-1), (j-1,k), (j-1, k+1), 
            (j, k-1), (j, k+1), 
            (j+1, k-1), (j, k+1), (j+1, k+1)]
    eps = 1e-4
    for jp,kp in pairs:
        #if(jp > max_j or kp < min_k): continue
        val = min(h.GetBinContent(i,jp,kp), 5.0)
        err = h.GetBinError(i,jp,kp)
        if(err < eps): continue
        avg += (1./err**2) * val
        w_sum += (1./err**2)

    x_out = avg / (w_sum+eps)
    return x_out
        



print(options)
if(not os.path.exists(options.outdir)): os.system("mkdir " + options.outdir)

f_ratio = ROOT.TFile.Open(options.fin, "READ")

h = ROOT.TH3F(f_ratio.Get("ratio_nom"))

max_j = h.GetYaxis().FindBin(np.log(0.8/0.005))
min_k = h.GetZaxis().FindBin(np.log(0.02))
print('maxj, min_k', max_j, min_k)
eps = 1e-4
h_pulls = ROOT.TH1F("h_pulls","Distribution of Pulls", 50, -4, 4)
ROOT.gStyle.SetOptStat(False)

for i in range(1, h.GetNbinsX()+1):
    clone1 = h.Clone("dummy")
    clone1.GetXaxis().SetRange(i,i)
    h_out = clone1.Project3D("zy")
    h_out.Reset()
    for j in range(h.GetNbinsY()+1):
        for k in range(h.GetNbinsZ()+1):
            if(j > max_j or k < min_k): continue
            c_nom = min(h.GetBinContent(i,j,k), 5.0)
            e_nom = h.GetBinError(i,j,k)
            if(e_nom < eps): continue
            e_nom = max(e_nom, 0.1)

            avg_neighs = get_neigh_avg(h, i, j, k, max_j, min_k)
            pull = (avg_neighs - c_nom) / e_nom
            h_out.SetBinContent(j,k, abs(pull))
            if(pull > 1.5): print(i,k,j, pull, c_nom, e_nom, avg_neighs)
            h_pulls.Fill(pull)

    h_out.SetTitle("Outliers PT Bin %i" % i)
    c_mc = ROOT.TCanvas("c", "", 1000, 1000)
    h_out.Draw("colz")
    c_mc.SetRightMargin(0.2)
    c_mc.Print(options.outdir + "lundPlane_outliers_bin%i.png" % i)

ROOT.gStyle.SetOptStat(True)
ROOT.gStyle.SetOptFit(1111)
c_mc = ROOT.TCanvas("c", "", 1000, 1000)
fitRes = h_pulls.Fit("gaus", "S")
fitRes.Print()
print("Mean", h_pulls.GetMean())
print("Std", h_pulls.GetStdDev())
h_pulls.GetXaxis().SetTitle("Bin Pulls")
c_mc.Print(options.outdir + "overall_pulls.png" )

