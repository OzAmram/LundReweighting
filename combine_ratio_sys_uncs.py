from Utils import *
import os

parser = input_options()
options = parser.parse_args()

print(options)

#f_ratio = ROOT.TFile.Open("ttbar_UL_top_rw_feb13/ratio.root", "UPDATE")
#f_ratio = ROOT.TFile.Open("ttbar_UL_feb14_W_rw/ratio.root", "UPDATE")
#f_ratio = ROOT.TFile.Open("ttbar_UL_feb14_W_rw_chg_only/ratio.root", "UPDATE")
f_ratio = ROOT.TFile.Open(options.fin, "UPDATE")

sys_list = list(sys_weights_map.keys())
sys_list.remove('nom_weight')
sys_ratios = []

h_nom = ROOT.TH3F(f_ratio.Get("ratio_nom"))
h_sys_up = ROOT.TH3F(h_nom.Clone("ratio_sys_tot_up"))
h_sys_down = ROOT.TH3F(h_nom.Clone("ratio_sys_tot_down"))
h_sys_up.Reset()
h_sys_down.Reset()

h_nom.Print()
h_sys_up.Print()


for sys in sys_list:
    h = f_ratio.Get("ratio_" + sys)
    sys_ratios.append(h)

for i in range(h.GetNbinsX()+1):
    for j in range(h.GetNbinsY()+1):
        for k in range(h.GetNbinsZ()+1):
            c_nom = h_nom.GetBinContent(i,j,k)
            c_err_up = 0.
            c_err_down = 0.
            #sum diffs (aka sys unc) in quadrature
            for h_sys in sys_ratios:
                c = h_sys.GetBinContent(i,j,k)
                if(c > c_nom): c_err_up += (c - c_nom)**2
                if(c < c_nom): c_err_down += (c - c_nom)**2

            c_err_up = c_err_up**0.5
            c_err_down = c_err_down**0.5

            eps = 1e-8
            h_sys_up.SetBinContent(i,j,k, c_nom + c_err_up)
            h_sys_down.SetBinContent(i,j,k, max(c_nom - c_err_down, eps))
            print(i,j,k, c_nom, c_nom + c_err_up, c_nom - c_err_down)

print("MEANs are:", h_nom.GetMean(), h_sys_up.GetMean(), h_sys_down.GetMean())
f_ratio.cd()
h_sys_up.Write()
h_sys_down.Write()
f_ratio.Close()

