import ROOT
import numpy as np

fname  = "W_RW_june9/ratio.root"
f = ROOT.TFile.Open(fname)
h = f.Get("ratio_nom")

ratios = []
errs = []
signifs = []
ref_val = 3.0

for i in range(h.GetNbinsX()+1):
    for j in range(h.GetNbinsY()+1):
        for k in range(h.GetNbinsZ()+1):
            c = h.GetBinContent(i,j,k)
            e = h.GetBinError(i,j,k)
            if( e > 0.):
                ratios.append(c)
                errs.append(e)
                signifs.append((c - ref_val)/e)

max_idx = np.argmax(signifs)
max_ridx = np.argmax(ratios)


print('max signif', signifs[max_idx], ratios[max_idx], errs[max_idx])
print('max ratio', signifs[max_ridx], ratios[max_ridx], errs[max_ridx])

ratios = np.array(ratios)
errs = np.array(errs)
signifs = np.array(signifs)

one_sig_aboves = signifs > 1.0
two_sig_aboves = signifs > 2.0

print("frac 1sigma above %.0f %.3f" % (ref_val, np.mean(one_sig_aboves)))
print("frac 2sigma above %.0f %.3f" % (ref_val, np.mean(two_sig_aboves)))

