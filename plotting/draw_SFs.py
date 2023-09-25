import mplhep as hep
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
import numpy as np

fout = "SF_plot.pdf"

plt.style.use(hep.style.CMS)
plt.figure(figsize=(12,9))


legend_elems = []
y_compare = [0.85, 1.03, 0.89, 0.91]
y_compare_errs = [0.04, 0.06, 0.04, 0.02]
y_lund = [0.8, 0.91, 0.92, 0.82, 0.94]
y_lund_errs = [0.095, 0.08, 0.08, 0.15, 0.29]
y_triboson = [0.66]
y_triboson_err =  [0.35]
xs = np.array([1,2,3,4,5])
offset = 0.1

plt.errorbar(xs[:4] - offset, y_compare, yerr = y_compare_errs, fmt = 's', color = 'blue', label = "Standard Calibration Technique", capsize = 2.0) 
plt.errorbar([5 - offset], y_triboson, yerr = y_triboson_err, fmt = 'x', color = 'red', label = "CMS Triboson Search", capsize = 2.0) 
plt.errorbar(xs+offset, y_lund, yerr = y_lund_errs, fmt = 'o', color = 'green', label = "Lund Jet Plane Reweighting", capsize = 2.0) 
               

#plt.xlabel("Jet Type (Tagging Variable)", labelpad =20)
plt.ylabel("Correction Factor")
plt.gca().minorticks_off()
y_minor = matplotlib.ticker.MultipleLocator(1)
plt.gca().yaxis.set_minor_locator(y_minor)
plt.ylim(0.3, 1.5)
plt.xlim(0.5, 5.5)
plt.xticks([1,2,3,4,5], ["W \n" r"($\tau_{21}$)", "W\n(DeepAK8)", "W\n(DeepAK8-MD)", "top\n" r"($\tau_{32})$",
    r"R $\rightarrow$ WW" "\n(DeepAK8)"])
hep.cms.label( data = True, lumi = 138)

leg1 = plt.legend(loc = 'upper left')

plt.savefig(fout , bbox_inches="tight")
plt.close()
