import mplhep as hep
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.PlotUtils import *



fout = "plots/SF_plot.png"

plt.style.use(hep.style.CMS)
plt.figure(figsize=(12,9))


legend_elems = []
labels = ["W \n" r"($\tau_{21}$)", "W\n(DeepAK8-MD)", "W\n(ParticleNet)", "top\n" r"($\tau_{32})$", r"R $\rightarrow$ WW" "\n(DeepAK8)"]
y_compare =         [0.85, 0.86, 0.99, 0.91]
y_compare_errs =    [0.04, 0.06, 0.06, 0.02 ] 
y_lund =           [0.85, 0.88, 0.94, 0.86,  0.99, ]
y_lund_errs_up =   [0.06, 0.05, 0.03, 0.13,  0.17, ]
y_lund_errs_down = [0.12, 0.07, 0.08, 0.21,  0.17, ]
y_triboson = [0.66]
y_triboson_err =  [0.35]
xs = np.array([1,2,3,4,5])
offset = 0.1

lund_errs = np.array(list(zip(y_lund_errs_up, y_lund_errs_down))).T

plt.errorbar(xs[:4] - offset, y_compare, yerr = y_compare_errs, fmt = 's', color = c_lightblue, label = "Standard Calibration Technique", capsize = 2.0) 
plt.errorbar([5 - offset], y_triboson, yerr = y_triboson_err, fmt = 'x', color = c_red, label = "CMS Triboson Search", capsize = 2.0) 
plt.errorbar(xs+offset, y_lund, yerr = lund_errs, fmt = 'o', color = c_purple, label = "Lund Jet Plane Reweighting", capsize = 2.0) 
               

plt.xlabel("Jet Type (Tagging Variable)", labelpad =20)
plt.ylabel("Correction Factor")
plt.gca().minorticks_off()
y_minor = matplotlib.ticker.MultipleLocator(1)
plt.gca().yaxis.set_minor_locator(y_minor)
plt.ylim(0.3, 1.5)
plt.xlim(0.5, 5.5)
plt.xticks([1,2,3,4,5], labels)
hep.cms.label( data = True, lumi = 138)

leg1 = plt.legend(loc = 'upper left')

plt.savefig(fout , bbox_inches="tight")
plt.close()
