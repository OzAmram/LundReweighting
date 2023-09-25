import numpy as np
from Utils import *

def pgq(z):
    return (1. + (1- z)**2)/z

#jet_pt = [50,100,200,300,400,500,600, 700]
jet_pts = np.linspace(50, 1000, 100)

log_kts = np.array([-1., -0.5, 0., 0.5, 1., 1.5, 2.0])
colors = ['g', 'b', 'purple', 'r', 'orange', 'yellow', 'pink']

outdir = "analytic_densities/"
log_1odeltas = [1,2,3,4]

deltas = 0.8 / np.exp(log_1odeltas)

kts = np.exp(log_kts)

ys = []


size = 0.4

for j,delta in enumerate(deltas):
    log_1odelta = log_1odeltas[j]
    fig, ax = plt.subplots()

    for i, kt in enumerate(kts):
        if(log_1odelta > 3 and log_kts[i] > 1): continue
        y_kt = []
        for jet_pt in jet_pts:
            zbar = kt / (delta * jet_pt)
            y = zbar * ( pgq(zbar) + pgq ( 1 - zbar))
            y_kt.append(y)
        ys.append(y_kt)
        ax.plot(jet_pts, y_kt, label = "log(kt/GeV) = %.2f" % log_kts[i], color = colors[i])
        




    ax.set_xlabel("Jet pT", fontsize=14)
    ax.set_ylabel("Lund Plane Rel. Density @ LO", fontsize=14)
    plt.ylim(1.0, 3.0)
    ax.legend(loc='upper right')
    plt.title(r"log(0.8/$\Delta$) = " + str(log_1odelta))

    plt.savefig(outdir + "density_logdelta%i.png" % log_1odelta)

