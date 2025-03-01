import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.Utils import *
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches




parser = input_options()
options = parser.parse_args()

def add_quad(v1,v2):
    return (v1**2 + v2**2)**0.5

# Function to plot error boxes
def makeErrorBoxes(xdata,ydata,xerror,yerror,fc='r',ec='None',alpha=0.5, label=''):

    # Create list for all the error patches
    errorboxes = []

    # Loop over data points; create box from errors at each point
    for xc,yc,xe,ye in zip(xdata,ydata,xerror,yerror):
        rect = Rectangle((xc-xe/2,yc-ye[0]), xe, (ye[0] + ye[1]))
        errorboxes.append(rect)

    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes,facecolor=fc,alpha=alpha,edgecolor=ec)

    return pc


def make_graph(xs, xerrs, ys, yerrs_stat, yerrs_sys, xlabel = "", title= "", fname = ""):

    hep.style.use("CMS")

    fig_size = (10,10)
    fig, ax = plt.subplots(figsize=fig_size)

    fontsize = 28

    yerrs_tot = [ [add_quad(yerrs_sys[i][0], yerrs_stat[i]), add_quad(yerrs_sys[i][1], yerrs_stat[i])] for i in range(len(yerrs_sys))]

    yerrs_tot = np.transpose(yerrs_tot, axes=[1,0])
    h3 = plt.errorbar(xs, ys, xerrs=None, yerr=yerrs_tot, fmt='o', color = 'black', ecolor = 'blue', label = 'Total error', linewidth = 5, markersize = 0, capsize = 15)

    #Patch version
    #alpha = 0.3
    #pc = makeErrorBoxes(xs,ys, xerrs, yerrs_tot, fc = 'blue', ec = 'None', alpha = alpha, label = 'Total Error')
    #ax.add_collection(pc)

    #h3 = plt.plot([1.0], [-999], color = 'blue', alpha = alpha)
    #h3 = mpatches.Patch(color='blue', alpha = alpha, label='Total error')




    h1 = plt.scatter(xs, ys, marker='o', color = 'black', s = 100, label = "Data/Sim. ratio" )
    
    h2 = plt.errorbar(xs, ys, xerr=None, yerr=yerrs_stat, fmt='o', color = 'black', ecolor = 'black', label = 'Statistical error', linewidth = 5, markersize = 0, capsize=10)




    handles = [h1,h2, h3]


    leg = plt.legend(handles=handles, loc='upper center', title=title, fontsize = fontsize, title_fontsize=fontsize)

    y_max = ax.get_ylim()[1] * 1.5
    plt.ylim([0, 2.5])

    #y_val = 0.82 * (ax.get_ylim()[1] - ax.get_ylim()[0]) + ax.get_ylim()[0]
    #x_val = 0.44 * (ax.get_xlim()[1] - ax.get_xlim()[0]) + ax.get_xlim()[0]
    #plt.text(x_val, y_val, title, horizontalalignment = 'left', fontweight = 'bold', fontsize = 18)

    plt.xlabel(xlabel, fontsize = fontsize*1.2)
    plt.ylabel("Lund plane data/sim. ratio", fontsize = fontsize*1.2)

    ax.tick_params(axis='both', which='major', labelsize=fontsize*0.9)


    hep.cms.lumitext(ax=ax, text=r"138 fb$^{-1}$                ")
    text = ""
    #text = "Preliminary"
    
    hep.cms.label(text, ax=ax, loc=0, data = True)

    if(fname != ""): 
        plt.savefig(fname, bbox_inches = 'tight')
        print("saving fig %s" %fname)





print(options)
if(not os.path.exists(options.outdir)): os.system("mkdir " + options.outdir)


fnames = ["data/ratio_2018.root", "data/ratio_2017.root", "data/ratio_2016.root"]
weights = [59.74, 41.4, 35.9]

h_nom = h_sys_up = h_sys_down = None
for i,fname in enumerate(fnames):
    f = ROOT.TFile.Open(fname)
    h_3d_ = f.Get("ratio_nom")
    h_sys_up_ = f.Get("ratio_sys_tot_up")
    h_sys_down_ = f.Get("ratio_sys_tot_down")

    if(h_nom is None):
        h_nom = h_3d_
        h_sys_up = h_sys_up_
        h_sys_down = h_sys_down_

        h_nom.SetDirectory(0)
        h_sys_up.SetDirectory(0)
        h_sys_down.SetDirectory(0)


        h_nom.Scale(weights[i])
        h_sys_up.Scale(weights[i])
        h_sys_down.Scale(weights[i])
    else:
        h_nom.Add(h_3d_, weights[i])
        h_sys_up.Add(h_sys_up_, weights[i])
        h_sys_down.Add(h_sys_down_, weights[i])
    f.Close()

weight_sum = np.sum(weights)
h_nom.Scale(1./weight_sum)
h_sys_up.Scale(1./weight_sum)
h_sys_down.Scale(1./weight_sum)

pt_bins_reg = array('f', [15., 65., 110., 175., 240., 300., 99999.])

max_j = h_nom.GetYaxis().FindBin(np.log(0.8/0.005))
min_k = h_nom.GetZaxis().FindBin(np.log(0.02))
print('maxj, min_k', max_j, min_k)
eps = 1e-4
h_pulls = ROOT.TH1F("h_pulls","Distribution of Pulls", 50, -4, 4)
ROOT.gStyle.SetOptStat(False)

#pt_bin_start = 4
#pt_bin_end = 5
pt_bin_start = 1
pt_bin_end = h_nom.GetNbinsX()+1

for i in range(pt_bin_start, pt_bin_end):
    clone1 = h_nom.Clone("dummy")
    clone1.GetXaxis().SetRange(i,i)
    pt_low, pt_high = pt_bins_reg[i-1], pt_bins_reg[i]
    for k in range(h_nom.GetNbinsZ()+1):
        xs = []
        x_errs = []
        vals = []
        errs_stat = []
        errs_sys = []
        ylow = h_nom.GetZaxis().GetBinLowEdge(k)
        yhigh = h_nom.GetZaxis().GetBinUpEdge(k)
        title = r"Subjet p$_T$ %.0f $-$ %.0f, ln(k$_T$) %.2f $-$ %.2f" % (pt_low, pt_high, ylow, yhigh)
        for j in range(h_nom.GetNbinsY()+1):
            if(j > max_j or k < min_k): continue
            c_nom = min(h_nom.GetBinContent(i,j,k), 5.0)
            e_nom = h_nom.GetBinError(i,j,k)
            if(e_nom < eps): continue
            sys1 = h_sys_up.GetBinContent(i,j,k) - c_nom
            sys2 = h_sys_down.GetBinContent(i,j,k) - c_nom
            sys_up = abs(max(sys1,sys2))
            sys_down = abs(max(sys1,sys2))

            xs.append(h_nom.GetYaxis().GetBinCenter(j))
            x_errs.append(h_nom.GetYaxis().GetBinWidth(j)/2.0)
            vals.append(c_nom)
            errs_stat.append(e_nom)
            errs_sys.append([sys_up, sys_down])

        if(len(xs) > 0):
            make_graph(xs, x_errs, vals, errs_stat, errs_sys, title = title, fname = options.outdir + "pt_bin%i_bin%i.png"  % (i,k), xlabel = r"ln(0.8/$\Delta$)")
            make_graph(xs, x_errs, vals, errs_stat, errs_sys, title = title, fname = options.outdir + "pt_bin%i_bin%i.pdf"  % (i,k), xlabel = r"ln(0.8/$\Delta$)")


