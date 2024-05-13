import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import gridspec
from sklearn.metrics import roc_curve,auc
import scipy.stats
import numpy as np
import ROOT
import CMS_lumi, tdrstyle
import copy
from array import array

#colors from CAT
#https://gitlab.cern.ch/cms-analysis/analysisexamples/plotting-demo/-/blob/master/1-tutorial_CAT_recommendations.ipynb?ref_type=heads

c_lightblue = "#5790fc"
c_orange = "#f89c20"
c_red = "#e42536"
c_purple = "#964a8b"
c_grey = "#9c9ca1"
c_indigo = "#7a21dd"

fig_size = (12,9)

def add_patch(legend, patch, name):
    from matplotlib.patches import Patch
    ax = legend.axes

    handles, labels = ax.get_legend_handles_labels()
    handles.append(patch)
    labels.append(name)

    legend._legend_box = None
    legend._init_legend_box(handles, labels)
    legend._set_loc(legend._loc)
    legend.set_title(legend.get_title().get_text())


def make_histogram(entries, labels, colors, xaxis_label="", title ="", num_bins = 10, logy = False, normalize = False, stacked = False, h_type = 'step', 
        h_range = None, fontsize = 24, fname="", yaxis_label = "", ymax = -1, mean_std = False, weights = None):
    alpha = 1.
    if(stacked): 
        h_type = 'barstacked'
        alpha = 0.2
    fig, ax = plt.subplots(figsize=fig_size)
    ns, bins, patches = plt.hist(entries, bins=num_bins, range=h_range, weights = weights, color=colors, alpha=alpha,label=labels, density = normalize, histtype=h_type)
    plt.xlabel(xaxis_label, fontsize =fontsize)
    plt.tick_params(axis='x', labelsize=fontsize)
    plt.tick_params(axis='y', labelsize=fontsize)


    if(logy): plt.yscale('log')
    elif(ymax > 0):
        plt.ylim([0,ymax])
    else:
        ymax = 1.3 * np.amax(ns)
        plt.ylim([0,ymax])

    if(yaxis_label != ""):
        plt.ylabel(yaxis_label, fontsize=fontsize)
        plt.tick_params(axis='y', labelsize=fontsize)
    plt.title(title, fontsize=fontsize)
    plt.legend(loc='upper right', fontsize = fontsize)

    if(mean_std):
        mean = np.mean(entries)
        std = np.std(entries)
        plt.text(0.6, 0.83, "Mean %.2f, Std. Dev. %.2f" % (mean, std), fontsize = 16, transform = ax.transAxes)

    if(fname != ""): 
        plt.savefig(fname, bbox_inches = 'tight')
        print("saving fig %s" %fname)
    #else: plt.show(block=False)
    return fig

def make_profile_hist(x,y, x_bins, xaxis_label="", yaxis_label="", fname = "", fontsize = 16, logy=False):

    x_bins_ids = np.digitize(x, bins = x_bins)
    bin_centers = 0.5 * (x_bins[:-1] + x_bins[1:])
    bin_width = x_bins[1] - x_bins[0]
    y_means = scipy.stats.binned_statistic(x, y, bins = x_bins, statistic = "mean").statistic
    y_sem = scipy.stats.binned_statistic(x, y, bins = x_bins, statistic=scipy.stats.sem).statistic

    fig = plt.figure(figsize=fig_size)
    plt.errorbar(x=bin_centers, y=y_means, xerr=bin_width, yerr=y_sem, linestyle='none')

    plt.xlabel(xaxis_label, fontsize =fontsize)
    plt.tick_params(axis='x', labelsize=fontsize)
    if(logy): plt.yscale('log')
    if(yaxis_label != ""):
        plt.ylabel(yaxis_label, fontsize=fontsize)
        plt.tick_params(axis='y', labelsize=fontsize)
    if(fname != ""): 
        print("saving fig %s" %fname)
        plt.savefig(fname)
    #else: plt.show(block=False)
    return fig



def make_outline_hist(stacks,outlines, labels, colors, xaxis_label, title, num_bins, logy = False,  normalize = False,  h_type = 'step', 
        h_range = None, fontsize = 16, fname="", yaxis_label = ""):
    alpha = 1.
    n_stacks = len(stacks)
    fig = plt.figure(figsize=fig_size)
    if(n_stacks > 0):
        plt.hist(stacks, bins=num_bins, range=h_range, color=colors[:n_stacks], alpha=0.2,label=labels[:n_stacks], density = normalize, histtype='barstacked')
    if(len(outlines) > 0):
        plt.hist(outlines, bins=num_bins, range=h_range, color=colors[n_stacks:], alpha=1.,label=labels[n_stacks:], density = normalize, histtype='step')
    plt.xlabel(xaxis_label, fontsize =fontsize)
    plt.tick_params(axis='x', labelsize=fontsize)
    if(logy): plt.yscale('log')
    if(yaxis_label != ""):
        plt.ylabel(yaxis_label, fontsize=fontsize)
        plt.tick_params(axis='y', labelsize=fontsize)
    plt.title(title, fontsize=fontsize)
    plt.legend(loc='upper right', fontsize = fontsize)
    if(fname != ""): 
        print("saving fig %s" %fname)
        plt.savefig(fname)
    #else: plt.show(block=False)
    return fig



#Mpl
def make_multi_ratio_histogram(entries, labels = None, colors = None, axis_label = None, title = None, num_bins=10, normalize = False, h_range = None, first_like_data = False,
        weights = None, fname="", ratio_range = -1, errors = False, logy = False, max_rw = 5, sys_weights = None):
    h_type= 'step'
    alpha = 1.
    fontsize = 22
    lw = 3
    fig = plt.figure(figsize = fig_size)
    gs = gridspec.GridSpec(2,1, height_ratios = [3,1])
    ax0 =  plt.subplot(gs[0])

    if(h_range is None):
        low = np.amin(entries[0])
        high = np.amax(entries[0])
    else:
        low,high = h_range

    if(first_like_data):
        data,bins = np.histogram(entries[0], bins = num_bins, range=(low,high), weights = weights[0])
        data_uncs = np.sqrt(data)

        if(normalize):
            norm = np.sum(data) * ((high - low)/num_bins) #counts * bin width
            data /= norm
            data_uncs /=norm

        ns, bins, patches  = ax0.hist(entries[1:], bins=num_bins, range=(low,high), color=colors[1:], alpha=alpha,label=labels[1:len(entries)], 
                density = normalize, weights = weights[1:], histtype=h_type, linewidth = lw)

        bincenters = 0.5*(bins[1:]+bins[:-1]) 
        ax0.errorbar(bincenters, data, yerr=data_uncs, fmt='o', ecolor = colors[0], color = colors[0], label = labels[0])

        ns.insert(0,data)

    else:
        ns, bins, patches  = ax0.hist(entries, bins=num_bins, range=(low,high), color=colors, alpha=alpha,label=labels[:len(entries)], 
                density = normalize, weights = weights, histtype=h_type, linewidth = lw)
    plt.xlim(h_range)


    plt.xlim([low, high])
    if(logy): plt.yscale("log")
    plt.title(title, fontsize=fontsize)

    bin_size = bins[1] - bins[0]
    bincenters = 0.5*(bins[1:]+bins[:-1]) 
    ax1 = plt.subplot(gs[1])

    ratios = []
    errs = []

    leg = ax0.legend(loc='best', fontsize = 14)
    frac_unc = None
    if(sys_weights is not None):
        for j in range(len(sys_weights)): # separate set for each observable
            sys_w = sys_weights[j]
            if(len(sys_w) == 0): continue
            nom = ns[j]

            uncs_up = np.zeros(nom.shape)
            uncs_down = np.zeros(nom.shape)
            for i in range(len(sys_w)):
                weights_up = sys_w[i][0]
                weights_down = sys_w[i][1]

                ns_sys_up, _ = np.histogram(entries[j], bins = bins, weights = weights_up, density =normalize )
                ns_sys_down, _ = np.histogram(entries[j], bins = bins, weights = weights_down, density =normalize )

                ns_up, ns_down = np.maximum(ns_sys_up, ns_sys_down), np.minimum(ns_sys_up, ns_sys_down)

                uncs_up += (ns_up - nom)**2
                uncs_down += (ns_down - nom)**2

            uncs_up = np.sqrt(uncs_up)
            uncs_down = np.sqrt(uncs_down)

            vals_up = nom + uncs_up
            vals_down = nom - uncs_down
            
            vals_up = np.append(vals_up, vals_up[-1])
            vals_down = np.append(vals_down, vals_down[-1])

            h_unc = ax0.fill_between(bins, vals_down, vals_up, color = colors[j], alpha = 0.5, step = 'post')

            #ratio panel
            ratios_unc_up = 1. + (uncs_up / nom)
            ratios_unc_down = 1. - (uncs_down / nom)

            #Add dummy value at end
            ratios_unc_up = np.append(ratios_unc_up, ratios_unc_up[-1])
            ratios_unc_down = np.append(ratios_unc_down, ratios_unc_down[-1])

            h_unc = ax1.fill_between(bins, ratios_unc_down, ratios_unc_up, color = colors[j], alpha = 0.5)

    for i in range(1, len(ns)):
        ratio = np.clip(ns[i], 1e-8, None) / np.clip(ns[0], 1e-8, None)
        ratios.append(ratio)

        ratio_err = None
        #if(errors):
        #    if(weights != None):
        #        w0 = weights[0]**2
        #        w1 = weights[i]**2
        #        norm0 = np.sum(weights[0])*bin_size
        #        norm1 = np.sum(weights[i])*bin_size
        #    else:
        #        w0 = w1 = None
        #        norm0 = entries[0].shape[0]*bin_size
        #        norm1 = entries[i].shape[0]*bin_size

        #    err0 = np.sqrt(np.histogram(entries[0], bins=bins, weights=w0)[0])/norm0
        #    err1 = np.sqrt(np.histogram(entries[i], bins=bins, weights=w1)[0])/norm1
        #    err0_alt  = np.sqrt(norm0*n0)/norm0
        #    err_alt1  = np.sqrt(norm1*n1)/norm1
        #    ratio_err = ratio * np.sqrt((err0/n0)**2 + (err1/n1)**2)
        #errs.append(ratio_err)

        ax1.errorbar(bincenters, ratio, yerr = ratio_err, alpha=alpha, markerfacecolor = colors[i], markeredgecolor = colors[i], fmt='ko')



    

    label_size = 18
    ax1.set_ylabel("Ratio", fontsize= label_size)
    ax1.set_xlabel(axis_label, fontsize = label_size)

    plt.xlim([low, high])

    if(type(ratio_range) == list or type(ratio_range) == tuple):
        plt.ylim(ratio_range[0], ratio_range[1])
    else:
        if(ratio_range > 0):
            plt.ylim([1-ratio_range, 1+ratio_range])

    plt.grid(axis='y')


    if(fname != ""): 
        plt.savefig(fname)
        print("saving fig %s" %fname)

    return ns, bins, ratios, frac_unc


def make_ratio_histogram(entries, labels, colors, axis_label, title, num_bins, normalize = False, h_range = None, 
        weights = None, fname="", ratio_range = -1, errors = False, extras = None, logy = False, max_rw = 5):


    h_type= 'step'
    alpha = 1.
    fontsize = 16
    fig = plt.figure(figsize = fig_size)
    gs = gridspec.GridSpec(2,1, height_ratios = [3,1])
    ax0 =  plt.subplot(gs[0])

    if(h_range is None):
        low = np.amin([np.amin(entries[0]), np.amin(entries[1])])
        high = np.amax([np.amax(entries[0]), np.amax(entries[1])])
        h_range = (low, high)
    else:
        low,high = h_range

    print(title, low,high)


    ns, bins, patches  = ax0.hist(entries, bins=num_bins, range=h_range, color=colors, alpha=alpha,label=labels, 
            density = normalize, weights = weights, histtype=h_type)

    if(extras is not None):
        ecolors = ['red', 'orange', 'cyan']
        for e_i, extra in enumerate(extras):
            ax0.hist(extra[0], bins= num_bins, range = h_range, color = ecolors[e_i], label = extra[2], density = normalize, weights = extra[1], histtype=h_type)


    plt.xlim([low, high])
    if(logy): plt.yscale("log")
    ax0.legend(loc='upper right')
    plt.title(title, fontsize=fontsize)
    n0 = np.clip(ns[0], 1e-8, None)
    n1 = np.clip(ns[1], 1e-8, None)
    ratio =  n0/ n1

    #low outliers more ok than high ones
    if(max_rw > 0):
        ratio = np.clip(ratio, 1./(2*max_rw), max_rw)

    ratio_err = None

    bin_size = bins[1] - bins[0]

    if(errors):
        if(weights != None):
            w0 = weights[0]**2
            w1 = weights[1]**2
            norm0 = np.sum(weights[0])*bin_size
            norm1 = np.sum(weights[1])*bin_size
        else:
            w0 = w1 = None
            norm0 = entries[0].shape[0]*bin_size
            norm1 = entries[1].shape[0]*bin_size

        err0 = np.sqrt(np.histogram(entries[0], bins=bins, weights=w0)[0])/norm0
        err1 = np.sqrt(np.histogram(entries[1], bins=bins, weights=w1)[0])/norm1
        err0_alt  = np.sqrt(norm0*n0)/norm0
        err_alt1  = np.sqrt(norm1*n1)/norm1


        ratio_err = ratio * np.sqrt((err0/n0)**2 + (err1/n1)**2)



    bincenters = 0.5*(bins[1:]+bins[:-1]) 
    ax1 = plt.subplot(gs[1])

    ax1.errorbar(bincenters, ratio, yerr = ratio_err, alpha=alpha, fmt='ko')

    plt.xlim([np.amin(entries[0]), np.amax(entries[0])])
    ax1.set_ylabel("Ratio")
    ax1.set_xlabel(axis_label)


    if(type(ratio_range) == list):
        plt.ylim(ratio_range[0], ratio_range[1])
    else:
        if(ratio_range > 0):
            plt.ylim([1-ratio_range, 1+ratio_range])

    plt.grid(axis='y')

    if(fname != ""):
        plt.savefig(fname)
        print("saving fig %s" %fname)

    return bins, ratio


def get_chi2(ratio):
    ys = ratio.GetY()
    xs = ratio.GetX()
    e_low = ratio.GetEYlow()
    e_high = ratio.GetEYhigh()
    ndof = len(ys)

    chi2 = 0.
    for i,y in enumerate(ys):
        e = e_low[i] if y > 1 else e_high[i]
        chi2 += ((1. - y)/ e)**2;
    return chi2

def makeCan(name, fname, histlist, bkglist=[],signals=[],totlist = [], colors=[],titles=[],dataName='Data',bkgNames=[],signalNames=[], drawSys = False, data_label = "data",
        draw_chi2 = False, stack = True, outlines = [], logy=False,rootfile=False,xtitle='',ytitle='',dataOff=False,datastyle='pe',year=-1, ratio_range = None, NDiv = 205, prelim = False):  

    sig_color = ROOT.kOrange - 7

    # histlist is just the generic list but if bkglist is specified (non-empty)
    # then this function will stack the backgrounds and compare against histlist as if 
    # it is data. The imporant bit is that bkglist is a list of lists. The first index
    # of bkglist corresponds to the index in histlist (the corresponding data). 
    # For example you could have:
    #   histlist = [data1, data2]
    #   bkglist = [[bkg1_1,bkg2_1],[bkg1_2,bkg2_2]]


    if len(histlist) == 1:
        width = 800
        height = 700
        padx = 1
        pady = 1
    elif len(histlist) == 2:
        width = 1200
        height = 700
        padx = 2
        pady = 1
    elif len(histlist) == 3:
        width = 1600
        height = 700
        padx = 3
        pady = 1
    elif len(histlist) == 4:
        width = 1200
        height = 1000
        padx = 2
        pady = 2
    elif len(histlist) == 6 or len(histlist) == 5:
        width = 1600
        height = 1000
        padx = 3
        pady = 2
    else:
        print('histlist of size ' + str(len(histlist)) + ' not currently supported')
        print(histlist)
        return 0

    tdrstyle.setTDRStyle()

    myCan = ROOT.TCanvas(name,name,width,height)
    myCan.Divide(padx,pady)

    # Just some colors that I think work well together and a bunch of empty lists for storage if needed
    default_colors = [ROOT.kRed,ROOT.kMagenta,ROOT.kGreen,ROOT.kCyan,ROOT.kBlue]
    if len(colors) == 0:   
        colors = default_colors
    stacks = []
    legends = []
    legends_list = []
    mains = []
    subs = []
    pulls = []
    logString = ''
    leg_align_right = True
    CMS_align_right = False

    # For each hist/data distribution
    for hist_index, hist in enumerate(histlist):
        # Grab the pad we want to draw in
        myCan.cd(hist_index+1)
        # if len(histlist) > 1:
        thisPad = myCan.GetPrimitive(name+'_'+str(hist_index+1))
        thisPad.cd()        

        # If this is a TH2, just draw the lego
        if hist.ClassName().find('TH2') != -1:
            print('ERROR: It seems you are trying to plot backgrounds with data on a 2D plot. This is not supported since there is no good way to view this type of distribution.')
        
        # Otherwise it's a TH1 hopefully
        else:
            titleSize = 0.09
            alpha = 1
            if dataOff:
                alpha = 0
            hist.SetLineColorAlpha(ROOT.kBlack,alpha)
            if 'pe' in datastyle.lower():
                hist.SetMarkerColorAlpha(ROOT.kBlack,alpha)
                hist.SetMarkerStyle(8)
            if 'hist' in datastyle.lower():
                hist.SetFillColorAlpha(0,0)

            hist.GetXaxis().SetTitle(xtitle)
            hist.GetYaxis().SetTitle(ytitle)
            
            # If there are no backgrounds, only plot the data (semilog if desired)
            if len(bkglist) == 0:
                hist.SetMaximum(1.13*hist.GetMaximum())
                if len(titles) > 0:
                    hist.SetTitle(titles[hist_index])
                    hist.SetTitleOffset(1.1)
                hist.Draw(datastyle)
            
            # Otherwise...
            else:
                # Create some subpads, a legend, a stack, and a total bkg hist that we'll use for the error bars
                if not dataOff:
                    mains.append(ROOT.TPad(hist.GetName()+'_main',hist.GetName()+'_main',0, 0.3, 1, 1))
                    subs.append(ROOT.TPad(hist.GetName()+'_sub',hist.GetName()+'_sub',0, 0, 1, 0.3))

                else:
                    mains.append(ROOT.TPad(hist.GetName()+'_main',hist.GetName()+'_main',0, 0.1, 1, 1))
                    subs.append(ROOT.TPad(hist.GetName()+'_sub',hist.GetName()+'_sub',0, 0, 0, 0))

                leg_align_right = True
                CMS_align_right = False
                x_max = totlist[hist_index].GetMaximumBin()
                nbins = totlist[hist_index].GetXaxis().GetNbins()
                #if(2 *x_max > nbins):
                #    print("Found max val in bin %i, aligning legend on the left" % x_max)
                #    leg_align_right = False
                #    CMS_align_right = True
                if not logy: 
                    y_end = 0.88
                    y_size = 0.2 + 0.025*(len(bkglist[0])+len(signals))
                    x_size = 0.53
                    if(leg_align_right):
                        x_start = 0.4
                    else:
                        x_start = 0.2

                    legends.append(ROOT.TLegend(x_start,y_end - y_size,x_start + x_size,y_end))
                else: 
                    legends.append(ROOT.TLegend(0.2,0.11,0.45,0.2+0.02*(len(bkglist[0])+len(signals))))

                legends_list.append([])


                # Set margins and make these two pads primitives of the division, thisPad
                mains[hist_index].SetBottomMargin(0.04)
                mains[hist_index].SetLeftMargin(0.2)
                mains[hist_index].SetRightMargin(0.05)
                mains[hist_index].SetTopMargin(0.08)

                subs[hist_index].SetLeftMargin(0.2)
                subs[hist_index].SetRightMargin(0.05)
                subs[hist_index].SetTopMargin(0.01)
                subs[hist_index].SetBottomMargin(0.5)
                mains[hist_index].SetFillColorAlpha(0,1)
                subs[hist_index].SetFillColorAlpha(0,1)
                mains[hist_index].SetLineColorAlpha(0,1)
                subs[hist_index].SetLineColorAlpha(0,1)
                mains[hist_index].Draw()
                subs[hist_index].Draw()

                # Build the stack
                if(stack):
                    stacks.append(ROOT.THStack(hist.GetName()+'_stack',hist.GetName()+'_stack'))
                    for bkg_index,bkg in enumerate(bkglist[hist_index]):     # Won't loop if bkglist is empty
                        # bkg.Sumw2()
                        bkg.SetLineColor(ROOT.kBlack)
                        if logy:
                            bkg.SetMinimum(1e-3)

                        if colors[bkg_index] != None:
                            bkg.SetFillColor(colors[bkg_index])
                            bkg.SetLineColor(colors[bkg_index])
                        else:
                            bkg.SetFillColor(default_colors[bkg_index])
                            bkg.Print()

                        stacks[hist_index].Add(bkg)
                        if bkgNames == []: this_bkg_name = bkg.GetName().split('_')[0]
                        elif type(bkgNames[0]) != list: this_bkg_name = bkgNames[bkg_index]
                        else: this_bkg_name = bkgNames[hist_index][bkg_index]
                        legends_list[hist_index].append((bkg,this_bkg_name,'f'))

                    histList = [stacks[hist_index],totlist[hist_index],hist]
                else:
                    histList = copy.copy(bkglist[hist_index])
                    histList += [totlist[hist_index], hist]

                    
                # Go to main pad, set logy if needed
                mains[hist_index].cd()


                # Set y max of all hists to be the same to accomodate the tallest
                max_scaling = 2.5

                yMax = histList[0].GetMaximum()
                maxHist = histList[0]
                for h in range(1,len(histList)):
                    if histList[h].GetMaximum() > yMax:
                        yMax = histList[h].GetMaximum()
                        maxHist = histList[h]
                for h in histList:
                    h.SetMaximum(yMax*max_scaling)
                    if logy == True:
                        h.SetMaximum(yMax*10)
                    else:
                        h.SetMinimum(0.)

                
                mLS = 0.07
                mTS = 0.1
                TOffset = 0.95
                # Now draw the main pad
                data_leg_title = hist.GetTitle()
                if len(titles) > 0:
                    hist.SetTitle(titles[hist_index])
                hist.GetYaxis().SetTitleOffset(TOffset)
                hist.GetXaxis().SetTitleOffset(TOffset)
                if(len(ytitle) == 0): ytitle = "Events / bin"
                hist.GetYaxis().SetTitle(ytitle)
                hist.GetYaxis().SetLabelSize(mLS)
                hist.GetYaxis().SetTitleSize(mTS)
                hist.GetYaxis().SetNdivisions(505)
                hist.GetXaxis().SetLabelOffset(999)
                hist.SetLineWidth(2)


                if logy == True:
                    hist.SetMinimum(1e-3)

                hist.SetBinErrorOption(ROOT.TH1.kPoisson)
                hist.Draw(datastyle)
                #print("Drawing %s %s \n" hist.GetName(), datastyle)

                if(stack): stacks[hist_index].Draw('same hist')
                else:
                    for i,h in enumerate(bkglist[hist_index]): 
                        h.SetLineWidth(2)
                        h.Draw('same hist')
                        legends_list[hist_index].append((h,bkgNames[i], 'L') )
                #print("Drawing %s same hist \n" stacks[hist_index].GetName())

                # Do the signals
                if len(signals) > 0: 
                    signals[hist_index].SetLineColor(sig_color)
                    signals[hist_index].SetLineWidth(4)
                    if logy == True:
                        signals[hist_index].SetMinimum(1e-3)
                    if signalNames == []: this_sig_name = signals[hist_index].GetName().split('_')[0]
                    legends_list[hist_index].append((signals[hist_index],this_sig_name,'L'))
                    signals[hist_index].Draw('hist same')

                if(stack):
                    ROOT.gStyle.SetHatchesLineWidth(3)
                    totlist[hist_index].SetLineColor(ROOT.kWhite)
                    totlist[hist_index].SetFillColor(ROOT.kBlack)
                    totlist[hist_index].SetFillStyle(3354)
                    totlist[hist_index].SetMarkerStyle(20)
                    totlist[hist_index].SetMarkerSize(0.01)

                    #if(drawSys):
                        #totlist[hist_index].Draw('e2 same')
                else:
                    totlist[hist_index].SetLineWidth(4)
                    totlist[hist_index].Draw('hist same')
                    legends_list[hist_index].append((totlist[hist_index], "Total", 'L') )


                if not dataOff:
                    legends_list[hist_index].append((hist,dataName,datastyle))
                    hist.Draw(datastyle+' same')


                #Draw helpful lines

                #legends[hist_index].SetHeader(titles[0], "c")
                legends[hist_index].SetNColumns(2)
                legends[hist_index].SetTextSize(0.04)

                for entry in legends_list[hist_index][::-1]:
                    legends[hist_index].AddEntry(entry[0], entry[1], entry[2])


                #if(drawSys):
                    #legends[hist_index].AddEntry(totlist[hist_index], "Sys. unc.", "f")


                ratio, ratio_sys_unc = makeRatio(hist,totlist[hist_index])
                ratio_signal = False
                if(len(signals) > 0):
                    ratio_signal = True
                    ratio_sig, _ = makeRatio(hist, signals[hist_index])
                    ratio_sig.SetLineColor(sig_color)
                    ratio_sig.SetMarkerColor(sig_color)

                chi2 = chi2_sig = 0.

                chi2 = get_chi2(ratio)
                if(len(signals) > 0): chi2_sig = get_chi2(ratio_sig)

                #print("Chi2/nbin for chan %s is %.1f/%i" % (titles[hist_index], chi2, pull.GetNbinsX()))

                #if(drawSys):
                    #legends[hist_index].AddEntry(ratio_sys_unc, "Sys. unc.", "f")




                legends[hist_index].SetBorderSize(0)
                legends[hist_index].Draw()
                ROOT.gPad.RedrawAxis()

                # Draw the pull
                subs[hist_index].cd()
                # Build the pull

                LS = mLS * 0.7/0.3
                #title size given as fraction of pad width, scale up to have same size as main pad
                YTS =  0.8 * mTS * 0.7/0.3
                XTS =  mTS * 0.7/0.3
                lTOffset = TOffset * 0.3 / 0.7


                ratio.SetMarkerStyle(8)


                r_axis_hist = ratio
                if(drawSys):
                    r_axis_hist = ratio_sys_unc

                r_axis_hist.GetYaxis().SetRangeUser(ratio_range[0], ratio_range[1])
                r_axis_hist.GetYaxis().SetTitleOffset(lTOffset)
                r_axis_hist.GetYaxis().SetTickLength(0.04)
                             
                r_axis_hist.GetYaxis().SetLabelSize(LS)
                r_axis_hist.GetYaxis().SetTitleSize(YTS)
                r_axis_hist.GetYaxis().SetNdivisions(NDiv)
                r_axis_hist.GetYaxis().SetTitle("Data / Sim.")

                r_axis_hist.GetXaxis().SetRangeUser(hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
                r_axis_hist.GetXaxis().SetTitleOffset(1.)
                r_axis_hist.GetXaxis().SetLabelOffset(0.05)
                r_axis_hist.GetXaxis().SetLabelSize(LS)
                r_axis_hist.GetXaxis().SetTitleSize(XTS)
                r_axis_hist.GetXaxis().SetTitle(xtitle)
                r_axis_hist.GetXaxis().SetTickLength(0.06)

                r_axis_hist.SetTitle("")

                #ratio_sys_unc.SetFillColor(ROOT.kBlack)
                #ratio_sys_unc.SetFillStyle(3015)
                ratio_sys_unc.SetLineColor(ROOT.kGray)
                ratio_sys_unc.SetFillColor(ROOT.kGray)
                #ratio_sys_unc.SetFillStyle(3015)



                ratio.SetLineStyle(1)
                ratio.SetLineWidth(2)
                if(drawSys):
                    ratio_sys_unc.Draw("A3 same")
                    ratio.Draw('p0e0Z same')
                else:
                    ratio.Draw('Ap0e0Z same')

                if(ratio_signal):
                    ratio_sig.Draw("p0e0Z same")

                
                line = ROOT.TLine(hist.GetXaxis().GetXmin(), 1.0, hist.GetXaxis().GetXmax(), 1.0)
                line.SetLineStyle(9)
                line.Draw()

                if(draw_chi2):
                    subs[hist_index].cd()
                    xs = ratio.GetX()

                    xscale = (np.amax(xs) - np.amin(xs))
                    x_start = np.amin(xs) + 0.22 * xscale
                    x_stop = x_start + 0.25* xscale

                    ndof = len(xs)


                    pave = ROOT.TPaveText(x_start, 0.75 * ratio_range[1], x_stop, 0.98*ratio_range[1])
                    pave.AddText("#chi^{2} / ndof = %.1f / %i" % (chi2, ndof))
                    pave.SetFillColor(ROOT.kWhite)
                    pave.SetFillColor(ROOT.kWhite)
                    #pave.SetBorderSize(1)
                    pave.SetBorderSize(0)
                    pave.Draw()

                    if(ratio_signal):
                        x_start = np.amin(xs) + 0.53 * xscale
                        x_stop = x_start + 0.25 * xscale

                        pave2 = ROOT.TPaveText(x_start, 0.75 * ratio_range[1], x_stop, 0.98*ratio_range[1])
                        pave2.AddText("#chi^{2} / ndof = %.1f / %i" % (chi2_sig, ndof))
                        pave2.SetTextColor(sig_color)
                        pave2.SetFillColor(ROOT.kWhite)
                        #pav2e.SetBorderSize(1)
                        pave2.SetBorderSize(0)
                        pave2.Draw()


                if logy == True:
                    mains[hist_index].SetLogy()

                if(prelim): 
                    print("Prelim")
                    CMS_lumi.writeExtraText = True
                else: CMS_lumi.writeExtraText = False
                if(CMS_align_right): CMS_loc = 33
                else: CMS_loc = 11
                CMS_lumi.CMS_lumi(mains[hist_index], year, CMS_loc)

    print("Creating " + fname)
    myCan.Print(fname)





        
def makeRatio( DATA,BKG):

    nbins = DATA.GetNbinsX()
    x = array('d')
    ratio = array('d')
    y_err_up = array('d')
    y_err_down = array('d')
    sys_err_up = array('d')
    sys_err_down = array('d')
    x_err_low = array('d')
    x_err_high = array('d')
    x_err_zero = array('d', [0.]* nbins)
    y_val1 = array('d', [1.]* nbins)
    
    for ibin in range(1,DATA.GetNbinsX()+1):
        DATAcont = DATA.GetBinContent(ibin)
        BKGcont = max(1e-6, BKG.GetBinContent(ibin))
        BKG_errup = BKG.GetBinErrorUp(ibin)
        BKG_errdown = BKG.GetBinErrorLow(ibin)
        DATA_errup = DATA.GetBinErrorUp(ibin)
        DATA_errdown = DATA.GetBinErrorLow(ibin)


        x.append(DATA.GetXaxis().GetBinCenter(ibin))
        x_err_low.append(DATA.GetXaxis().GetBinCenter(ibin) - DATA.GetXaxis().GetBinLowEdge(ibin))
        x_err_high.append(-DATA.GetXaxis().GetBinCenter(ibin) + DATA.GetXaxis().GetBinUpEdge(ibin))

    

        ratio.append(DATAcont/BKGcont)
        y_err_up.append(DATA_errup / BKGcont)
        y_err_down.append(DATA_errdown / BKGcont)
        sys_err_up.append(BKG_errup/BKGcont)
        sys_err_down.append(BKG_errdown/BKGcont)
        
    pull = ROOT.TGraphAsymmErrors(nbins,x,ratio, x_err_zero, x_err_zero, y_err_down, y_err_up)
    #add extra at high-x for ratio plot
    x.append(x[-1] + x_err_high[-1])
    y_val1.append(1)
    x_err_low.append(x_err_low[-1])
    x_err_high.append(0)
    sys_err_up.append(sys_err_up[-1])
    sys_err_down.append(sys_err_down[-1])
    #add extra at low-x for ratio plot
    x.insert(0, x[0] - x_err_low[0])
    y_val1.insert(0,1)
    x_err_low.insert(0,0)
    x_err_high.insert(0, x_err_high[0])
    sys_err_up.insert(0,sys_err_up[0])
    sys_err_down.insert(0,sys_err_down[0])



    sys_unc = ROOT.TGraphAsymmErrors(nbins+2, x, y_val1, x_err_low, x_err_high, sys_err_down,  sys_err_up)
    return pull, sys_unc

def Make_up_down(hist):
    hist_up = hist.Clone(hist.GetName()+'_up')
    hist_down = hist.Clone(hist.GetName()+'_down')

    for xbin in range(1,hist.GetNbinsX()+1):
        errup = hist.GetBinErrorUp(xbin)
        errdown = hist.GetBinErrorLow(xbin)
        nom = hist.GetBinContent(xbin)

        hist_up.SetBinContent(xbin,nom+errup)
        hist_down.SetBinContent(xbin,nom-errdown)

    return hist_up,hist_down


def get_eff_unc(data = None, weight = None, name = "h", num_bins = 1000, unc = None, cut_val = 0.):
    #assumes selecting vals greater than cut_val
    if(unc is None or weight is None): 
        print("Missing args")
        return None

    bin_low = np.amin(data) - 1e-8
    bin_high = np.amax(data) + 1e-8

    h = ROOT.TH1F(name, "", num_bins, bin_low, bin_high)
    h_up = ROOT.TH1F(name+"_up", "", num_bins, bin_low, bin_high)
    h_down = ROOT.TH1F(name+"_down", "", num_bins, bin_low, bin_high)

    weights_up = np.clip(weight + unc, 0, 9999)
    weights_down = np.clip(weight - unc, 0, 9999)

    for i,e in enumerate(data):
        h.Fill(e, weight[i])
        h_up.Fill(e, weights_up[i])
        h_down.Fill(e, weights_down[i])


    h.Print()
    cut_bin = h.GetXaxis().FindBin(cut_val)
    tot_nom = h.Integral()
    tot_up = h_up.Integral()
    tot_down = h_down.Integral()
    cut_nom = cut_up = cut_down = 0.
    #add difference as uncertainty
    for ibin in range(cut_bin,num_bins +1):
        cut_nom += h.GetBinContent(ibin)
        cut_up += h_up.GetBinContent(ibin)
        cut_down += h_down.GetBinContent(ibin)

    eff_nom = cut_nom / tot_nom
    #eff_unc_up = (cut_up - cut_nom) / tot_nom 
    #eff_unc_down = (cut_down - cut_nom) / tot_nom
    eff_unc_up = eff_nom - cut_up/ tot_up
    eff_unc_down = eff_nom - cut_down / tot_down
    return eff_nom, eff_unc_up, eff_unc_down

def fill_hist(h, data, weight):

    for i,e in enumerate(data):
        if(weight is not None): w = weight[i]
        else: w = 1.

        if(hasattr(e, "__len__")): 
            for en in e: h.Fill(en, w)
        else: h.Fill(e, w)


def make_root_hist(data = None, weight = None, name = "h", bins = None, num_bins = 1, bin_low = 0, bin_high = 1, unc = None):

    if(bins is None):
        bins = array('f', np.linspace(bin_low, bin_high, num_bins+1))

    h = ROOT.TH1F(name, "", num_bins, bins)
    h.Sumw2()

    fill_hist(h, data, weight)


    if(unc is not None and type(unc) == float):
        #single fractional unc
        for ibin in range(1,num_bins +1):

            cont = h.GetBinContent(ibin)
            h.SetBinError(ibin, unc * cont)


    elif(unc is not None and np.sum(unc) > 1e-4):
        #unc weights
        weights_up = np.clip(weight + unc, 0, 9999)
        weights_down = np.clip(weight - unc, 0, 9999)

        h_up = ROOT.TH1F(name+"_up", "", num_bins, bin_low, bin_high)
        h_down = ROOT.TH1F(name+"_down", "", num_bins, bin_low, bin_high)
        for i,e in enumerate(data):
            h_up.Fill(e, weights_up[i])
            h_down.Fill(e, weights_down[i])

        #add difference as uncertainty
        for ibin in range(1,num_bins +1):

            sys_err = abs(h_up.GetBinContent(ibin) - h_down.GetBinContent(ibin)) / 2.0
            old_err = h.GetBinError(ibin)
            new_err = (sys_err**2 + old_err**2)**(0.5)
            h.SetBinError(ibin, new_err)




    return h




#ROOT version
def make_multi_sum_ratio_histogram(data = None, entries = None, labels = None, uncs = None, colors = None, axis_label = None, title = None, num_bins = None, drawSys = False, stack = True,
    draw_chi2 = False, normalize = False, h_range = None, weights = None, fname="", ratio_range = -1, errors = False, extras = None, logy = False, max_rw = 5, data_label = "data", year = -1):
    if(h_range is None):
        low = np.amin(entries[0])
        high = np.amax(entries[0])
        h_range = (low, high)
    else:
        low,high = h_range

    hists = []
    for i,e in enumerate(entries):
        unc = weight = None
        if(weights is not None): weight = weights[i]
        if(uncs is not None): unc = uncs[i]
        hist = make_root_hist(data = e, weight = weight, unc = unc, name = "h" + labels[i], num_bins = num_bins, bin_low = low, bin_high = high)
        if(stack): hist.SetFillColor(colors[i])
        hist.SetLineColor(colors[i])
        hists.append(hist)
    
    h_data = make_root_hist(data = data, weight = None, unc = None, name = data_label, num_bins = num_bins, bin_low = low, bin_high = high)
    h_tot = hists[0].Clone("h_tot")
    h_tot.Reset()
    for h in hists: h_tot.Add(h)
    if(not stack): 
        h_tot.SetLineColor(ROOT.kBlue)


    fname_root = fname.replace(".png", "").replace(".pdf", "") + ".root"
    print("Saving ROOT: " + fname_root)
    f = ROOT.TFile.Open(fname_root, "RECREATE")
    for h in hists:
        h.Write()
        h_tot.Write()
        h_data.Write()
    f.Close()
        

    return makeCan("temp", fname, [h_data], bkglist = [hists], totlist = [h_tot], colors = colors, bkgNames = labels, titles = [title], logy = logy, xtitle = axis_label,
        drawSys = drawSys, ratio_range = ratio_range, stack = stack, draw_chi2 = draw_chi2, prelim = True, year = year)



def draw_jet_image(image, title, fname = "", do_log = False):
    fontsize = 20
    image = np.clip(np.squeeze(image).astype('float'), 1e-8, None)
    if(do_log): image = np.log(image)
    fig = plt.figure(figsize=fig_size)
    plt.imshow(image, cmap = 'Blues', interpolation = 'nearest')
    plt.title(title, fontsize = fontsize)
    if(fname != ""):
        plt.savefig(fname)

def make_scatter_plot(x, y, color, axis_names, fname= ""  ):

    fig, ax = plt.subplots()
    alpha = 0.5
    size = 0.4
    ax.scatter(x,y, alpha = alpha, c = color, s=size)

    correlation = np.corrcoef(x,y)[0,1]
    text_str = r'$\rho_{x,y} $ = %.3f' % correlation
    plt.annotate(text_str, xy = (0.05, 0.95), xycoords = 'axes fraction', fontsize=14)

    ax.set_xlabel(axis_names[0], fontsize=14)
    ax.set_ylabel(axis_names[1], fontsize=14)
    plt.tick_params(axis='y', labelsize=12)
    plt.tick_params(axis='x', labelsize=12)
    if(fname != ""):
        print("saving %s" % fname)
        plt.savefig(fname)

def horizontal_bar_chart(vals, labels, fname = "", xaxis_label = ""):
    fig, ax = plt.subplots()

    order = np.flip(vals.argsort())
    labels = np.array(labels)
    y_pos = np.arange(len(vals))
    ax.barh(y_pos, vals[order], align = 'center')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels[order])
    ax.invert_yaxis()
    ax.set_xlabel(xaxis_label)

    if(fname != ""):
        print("saving %s" % fname)
        plt.savefig(fname)

def compute_chi2(data, data_unc, y):
    chi2 = 0.
    eps = 1e-8
    for i in range(len(data)):
        if(data_unc[i] <= 0.): data_unc[i] = max(y[i], data[i]) + eps
        chi2 += (data[i] - y[i])**2/data_unc[i]**2
    return chi2

def make_herwig_ratio_histogram(entries = None, labels = None, colors = None, axis_label = None, title = None, num_bins=10, normalize = False, h_range = None, first_like_data = True,
        weights = None, fname="", ratio_range = -1, errors = False, logy = False, max_rw = 5, sys_weights = None, stat_weights = None, draw_chi2 = True):
    h_type= 'step'
    alpha = 1.
    fontsize = 22
    label_size = 18
    lw = 3
    fig = plt.figure(figsize = fig_size)
    gs = gridspec.GridSpec(2,1, height_ratios = [3,1])
    ax0 =  plt.subplot(gs[0])

    if(h_range is None):
        low = np.amin(entries[0])
        high = np.amax(entries[0])
    else:
        low,high = h_range

    data,bins = np.histogram(entries[0], bins = num_bins, range=(low,high), weights = weights[0])
    data_uncs = np.sqrt(data)

    if(normalize):
        norm = np.sum(data) * ((high - low)/num_bins) #counts * bin width
        data /= norm
        data_uncs /=norm

    ns, bins, patches  = ax0.hist(entries[1:], bins=num_bins, range=(low,high), color=colors[1:], alpha=alpha,label=labels[1:len(entries)], 
            density = normalize, weights = weights[1:], histtype=h_type, linewidth = lw)

    bincenters = 0.5*(bins[1:]+bins[:-1]) 
    ax0.errorbar(bincenters, data, yerr=data_uncs, fmt='ko', markerfacecolor = colors[0], ecolor = colors[0], markeredgecolor = colors[0], label = labels[0])

    plt.xlim(h_range)

    if(logy): plt.yscale("log")
    plt.title(title, fontsize=fontsize)

    bin_size = bins[1] - bins[0]
    bincenters = 0.5*(bins[1:]+bins[:-1]) 

    ax1 = plt.subplot(gs[1])
    ax1.errorbar(bincenters, np.ones_like(bincenters), yerr = data_uncs/data, markerfacecolor = colors[0], ecolor = colors[0], markeredgecolor = colors[0], fmt='ko')

    ratios = []
    chi2s = []

    for i in range(len(ns)):
        ratio =  np.clip(data, 1e-8, None)/ np.clip(ns[i], 1e-8, None)
        ratio = np.append(ratio, ratio[-1])
        ratios.append(ratio)

        plt.step(bins, ratio, color = colors[i+1], linewidth = lw, where = 'post')

        chi2s.append(compute_chi2(data, data_uncs, ns[i]))

    print(chi2s)

    leg = ax0.legend(loc='best', fontsize = 14)

    #draw sys unc band
    if(sys_weights is not None or stat_weights is not None):
        for j in range(len(weights)-1): # separate set for each observable
            nom = ns[j]
            uncs_up = np.zeros(nom.shape)
            uncs_down = np.zeros(nom.shape)

            if(sys_weights is not None):
                sys_w = sys_weights[j]
                if(len(sys_w) == 0): continue

                for i in range(len(sys_w)):
                    weights_up = sys_w[i][0]
                    weights_down = sys_w[i][1]

                    ns_sys_up, _ = np.histogram(entries[j+1], bins = bins, weights = weights_up, density =normalize )
                    ns_sys_down, _ = np.histogram(entries[j+1], bins = bins, weights = weights_down, density =normalize )

                    ns_up, ns_down = np.maximum(ns_sys_up, ns_sys_down), np.minimum(ns_sys_up, ns_sys_down)

                    uncs_up += (ns_up - nom)**2
                    uncs_down += (ns_down - nom)**2

            if(stat_weights is not None):
                stat_w = stat_weights[j]
                if(len(stat_w) == 0): continue
                for i in range(len(stat_w)):
                    ns_var = []
                    for k in range(stat_w[i].shape[1]):
                        ns_, _ = np.histogram(entries[j+1], bins = bins, weights = stat_w[i][:,k], density =normalize )
                        ns_var.append(ns_)

                    ns_var = np.array(ns_var)
                    std_devs = np.std(ns_var, axis=0)
                    means = np.mean(ns_var, axis=0)

                    uncs_up +=  (nom - means)**2
                    uncs_down +=  (nom - means)**2

                    uncs_up +=  std_devs**2
                    uncs_down +=  std_devs**2

            uncs_up = np.sqrt(uncs_up)
            uncs_down = np.sqrt(uncs_down)

            vals_up = nom + uncs_up
            vals_down = nom - uncs_down

            #compute ratio unc
            ratio_up =  np.clip(data, 1e-8, None)/ np.clip(vals_up, 1e-8, None) 
            ratio_down =  np.clip(data, 1e-8, None)/ np.clip(vals_down, 1e-8, None) 
            
            #Pad with dummy val
            vals_up = np.append(vals_up, vals_up[-1])
            vals_down = np.append(vals_down, vals_down[-1])

            ratio_up = np.append(ratio_up, ratio_up[-1])
            ratio_down = np.append(ratio_down, ratio_down[-1])

            ax0.fill_between(bins, vals_down, vals_up, color = colors[j+1], alpha = 0.5, step = 'post')
            ax1.fill_between(bins, ratio_down, ratio_up, color = colors[j+1], alpha = 0.5, step = 'post')





    ax1.set_ylabel("Herwig/Pythia", fontsize= label_size)
    ax1.set_xlabel(axis_label, fontsize = label_size)

    plt.xlim([low, high])

    if(type(ratio_range) == list or type(ratio_range) == tuple):
        plt.ylim(ratio_range[0], ratio_range[1])
    else:
        if(ratio_range > 0):
            plt.ylim([1-ratio_range, 1+ratio_range])

    plt.grid(axis='y')

    if(draw_chi2):
        plt.sca(ax0)
        ndof = len(bins)-1
        y_val = ax0.get_ylim()[1] * 0.7
        x_val = ax0.get_xlim()[1] * 0.1
        for j,chi2 in enumerate(chi2s):
            print(chi2, ndof)
            txt = r"$\chi^2$ / ndof = %.1f / %i" % (chi2, ndof)
            plt.text(bins[2], y_val, txt, color = colors[j+1], horizontalalignment = 'left', fontweight = 'bold')
            y_val -= y_val * 0.05



    if(fname != ""): 
        plt.savefig(fname)
        print("saving fig %s" %fname)

    return 

