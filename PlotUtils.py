import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import gridspec
from sklearn.metrics import roc_curve,auc
import scipy.stats
import numpy as np
import ROOT
import tdrstyle
from array import array

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





def plot_training(hist, fname =""):
    #plot trianing and validation loss

    loss = hist['loss']

    epochs = range(1, len(loss) + 1)
    colors = ['b', 'g', 'grey', 'r']
    idx=0

    plt.figure(figsize=fig_size)
    for label,loss_hist in hist.items():
        if(len(loss_hist) > 0): plt.plot(epochs, loss_hist, colors[idx], label=label)
        idx +=1
    plt.title('Training and validation loss')
    plt.yscale("log")
    plt.xlabel('Steps')
    plt.ylabel('Loss')
    plt.legend()
    if(fname != ""): 
        plt.savefig(fname)
        print("saving fig %s" %fname)
    #else: 
        #plt.show(block=False)


def make_roc_curve(classifiers, y_true, colors = None, logy=True, labels = None, fname=""):
    plt.figure(figsize=fig_size)

    y_true = np.clip(y_true, 0,1)

    fs = 18
    fs_leg = 16
    sig_effs = []
    bkg_effs = []
    aucs = []
    for idx,scores in enumerate(classifiers):

        fpr, tpr, thresholds = roc_curve(y_true, scores)
        roc_auc = auc(fpr, tpr)
        sig_effs.append(tpr)
        bkg_effs.append(fpr)

    make_roc_plot(sig_effs, bkg_effs, colors = colors, labels = labels, logy = logy, fname = fname)
        
def make_roc_plot(sig_effs, bkg_effs, colors = None, logy = True, labels = None, fname = ""):

    plt.figure(figsize=fig_size)
    fs = 18
    fs_leg = 16
    sic_max = 10.

    for idx in range(len(sig_effs)):
        tpr = sig_effs[idx]
        fpr = bkg_effs[idx]
        auc_ = auc(fpr, tpr)

        fpr = np.clip(fpr, 1e-8, 1.)
        #guard against division by 0
        ys = 1./fpr

        lbl = 'auc %.3f' % auc_
        clr = 'navy'
        if(labels is not None): lbl = labels[idx] + " = %.3f" % auc_
        if(colors is not None): clr = colors[idx]

        print(lbl)

        plt.plot(tpr, ys, lw=2, color=clr, label=lbl)



    plt.xlim([0, 1.0])
    plt.xlabel('Signal Efficiency', fontsize=fs)
    if(logy): 
        plt.yscale('log')
    plt.ylim([1., 1e4])
    plt.ylabel('QCD Rejection Rate', fontsize=fs)

    plt.legend(loc="upper right", fontsize = fs_leg)
    if(fname != ""): 
        print("Saving roc plot to %s" % fname)
        plt.savefig(fname)
    #else: 
        #plt.show(block=False)


def make_sic_plot(sig_effs, bkg_effs, colors = None, labels = None, eff_min = 1e-3, ymax = -1, fname=""):


    plt.figure(figsize=fig_size)

    fs = 18
    fs_leg = 16
    sic_max = 10.



    for idx in range(len(sig_effs)):
        sic = sig_effs[idx] / np.sqrt(bkg_effs[idx])

        lbl = 'auc'
        clr = 'navy'
        if(labels != None): lbl = labels[idx]
        if(colors != None): clr = colors[idx]
        print(lbl, "max sic: ", np.amax(sic))
        plt.plot(bkg_effs[idx], sic, lw=2, color=clr, label='%s' % lbl)


    
    plt.xlim([eff_min, 1.0])
    if(ymax < 0):
        ymax = sic_max
        
    plt.ylim([0,ymax])
    plt.xscale('log')
    plt.xlabel('Background Efficiency', fontsize = fs)
    plt.ylabel('Significance Improvement', fontsize = fs)
    plt.tick_params(axis='x', labelsize=fs_leg)
    plt.tick_params(axis='y', labelsize=fs_leg)
    plt.grid(axis = 'y', linestyle='--', linewidth = 0.5)
    plt.legend(loc="best", fontsize= fs_leg)
    if(fname != ""):
        plt.savefig(fname)
        print("Saving file to %s " % fname)


def make_sic_curve(classifiers, y_true, colors = None, logy=False, labels = None, eff_min = 1e-3, ymax = -1, fname=""):

    y_true = np.clip(y_true, 0,1)
    plt.figure(figsize=fig_size)

    fs = 18
    fs_leg = 16

    sic_max = 10.



    for idx,scores in enumerate(classifiers):
        fpr, tpr, thresholds = roc_curve(y_true, scores)
        fpr= np.clip(fpr, 1e-8, 1.)
        tpr = np.clip(tpr, 1e-8, 1.)

        mask = tpr > eff_min

        xs = fpr[mask]
        ys = tpr[mask]/np.sqrt(xs)

        sic_max = max(np.amax(ys), sic_max)


        


        lbl = 'auc'
        clr = 'navy'
        if(labels != None): lbl = labels[idx]
        if(colors != None): clr = colors[idx]
        print(lbl, "max sic: ", np.amax(ys))
        plt.plot(xs, ys, lw=2, color=clr, label='%s' % lbl)


    
    plt.xlim([eff_min, 1.0])
    if(ymax < 0):
        ymax = sic_max
        
    plt.ylim([0,ymax])
    plt.xscale('log')
    plt.xlabel('Background Efficiency', fontsize = fs)
    plt.ylabel('Significance Improvement', fontsize = fs)
    plt.tick_params(axis='x', labelsize=fs_leg)
    plt.tick_params(axis='y', labelsize=fs_leg)
    plt.grid(axis = 'y', linestyle='--', linewidth = 0.5)
    plt.legend(loc="best", fontsize= fs_leg)
    if(fname != ""):
        plt.savefig(fname)
        print("Saving file to %s " % fname)

def make_histogram(entries, labels, colors, xaxis_label="", title ="", num_bins = 10, logy = False, normalize = False, stacked = False, h_type = 'step', 
        h_range = None, fontsize = 16, fname="", yaxis_label = "", ymax = -1):
    alpha = 1.
    if(stacked): 
        h_type = 'barstacked'
        alpha = 0.2
    fig = plt.figure(figsize=fig_size)
    ns, bins, patches = plt.hist(entries, bins=num_bins, range=h_range, color=colors, alpha=alpha,label=labels, density = normalize, histtype=h_type)
    plt.xlabel(xaxis_label, fontsize =fontsize)
    plt.tick_params(axis='x', labelsize=fontsize)

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
    if(fname != ""): 
        plt.savefig(fname)
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



def make_outline_hist(stacks,outlines, labels, colors, xaxis_label, title, num_bins, logy = False,  normalize = False, save=False, h_type = 'step', 
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
    if(save): 
        print("saving fig %s" %fname)
        plt.savefig(fname)
    #else: plt.show(block=False)
    return fig



def make_multi_ratio_histogram(entries, labels, colors, axis_label, title, num_bins, normalize = False, save=False, h_range = None, 
        weights = None, fname="", ratio_range = -1, errors = False, logy = False, max_rw = 5, unc_band_norm = -1):
    h_type= 'step'
    alpha = 1.
    fontsize = 22
    fig = plt.figure(figsize = fig_size)
    gs = gridspec.GridSpec(2,1, height_ratios = [3,1])
    ax0 =  plt.subplot(gs[0])

    low = np.amin(entries[0])
    high = np.amax(entries[0])
    if(h_range == None):
        h_range = (low, high)

    ns, bins, patches  = ax0.hist(entries, bins=num_bins, range=h_range, color=colors, alpha=alpha,label=labels[:len(entries)], 
            density = normalize, weights = weights, histtype=h_type)


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
    if(unc_band_norm > 0):
        raw_dist = np.copy(ns[0])
        if(not normalize): norm_dist = raw_dist * unc_band_norm / np.sum(raw_dist)
        else: normed_dist = raw_dist * bin_size * unc_band_norm
        #print(raw_dist)
        unc_dist = np.sqrt(normed_dist)
        frac_unc = unc_dist / normed_dist
        bincenters_band = np.copy(bincenters)
        bincenters_band[0] = low
        bincenters_band[-1] = high

        if(len(labels) <= len(entries)): label = "Signal Injection Stat. Unc."
        else: label = labels[-1]

        l_unc = ax1.fill_between(bincenters_band, 1. - frac_unc, 1. + frac_unc, color = 'gray', alpha = 0.5, label = labels)
        add_patch(leg, l_unc, label)

    for i in range(1, len(ns)):
        ratio = np.clip(ns[i], 1e-8, None) / np.clip(ns[0], 1e-8, None)
        ratios.append(ratio)

        ratio_err = None
        if(errors):
            if(weights != None):
                w0 = weights[0]**2
                w1 = weights[i]**2
                norm0 = np.sum(weights[0])*bin_size
                norm1 = np.sum(weights[i])*bin_size
            else:
                w0 = w1 = None
                norm0 = entries[0].shape[0]*bin_size
                norm1 = entries[i].shape[0]*bin_size

            err0 = np.sqrt(np.histogram(entries[0], bins=bins, weights=w0)[0])/norm0
            err1 = np.sqrt(np.histogram(entries[i], bins=bins, weights=w1)[0])/norm1
            err0_alt  = np.sqrt(norm0*n0)/norm0
            err_alt1  = np.sqrt(norm1*n1)/norm1
            ratio_err = ratio * np.sqrt((err0/n0)**2 + (err1/n1)**2)
        errs.append(ratio_err)

        ax1.errorbar(bincenters, ratio, yerr = ratio_err, alpha=alpha, markerfacecolor = colors[i], markeredgecolor = colors[i], fmt='ko')



    

    label_size = 18
    ax1.set_ylabel("Ratio", fontsize= label_size)
    ax1.set_xlabel(axis_label, fontsize = label_size)

    plt.xlim([low, high])

    if(type(ratio_range) == list):
        plt.ylim(ratio_range[0], ratio_range[1])
    else:
        if(ratio_range > 0):
            plt.ylim([1-ratio_range, 1+ratio_range])

    plt.grid(axis='y')

    if(save): 
        plt.savefig(fname)
        print("saving fig %s" %fname)

    return ns, bins, ratios, frac_unc

def makeCan(name, fname, histlist, bkglist=[],signals=[],totlist = [], colors=[],titles=[],dataName='Data',bkgNames=[],signalNames=[],
        logy=False,rootfile=False,xtitle='',ytitle='',dataOff=False,datastyle='pe',year=1, ratio_range = None, NDiv = 205, prelim = False):  

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
        print 'histlist of size ' + str(len(histlist)) + ' not currently supported'
        print histlist
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
            print 'ERROR: It seems you are trying to plot backgrounds with data on a 2D plot. This is not supported since there is no good way to view this type of distribution.'
        
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
                    y_size = 0.2 + 0.03*(len(bkglist[0])+len(signals))
                    x_size = 0.5
                    if(leg_align_right):
                        x_start = 0.42
                    else:
                        x_start = 0.2

                    legends.append(ROOT.TLegend(x_start,y_end - y_size,x_start + x_size,y_end))
                else: 
                    legends.append(ROOT.TLegend(0.2,0.11,0.45,0.2+0.02*(len(bkglist[0])+len(signals))))

                stacks.append(ROOT.THStack(hist.GetName()+'_stack',hist.GetName()+'_stack'))
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
                    
                # Go to main pad, set logy if needed
                mains[hist_index].cd()


                # Set y max of all hists to be the same to accomodate the tallest
                max_scaling = 3.0
                histList = [stacks[hist_index],totlist[hist_index],hist]

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
                TOffset = 1.
                # Now draw the main pad
                data_leg_title = hist.GetTitle()
                if len(titles) > 0:
                    hist.SetTitle(titles[hist_index])
                hist.GetYaxis().SetTitleOffset(TOffset)
                hist.GetXaxis().SetTitleOffset(1.2)
                hist.GetYaxis().SetTitle('Events / bin')
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

                stacks[hist_index].Draw('same hist')
                #print("Drawing %s same hist \n" stacks[hist_index].GetName())

                # Do the signals
                if len(signals) > 0: 
                    signals[hist_index].SetLineColor(kBlue)
                    signals[hist_index].SetLineWidth(2)
                    if logy == True:
                        signals[hist_index].SetMinimum(1e-3)
                    if signalNames == []: this_sig_name = signals[hist_index].GetName().split('_')[0]
                    legends_list[hist_index].append((signals[hist_index],this_sig_name,'L'))
                    signals[hist_index].Draw('hist same')

                ROOT.gStyle.SetHatchesLineWidth(2)
                totlist[hist_index].SetLineColor(ROOT.kWhite)
                totlist[hist_index].SetFillColor(ROOT.kBlack)
                totlist[hist_index].SetFillStyle(3354)
                totlist[hist_index].SetMarkerStyle(20)
                totlist[hist_index].SetMarkerSize(0.01)

                totlist[hist_index].Draw('e2 same')

                if not dataOff:
                    legends_list[hist_index].append((hist,dataName,datastyle))
                    hist.Draw(datastyle+' same')


                #Draw helpful lines

                #legends[hist_index].SetHeader(titles[0], "c")
                legends[hist_index].SetNColumns(2)
                legends[hist_index].SetTextSize(0.05)
                
                for entry in legends_list[hist_index][::-1]:
                    legends[hist_index].AddEntry(entry[0], entry[1], entry[2])

                legends[hist_index].AddEntry(totlist[hist_index], "Sys. unc.", "f")


                ratio, ratio_sys_unc = makeRatio(hist,totlist[hist_index])
                chi2 = 0.
                #for i in range(1, pull.GetNbinsX()+1):
                    #chi2 += pull.GetBinContent(i)**2;
                #print("Chi2/nbin for chan %s is %.1f/%i" % (titles[hist_index], chi2, pull.GetNbinsX()))

                legends[hist_index].AddEntry(ratio_sys_unc, "Total fit. unc.", "f")




                legends[hist_index].SetBorderSize(0)
                legends[hist_index].Draw()
                ROOT.gPad.RedrawAxis()

                # Draw the pull
                subs[hist_index].cd()
                # Build the pull

                LS = mLS * 0.7/0.3
                #title size given as fraction of pad width, scale up to have same size as main pad
                YTS =  mTS * 0.7/0.3
                XTS =  mTS * 0.7/0.3
                lTOffset = TOffset * 0.27 / 0.7


                ratio.SetMarkerStyle(8)


                ratio_sys_unc.GetYaxis().SetRangeUser(ratio_range[0], ratio_range[1])
                ratio_sys_unc.GetYaxis().SetTitleOffset(lTOffset)
                ratio_sys_unc.GetYaxis().SetTickLength(0.04)
                             
                ratio_sys_unc.GetYaxis().SetLabelSize(LS)
                ratio_sys_unc.GetYaxis().SetTitleSize(YTS)
                ratio_sys_unc.GetYaxis().SetNdivisions(NDiv)
                ratio_sys_unc.GetYaxis().SetTitle("Data / fit")

                ratio_sys_unc.GetXaxis().SetRangeUser(hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax())
                ratio_sys_unc.GetXaxis().SetTitleOffset(1.)
                ratio_sys_unc.GetXaxis().SetLabelOffset(0.05)
                ratio_sys_unc.GetXaxis().SetLabelSize(LS)
                ratio_sys_unc.GetXaxis().SetTitleSize(XTS)
                ratio_sys_unc.GetXaxis().SetTitle(xtitle)
                ratio_sys_unc.GetXaxis().SetTickLength(0.06)


                #ratio_sys_unc.SetFillColor(ROOT.kBlack)
                #ratio_sys_unc.SetFillStyle(3015)
                ratio_sys_unc.SetLineColor(ROOT.kGray)
                ratio_sys_unc.SetFillColor(ROOT.kGray)
                #ratio_sys_unc.SetFillStyle(3015)

                ratio_sys_unc.SetTitle("")


                ratio.SetLineStyle(1)
                ratio.SetLineWidth(2)
                ratio_sys_unc.Draw("A3 same")
                ratio.Draw('p0e0Z same')
                
                line = ROOT.TLine(0, 1.0, hist.GetNbinsX() - 0.08, 1.0)
                line.SetLineStyle(9)
                line.Draw()

                if logy == True:
                    mains[hist_index].SetLogy()

                #if(prelim): 
                #    print("Prelim")
                #    CMS_lumi.writeExtraText = True
                #else: CMS_lumi.writeExtraText = False
                #if(CMS_align_right): CMS_loc = 33
                #else: CMS_loc = 11
                #CMS_lumi.CMS_lumi(mains[hist_index], year, CMS_loc)

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


def make_root_hist(data = None, weight = None, name = "h", num_bins = 1, bin_low = 0, bin_high = 1):
    h = ROOT.TH1F(name, "", num_bins, bin_low, bin_high)
    for i,e in enumerate(data):
        if(weight is not None): w = weight[i]
        else: w = 1.
        h.Fill(e, w)
    return h




def make_stack_ratio_histogram(data = None, entries = None, labels = None, colors = None, axis_label = None, title = None, num_bins = None, 
    normalize = False, h_range = None, weights = None, fname="", ratio_range = -1, errors = False, extras = None, logy = False, max_rw = 5):
    if(h_range == None):
        low = np.amin(entries[0])
        high = np.amax(entries[0])
        h_range = (low, high)
    else:
        low,high = h_range

    hists = []
    for i,e in enumerate(entries):
        if(weights is not None): weight = weights[i]
        hist = make_root_hist(data = e, weight = weight, name = "h" + labels[i], num_bins = num_bins, bin_low = low, bin_high = high)
        hist.SetFillColor(colors[i])
        hist.SetLineColor(colors[i])
        hists.append(hist)

    h_data = make_root_hist(data = data, weight = None, name = 'data', num_bins = num_bins, bin_low = low, bin_high = high)
    h_tot = hists[0].Clone("h_tot")
    h_tot.Reset()
    for h in hists: h_tot.Add(h)

    return makeCan("temp", fname, [h_data], bkglist = [hists], totlist = [h_tot], colors = colors, bkgNames = labels, titles = [title], logy = logy, xtitle = axis_label,
        ratio_range = ratio_range)


def make_ratio_histogram(entries, labels, colors, axis_label, title, num_bins, normalize = False, save=False, h_range = None, 
        weights = None, fname="", ratio_range = -1, errors = False, extras = None, logy = False, max_rw = 5):
    h_type= 'step'
    alpha = 1.
    fontsize = 16
    fig = plt.figure(figsize = fig_size)
    gs = gridspec.GridSpec(2,1, height_ratios = [3,1])
    ax0 =  plt.subplot(gs[0])

    if(h_range == None):
        low = np.amin([np.amin(entries[0]), np.amin(entries[1])])
        high = np.amax([np.amax(entries[0]), np.amax(entries[1])])
        h_range = (low, high)
    else:
        low,high = h_range


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

    if(save): 
        plt.savefig(fname)
        print("saving fig %s" %fname)

    return bins, ratio

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

