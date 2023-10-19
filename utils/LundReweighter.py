import h5py
import argparse
import fastjet as fj
import numpy as np
import sys
from .PlotUtils import *
import ROOT
from array import array
import copy

def cleanup_hist(h):
    if(type(h) == ROOT.TH3F):
        for i in range(h.GetNbinsX()+1):
            for j in range(h.GetNbinsY()+1):
                for k in range(h.GetNbinsZ()+1):
                    c = h.GetBinContent(i,j,k)
                    if( c < 0): 
                        h.SetBinContent(i,j,k, 0.)
                        h.SetBinError(i,j,k, 0.)
    elif(type(h) == ROOT.TH2F):
        for i in range(h.GetNbinsX()+1):
            for j in range(h.GetNbinsY()+1):
                c = h.GetBinContent(i,j)
                if( c < 0): 
                    h.SetBinContent(i,j, 0.)
                    h.SetBinError(i,j, 0.)

def cleanup_ratio(h, h_min = 0., h_max = 2.):
    for i in range(0, h.GetNbinsX() + 2):
        for j in range(0, h.GetNbinsY() + 2):
            cont = h.GetBinContent(i,j)
            cont = max(h_min, min(cont, h_max))
            h.SetBinContent(i,j,cont)
    #h.GetZAxis().SetRangeUser(h_min, h_max);

def copy_proj(bin_i, h_ratio_proj, h_ratio):
    for j in range(h_ratio.GetNbinsY()):
        for k in range(h_ratio.GetNbinsZ()):
            h_ratio.SetBinContent(bin_i,j,k, h_ratio_proj.GetBinContent(j,k))
            h_ratio.SetBinError(bin_i,j,k, h_ratio_proj.GetBinError(j,k))
    return

def get_unc_hist(h):
    h_unc = h.Clone(h.GetName() + "_unc")
    for i in range(1, h.GetNbinsX() + 1):
        for j in range(1, h.GetNbinsY() + 1):
            err = h.GetBinError(i,j)
            cont = h.GetBinContent(i,j)

            if(cont > 0):
                h_unc.SetBinContent(i,j, err/cont)
                h_unc.SetBinError(i,j, 0.)
            else:
                h_unc.SetBinContent(i,j, 0.)
    return h_unc

def ang_dist(phi1, phi2):
    phi1 = phi1 % (2. * np.pi)
    phi2 = phi2 % (2. * np.pi)
    dphi = phi1 - phi2
    if(len(dphi.shape) > 0):
        dphi[dphi < -np.pi] += 2.*np.pi
        dphi[dphi > np.pi] -= 2.*np.pi
    else:
        if(dphi < -np.pi): dphi += 2.*np.pi
        if(dphi > np.pi): dphi -= 2.*np.pi

    return dphi

def get_dRs(gen_eta_phi, j_4vec):
    dR = np.sqrt(np.square(gen_eta_phi[:,0] - j_4vec[1]) + 
            np.square(ang_dist(gen_eta_phi[:,1], j_4vec[2] )))
    return dR

def get_subjet_dist(q_eta_phis, subjets_eta_phis):
    q_eta_phis = np.expand_dims(q_eta_phis, 0)
    subjets_eta_phis = np.expand_dims(subjets_eta_phis, 1)
    return np.sqrt(np.square(subjets_eta_phis[:,:,0] - q_eta_phis[:,:,0]) + 
            np.square(ang_dist(subjets_eta_phis[:,:,1], q_eta_phis[:,:,1] )))


class LundReweighter():

    def __init__(self, jetR = -1, maxJets = -1, pt_extrap_dir = None, pt_extrap_val = 350., pf_pt_min = 1.0, charge_only = False, 
             min_kt = 0.002, max_kt = 99999., min_delta = 0.005, max_delta = 99999.) :

        self.jetR = jetR
        self.maxJets = maxJets
        self.charge_only = charge_only
        self.dR = 0.8
        self.pt_extrap_dir = pt_extrap_dir
        self.pt_extrap_val = pt_extrap_val
        self.pf_pt_min = pf_pt_min
        self.charge_only = charge_only
        self.max_rw = 5.
        self.min_rw = 0.2
        self.min_kt, self.max_kt = min_kt, max_kt
        self.min_delta, self.max_delta = min_delta, max_delta
        self.func_dict = {}



    def check_bad_subjet_matching(self, gen_parts_eta_phi, subjets):
        # check if subjets fail matching criteria
        if(gen_parts_eta_phi is None): return None,None,None

        deltaR_cut = 0.2
        #Matrix of quark-subjet dRs
        dists = get_subjet_dist(gen_parts_eta_phi, np.array(subjets)[:,1:3])
        #DR to closest subjet for each quark
        j_closest = np.amin(dists, axis = 0)
        #index of closest subjet for each quark
        j_which = np.argmin(dists, axis = 0)
        #Which subjets are well matched
        matches = j_which[j_closest < deltaR_cut]

        #check if each subjet within Delta < 0.2 of a quark 
        j_matched = [i in matches for i in range(len(subjets))]
        #check if two quarks matched to same subjet

        repeats = matches.shape[0] != np.unique(matches).shape[0]

        return j_matched, repeats, j_closest

    def get_splittings_and_matching(self, pf_cands, gen_particles_eta_phi, ak8_jet, rescale_subjets = "", rescale_val = 1.0):
        """Given a list of pf_candidates (px, py,pz,E), and gen_particles (eta, phi), and an AK8 jet 4 vector (pt, eta,phi, M) 
        Recluster into a number of subjets based on the number of gen-level quarks inside the AK8 jet
        Also returns the fraction of bad matches.
        The momentum of these subjets is scaled based on the rescale_subjets and rescale_val args.

        rescale_subjets (optional): Method to rescale the momentum of the subjets ('jec' or 'vec'). 
                                    'vec' ensures the pt vector sum of the subjets adds up to rescale_val (ie total AK8 jet pt).  
                                    'jec' multiplies each subjet by the value of rescale_val (ie a jec value). 

        rescale_val (optional): Value used in subjet scaling.

        """


        RO = ReclusterObj()
        RO_prongsUp = None
        RO_prongsDown = None

        #deltaR's between AK8 jet and gen quarks
        AK8_dRs = get_dRs(gen_particles_eta_phi, ak8_jet)
        #ensure at least 1 prong or reclustering will crash
        RO.n_prongs = max(1, np.sum(AK8_dRs < 0.8))

        #check for quarks near boundary of jet
        prongs_up = np.any((0.8 < AK8_dRs) & (AK8_dRs < 0.9))
        prongs_down = np.any((0.7 < AK8_dRs) & (AK8_dRs < 0.8))


        RO.subjet, RO.split = self.get_splittings(pf_cands, num_excjets = RO.n_prongs, rescale_subjets = rescale_subjets, rescale_val = rescale_val)

        #check subjets matched to quarks
        RO.subjet_match, RO.double_match, RO.subjet_dRs = self.check_bad_subjet_matching(gen_particles_eta_phi, RO.subjet)

        #If any subjets not matched, or double matched, perform N+1/N-1 subjets reclustering
        badmatch_reclust = (not np.all(RO.subjet_match)) or RO.double_match

        #Recluster with one more subjet
        if(badmatch_reclust or prongs_up):
            RO_prongsUp = ReclusterObj()
            RO_prongsUp.n_prongs = RO.n_prongs + 1
            RO_prongsUp.badmatch_reclust = badmatch_reclust
            RO_prongsUp.prongs_up = prongs_up
            RO_prongsUp.subjet, RO_prongsUp.split = self.get_splittings(pf_cands, num_excjets = RO.n_prongs+1, rescale_subjets = rescale_subjets, rescale_val = rescale_val)
            #check subjets matched to quarks
            RO_prongsUp.subjet_match, _, RO_prongsUp.subjet_dRs = self.check_bad_subjet_matching(gen_particles_eta_phi, RO_prongsUp.subjet)

        #Recluster with one less subjet
        if(RO.badmatch_reclust or prongs_down):
            RO_prongsDown = ReclusterObj()
            RO_prongsDown.n_prongs = RO.n_prongs - 1
            RO_prongsDown.badmatch_reclust = badmatch_reclust
            RO_prongsDown.prongs_down = prongs_down
            RO_prongsDown.subjet, RO_prongsDown.split = self.get_splittings(pf_cands, num_excjets = RO.n_prongs-1, rescale_subjets = rescale_subjets, rescale_val = rescale_val)
            #check subjets matched to quarks
            RO_prongsDown.subjet_match, _, RO_prongsDown.subjet_dRs = self.check_bad_subjet_matching(gen_particles_eta_phi, RO_prongsDown.subjet)


        return RO, RO_prongsUp, RO_prongsDown




    def get_splittings(self, pf_cands, num_excjets = -1, rescale_subjets = "", rescale_val = 1.0):
        """Given a list of pf_candidates (px, py,pz,E), recluster into a given (num_excjets) number of subjets (-1 for variable number, not recommended). 
        the momentum of these subjets is scaled based on the rescale_subjets and rescale_val args.

        rescale_subjets (optional): Method to rescale the momentum of the subjets ('jec' or 'vec'). 
                                    'vec' ensures the pt vector sum of the subjets adds up to rescale_val (ie total AK8 jet pt).  
                                    'jec' multiplies each subjet by the value of rescale_val (ie a jec value). 

        rescale_val (optional): Value used in subjet scaling.
        """

        pjs = []
        pfs_cut = []
        for i,c in enumerate(pf_cands):
            if(c[3] > 0.0001):
                pj = fj.PseudoJet(c[0], c[1], c[2], c[3])

                if(pj.pt() > 1.0):
                    pfs_cut.append(c)

                pjs.append(pj)

        if(self.jetR < 0): R = 1000.0
        else: R = self.jetR
        #jet_algo = fj.cambridge_algorithm
        jet_algo = fj.kt_algorithm
        jet_def = fj.JetDefinition(jet_algo, R)
        cs = fj.ClusterSequence(pjs, jet_def)
        if(num_excjets < 0):
            js = fj.sorted_by_pt(cs.inclusive_jets())
            if(self.maxJets > 0):
                nMax = min(len(js), self.maxJets)
                js = js[:nMax]
        else:
            js = list(fj.sorted_by_pt(cs.exclusive_jets_up_to(int(num_excjets))))

            #for kt jets, recluster to get CA splittings
            if (jet_algo is fj.kt_algorithm):
                CA_R = 1000.
                js_new = []
                clust_seqs = []
                for i, j in enumerate(js):
                    CA_jet_def = fj.JetDefinition(fj.cambridge_algorithm, CA_R)
                    constituents = j.validated_cs().constituents(j)

                    cs_CA = []
                    for c in constituents:
                        #apply a cut on the constituents
                        if(c.pt() > self.pf_pt_min):
                            if(self.charge_only):
                                pf = find_matching_pf(pfs_cut, c)
                                if(pf is None):
                                    print("NO match!")
                                    print(c)
                                    print(pfs)
                                    exit(1)
                                #4th entry is PUPPI weight, 5th entry is charge of PFCand
                                eps = 1e-4
                                if(pf is not None and abs(pf[5] ) > eps):
                                    cs_CA.append(c)
                            else:
                                cs_CA.append(c)

                    if(len(cs_CA) > 0):
                        CA_cs = fj.ClusterSequence(cs_CA, CA_jet_def)
                        CA_jet = fj.sorted_by_pt(CA_cs.inclusive_jets())
                        js_new.append(j)
                        clust_seqs.append(CA_cs) #prevent from going out of scope

                js = js_new

        subjets = []
        splittings = []
        #print("%i subjets " % len(js))
        total_jet = fj.PseudoJet()
        for i, j in enumerate(js):
            #print("sj %i" % i)
            pseudojet = j
            jet_pt = j.pt()
            subjets.append([j.pt(), j.eta(), j.phi(), j.m()])
            total_jet += j
            while True:
                j1 = fj.PseudoJet()
                j2 = fj.PseudoJet()
                if pseudojet and pseudojet.has_parents(j1, j2):
                    # order the parents in pt
                    if (j2.pt() > j1.pt()):
                        j1, j2 = j2, j1
                    delta = j1.delta_R(j2)
                    kt = j2.pt() * delta
                    splittings.append([i, delta, kt])
                    pseudojet = j1
                else:
                    break
    

        #Rescale subjet momenta
        if(rescale_subjets == "jec"):
            for i in range(len(subjets)):
                subjets[i][0] *= rescale_val
        elif(rescale_subjets == "vec"):
            rescale_val = rescale_val / total_jet.pt()
            for i in range(len(subjets)):
                subjets[i][0] *= rescale_val

        return subjets, splittings


    def fill_lund_plane(self, h, pf_cands = None, subjets = None,  splittings = None, num_excjets = -1, weight = 1., subjet_idx = -1,
            rescale_subjets = "vec", rescale_val = 1.0):
        """Fill Lund Plane based on splittings"""

        if(type(h) != list):
            hists = [h]
            weights = [weight]
        else:
            hists = h
            weights = weight

        if(subjets is None or splittings is None):
            subjets, splittings = self.get_splittings(pf_cands, num_excjets = num_excjets, rescale_subjets = rescale_subjets, rescale_val = rescale_val)
            if(len(subjets) == 0): subjets = [[0,0,0,0]]



        no_idx = (len(subjets) == 1)
        subjets_reshape = np.array(subjets).reshape(-1)

        for jet_i, delta, kt in splittings:
            if(subjet_idx >= 0 and jet_i != subjet_idx): continue
            jet_int = int(np.round(jet_i))
            jet_pt = subjets_reshape[0] if no_idx else subjets_reshape[jet_int*4]
            if(delta > 0. and kt > 0.):
                for h_idx, h in enumerate(hists):
                    if(type(h) == ROOT.TH3F): h.Fill(jet_pt, np.log(self.dR/delta), np.log(kt), weights[h_idx])
                    else: h.Fill(np.log(self.dR/delta), np.log(kt), weights[h_idx])
        return subjets, splittings
    
    def get_lund_plane_idxs(self, h,  subjets = None,  splittings = None, subjet_idx = -1):
        """Get LP bin indices  for some splittings"""
        no_idx = (len(subjets) == 1)
        subjets_reshape = np.array(subjets).reshape(-1)

        idxs = []

        binx,biny,binz = array('i', [0]), array('i', [0]), array('i', [0])
        xmax, ymax, zmax = h.GetNbinsX(), h.GetNbinsY(),h.GetNbinsZ()
        for jet_i, delta, kt in splittings:
            if(subjet_idx >= 0 and jet_i != subjet_idx): continue
            jet_int = int(np.round(jet_i))
            jet_pt = subjets_reshape[0] if no_idx else subjets_reshape[jet_int*4]
            if(delta > 0. and kt > 0. and delta > self.min_delta and delta < self.max_delta 
                    and kt > self.min_kt and kt < self.max_kt):
                bin_idx = h.FindBin(jet_pt, np.log(self.dR/delta), np.log(kt))

                h.GetBinXYZ(bin_idx, binx, biny, binz)
                idxs.append((int(np.clip(binx[0], 1, xmax)), int(np.clip(biny[0], 1, ymax)), int(np.clip(binz[0], 1, zmax))))
        return idxs



    def reweight_pt_extrap(self,  subjet_pt, lp_idxs, rw, smeared_rw, pt_smeared_rw, pt_rand_noise = None, sys_str = ""):
        """Reweight based on pt extrapolated functions"""


        for (i,j,k) in lp_idxs:

            f_str = "func_%s%i_%i" % (sys_str, j,k)
            if(f_str in self.func_dict.keys()):
                f = self.func_dict[f_str]
            else:
                f = self.pt_extrap_dir.Get(f_str)
                self.func_dict[f_str] = f
            #val = f.Eval(subjet_pt)
            val = f.Eval(1./subjet_pt)
            val = np.clip(val, self.min_rw, self.max_rw)

            rw *= val 
            #keep noise smeared vals consistent
            if(smeared_rw is not None): smeared_rw *= val

            if(pt_rand_noise is not None):
                for n in range(pt_rand_noise.shape[0]):
                #up_down = [-1., 1.]
                #for n in range(2):
                    pars = array('d')
                    for p in range(f.GetNpar()):
                        pnom = f.GetParameter(p)
                        perr = f.GetParError(p)
                        pnew = pnom + perr * pt_rand_noise[n, j-1, k-1, p]
                        pars.append(pnew)

                    smeared_val = f.EvalPar(array('d', [1./subjet_pt]), pars)
                    smeared_val = np.clip(smeared_val, self.min_rw, self.max_rw)
                    pt_smeared_rw[n] *= smeared_val



        return rw, smeared_rw, pt_smeared_rw


    def reweight(self, h_rw, lp_idxs, rw, smeared_rw, pt_smeared_rw, rand_noise = None):
        """Reweight based on directly measured data/MC LP ratio"""

        for (i,j,k) in lp_idxs:
            #print("Rw %.3f, cont %.3f, i %i j %i k %i n %i" % (rw, h_rw.GetBinContent(i,j,k), i,j,k, n_cands))
            val = h_rw.GetBinContent(i,j,k)
            err = h_rw.GetBinError(i,j,k)

            if(val <= 1e-4 and err <= 1e-4):
                val = 1.0
                err = 1.0
                #print("EMPTY BIN")

            val = np.clip(val, self.min_rw, self.max_rw)
            rw *= val
            #keep pt smearing vals consistent
            if(pt_smeared_rw is not None): pt_smeared_rw *= val

            if(rand_noise is not None):
                smeared_vals = val + rand_noise[:,i-1,j-1,k-1] * err
                smeared_vals = np.clip(smeared_vals, self.min_rw, self.max_rw)
                smeared_rw *= smeared_vals


        return rw, smeared_rw, pt_smeared_rw



    def reweight_lund_plane(self, h_rw, pf_cands = None, reclust_obj = None, splittings = None, subjets = None, num_excjets = -1, 
                            rand_noise = None, pt_rand_noise = None,  sys_str = "", rescale_subjets = "", rescale_val = 1.0):
        """ Main function to compute reweighting factors. Can take in already computed subjets + splittings or recluster 
        itself using the PF candidates and the num_excjets args

        Args:
        h_rw : 3D histogram of data/MC ratio
        pf_cands (optional): List of PF candidates 
        num_excjets (optional): Number of subjets to recluster to
        splittings (optional): List of splittings of the subjets (subjet_idx, deltaR, kt)
        subjets (optional): List of subjets (pt, eta, phi, m)
        rand_noise (optional): Used for the computation of the statistical uncertainty on the weights. 
                              A 3D list of random numbers (taken from a std normal distribution) of size (nToys, n_bins_X, n_bins_Y, n_bins_Z)
                               where the latter three numbers are the numbers of bins in h_rw.
        pt_rand_noise (optional): Used for the computation of the pt extrapolation uncertainty on the weights. 
                              A 3D list of random numbers (taken from a std normal distribution) of size (nToys, n_bins_X, n_bins_Y, n_bins_Z)
                               where the latter three numbers are the numbers of bins in h_rw.

        rescale_subjets (optional): Method to rescale the momentum of the subjets ('jec' or 'vec'). 
                                    'vec' ensures the pt vector sum of the subjets adds up to rescale_val (ie total AK8 jet pt).  
                                    'jec' multiplies each subjet by the value of rescale_val (ie a jec value). 

        rescale_val (optional): Value used in subjet scaling.


        Returns : 
        A tuple
        (Event reweighting factor, reweighting factors from statistical variation toys,  reweighting factor from pt extrapolation toys)



         """

        #assume matched if we don't have other info
        subjet_match = [True,]*100
        if(reclust_obj is not None):
            subjets, splittings = reclust_obj.subjet, reclust_obj.split
            subjet_match = reclust_obj.subjet_match

        if(subjets is None or splittings is None):
            subjets, splittings = self.get_splittings(pf_cands, num_excjets = num_excjets, rescale_subjets = rescale_subjets, rescale_val = rescale_val )

        rw = 1.0
        pt_smeared_rw = smeared_rw = None

        if(rand_noise is not None): smeared_rw = np.array([1.0]*rand_noise.shape[0])
        if(pt_rand_noise is not None): pt_smeared_rw = np.array([1.0]*pt_rand_noise.shape[0])

        #splittings save idx of associated subjet

        for i in range(len(subjets)):

            #ignore subjets that aren't matched to anything
            if(not subjet_match[i]): continue

            if(len(splittings) > 0):
                lp_idxs = self.get_lund_plane_idxs(h_rw, subjet_idx = i, splittings = splittings, subjets = [subjets[i]])

                if(self.pt_extrap_dir is None or subjets[i][0] < self.pt_extrap_val or 'bquark' in sys_str):
                    rw, smeared_rw, pt_smeared_rw = self.reweight(h_rw, lp_idxs, rw, smeared_rw, pt_smeared_rw, rand_noise = rand_noise)
                else:
                    rw, smeared_rw, pt_smeared_rw = self.reweight_pt_extrap(subjets[i][0], lp_idxs, rw, smeared_rw, pt_smeared_rw, pt_rand_noise = pt_rand_noise, sys_str = sys_str)

        
        return rw, smeared_rw, pt_smeared_rw


    def make_LP_ratio(self, h_data, h_bkg, h_mc,  h_data_subjet_pt = None, h_bkg_subjet_pt = None, h_mc_subjet_pt = None, pt_bins = None, outdir = "", save_plots = False):
        """ Function to construct data/MC LP ratio"""


        h_data.Print()
        h_bkg.Print()
        h_mc.Print()

        do_jet_pt_norm  = (h_data_subjet_pt is not None) and (h_mc_subjet_pt is not None) and (h_bkg_subjet_pt is not None)

        cleanup_hist(h_mc)
        cleanup_hist(h_bkg)

        h_bkg_clone = h_bkg.Clone(h_bkg.GetName() + "_clone")
        h_mc_clone = h_mc.Clone(h_mc.GetName() + "_clone")



        h_ratio = h_data.Clone(h_mc_clone.GetName() + "_ratio")
        h_ratio.SetTitle("(Data - Bkg ) / TTbar MC")

        data_norm = h_data.Integral()
        est = h_bkg_clone.Integral() + h_mc_clone.Integral()


        if(do_jet_pt_norm):

            h_data_subjet_pt.Print()
            h_bkg_subjet_pt.Print()
            h_mc_subjet_pt.Print()

            h_data_subjet_pt_clone = h_data_subjet_pt.Clone(h_data_subjet_pt.GetName() + "_clone")
            h_bkg_subjet_pt_clone = h_bkg_subjet_pt.Clone(h_bkg_subjet_pt.GetName() + "_clone")
            h_mc_subjet_pt_clone = h_mc_subjet_pt.Clone(h_mc_subjet_pt.GetName() + "_clone")
            
            data_norm = h_data_subjet_pt.Integral()
            est = h_bkg_subjet_pt_clone.Integral() + h_mc_subjet_pt_clone.Integral()

            h_mc_subjet_pt_clone.Scale(data_norm / est)
            h_bkg_subjet_pt_clone.Scale(data_norm / est)


        h_bkg_clone.Scale(data_norm / est)
        h_mc_clone.Scale(data_norm / est)
            
        h_data_sub = h_data.Clone("h_data_sub")
        h_data_sub.Add(h_bkg_clone, -1.)
        h_data_sub.Print()



        cleanup_hist(h_data_sub)


        for i in range(1, h_data.GetNbinsX() + 1):
            h_bkg_clone1 = h_bkg_clone.Clone("h_bkg_clone%i" %i)
            h_mc_clone1 = h_mc_clone.Clone("h_mc_clone%i"% i)
            h_data_clone1 = h_data_sub.Clone("h_data_clone%i" %i)

            h_mc_clone1.GetXaxis().SetRange(i,i)
            h_bkg_clone1.GetXaxis().SetRange(i,i)
            h_data_clone1.GetXaxis().SetRange(i,i)


            h_mc_proj = h_mc_clone1.Project3D("zy")
            h_bkg_proj = h_bkg_clone1.Project3D("zy")
            h_data_proj = h_data_clone1.Project3D("zy")

            #normalize by number of subjets rather than number of splittings
            if(do_jet_pt_norm):
                mc_norm = h_mc_subjet_pt_clone.GetBinContent(i)
                bkg_norm = h_bkg_subjet_pt_clone.GetBinContent(i)
                data_norm = h_data_subjet_pt_clone.GetBinContent(i) - bkg_norm
            else:
                data_norm = h_data_proj.Integral()
                bkg_norm = h_bkg_proj.Integral()
                mc_norm = h_mc_proj.Integral()

            print(mc_norm, bkg_norm, data_norm)


            h_bkg_proj.Scale(1./bkg_norm)

            h_data_proj.Scale(1./data_norm)
            h_mc_proj.Scale(1./mc_norm)


            h_ratio_proj = h_data_proj.Clone("h_ratio_proj%i" %i)
            h_ratio_proj.Divide(h_mc_proj)

            #if(i == 1): 
            #    h_mc_proj.Print("range")
            #    h_data_proj.Print("range")
            #    h_ratio_proj.Print("range")

            copy_proj(i, h_ratio_proj, h_ratio)



            if(save_plots): 

                h_bkg_proj.SetTitle("Bkg MC pT %.0f - %.0f" % (pt_bins[i-1], pt_bins[i]))
                h_mc_proj.SetTitle("TTbar MC pT %.0f - %.0f" % (pt_bins[i-1], pt_bins[i]))
                h_data_proj.SetTitle("Data - Bkg pT %.0f - %.0f (N = %.0f)" % (pt_bins[i-1], pt_bins[i], data_norm))
                h_ratio_proj.SetTitle("Ratio pT %.0f - %.0f (N = %.0f)" % (pt_bins[i-1], pt_bins[i], data_norm))


                c_mc = ROOT.TCanvas("c", "", 1000, 1000)
                h_mc_proj.Draw("colz")
                c_mc.SetRightMargin(0.2)
                c_mc.Print(outdir + "lundPlane_bin%i_MC.png" % i)


                c_bkg = ROOT.TCanvas("c", "", 1000, 800)
                h_bkg_proj.Draw("colz")
                c_bkg.SetRightMargin(0.2)
                c_bkg.Print(outdir + "lundPlane_bin%i_bkg.png" % i)

                c_data = ROOT.TCanvas("c", "", 1000, 800)
                h_data_proj.Draw("colz")
                c_data.SetRightMargin(0.2)
                c_data.Print(outdir + "lundPlane_bin%i_data.png" %i )



                c_ratio = ROOT.TCanvas("c", "", 1000, 800)
                cleanup_ratio(h_ratio_proj, h_min =0., h_max = 2.0)
                h_ratio_proj.Draw("colz")
                c_ratio.SetRightMargin(0.2)
                c_ratio.Print(outdir + "lundPlane_bin%i_ratio.png" % i)

                h_ratio_unc = get_unc_hist(h_ratio_proj)
                cleanup_ratio(h_ratio_unc, h_min = 0., h_max = 1.0)
                c_ratio_unc = ROOT.TCanvas("c_unc", "", 800, 800)
                h_ratio_unc.SetTitle("Ratio pT %.0f - %.0f (N = %.0f) Relative Unc." % (pt_bins[i-1], pt_bins[i], data_norm))
                h_ratio_unc.Draw("colz")
                c_ratio_unc.SetRightMargin(0.2)
                c_ratio_unc.Print(outdir + "lundPlane_bin%i_ratio_unc.png" % i)
                h_ratio_unc.Reset()

        return h_ratio

    def normalize_weights(self, lund_weights, w_min = 0.1, w_max = 10.):
        """Normalize lund plane weights so average weight is 1 (necessary to preserve normalization of MC.
        Also clip outlier weights so to not be dominated by statistical fluctuations. """

        lund_weights = np.clip(lund_weights, 0., w_max)
        if(len(lund_weights.shape) > 1): lund_weights /= np.mean(lund_weights, axis = 0, keepdims=True)
        else: lund_weights /= np.mean(lund_weights)

        lund_weights  = np.clip(lund_weights, w_min, w_max)
        if(len(lund_weights.shape) > 1): lund_weights /= np.mean(lund_weights, axis = 0, keepdims=True)
        else: lund_weights /= np.mean(lund_weights)

        return lund_weights

class ReclusterObj():
    def __init__(self,subjets = None, split = None, dRs = None):
        self.subjet = None
        self.split = None
        self.dRs = None
        self.subjet_match = None
        self.badmatch_reclust = False
        self.prongs_up = False
        self.prongs_down = False





def matched(c,cj):
    eps = 1e-4
    return (abs(c[0] - cj.px()) < eps) and (abs(c[1] - cj.py()) < eps) and (abs(c[2] - cj.pz()) < eps)


def find_matching_pf(cj_list, cj):

    for c in cj_list:
        if(matched(c, cj)): 
            return c
    return None



