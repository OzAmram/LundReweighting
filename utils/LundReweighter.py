import h5py
import argparse
import fastjet as fj
import numpy as np
import sys
from .PlotUtils import *
from .Consts import *
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

def symmetrize(up, down, nom):
    ratio_up = up/nom
    ratio_down = down/nom
    if(abs(ratio_up -1.0) > abs(ratio_down -1.0)):
        down = nom * (nom/up)
    else:
        up = nom * (nom/down)
    return up,down


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

    def __init__(self, f_ratio = None, jetR = -1, maxJets = -1, pt_extrap_dir = None, pt_extrap_val = 350., pf_pt_min = 1.0, charge_only = False, 
             min_kt = 0.002, max_kt = 99999., min_delta = 0.005, max_delta = 99999., LP_order = 1) :

        self.jetR = jetR
        self.maxJets = maxJets
        self.charge_only = charge_only
        self.dR = 0.8
        self.pt_extrap_val = pt_extrap_val
        self.pf_pt_min = pf_pt_min
        self.charge_only = charge_only
        self.max_rw = 5.
        self.min_rw = 0.2
        self.min_kt, self.max_kt = min_kt, max_kt
        self.min_delta, self.max_delta = min_delta, max_delta
        self.LP_order = LP_order
        self.func_dict = {}

        #nominal data/MC Lund plane ratio (3d histogram)
        self.f_ratio = f_ratio

        if(f_ratio is not None):
            self.h_ratio = f_ratio.Get("ratio_nom")
            self.h_mc = f_ratio.Get("mc_nom")
            #systematic variations
            self.h_ratio_sys_up = f_ratio.Get("ratio_sys_tot_up")
            self.h_ratio_sys_down = f_ratio.Get("ratio_sys_tot_down")
            if(not isinstance(self.h_ratio_sys_up, ROOT.TH3) or not isinstance(self.h_ratio_sys_down, ROOT.TH3)):
                print("\nMissing LP sys up/down variations! Will ignore this unc.\n")
                self.h_ratio_sys_up = self.h_ratio_sys_down = self.h_ratio
            #MC ratio of b to light quarks
            self.b_light_ratio = f_ratio.Get("h_bl_ratio")


            #directory of pt extrapolation fits
            f_ratio.cd('pt_extrap')
            self.pt_extrap_dir = ROOT.gDirectory




    def check_bad_subjet_matching(self, gen_parts_eta_phi, subjets):
        # check if subjets fail matching criteria
        if(gen_parts_eta_phi is None or len(subjets) == 0): return None,None,None

        deltaR_cut = 0.2
        #Matrix of quark-subjet dRs
        dists = get_subjet_dist(gen_parts_eta_phi, np.array(subjets)[:,1:3])
        #DR to closest subjet for each quark
        subjet_closest_dR = np.amin(dists, axis = 0)
        #index of closest subjet for each quark
        subjet_closest_idx = np.argmin(dists, axis = 0)
        #Which subjets are matched to a quark
        matches = subjet_closest_idx[subjet_closest_dR < deltaR_cut]

        #check if each subjet within Delta < 0.2 of a quark 
        is_subjet_matched = [i in matches for i in range(len(subjets))]
        #check if two quarks matched to a given subjet
        is_subjet_double_matched = [np.sum(matches == i)>=2 for i in range(len(subjets))]

        return is_subjet_matched, is_subjet_double_matched, subjet_closest_dR

    def get_splittings_and_matching(self, pf_cands, gen_particles_eta_phi, ak8_jet, rescale_subjets = "", rescale_val = 1.0, pf_cands_PtEtaPhiE_format = False):
        """Given a list of pf_candidates (px, py,pz,E), and gen_particles (eta, phi), and an AK8 jet 4 vector (pt, eta,phi, M) 
        Recluster into a number of subjets based on the number of gen-level quarks inside the AK8 jet
        Also returns the fraction of bad matches.
        The momentum of these subjets is scaled based on the rescale_subjets and rescale_val args.

        rescale_subjets (optional): Method to rescale the momentum of the subjets ('jec' or 'vec'). 
                                    'vec' ensures the pt vector sum of the subjets adds up to rescale_val (ie total AK8 jet pt).  
                                    'jec' multiplies each subjet by the value of rescale_val (ie a jec value). 

        rescale_val (optional): Value used in subjet scaling.
        pf_cands_PtEtaPhiE_format (optional): Alternate representation of pf candidates (default is px,py,pz,E)

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


        RO.subjet, RO.split = self.get_splittings(pf_cands, num_excjets = RO.n_prongs, rescale_subjets = rescale_subjets, 
                                rescale_val = rescale_val, pf_cands_PtEtaPhiE_format = pf_cands_PtEtaPhiE_format)

        #check subjets matched to quarks
        RO.subjet_match, RO.subjet_double_matched, RO.subjet_dRs = self.check_bad_subjet_matching(gen_particles_eta_phi, RO.subjet)

        #If any subjets not matched, or double matched, perform N+1/N-1 subjets reclustering
        RO.badmatch = (not np.all(RO.subjet_match)) or np.any(RO.subjet_double_matched)

        #Recluster with one more subjet
        if(RO.badmatch or prongs_up):
            RO_prongsUp = ReclusterObj()
            RO_prongsUp.subjet, RO_prongsUp.split = self.get_splittings(pf_cands, num_excjets = RO.n_prongs+1, rescale_subjets = rescale_subjets, 
                                            rescale_val = rescale_val, pf_cands_PtEtaPhiE_format = pf_cands_PtEtaPhiE_format)
            RO_prongsUp.n_prongs = RO.n_prongs + 1
            RO_prongsUp.from_badmatch = RO.badmatch
            RO_prongsUp.from_prongs_up = prongs_up
            #check subjets matched to quarks
            RO_prongsUp.subjet_match, RO_prongsUp.subjet_double_matched, RO_prongsUp.subjet_dRs = self.check_bad_subjet_matching(gen_particles_eta_phi, RO_prongsUp.subjet)

        #Recluster with one less subjet
        if((RO.badmatch or prongs_down) and RO.n_prongs > 1):
            RO_prongsDown = ReclusterObj()
            RO_prongsDown.subjet, RO_prongsDown.split = self.get_splittings(pf_cands, num_excjets = RO.n_prongs-1, rescale_subjets = rescale_subjets, 
                    rescale_val = rescale_val, pf_cands_PtEtaPhiE_format = pf_cands_PtEtaPhiE_format)
            RO_prongsDown.n_prongs = RO.n_prongs - 1
            RO_prongsDown.from_badmatch = RO.badmatch
            RO_prongsDown.from_prongs_down = prongs_down
            #check subjets matched to quarks
            RO_prongsDown.subjet_match, RO_prongsDown.subjet_double_matched, RO_prongsDown.subjet_dRs = self.check_bad_subjet_matching(gen_particles_eta_phi, RO_prongsDown.subjet)


        return RO, RO_prongsUp, RO_prongsDown




    def get_splittings(self, pf_cands, num_excjets = -1, rescale_subjets = "", rescale_val = 1.0, pf_cands_PtEtaPhiE_format = False):
        """Given a list of pf_candidates (px, py,pz,E), recluster into a given (num_excjets) number of subjets (-1 for variable number, not recommended). 
        the momentum of these subjets is scaled based on the rescale_subjets and rescale_val args.

        rescale_subjets (optional): Method to rescale the momentum of the subjets ('jec' or 'vec'). 
                                    'vec' ensures the pt vector sum of the subjets adds up to rescale_val (ie total AK8 jet pt).  
                                    'jec' multiplies each subjet by the value of rescale_val (ie a jec value). 

        rescale_val (optional): Value used in subjet scaling.
        pf_cands_PtEtaPhiE_format (optional): Alternate representation of pf candidates (default is px,py,pz,E)
        """

        if(num_excjets == 0):
            print("Reclustering into 0 subjets?! Something went wrong")
            exit(1)
            return [],[]
        pjs = []
        pfs_cut = []
        for i,c in enumerate(pf_cands):
            if(pf_cands_PtEtaPhiE_format):
                v = ROOT.Math.PtEtaPhiEVector(c[0],c[1],c[2],c[3])
                px,py,pz,E = v.px(), v.py(), v.pz(), v.E()

            else: px,py,pz,E = c[0], c[1], c[2], c[3]

            if(E > 0.0001):
                pj = fj.PseudoJet(px,py,pz,E)

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
        total_jet = fj.PseudoJet()
        for i, j in enumerate(js):
            pseudojet = j
            subjets.append([j.pt(), j.eta(), j.phi(), j.m()])
            total_jet += j
            j.subjet_pt = j.pt()
            pj_cands = [j]
            j.order = 1
            while True:
                if(len(pj_cands) == 0): break
                pseudojet = pj_cands.pop(0)

                j1 = fj.PseudoJet()
                j2 = fj.PseudoJet()
                if pseudojet and pseudojet.has_parents(j1, j2):
                    # order the parents in pt
                    if (j2.pt() > j1.pt()):
                        j1, j2 = j2, j1
                    delta = j1.delta_R(j2)
                    kt = j2.pt() * delta
                    splittings.append([i, pseudojet.subjet_pt, pseudojet.order, delta, kt])
                    
                    #'harder' branch has same order
                    j1.order = pseudojet.order
                    j1.subjet_pt = pseudojet.subjet_pt
                    pj_cands.append(j1)

                    #softer branch is one order higher, effective a new subjet
                    j2.order = pseudojet.order + 1
                    if(j2.order <= self.LP_order):
                        j2.subjet_pt = j2.pt()
                        pj_cands.append(j2)
                else:
                    continue
    

        #Rescale subjet momenta
        if(rescale_subjets == "jec"):
            for i in range(len(subjets)):
                subjets[i][0] *= rescale_val
        elif(rescale_subjets == "vec"):
            rescale_val = rescale_val / total_jet.pt()
            for i in range(len(subjets)):
                subjets[i][0] *= rescale_val


        return subjets, splittings


    def fill_lund_plane(self, h, h_subjet = None, pf_cands = None, subjets = None,  splittings = None, reclust_obj = None, num_excjets = -1, weight = 1., subjet_idx = -1,
            rescale_subjets = "vec", rescale_val = 1.0):
        """Fill Lund Plane based on splittings"""

        if(type(h) != list):
            hists = [h]
            hists_subjets = [h_subjet]
            weights = [weight]
        else:
            hists = h
            weights = weight
            hists_subjets = h_subjet

        #assume matched if we don't have other info
        subjet_match = [True,]*100
        if(reclust_obj is not None):
            subjets, splittings = reclust_obj.subjet, reclust_obj.split
            subjet_match = reclust_obj.subjet_match
        if(subjets is None or splittings is None):
            subjets, splittings = self.get_splittings(pf_cands, num_excjets = num_excjets, rescale_subjets = rescale_subjets, rescale_val = rescale_val)
            if(len(subjets) == 0): subjets = [[0,0,0,0]]



        no_idx = (len(subjets) == 1)
        subjets_reshape = np.array(subjets).reshape(-1)

        filled = []
        for subjet_i, subjet_pt, order, delta, kt in splittings:
            if((subjet_idx >= 0 and subjet_i != subjet_idx) or not subjet_match[subjet_i] ): continue
            if(delta > 0. and kt > 0.):
                for h_idx, h in enumerate(hists):
                    if(type(h) == ROOT.TH3F): h.Fill(subjet_pt, np.log(self.dR/delta), np.log(kt), weights[h_idx])
                    else: h.Fill(np.log(self.dR/delta), np.log(kt), weights[h_idx])

                #fill subjet pt once per subjet / order 
                if((subjet_i, order) not in filled and  hists_subjets[0] is not None): 
                    filled.append((subjet_i, order))
                    for h_idx, h_sj in enumerate(hists_subjets): h_sj.Fill(subjet_pt, weights[h_idx])


        return subjets, splittings
    
    def get_lund_plane_idxs(self, h,  subjets = None,  splittings = None, subjet_idx = -1, LP_order = -1):
        """Get LP bin indices  for some splittings"""
        no_idx = (len(subjets) == 1)
        subjets_reshape = np.array(subjets).reshape(-1)

        idxs = []

        binx,biny,binz = array('i', [0]), array('i', [0]), array('i', [0])
        xmax, ymax, zmax = h.GetNbinsX(), h.GetNbinsY(),h.GetNbinsZ()
        for jet_i, subjet_pt, order, delta, kt in splittings:
            if(subjet_idx >= 0 and jet_i != subjet_idx): continue
            if(LP_order > 0 and order != LP_order): continue

            jet_int = int(np.round(jet_i))
            jet_pt = subjets_reshape[0] if no_idx else subjets_reshape[jet_int*4]
            if(delta > 0. and kt > 0. and delta > self.min_delta and delta < self.max_delta 
                    and kt > self.min_kt and kt < self.max_kt):
                bin_idx = h.FindBin(subjet_pt, np.log(self.dR/delta), np.log(kt))

                h.GetBinXYZ(bin_idx, binx, biny, binz)
                idxs.append((int(np.clip(binx[0], 1, xmax)), int(np.clip(biny[0], 1, ymax)), int(np.clip(binz[0], 1, zmax))))
        return idxs



    def reweight_pt_extrap(self,  subjet_pt, lp_idxs, rw, smeared_rw, pt_smeared_rw, pt_rand_noise = None, sys_str = ""):
        """Reweight based on pt extrapolated functions"""


        for (i,j,k) in lp_idxs:

            f_str = "func_%s%i_%i" % (sys_str, j,k)
            covar_str = "func_%scovar_%i_%i" % (sys_str, j,k)

            if(f_str in self.func_dict.keys()):
                f,covar = self.func_dict[f_str]
            else:
                f = self.pt_extrap_dir.Get(f_str)
                if(not self.pt_extrap_dir.GetListOfKeys().Contains(covar_str)):
                    print("Missing covariance for pt extrap in LP bin %i, %i! (setting to zero)" % (j,k))
                    covar = 0.
                else:
                    covar_o = self.pt_extrap_dir.Get(covar_str)
                    covar = covar_o.GetVal()
                self.func_dict[f_str] = (f,covar)


            val = f.Eval(1./subjet_pt)
            val = np.clip(val, self.min_rw, self.max_rw)

            #nominal
            rw *= val 
            if(smeared_rw is not None): smeared_rw *= val

            if(pt_rand_noise is not None):

                #single parameter, vary by unc
                if(f.GetNpar() <= 1):
                    pnom = f.GetParameter(0)
                    perr = f.GetParError(0)

                    for n in range(pt_rand_noise.shape[0]):
                        pnew = pnom + perr * pt_rand_noise[n, j-1, k-1, 0]
                        pars = array('d', [pnew])

                        smeared_val = f.EvalPar(array('d', [1./subjet_pt]), pars)
                        smeared_val = np.clip(smeared_val, self.min_rw, self.max_rw)
                        pt_smeared_rw[n] *= smeared_val


                #two parameters, need to account for covariance
                elif(f.GetNpar() == 2):
                    p0 = f.GetParameter(0)
                    e0 = f.GetParError(0)

                    p1 = f.GetParameter(1)
                    e1 = f.GetParError(1)

                    cov_mat = np.array([[e0**2, covar], [covar, e1**2]])

                    #Use cholesky decomp to find diagonal basis 
                    #https://github.com/numpy/numpy/blob/main/numpy/random/_generator.pyx#L3930
                    l = np.linalg.cholesky(cov_mat)
                    p_sampled = np.array([p0,p1]) + pt_rand_noise[:,j-1,k-1,:2].dot(l.T)

                    #p_sampled = np.random.multivariate_normal([p0,p1], cov = cov_mat, size = pt_rand_noise.shape[0])
                    for n in range(pt_rand_noise.shape[0]):
                        pars = array('d', p_sampled[n])

                        smeared_val = f.EvalPar(array('d', [1./subjet_pt]), pars)
                        smeared_val = np.clip(smeared_val, self.min_rw, self.max_rw)
                        pt_smeared_rw[n] *= smeared_val
                else:
                    print("Pt extrap function order %i not implemented!" % f.GetNpar())
                    exit(1)





        return rw, smeared_rw, pt_smeared_rw


    def reweight(self, h_rw, lp_idxs, rw, smeared_rw, pt_smeared_rw, rand_noise = None):
        """Reweight based on directly measured data/MC LP ratio"""

        for (i,j,k) in lp_idxs:
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



    def get_up_down_prongs_weights(self, h_rw, reclust_prongs_up = None, reclust_prongs_down = None, nom_weight = None, do_symmetrize = False):
    
        prongs_up_weight = prongs_down_weight =  nom_weight

        #Compute weight for prong variations if needed
        #Separate prongs up from non-fully prongd decays and from bad matching
        if(reclust_prongs_up is not None): 
            prongs_up_weight, _, _ = self.reweight_lund_plane(h_rw = h_rw, reclust_obj = reclust_prongs_up)

        if(reclust_prongs_down is not None): 
            prongs_down_weight, _, _ = self.reweight_lund_plane(h_rw = h_rw, reclust_obj = reclust_prongs_down)

        return prongs_up_weight, prongs_down_weight

    def check_reclust_still_bad(self, reclust_prongs_up, reclust_prongs_down):
        still_bad = False
        if(reclust_prongs_up is not None):

            if( reclust_prongs_up.from_badmatch and (np.sum(reclust_prongs_up.subjet_match) != (reclust_prongs_up.n_prongs-1) or 
                np.sum(reclust_prongs_up.subjet_double_matched) != 0)):
                   still_bad= True

        if(reclust_prongs_down is not None):

            if( reclust_prongs_down.from_badmatch and (np.sum(reclust_prongs_down.subjet_match) != (reclust_prongs_down.n_prongs) or 
                np.sum(reclust_prongs_down.subjet_double_matched) != 0)):
                   still_bad = True

        return still_bad

    def init_weight_dict(self, nEvts, nToys):

        out = {
                'nom': np.zeros((nEvts)),
                'stat_vars': np.zeros((nEvts, nToys)),
                'pt_vars': np.zeros((nEvts, nToys)),
                'sys_up': np.zeros((nEvts)),
                'sys_down': np.zeros((nEvts)),
                'bquark_up': np.zeros((nEvts)),
                'bquark_down': np.zeros((nEvts)),
                'prongs_up': np.zeros((nEvts)),
                'prongs_down': np.zeros((nEvts)),
                'unclust_up': np.zeros((nEvts)),
                'unclust_down': np.zeros((nEvts)),
                'distortion_up': np.zeros((nEvts)),
                'distortion_down': np.zeros((nEvts)),
                'n_prongs': np.zeros((nEvts), dtype=np.int32),
                'subjet_pts': [],
                'bad_match': [False,]*nEvts,
                'reclust_still_bad_match': [False,]*nEvts,
                'reclust_nom': [],
                'reclust_prongs_up': [],
                'reclust_prongs_down': [],
        }
        return out


    def get_all_weights(self, pf_cands, gen_parts_eta_phi, ak8_jets, gen_parts_pdg_ids = None, do_sys_weights = True, distortion_sys = True, nToys = 100, pf_cands_PtEtaPhiE_format = False):
        """ Master function for the lund plane reweighting method. Takes in collection of events and computes nominal set of weights and variations from uncertainties
            All weights are normalized to average to one, so that the sample normalization is preserved
        Inputs:
        pf_cands : List of PF candidates (px,py,pz, E) for each event
        gen_parts_eta_phi : List of (eta, phi) of generator level quarks from the hard process (defines number of prongs for reclustering) for each event
        ak8_jets : 4 vector (pt, eta, phi, M) of AK8 jet for each event
        gen_parts_pdg_ids (optional): Pdg ids of the generator level quarks, used to check for the presence of b quarks which get a special uncertainty (assumes no b quarks if not given)
        do_sys_weights (optional): Compute all the systematic weight variations (default is true)
        distortion_sys (optional): Compute the systematic due to distortion of the Lund plane for this sample 
                    as compared to the sample the correction was derived from (stems from imperfections in reclustering)
                    Requires a sufficiently large sample events to get a good estimate, minimum is 1k but at least ~5k is recommended 
        nToys (optional): Number of toys to use for pt and stat variations of weights (default is 100)
        pf_cands_PtEtaPhiE_format (optional): Alternate representation of pf candidates (default is px,py,pz,E)


        Output:
            Dictionary with computed weights
            Nominal weights are under 'nom' and variations given under other keys. Naming is relatively explanatory
            Some additional info is also saved like number of the prongs for each event, and the reclustering object which saves the subjets and the splittings
        """


        nEvts = len(pf_cands)

        #dict for all the outputs
        out = self.init_weight_dict(nEvts, nToys)


        rand_noise = np.random.normal(size = (nToys, self.h_ratio.GetNbinsX(), self.h_ratio.GetNbinsY(), self.h_ratio.GetNbinsZ()))
        pt_rand_noise = np.random.normal(size = (nToys, self.h_ratio.GetNbinsY(), self.h_ratio.GetNbinsZ(), 3))

        unclust_factor = 5.0

        h_lp_signal = self.h_mc.Clone("h_lp_signal")
        h_lp_signal.Reset()

        if(distortion_sys and len(pf_cands) < 1000):
            print("Only %i jets given, will not include LP distortion systematic" % len(pf_cands))
            distortion_sys = False


        n_badmatch = 0
        #Recluster into subjets
        for i,cands in enumerate(pf_cands):
            reclust_nom, reclust_prongs_up, reclust_prongs_down = self.get_splittings_and_matching(cands, gen_parts_eta_phi[i], ak8_jets[i], pf_cands_PtEtaPhiE_format = pf_cands_PtEtaPhiE_format)

            out['bad_match'][i] = reclust_nom.badmatch
            out['n_prongs'][i] = reclust_nom.n_prongs

            out['reclust_nom'].append(reclust_nom)
            out['reclust_prongs_up'].append(reclust_prongs_up)
            out['reclust_prongs_down'].append(reclust_prongs_down)

            if(distortion_sys and not out['bad_match'][i]): self.fill_lund_plane(h_lp_signal, reclust_obj = reclust_nom)
            #print(i, reclust_nom.n_prongs, reclust_nom.subjet, reclust_nom.split[:5])
            #exit(1)


        #compute how much signal LP differs from W's MC used to derive correction
        if(distortion_sys):
            h_dummy = self.h_mc.Clone("h_dummy")
            h_dummy.Reset()
            #h_distortion_ratio = self.make_LP_ratio(self.h_mc, h_dummy, h_lp_signal)
            h_distortion_ratio = self.make_LP_ratio(self.h_mc, h_dummy, h_lp_signal)

        for i in range(len(out['reclust_nom'])):
            reclust_nom, reclust_prongs_up, reclust_prongs_down = out['reclust_nom'][i], out['reclust_prongs_up'][i], out['reclust_prongs_down'][i]


            #Gets the nominal LP reweighting factor for this event and statistical + pt extrapolation toys
            out['nom'][i], out['stat_vars'][i], out['pt_vars'][i] = self.reweight_lund_plane(h_rw = self.h_ratio,
                                    reclust_obj = reclust_nom, rand_noise = rand_noise, pt_rand_noise = pt_rand_noise, )


            
            out['prongs_up'][i], out['prongs_down'][i]  = self.get_up_down_prongs_weights(h_rw = self.h_ratio, 
                    reclust_prongs_up = reclust_prongs_up, reclust_prongs_down = reclust_prongs_down, nom_weight = out['nom'][i])

            #quarks which are still not matched despite varying prongs up/down
            out['reclust_still_bad_match'][i] = self.check_reclust_still_bad(reclust_prongs_up, reclust_prongs_down)
            out['unclust_up'][i] = out['nom'][i]
            out['unclust_down'][i] = out['nom'][i]

            #un-matched quarks are not calibrated by the procedure, varying their weight up/down by conservative factor
            if(out['reclust_still_bad_match'][i]):
                out['unclust_up'][i] *= unclust_factor
                out['unclust_down'][i] /= unclust_factor



            for sj in reclust_nom.subjet: out['subjet_pts'].append(sj[0])

            #compute systematic due to distorted LP 
            if(distortion_sys):
                distortion_weight,_,_ = self.reweight_lund_plane(h_rw = h_distortion_ratio, reclust_obj = reclust_nom, sys_str = 'distortion')

                out['distortion_up'][i] = out['nom'][i] * distortion_weight
                out['distortion_down'][i] = out['nom'][i] / distortion_weight


            if(do_sys_weights):
                #Now get systematic variations due to systemtatic uncertainties on LP
                out['sys_up'][i],_,_ = self.reweight_lund_plane(h_rw = self.h_ratio_sys_up, reclust_obj = reclust_nom)
                out['sys_down'][i],_,_ = self.reweight_lund_plane(h_rw = self.h_ratio_sys_down, reclust_obj = reclust_nom)

                #compute special systematic for subjets matched to b quarks
                #not needed if signal does not specifically produce b quark subjets
                if(gen_parts_pdg_ids is None): gen_bs = []
                else: gen_bs = [j for j in range(len(gen_parts_pdg_ids[i])) if abs(gen_parts_pdg_ids[i][j]) == B_PDG_ID]

                if(len(gen_bs) == 0): b_rw = 1.0
                else:
                    eta_phi = np.array(gen_parts_eta_phi[i]) #this event only (ignore shape irregularities)
                    dists = get_subjet_dist(eta_phi[gen_bs,:], np.array(reclust_nom.subjet)[:,1:3])

                    deltaR_cut = 0.2
                    b_matches = []
                    #which subjet is each quark closest to
                    j_closest = np.amin(dists, axis = 0)
                    j_which = np.argmin(dists, axis = 0)
                    b_matches = np.unique(j_which[j_closest < deltaR_cut])

                    #reweight only subjets matched to b quarks
                    if(len(b_matches) > 0):
                        b_subjet = [reclust_nom.subjet[j] for j in range(len(reclust_nom.subjet)) if j in b_matches]
                        b_split  = [reclust_nom.split[j]  for j in range(len(reclust_nom.split)) if reclust_nom.split[j][0] in b_matches]

                        b_rw, _,_   = self.reweight_lund_plane(self.b_light_ratio, subjets = b_subjet, splittings = b_split, sys_str = 'bquark')

                    else: b_rw = 1.0
                out['bquark_up'][i] = out['nom'][i]* b_rw
                out['bquark_down'][i] = out['nom'][i]/b_rw


        for key in out.keys(): 
            if(('nom' in key) or ('up' in key) or ('down' in key) or ('vars' in key)):
                if(isinstance(out[key], np.ndarray)): out[key] = self.normalize_weights(out[key], n_prongs = out['n_prongs'])


        return out





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

                if(self.pt_extrap_dir is None or subjets[i][0] < self.pt_extrap_val or ('bquark' in sys_str) or ('distortion' in sys_str)):
                    rw, smeared_rw, pt_smeared_rw = self.reweight(h_rw, lp_idxs, rw, smeared_rw, pt_smeared_rw, rand_noise = rand_noise)
                else:
                    rw, smeared_rw, pt_smeared_rw = self.reweight_pt_extrap(subjets[i][0], lp_idxs, rw, smeared_rw, pt_smeared_rw, pt_rand_noise = pt_rand_noise, sys_str = sys_str)

        
        return rw, smeared_rw, pt_smeared_rw


    def make_LP_ratio(self, h_data, h_bkg, h_mc,  h_data_subjet_pt = None, h_bkg_subjet_pt = None, h_mc_subjet_pt = None, pt_bins = None, outdir = "", save_plots = False):
        """ Function to construct data/MC LP ratio"""


        #h_data.Print()
        #h_bkg.Print()
        #h_mc.Print()

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

            #h_data_subjet_pt.Print()
            #h_bkg_subjet_pt.Print()
            #h_mc_subjet_pt.Print()

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
        #h_data_sub.Print()



        cleanup_hist(h_data_sub)

        if (pt_bins is None):
            pt_bins = np.arange(0, h_data.GetNbinsX()+1)


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



            if(bkg_norm > 0): h_bkg_proj.Scale(1./bkg_norm)

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

                h_mc_proj.SetTitle("TTbar MC pT %.0f - %.0f" % (pt_bins[i-1], pt_bins[i]))
                h_data_proj.SetTitle("Data - Bkg pT %.0f - %.0f (N = %.0f)" % (pt_bins[i-1], pt_bins[i], data_norm))
                h_ratio_proj.SetTitle("Ratio pT %.0f - %.0f (N = %.0f)" % (pt_bins[i-1], pt_bins[i], data_norm))


                c_mc = ROOT.TCanvas("c", "", 1000, 1000)
                h_mc_proj.Draw("colz")
                c_mc.SetRightMargin(0.2)
                c_mc.Print(outdir + "lundPlane_bin%i_MC.png" % i)


                if(bkg_norm > 0):
                    h_bkg_proj.SetTitle("Bkg MC pT %.0f - %.0f" % (pt_bins[i-1], pt_bins[i]))
                    c_bkg = ROOT.TCanvas("c", "", 1000,1000)
                    h_bkg_proj.Draw("colz")
                    c_bkg.SetRightMargin(0.2)
                    c_bkg.Print(outdir + "lundPlane_bin%i_bkg.png" % i)

                c_data = ROOT.TCanvas("c", "", 1000, 1000)
                h_data_proj.Draw("colz")
                c_data.SetRightMargin(0.2)
                c_data.Print(outdir + "lundPlane_bin%i_data.png" %i )



                c_ratio = ROOT.TCanvas("c", "", 1000, 1000)
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

    def normalize_weights(self, lund_weights, w_min = 0.1, w_max = 10., n_prongs = None):
        """Normalize lund plane weights so average weight is 1 (necessary to preserve normalization of MC.)
        Done separately for jets of different number of prongs, so not biased
        Also clip outlier weights so to not be dominated by statistical fluctuations. """

        #separate norm per each number of prongs (so dist is not biased)
        if(n_prongs is None): n_prongs = np.ones_like(lund_weights, dtype=np.int32)
        max_prongs = int(round(np.amax(n_prongs)))

        for n in range(1, max_prongs+1):
            mask = (n_prongs == n)
            weights = lund_weights[mask]

            weights = np.clip(weights, 0., w_max)
            if(len(weights.shape) > 1): weights /= np.mean(weights, axis = 0, keepdims=True)
            else: weights /= np.mean(weights)

            weights  = np.clip(weights, w_min, w_max)
            if(len(weights.shape) > 1): weights /= np.mean(weights, axis = 0, keepdims=True)
            else: weights /= np.mean(weights)

            lund_weights[mask] = weights

        return lund_weights

class ReclusterObj():
    def __init__(self,subjets = None, split = None, dRs = None):
        self.subjet = None
        self.split = None
        self.dRs = None
        self.subjet_match = None #is subjet matched to a quark
        self.subjet_double_matched = None #is subjet matched to multiple quarks
        self.subjet_dRs = None #deltaRs between each quark and closest subjet
        self.from_badmatch = False
        self.from_prongs_up = False
        self.from_prongs_down = False
        self.n_prongs = 0





def matched(c,cj):
    eps = 1e-4
    return (abs(c[0] - cj.px()) < eps) and (abs(c[1] - cj.py()) < eps) and (abs(c[2] - cj.pz()) < eps)


def find_matching_pf(cj_list, cj):

    for c in cj_list:
        if(matched(c, cj)): 
            return c
    return None



