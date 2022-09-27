import h5py
import fastjet as fj
import numpy as np
import sys
from PlotUtils import *
import ROOT
from array import array
import copy

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(False)


        


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



def convert_4vec(vec):
    rvec = ROOT.Math.PtEtaPhiMVector(vec[0], vec[1], vec[2], vec[3])
    return [rvec.Px(), rvec.Py(), rvec.Pz(), rvec.E()]

def get_splittings(pf_cands, jetR = -1, boost_vec = None, maxJets = -1, num_excjets = -1):
    pjs = []
    for c in pf_cands:
        if(c[3] > 0.0001):
            pj = fj.PseudoJet(c[0], c[1], c[2], c[3])
            if(boost_vec is not None):
                pj_boost = pj.unboost(boost_vec)
                #idk why you have to do this...
                pjp = fj.PseudoJet(pj_boost.px(), pj_boost.py(), pj_boost.pz(), pj_boost.E())
                pjs.append(pjp)
            else:
                pjs.append(pj)

    if(jetR < 0): R = 1000.0
    else: R = jetR
    jet_def = fj.JetDefinition(fj.cambridge_algorithm, R)
    cs = fj.ClusterSequence(pjs, jet_def)
    if(num_excjets < 0):
        js = fj.sorted_by_pt(cs.inclusive_jets())
        if(boost_vec is None): #only do first vec
            js = [js[0]]
        elif(maxJets > 0):
            nMax = min(len(js), maxJets)
            js = js[:nMax]
    else:
        js = fj.sorted_by_pt(cs.exclusive_jets_up_to(num_excjets))

    subjets = []
    splittings = []
    #print("%i subjets " % len(js))
    for i, j in enumerate(js):
        #print("sj %i" % i)
        pseudojet = j
        jet_pt = j.pt()
        subjets.append([j.pt(), j.eta(), j.phi(), j.m()])
        while True:
            j1 = fj.PseudoJet()
            j2 = fj.PseudoJet()
            if pseudojet and pseudojet.has_parents(j1, j2):
                # order the parents in pt
                if (j2.pt() > j1.pt()):
                    j1, j2 = j2, j1
                # check if we satisfy cuts
                delta = j1.delta_R(j2)
                kt = j2.pt() * delta
                splittings.append([i, delta, kt])
                pseudojet = j1
            else:
                break
    return subjets, splittings


def fill_lund_plane(h, pf_cands = None, subjets = None,  splittings = None, fill_z = True, jetR = -1, dR = 0.8, boost_vec = None, maxJets = -1, num_excjets = -1, pt_min = 0., weight = 1.):

    if(subjets is None or splittings is None):
        subjets, splittings = get_splittings(pf_cands, jetR = jetR, boost_vec = boost_vec, maxJets = maxJets, num_excjets = num_excjets)
        subjets = np.array(subjets).reshape(-1)

    for jet_i, delta, kt in splittings:
        jet_pt = subjets[jet_i*4]
        if(delta > 0. and kt > 0.):
            if(fill_z): 
                print("FillZ DEPRECATED")
                break
            else: 
                if(type(h) == ROOT.TH3F): h.Fill(jet_pt, np.log(dR/delta), np.log(kt), weight)
                else: h.Fill(np.log(dR/delta), np.log(kt), weight)
    return



def reweight_lund_plane(h_rw, pf_cands, dR = 0.8, fill_z = True, jetR = -1, boost_vec = None, maxJets = -1, num_excjets = -1, pt_min = 0., weight = 1., uncs = False):
    pjs = []

    h_jet = h_rw.Clone("temp")
    h_jet.Reset()
    fill_lund_plane(h_jet, pf_cands, dR, fill_z = fill_z, jetR = jetR, boost_vec = boost_vec, maxJets = maxJets, num_excjets = num_excjets, weight = weight)
    #h_jet.Scale(1./h_jet.Integral())
    #h_jet.Multiply(h_rw)


    rw = 1.0
    unc = 0.0
    eps = 1e-4
    if(type(h_rw) == ROOT.TH3F):
        for i in range(1, h_jet.GetNbinsX() + 1):
            for j in range(1, h_jet.GetNbinsY() + 1):
                for k in range(1, h_jet.GetNbinsZ() + 1):
                    n_cands = h_jet.GetBinContent(i,j,k)
                    if(n_cands > 0): 
                        #print("Rw %.3f, cont %.3f, i %i j %i k %i n %i" % (rw, h_rw.GetBinContent(i,j,k), i,j,k, n_cands))
                        val = h_rw.GetBinContent(i,j,k)
                        if(uncs): 
                            err = h_rw.GetBinError(i,j,k)
                            #uncertainty propagation
                            unc = ( (n_cands * rw * val**(n_cands -1) * err)**2 + (val**n_cands * unc)**2) ** (0.5)
                        rw *= val ** n_cands

    else:
        for i in range(1, h_jet.GetNbinsX() + 1):
            for j in range(1, h_jet.GetNbinsY() + 1):
                n_cands = h_jet.GetBinContent(i,j)
                if(n_cands > 0): 
                    val = h_rw.GetBinContent(i,j)
                    if(uncs): 
                        err = h_rw.GetBinError(i,j)
                        #uncertainty propagation
                        unc = ( (n_cands * rw * val**(n_cands -1) * err)**2 + (val**n_cands * unc)**2) ** (0.5)
                    rw *= val ** n_cands

    #h_jet.Print("range")
    
    return rw, unc

#No longer used
def reweight_subjet_lund_plane(h_rw, pf_cands, boost_vec, dR = 0.4, fill_z = True, jetR = 0.4):
    pjs = []

    h_jet = h_rw.Clone("temp")
    h_jet.Reset()


    #fill_lund_plane(h_jet, pf_cands, dR, fill_z)

    pjs = []

    for c in pf_cands:
        if(c[3] > 0.0001):
            pj = fj.PseudoJet(c[0], c[1], c[2], c[3])
            pj_boost = pj.unboost(boost_vec)
            #idk why you have to do this...
            pjp = fj.PseudoJet(pj_boost.px(), pj_boost.py(), pj_boost.pz(), pj_boost.E())
            pjs.append(pjp)



    jetR = 0.4
    jet_def = fj.JetDefinition(fj.cambridge_algorithm, jetR)
    cs = fj.ClusterSequence(pjs, jet_def)

    jets = cs.inclusive_jets()

    for i,j in enumerate(jets):
        pseudojet = j
        while True:
            j1 = fj.PseudoJet()
            j2 = fj.PseudoJet()
            if pseudojet and pseudojet.has_parents(j1, j2):
                # order the parents in pt
                if (j2.pt() > j1.pt()):
                    j1, j2 = j2, j1
                delta = j1.delta_R(j2)
                kt = j2.pt() * delta
                z = j2.pt() / (j1.pt() + j2.pt())
                if(fill_z): h.Fill(np.log(dR/delta), np.log(1./z))
                else: h_jet.Fill(np.log(dR/delta), np.log(kt))
                pseudojet = j1
            else:
                break

    rw = 1.0
    eps = 1e-4
    for i in range(1, h_jet.GetNbinsX() + 1):
        for j in range(1, h_jet.GetNbinsY() + 1):
            n_cands = h_jet.GetBinContent(i,j)
            if(n_cands > 0): rw *= h_rw.GetBinContent(i,j)

    return rw



class Dataset():
    def __init__(self, f, is_data = False, label = "", color = ""):

        self.f = f
        self.is_data = is_data

        self.label = label
        self.color = color

        self.n_evts = f['event_info'].shape[-1]
        self.mask = f['jet_kinematics'][:,0] > 0.
        self.norm_factor = 1.0

    def n(self):
        return np.sum(self.mask)

    def apply_cut(self, cut):
        self.mask = self.mask & cut
    
    def get_masked(self, key):
        return self.f[key][()][self.mask]

    def get_weights(self):
        if(self.is_data): return np.ones(self.n())
        else: return self.get_masked('norm_weights') * self.norm_factor

    def compute_obs(self):
        eps = 1e-8
        feats = self.get_masked('jet1_extraInfo')
        kins = self.get_masked('jet_kinematics')
        self.tau21 = (feats[:,1] / (feats[:,0] + eps))
        self.tau32 = (feats[:,2] / (feats[:,1] + eps))
        self.tau43 = (feats[:,3] / (feats[:,2] + eps))
        self.mSoftDrop = kins[:,3]
        self.pt = kins[:,0]
        self.nPF= feats[:,6]



    def fill_LP(self, h, subjet_rw = False, fill_z = False, jetR = 0.8, num_excjets = 2):

        pf_cands = self.get_masked("jet1_PFCands").astype(np.float64)
        jet_kinematics = self.get_masked("jet_kinematics")
        weights = self.get_weights()
        for i,pf_cand in enumerate(pf_cands):
            weight =weights[i]
            if(subjet_rw):
                pt_eta_phi_m_vec = jet_kinematics[i]
                jet_4vec = convert_4vec(pt_eta_phi_m_vec)
                boost_vec = fj.PseudoJet(jet_4vec[0], jet_4vec[1], jet_4vec[2], jet_4vec[3])
                fill_lund_plane(h, pf_cand,  boost_vec = boost_vec, fill_z =fill_z, dR = jetR, jetR = jetR, weight = weight)
            else: fill_lund_plane(h, pf_cand, fill_z = fill_z, jetR = jetR, num_excjets = num_excjets, weight = weight)

    def reweight_LP(self, h_ratio, subjet_rw = False, fill_z = False, jetR = 0.8, num_excjets = 2, uncs = True, max_evts =-1):
    #always uncs for now

        LP_weights = []
        LP_uncs = []
        pf_cands = self.get_masked("jet1_PFCands").astype(np.float64)
        jet_kinematics = self.get_masked("jet_kinematics")

        if(max_evts > 0 and pf_cands.shape[0] > max_evts):
            pf_cands = pf_cands[:max_evts]
            jet_kinematics = jet_kinematics[:max_evts]

        for i,pf_cand in enumerate(pf_cands):
            if(subjet_rw):
                pt_eta_phi_m_vec = jet_kinematics[i]
                jet_4vec = convert_4vec(pt_eta_phi_m_vec)
                boost_vec = fj.PseudoJet(jet_4vec[0], jet_4vec[1], jet_4vec[2], jet_4vec[3])
                rw, unc = reweight_lund_plane(h_ratio, pf_cand,  boost_vec = boost_vec, fill_z =fill_z, dR = jetR, jetR = jetR, uncs = True)
            else: 
                rw, unc = reweight_lund_plane(h_ratio, pf_cand, fill_z = fill_z, jetR = jetR, num_excjets = num_excjets, uncs = True)

            LP_weights.append(rw)
            if(rw >= 1e-6):
                LP_uncs.append(unc/rw)
            else:
                LP_uncs.append(0.)


        LP_weights = np.clip(np.array(LP_weights), 0., 10.)
        LP_uncs = np.clip(np.array(LP_uncs), 0., 1.5)
        mean = np.mean(LP_weights)
        LP_weights /= mean
        return LP_weights, LP_uncs
