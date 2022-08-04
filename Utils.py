import h5py
import fastjet as fj
import numpy as np
import sys
from PlotUtils import *
import ROOT
from array import array

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(False)

def convert_4vec(vec):
    rvec = ROOT.Math.PtEtaPhiMVector(vec[0], vec[1], vec[2], vec[3])
    return [rvec.Px(), rvec.Py(), rvec.Pz(), rvec.E()]


def fill_lund_plane(h, pf_cands, dR = 0.8,  fill_z = True, jetR = -1, boost_vec = None, maxJets = -1, num_excjets = -1, pt_min = 0., weight = 1.):
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


    for i, j in enumerate(js):
        pseudojet = j
        jet_pt = j.pt()
        #print(i, j.pt())
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
                z = j2.pt() / (j1.pt() + j2.pt())
                if(delta > 0. and z > 0.):
                    if(fill_z): 
                        if(type(h) == ROOT.TH3F): h.Fill(jet_pt, np.log(dR/delta), np.log(1./z), weight)
                        else: h.Fill(np.log(dR/delta), np.log(1./z), weight)
                    else: 
                        if(type(h) == ROOT.TH3F): h.Fill(jet_pt, np.log(dR/delta), np.log(kt), weight)
                        else: h.Fill(np.log(dR/delta), np.log(kt), weight)
                pseudojet = j1
            else:
                break
    return



def reweight_lund_plane(h_rw, pf_cands, dR = 0.8, fill_z = True, jetR = -1, boost_vec = None, maxJets = -1, num_excjets = -1, pt_min = 0., weight = 1.):
    pjs = []

    h_jet = h_rw.Clone("temp")
    h_jet.Reset()
    fill_lund_plane(h_jet, pf_cands, dR, fill_z = fill_z, jetR = jetR, boost_vec = boost_vec, maxJets = maxJets, num_excjets = num_excjets, weight = weight)
    #h_jet.Scale(1./h_jet.Integral())
    #h_jet.Multiply(h_rw)

    rw = 1.0
    eps = 1e-4
    if(type(h_rw) == ROOT.TH3F):
        for i in range(1, h_jet.GetNbinsX() + 1):
            for j in range(1, h_jet.GetNbinsY() + 1):
                for k in range(1, h_jet.GetNbinsZ() + 1):
                    n_cands = h_jet.GetBinContent(i,j,k)
                    if(n_cands > 0): 
                        #print("Rw %.3f, cont %.3f, i %i j %i k %i n %i" % (rw, h_rw.GetBinContent(i,j,k), i,j,k, n_cands))
                        rw *= h_rw.GetBinContent(i,j,k)**n_cands

    else:
        for i in range(1, h_jet.GetNbinsX() + 1):
            for j in range(1, h_jet.GetNbinsY() + 1):
                n_cands = h_jet.GetBinContent(i,j)
                if(n_cands > 0): rw *= h_rw.GetBinContent(i,j)**n_cands

    #h_jet.Print("range")
    return rw

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
