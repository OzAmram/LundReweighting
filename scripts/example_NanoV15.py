#!/usr/bin/env python3
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import *

import sys, os

sys.path.insert(0, "")
sys.path.append("../")
from utils.Utils import *


""" 
An example how to use the Lund Plane reweighting code on a signal file in NanoV15 
Computes reweighting factors for an example signal. 
It then uses them to compute the efficiency and uncertainty of a given substructure cut

This example uses NanoAODTools to read the relevant branches from the Nano V15
This tool is included in CMSSW, but can be used standalone as well : 
 https://github.com/cms-sw/cmssw/tree/master/PhysicsTools/NanoAODTools

Note that the correction itself does not use NanoAODTools, the user is welcome to use whatever 
framework they like for reading the values from NanoAOD.

"""


def get_ttbar_gen_parts(event, verbose=True):
    """
    Function which finds parses the gen particles to find the correct gen-level quarks which define the prongs.
    Here this is done for ttbar
    Note that this function will have be custom for every signal to correctly identify the gen quarks
    """

    def isFinal(genPart):
        # utility function
        # check if isLastCopy flag is set (Pythia)
        mask = 1 << 13  # 13th bit of status flag
        return (genPart.statusFlags & mask) != 0

    # pdg IDs
    top_ID = 6
    H_ID = 25
    W_ID = 24
    Z_ID = 23
    B_ID = 5
    MAXLEP_ID = 16
    MAXLIGHTQUARK_ID = 5

    GenPartsColl = Collection(event, "GenPart")

    top = anti_top = W = anti_W = fermion1 = anti_fermion1 = b_quark1 = fermion2 = (
        anti_fermion2
    ) = b_quark2 = None

    # Find top and W's
    for genPart in GenPartsColl:
        # tops
        if abs(genPart.pdgId) == top_ID and isFinal(genPart):
            if genPart.pdgId > 0:
                if top is None:
                    top = genPart
                else:
                    print("WARNING : Extra top ? ")
            else:
                if anti_top is None:
                    anti_top = genPart
                else:
                    print("WARNING : Extra antitop ? ")
        m = genPart.genPartIdxMother
        # W's
        if abs(genPart.pdgId) == W_ID and isFinal(genPart):
            if genPart.pdgId > 0:
                if W is None:
                    W = genPart
                else:
                    print("WARNING : Extra W ? ")
            else:
                if anti_W is None:
                    anti_W = genPart
                else:
                    print("WARNING : Extra anti W ? ")

    if top is None or anti_top is None or W is None or anti_W is None:
        if verbose:
            print("Couldnt find top or W: ")
            print(top, anti_top, W, anti_W)
            # count = 0
            for genPart in GenPartsColl:
                print(count, genPart.pdgId, genPart.pt, genPart.genPartIdxMother)
                count += 1

    for genPart in GenPartsColl:
        # Find quarks or leptons from W decay
        m = genPart.genPartIdxMother
        mother = GenPartsColl[m] if m > 0 else None
        w_mother_match = mother is W
        anti_w_mother_match = mother is anti_W
        if abs(genPart.pdgId) <= MAXLEP_ID and m > 0 and w_mother_match:
            if genPart.pdgId > 0:
                if fermion1 is None:
                    fermion1 = genPart
                elif verbose:
                    print("WARNING : Extra quark ? ")
            else:
                if anti_fermion1 is None:
                    anti_fermion1 = genPart
                elif verbose:
                    print("WARNING : Extra anti quark ? ")

        elif abs(genPart.pdgId) <= MAXLEP_ID and m > 0 and anti_w_mother_match:
            if genPart.pdgId > 0:
                if fermion2 is None:
                    fermion2 = genPart
                elif verbose:
                    print("WARNING : Extra quark ? ")
            else:
                if anti_fermion2 is None:
                    anti_fermion2 = genPart
                elif verbose:
                    print("WARNING : Extra anti quark ? ")

        # find b quark from top
        top_mother_match = mother is top
        anti_top_mother_match = mother is anti_top
        if abs(genPart.pdgId) == B_ID and top_mother_match:
            if b_quark1 is None:
                b_quark1 = genPart
            elif verbose:
                print("WARNING : Extra quark ? ")

        elif abs(genPart.pdgId) == B_ID and anti_top_mother_match:
            if b_quark2 is None:
                b_quark2 = genPart
            elif verbose:
                print("WARNING : Extra quark ? ")

    return (
        top,
        anti_top,
        W,
        anti_W,
        fermion1,
        anti_fermion1,
        b_quark1,
        fermion2,
        anti_fermion2,
        b_quark2,
    )


# utility function
def ang_dist(phi1, phi2):
    dphi = phi1 - phi2
    if dphi < -math.pi:
        dphi += 2.0 * math.pi
    if dphi > math.pi:
        dphi -= 2.0 * math.pi
    return dphi


# utility function
def deltaR(o1, o2):
    return ((o1.eta - o2.eta) ** 2 + ang_dist(o1.phi, o2.phi) ** 2) ** (0.5)


def check_matching(jet, f1, f2, b_quark):
    # check if quarks are inside ak8 jet
    # 0 = no matching, 1 = W_matched, 2 = top_matched

    f1_in = f1 is not None and deltaR(jet, f1) < 0.8
    f2_in = f2 is not None and deltaR(jet, f2) < 0.8
    b_in = b_quark is not None and deltaR(jet, b_quark) < 0.8

    return (f1_in + f2_in + b_in) > 1

def get_inputs(inputFile, max_events=5000):
    """Function to read the NanoV15 file of the signal and get the necessary ingredients
    for the the correction (the jets, then gen-quarks forming the prongs, and the PF Candidates inside the jet.
    This is just an example, the user is welcome to read the Nano with whatever tools/framework they like
    """

    triggers = [
        "HLT_PFHT890",
        "HLT_PFHT1050",
        "HLT_PFJet450",
        "HLT_PFJet500",
    ]

    # Loop through the Nano file to select jets passing our preselection
    TTree = inputFile.Get("Events")

    nTotal = TTree.GetEntries()
    inTree = InputTree(TTree)

    print("Running over %i entries \n" % nTotal)

    nEvents = 0
    jets = []
    quarks = []
    cands = []

    for entry in range(inTree.entries):
        # Grab the event

        if entry % 10000 == 0:
            print(
                "--------- Processing Event %i. %i jets saved so far" % (entry, nEvents)
            )

        event = Event(inTree, entry)

        PFCands = Collection(event, "PFCand")
        AK8Jets = Collection(event, "FatJet")
        FatJetPFCands = Collection(
            event, "FatJetPFCand"
        )  # association between jets and PFCands

        # You would apply your basic preselection here
        # For this simple example we will just apply a trigger
        passTrigger = False
        for trig in triggers:
            passTrigger = passTrigger or inTree.readBranch(trig)
        if not passTrigger:
            continue



        # Select the leading jet with pt > 400 
        # We also require our fatjet to be a genuine mutli-prong jet, checking that it contains two or more gen level quarks from the top decay inside it
        # Note you don't have to require a truth-matching like this in your analysis, we just use it here for a purer sample of boosted top/W's
        jet_min_pt = 400
        my_jet = None

        if(len(AK8Jets) == 0 or AK8Jets[0].pt < jet_min_pt): 
            continue

        #Get gen-level quarks which define our prongs
        top, anti_top, W, anti_W, q1a, q1b, b1, q2a, q2b, b2 = (
            get_ttbar_gen_parts(event)
        )


        jet_index = 0
        for i, jet in enumerate(AK8Jets):
            jet.idx = i
            if (
                jet.pt > jet_min_pt
                and abs(jet.eta) < 2.4
                and (check_matching(jet, q1a, q1b, b1) or
                     check_matching(jet, q2a, q2b, b2))
            ):
                my_jet = jet
                break

        if my_jet is None:
            continue
        nEvents += 1

        # save our jet

        # Only need the 4 vector but we will save some tagging info as well
        eps = 1e-6
        tau21 = my_jet.tau2 / (my_jet.tau1 + eps)
        jets.append(
            [
                my_jet.pt,
                my_jet.eta,
                my_jet.phi,
                my_jet.msoftdrop,
                tau21,
                my_jet.particleNetLegacy_Xqq,
            ]
        )

        # Save the quarks from the gen-level process which will define our prongs
        # Save all 6 quarks, reweighting code will check how many actuall end up inside the AK8 jet
        quark_list = [
            [q1a.pt, q1a.eta, q1a.phi, q1a.pdgId],
            [q1b.pt, q1b.eta, q1b.phi, q1b.pdgId],
            [b1.pt, b1.eta, b1.phi, b1.pdgId],
            [q2a.pt, q2a.eta, q2a.phi, q2a.pdgId],
            [q2b.pt, q2b.eta, q2b.phi, q2b.pdgId],
            [b2.pt, b2.eta, b2.phi, b2.pdgId],
        ]
        quarks.append(quark_list)

        # Save the PF Candidates
        jet_PFCands = []
        for cand_map in FatJetPFCands:
            if cand_map.jetIdx == my_jet.idx:
                idx = cand_map.pfCandIdx
                cand = ROOT.Math.PtEtaPhiMVector(
                    PFCands[idx].pt,
                    PFCands[idx].eta,
                    PFCands[idx].phi,
                    PFCands[idx].mass,
                )
                jet_PFCands.append([cand.Px(), cand.Py(), cand.Pz(), cand.E()])

        cands.append(jet_PFCands)

        if nEvents >= max_events:
            break

    jets = np.array(jets)
    quarks = np.array(quarks)
    return jets, quarks, cands

# ----------- Setup --------------

# Input NanoV15 file  for signal, here we use a ttbar file for example
# NOTE you may want to change the file redirector depending on your region
# Alternatively you can use xrootd to copy the file locally which may be faster (its ~1.5 GB)
fname = "root://cms-xrd-global.cern.ch//store/mc/RunIII2024Summer24NanoAODv15/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/150X_mcRun3_2024_realistic_v1-v2/100000/013c4b44-92a8-42a8-ac27-c127158c4726.root"
#fname = "TTbar_Nano_test.root"

# File with the correction ingredients (centrally provided)
# NOTE This should usually match the corresponding era of your MC sample
# Run3 corrections are WIP, for now just use 2018 for example purposes
f_ratio_name = "data/ratio_2018.root"

print("Opening signal file %s" % fname)
print("Opening Lund Plane Correction file %s" % f_ratio_name)
f_sig = ROOT.TFile.Open(fname)
f_ratio = ROOT.TFile.Open(f_ratio_name)

# For this example just run over small sample
# Generally one should run over the entire signal MC
max_evts = 1000

# Get the needed inputs for the reweighting
print("Getting inputs from Nano")
jets, quarks, cands = get_inputs(f_sig, max_evts)

# 4 vector of AK8 jets we are reweighting
ak8_jets = jets[:, :4]


# eta phi of gen-level quarks from the hard process which define our prongs
gen_parts_eta_phi = quarks[:, :, 1:3]

# pdgID of the quarks
gen_parts_pdg_id = quarks[:, :, 3]

# List of PFCands for each jet
# (px, py, pz, E) format
pf_cands = cands

# Initialize reweighting tool
LP_rw = LundReweighter(f_ratio=f_ratio)

# Nominal event weights of the MC, (assume every event is weight '1' for this example)
nom_weights = np.ones(len(pf_cands))


# ----------- Do the reweighting ---------------------
# Use the tool to compute the weights

# The 'get_all_weights' 'master' function, it computes all the weights for you and the systematic variations
# Note that it takes the whole collection of jets in one go, it cannot be used on one jet at at time
# This is because the normalization of the weights needs to be computed on the whole sample
print("Computing reweighting factors")
LP_weights = LP_rw.get_all_weights(
    pf_cands, gen_parts_eta_phi, ak8_jets, gen_parts_pdg_ids=gen_parts_pdg_id
)
print("Done!")


# multiply Lund plane reweighting factor with nominal event weights to get total weight
for key in LP_weights.keys():
    if "nom" in key or "up" in key or "down" in key:
        if isinstance(LP_weights[key], np.ndarray):
            LP_weights[key] *= nom_weights

# The LP_weights object is a dictionary containing the baseline Lund weights ('nom')
# As well as the weights for many systematic variations
# And some diagnostic information

print("Keys inside Lund reweighting object: ")
print(list(LP_weights.keys()))

# The variations with 'up' and 'down' in their name define up/down systematics on the procedure and should be used in the normal way
# The statistical uncertainty and the pt extrapolation uncertainty are defined using 100 toy variations of the weights
# Later in the script we will see how to use them to get the final uncertainty

# Print out the fraction of prongs that are not well matched to our reclustered subjets (want this to be low)
print("Bad match frac %.2f" % np.mean(LP_weights["bad_match"]))
# Fraction of prongs that are still not well matched after reclustering with varied number of prongs
print(
    "Reclustered bad match frac %.2f" % np.mean(LP_weights["reclust_still_bad_match"])
)

# Thats it! If you use these new event weights and the systematic variations
# then you will be applying the correction and its uncertainties
# You may want to save them in whatever data format you use for your analysis


f_ratio.Close()
f_sig.Close()

# -------- Use weights to compute efficiency/SF of a cut -------

# This part shows an example of how these weights can be used to compute a 'scale factor' on the efficiency
# of a substrcture cut


# First we define the substructure cut we want to correct

tau21 = jets[:, 4]
tau21_cut = 0.4
score_cut = tau21 < tau21_cut

# The 'scale factor' of the Lund Plane reweighting is defined as the ratio
# between the efficiency of the cut in the nominal MC to the efficiency of the cut
# with the Lund weights applied

# Efficiency of the cut in nominal MC (without the Lund weights)
eff_nom = np.average(score_cut, weights=nom_weights)

# Efficiency of the cut after the Lund Plane reweighting
eff_rw = np.average(score_cut, weights=LP_weights["nom"])

# Nominal 'scale factor'
SF = eff_rw / eff_nom

print(
    "Nominal efficiency %.3f, Corrected efficiency %.3f, SF (corrected / nom) %.3f"
    % (eff_nom, eff_rw, SF)
)

# NOTE, because there is kinematic dependence to the correction, it is better to use corrected efficiency computed
# separately for each MC sample rather than a single 'SF'


# ----- Compute uncertainties on the efficiency from the various weight variations --------

# Statistical and pt extrapolation uncertainties are derived from 100 variations of the weights
# Take std dev of the observable under these 100 variations to determine unc
# Generally these uncertainties are subdominant

nToys = LP_weights["stat_vars"].shape[1]
eff_toys = []
pt_eff_toys = []
for i in range(nToys):
    eff = np.average(score_cut, weights=LP_weights["stat_vars"][:, i])
    eff_toys.append(eff)

    eff1 = np.average(score_cut, weights=LP_weights["pt_vars"][:, i])
    pt_eff_toys.append(eff1)

# Compute stat and pt uncertainty based on variation in the toys
toys_mean = np.mean(eff_toys)
toys_std = np.std(eff_toys)
pt_toys_mean = np.mean(pt_eff_toys)
pt_toys_std = np.std(pt_eff_toys)

# if mean of toys is biased, also include it as an unc (should be zero)
eff_stat_unc = abs(toys_mean - eff_rw) + toys_std
eff_pt_unc = abs(pt_toys_mean - eff_rw) + pt_toys_std

print("Stat variation toys eff. avg %.3f, std dev %.3f" % (toys_mean, toys_std))
print("Pt variation toys eff. avg %.3f, std dev %.3f" % (pt_toys_mean, pt_toys_std))

# Other systematics come from up/down variations of the weights
sys_keys = ["sys", "bquark", "prongs", "unclust", "distortion"]
sys_uncs = dict()

# 'unclust' and 'distortion' are generally the largest uncertainties but there is some variation
# The up/down variations may not be symmetric, it depends on the use case


for sys in sys_keys:
    sys_uncs[sys] = [0.0, 0.0]


# Compute difference in efficiency due to weight variations as uncertainty
def get_uncs(cut, weights_up, weights_down, eff_baseline):
    eff_up = np.average(cut, weights=weights_up)
    eff_down = np.average(cut, weights=weights_down)

    unc_up = eff_up - eff_baseline
    unc_down = eff_down - eff_baseline
    return unc_up, unc_down


for sys in sys_keys:
    unc_up, unc_down = get_uncs(
        score_cut, LP_weights[sys + "_up"], LP_weights[sys + "_down"], eff_rw
    )
    sys_uncs[sys] = [unc_up, unc_down]


# Print uncertainty breakdown
eff_str = "Calibrated efficiency  is %.2f +/- %.2f (stat) +/- %.2f (pt)" % (
    eff_rw,
    eff_stat_unc,
    eff_pt_unc,
)
tot_unc_up = tot_unc_down = eff_stat_unc**2 + eff_pt_unc**2

for sys in sys_keys:
    eff_str += " %.2f/%.2f (%s)" % (sys_uncs[sys][0], sys_uncs[sys][1], sys)
    up_var = max(sys_uncs[sys][0], sys_uncs[sys][1])
    down_var = min(sys_uncs[sys][0], sys_uncs[sys][1])
    tot_unc_up += up_var**2
    tot_unc_down += down_var**2


tot_unc_up = tot_unc_up**0.5
tot_unc_down = tot_unc_down**0.5

# Note to get the 'scale factor' we just divide by the nominal efficiency

# Print final calibrated efficiency and total uncertaintiy
eff_str += "\n Original Eff. %.2f, Lund Plane Corrected Eff. %.2f +%.2f/-%.2f. SF %.2f +%.2f/-%.2f\n" % (
    eff_nom,
    eff_rw,
    tot_unc_up,
    tot_unc_down,
    eff_rw / eff_nom,
    tot_unc_up / eff_nom,
    tot_unc_down / eff_nom,
)

print(eff_str)
