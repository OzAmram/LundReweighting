from .LundReweighter import *
from .Consts import *

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(False)
h_dummy = None

sys_weights_map = {
        'nom_weight' : 0,
        'pdf_up' : 1,
        'pdf_down': 2,
        'prefire_up': 3,
        'prefire_down' : 4,
        'pileup_up' : 5 ,
        'pileup_down' : 6,
        'btag_up' : 7,
        'btag_down' : 8,
        'PS_ISR_up' : 9,
        'PS_ISR_down' : 10,
        'PS_FSR_up' : 11,
        'PS_FSR_down' : 12,
        'F_up' : 13,
        'F_down' : 14,
        'R_up' : 15,
        'R_down' : 16,
        'RF_up' : 17,
        'RF_down' : 18,
        'top_ptrw_up' : 19,
        'top_ptrw_down' : 20,
        'mu_trigger_up': 21,
        'mu_trigger_down': 22,
        'mu_id_up': 23,
        'mu_id_down': 24,
        #'mu_iso_up': 25,
        #'mu_iso_down': 26,
        #'puID_up': 27,
        #'puID_down': 28,
        'QCD_norm_up': -999,
        'QCD_norm_down': -999,
        'Diboson_norm_up': -999,
        'Diboson_norm_down': -999,
        'Single_norm_up': -999, #single top
        'Single_norm_down': -999,
        'unmatched_norm_up': -999, #unmatched ttbar
        'unmatched_norm_down': -999,
        }

sig_sys = {
        'pdf_up' : 1,
        'pdf_down': 2,
        'prefire_up': 3,
        'prefire_down' : 4,
        'btag_up' : 7,
        'btag_down' : 8,
        'F_up' : 13,
        'F_down' : 14,
        'R_up' : 15,
        'R_down' : 16,
        'RF_up' : 17,
        'RF_down' : 18,
        'top_ptrw_up' : 19,
        'top_ptrw_down' : 20,
        'mu_trigger_up': 21,
        'mu_trigger_down': 22,
        'mu_id_up': 23,
        'mu_id_down': 24,
        'mu_iso_up': 25,
        'mu_iso_down': 26,
        'puID_up': 27,
        'puID_down': 28,
        }

        
def input_options():
#input options for all the of the different scripts. Not all options are used for everything

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--fin", default='', help="Input file")
    parser.add_argument("--fout", default='', help="Output file")
    parser.add_argument("-r", "--f_ratio", default='', help="Input ratio file")
    parser.add_argument("-o", "--outdir", default='test/', help="Output directory")
    parser.add_argument("--charge_only", default=False, action='store_true', help="Only charged particles in Lund Plane")
    parser.add_argument("--no_pt_extrap", default=False, action='store_true', help="Only charged particles in Lund Plane")
    parser.add_argument("--no_sys", default=False, action='store_true', help="No systematics")
    parser.add_argument("--max_evts", default=-1, type = int, help="Max number of evts to reweight")
    parser.add_argument("--num_jobs", default=1, type = int, help="Max number of evts to reweight")
    parser.add_argument("--job_idx", default=0, type = int, help="Max number of evts to reweight")
    parser.add_argument("-y", "--year", default=0, type = int, help="Year")
    parser.add_argument("--mode", default="",  help="Running mode")
    return parser




def convert_4vec(vec):
    rvec = ROOT.Math.PtEtaPhiMVector(vec[0], vec[1], vec[2], vec[3])
    return [rvec.Px(), rvec.Py(), rvec.Pz(), rvec.E()]




def get_dists(q_eta_phis, subjets_eta_phis):
    q_eta_phis = np.expand_dims(q_eta_phis, 2)
    subjets_eta_phis = np.expand_dims(subjets_eta_phis, 1)
    #print(q_eta_phis.shape)
    #print(subjets_eta_phis.shape)
    return np.sqrt(np.square(subjets_eta_phis[:,:,:,0] - q_eta_phis[:,:,:,0]) + 
            np.square(ang_dist(subjets_eta_phis[:,:,:,1], q_eta_phis[:,:,:,1] )))



def deltaR(v1, v2):
    dR = np.sqrt(np.square(v1[1] - v2[1]) + 
            np.square(ang_dist(v1[2], v2[2] )))
    return dR

class Dataset():
    def __init__(self, f, is_data = False, label = "", color = "", jms_corr = 1.0, dtype = 0):

        self.f = f
        self.is_data = is_data

        self.label = label
        self.color = color


        self.n_evts = f['event_info'].shape[-1]
        self.mask = f['jet_kinematics'][:,0] > 0.
        self.norm_factor = 1.0

        self.jms_corr = jms_corr
        self.sys_key = ""
        self.sys_power = 1.0

        self.norm_unc = 0.0
        self.dtype = dtype


    def n(self):
        return np.sum(self.mask)

    def apply_cut(self, cut):
        if(len(self.mask) == len(cut)): self.mask = self.mask & cut
        elif(len(self.mask[self.mask]) == len(cut)): self.mask[self.mask] = cut
        else:
            print("Mask length (%i, %i) and cut length (%i) incompatable! Skipping" % (len(self.mask), len(self.mask[self.mask]), len(cut)))
            exit(1)
    
    def get_masked(self, key):
        return self.f[key][()][self.mask]

    def apply_sys(self, sys_key):
        if(sys_key not in sys_weights_map.keys()): 
            print("Invalid sys %s!")
            exit(1)
        self.sys_key = sys_key

    def get_weights(self):
        max_weight = 50.
        if(self.is_data): return np.ones(self.n()) * self.norm_factor
        weights = self.get_masked('norm_weights') * self.norm_factor
        if(len(self.sys_key) > 0):
            sys_idx = sys_weights_map[self.sys_key]
            reweighting = np.power(self.get_masked('sys_weights')[:,sys_idx], self.sys_power)
            #reweighting = self.get_masked('sys_weights')[:,sys_idx]
            np.clip(reweighting, 0., 10.0)
            reweighting /= np.mean(reweighting)
            weights *=  reweighting
        weights = np.clip(weights, -max_weight, max_weight)
        return weights

    def compute_obs(self):
        eps = 1e-8
        feats = self.get_masked('jet1_extraInfo')
        kins = self.get_masked('jet_kinematics')
        self.tau21 = (feats[:,1] / (feats[:,0] + eps))
        self.tau32 = (feats[:,2] / (feats[:,1] + eps))
        self.tau43 = (feats[:,3] / (feats[:,2] + eps))
        if(feats.shape[1] > 8): self.DeepAK8_W_MD = feats[:,8]
        if(feats.shape[1] > 9): self.DeepAK8_W = feats[:,9]
        if(feats.shape[1] > 10): self.ParticleNet_W = feats[:,10]
        if(self.dtype == 1): self.mSoftDrop = kins[:,5] * self.jms_corr
        else: self.mSoftDrop = kins[:,3] * self.jms_corr
        self.pt = kins[:,0]
        self.nPF= feats[:,6]

    def compute_kinematics(self):
        mu = self.get_masked('mu_info')
        evt = self.get_masked('event_info')
        bjet = self.get_masked('btag_jet_info')
        ak8 = self.get_masked('jet_kinematics')

        self.mu_pt, self.mu_eta, self.mu_phi = mu[:,0], mu[:,1], mu[:,2]
        self.met_pt, self.met_phi = evt[:,1], evt[:,2]

        w_cand_px = self.mu_pt * np.cos(self.mu_phi) + self.met_pt * np.cos(self.met_phi)
        w_cand_py = self.mu_pt * np.sin(self.mu_phi) + self.met_pt * np.sin(self.met_phi)
        self.w_cand_pt = (w_cand_px**2 + w_cand_py**2)**0.5

        self.bjet_pt, self.bjet_eta, self.bjet_phi = bjet[:,0], bjet[:,1], bjet[:,2]

        self.ak8_pt, self.ak8_eta, self.ak8_phi = ak8[:,0], ak8[:,1], ak8[:,2]

        self.dphi_mu_bjet = ang_dist(self.mu_phi, self.bjet_phi)
        self.dphi_mu_ak8 = ang_dist(self.mu_phi, self.ak8_phi)

        self.dR_mu_bjet = (self.dphi_mu_bjet ** 2 + (self.mu_eta - self.bjet_eta)**2)**0.5
        self.dR_mu_ak8 = (self.dphi_mu_ak8 ** 2 + (self.mu_eta - self.ak8_eta)**2)**0.5



    def fill_LP(self, LP_rw, h, num_excjets = 2, prefix = "2prong", sys_variations = None, rescale_subjets = "vec"):


        pf_cands = self.get_masked("jet1_PFCands").astype(np.float64)
        splittings = subjets =  split = subjet = None



        jet_kinematics = self.get_masked("jet_kinematics")
        nom_weights = self.get_weights()

        rescale_vals = [1.] * len(pf_cands)
        if(rescale_subjets == "jec"):
            rescale_vals = self.get_masked("jet1_JME_vars")[:,-1]
        elif(rescale_subjets == "vec"):
            if(self.dtype ==1): 
                if(which_j == 1): j_4vec = self.get_masked('jet_kinematics')[:,2:6].astype(np.float64)
                else: j_4vec = self.get_masked('jet_kinematics')[:,6:10].astype(np.float64)

            else:
                j_4vec = self.get_masked('jet_kinematics')[:,:4].astype(np.float64)

            rescale_vals = j_4vec[:, 0]


        hists = [h]
        weights = [nom_weights]

        if(sys_variations is not None):

            all_sys_weights = self.get_masked('sys_weights')


            for sys in sys_variations.keys():
                if('norm' in sys): #normalization uncs split by process
                    process = sys.split("_")[0]
                    print(sys, self.label, process)
                    if(process in self.label and 'up' in sys): weights_sys = nom_weights * (1. + self.norm_unc)
                    elif(process in self.label and 'down' in sys): weights_sys = nom_weights  * (1. - self.norm_unc)
                    else: weights_sys = nom_weights
                else:
                    sys_idx = sys_weights_map[sys]
                    weights_sys = nom_weights * all_sys_weights[:, sys_idx]

                weights.append(weights_sys)
                hists.append(sys_variations[sys])


            #for idx,h in enumerate(hists):
            #    h.Print()
            #    print(weights[idx][:10])


        weights = np.array(weights, dtype = np.float32)
        subjets = []
        for i,pf_cand in enumerate(pf_cands):
            if(splittings is not None):
                split = splittings[i]
                subjet = subjets[i]


            subjet, _ = LP_rw.fill_lund_plane(hists, pf_cands = pf_cand, subjets = subjet, splittings = split, 
                    num_excjets = num_excjets, weight = weights[:,i], rescale_subjets = rescale_subjets, rescale_val = rescale_vals[i])
            subjets.append(subjet)

        return subjets
            

    def get_pt_response(self, gen_parts, subjets):
        deltaR_cut = 0.2
        gen_parts_eta_phi = gen_parts[:, 1:3]
        subjets = np.array(subjets)
        dists = get_subjet_dist(gen_parts_eta_phi, subjets[:,1:3])

        j_closest = np.amin(dists, axis = -1)
        j_which = np.argmin(dists, axis = -1)
        matches = j_which[j_closest < deltaR_cut]

        #check all quarks within 0.2 of subjet and no two quarks matched to same subjet
        no_match = np.sum(j_closest < deltaR_cut) != len(subjets)
        repeats = matches.shape[0] != np.unique(matches).shape[0]
        bad_match = no_match or repeats

        responses = subjets[:,0] / gen_parts[j_which,0]
        return responses





    def get_all_matched_splittings(self, LP_rw, num_excjets = 2, min_evts = None, max_evts = None, which_j =1, rescale_subjets = "vec"):


        pf_cands = self.get_masked("jet%i_PFCands" % which_j).astype(np.float64)[min_evts:max_evts]
        if(self.dtype ==1): 
            if(which_j == 1): j_4vec = self.get_masked('jet_kinematics')[min_evts:max_evts][:,2:6].astype(np.float64)
            else: j_4vec = self.get_masked('jet_kinematics')[min_evts:max_evts][:,6:10].astype(np.float64)

        else:
            j_4vec = self.get_masked('jet_kinematics')[min_evts:max_evts][:,:4].astype(np.float64)

        rescale_vals = [1.] * len(j_4vec)
        if(rescale_subjets == "jec"):
            rescale_vals = self.get_masked("jet%i_JME_vars" % which_j)[min_evts:max_evts, 12]
        elif(rescale_subjets == "vec"):
            rescale_vals = j_4vec[:,0]


        num_excjets_l = [num_excjets]*len(pf_cands)


        if(num_excjets > 0 and self.dtype < 0): 
            #No gen info
            subjets = []
            splittings = []

            for i in range(len(pf_cands)):
                subjet, split = LP_rw.get_splittings(pf_cands[i], num_excjets = num_excjets_l[i], rescale_subjets = rescale_subjets, rescale_val = rescale_vals[i])
                subjets.append(subjet)
                splittings.append(split)

            RO = ReclusterObj()
            RO.subjet, RO.split = subjets, splittings
            RO.n_prongs = num_excjets

            return RO, None, None

        if(self.dtype == 1): #CASE h5
            gen_parts = self.get_masked('gen_info')[min_evts:max_evts]
            n_evts = gen_parts.shape[0]
            gen_parts_eta_phi_raw = gen_parts[:,:,1:3]
            gen_pdg_id = np.abs(gen_parts[:,:,3])
            #neutrino pdg ids are 12,14,16
            is_lep = gen_pdg_id > 10
            not_neutrinos = ((~np.isclose(gen_pdg_id, 12)) & (~np.isclose(gen_pdg_id, 14)) & (~np.isclose(gen_pdg_id, 16)))
            
            gen_parts_eta_phi = [gen_parts_eta_phi_raw[i][not_neutrinos[i]] for i in range(n_evts)]
            #gen_parts_eta_phi = gen_parts_eta_phi[not_neutrinos].reshape(n_evts, -1, 2)
            is_lep = is_lep[not_neutrinos]


        else:#W or t matched MC
            gen_parts = self.get_masked('gen_parts')[min_evts:max_evts]
            q1_eta_phi = gen_parts[:,18:20]
            q2_eta_phi = gen_parts[:,22:24]
            b_eta_phi = gen_parts[:,26:28]
            if(self.dtype == 2): gen_parts_eta_phi = np.stack([q1_eta_phi,q2_eta_phi], axis = 1)
            else: gen_parts_eta_phi = np.stack([q1_eta_phi, q2_eta_phi, b_eta_phi], axis = 1)



        ro_nom = []
        ro_up = []
        ro_down = []
        for i in range(len(pf_cands)):
            ro, ro_prongs_up, ro_prongs_down = LP_rw.get_splittings_and_matching(pf_cands[i], gen_parts_eta_phi[i], j_4vec[i], rescale_subjets = rescale_subjets, rescale_val = rescale_vals[i])

            ro_nom.append(ro)
            ro_up.append(ro_prongs_up)
            ro_down.append(ro_prongs_down)

        return ro_nom, ro_up, ro_down


    def reweight_all_up_down_prongs(self, LP_rw, h_ratio = None, reclust_prongs_up = None, reclust_prongs_down = None, nom_weights = None):
        weights_match_up = np.copy(nom_weights)
        weights_match_down = np.copy(nom_weights)
        weights_prong_up = np.copy(nom_weights)
        weights_prong_down = np.copy(nom_weights)
        for i in range(len(reclust_prongs_up)):
            weights_prong_up[i], weights_prong_down[i], weights_match_up[i], weights_match_down[i] = LP_rw.get_up_down_prongs_weights(h_ratio, 
                                                                                    reclust_prongs_up[i], reclust_prongs_down[i], nom_weights[i])

        return weights_prong_up, weights_prong_down, weights_match_up, weights_match_down


    def reweight_all(self, LP_rw, h_ratio, num_excjets = 2, min_evts = None, max_evts =None, prefix = "", reclust_objs = None, 
            rand_noise = None,  pt_rand_noise = None, sys_str = "", subjets = None, splittings = None, norm = True):

        LP_weights = []
        LP_smeared_weights = []
        pt_smeared_weights = []
        pf_cands = self.get_masked("jet1_PFCands").astype(np.float64)[min_evts:max_evts]
        if('bquark' in sys_str): 
            if(self.dtype == 1):
                gen_parts = self.get_masked('gen_info')[min_evts:max_evts]
            else: 
                gen_parts = self.get_masked('gen_parts')[min_evts:max_evts]



        if(reclust_objs is None):
            reclust_objs, reclust_prongs_up, reclust_prongs_down = self.get_all_matched_splittings(LP_rw, num_excjets, min_evts = min_evts, max_evts =max_evts)

        n_b_matched = 0
        for i in range(len(pf_cands)):

            ro = reclust_objs[i]
            split, subjet = ro.split, ro.subjet

            if(sys_str == 'bquark'): # systematic only for bquarks
                deltaR_cut = 0.2
                if(self.dtype == 1): #CASE h5 saves all gen decays with pdg ID
                    B_ID = 5 
                    #pick out subjets matched to a b quark
                    gen_bs = [j for j in range(len(gen_parts[i])) if abs(gen_parts[i,j,3]) == B_ID]
                    dists = get_subjet_dist(gen_parts[i,gen_bs,1:3], np.array(subjet)[:,1:3])
                else: #ttbar MC saves in order
                    b_eta_phi = gen_parts[i,26:28]
                    dists = get_subjet_dist([b_eta_phi], np.array(subjet)[:,1:3])

                b_matches = []

                j_closest = np.amin(dists, axis = 0)
                j_which = np.argmin(dists, axis = 0)

                b_matches = np.unique(j_which[j_closest < deltaR_cut])

                #reweight only matched subjets
                if(len(b_matches) > 0):
                    b_subjet = [subjet[j] for j in range(len(subjet)) if j in b_matches]

                    #splittings save idx of associated subjet
                    b_split  = [split[j]  for j in range(len(split)) if split[j][0] in b_matches]

                    rw, smeared_rw, pt_smeared_rw  = LP_rw.reweight_lund_plane(h_ratio, subjets = b_subjet, splittings = b_split, sys_str = sys_str)
                    n_b_matched +=1
                else: 
                    rw = 1.0


            else:
                rw, smeared_rw, pt_smeared_rw  = LP_rw.reweight_lund_plane(h_ratio, reclust_obj = ro,
                        rand_noise = rand_noise, pt_rand_noise = pt_rand_noise, sys_str = sys_str)

            eps = 1e-6
            rw = max(rw, eps)
            LP_weights.append(rw)
            if(rand_noise is not None):
                LP_smeared_weights.append(smeared_rw)
            if(pt_rand_noise is not None):
                pt_smeared_weights.append(pt_smeared_rw)


        if(sys_str == 'bquark'):
            print("%i / %i jets were b-matched" % (n_b_matched, len(pf_cands)))
        if(norm):
            LP_weights = np.clip(np.array(LP_weights), 0., LP_rw.max_rw)
            LP_weights /= np.mean(LP_weights)
            LP_weights = np.clip(np.array(LP_weights), LP_rw.min_rw, LP_rw.max_rw)
            LP_weights /= np.mean(LP_weights)

        if(rand_noise is None):
            return LP_weights
        else:
            if(norm):
                LP_smeared_weights = np.clip(np.array(LP_smeared_weights), 0., LP_rw.max_rw)
                LP_smeared_weights /= np.mean(LP_smeared_weights, axis = 0)
                LP_smeared_weights = np.clip(np.array(LP_smeared_weights), LP_rw.min_rw, LP_rw.max_rw)
                LP_smeared_weights /= np.mean(LP_smeared_weights, axis = 0)

                pt_smeared_weights = np.clip(np.array(pt_smeared_weights), 0., LP_rw.max_rw)
                pt_smeared_weights /= np.mean(pt_smeared_weights, axis = 0)
                pt_smeared_weights = np.clip(np.array(pt_smeared_weights), LP_rw.min_rw, LP_rw.max_rw)
                pt_smeared_weights /= np.mean(pt_smeared_weights, axis = 0)

            return LP_weights, LP_smeared_weights, pt_smeared_weights
        

def add_dset(f, key, data):
    if(key in f.keys()):
        prev_size = f[key].shape[0]
        f[key].resize(( prev_size + data.shape[0]), axis=0)
        f[key][prev_size:] = data
    else:
        if(len(data) > 1):
            shape = list(data.shape)
            shape[0] = None
            f.create_dataset(key, data = data, chunks = True, maxshape = shape)
        else:
            f.create_dataset(key, data = data)

def fit_ratio(data, s_val, b_val):
    w = ROOT.RooWorkspace("w", "w")
    n_obs = ROOT.RooRealVar("n_obs", "", 10, 0, 1000000)

    s = ROOT.RooRealVar("s", "", s_val, 0, 1000000)
    s.setVal(float(s_val))
    s.setConstant(True)
    b = ROOT.RooRealVar("b", "", b_val, 0, 1000000)
    b.setVal(float(b_val))
    b.setConstant(True)
    r_guess = (data - b_val) / s_val
    r = ROOT.RooRealVar("r", "", r_guess, -10000, 10000)

    n_exp = ROOT.RooFormulaVar("n_exp", "@0 + @1 * @2", ROOT.RooArgList(b,s,r))
    model = ROOT.RooPoisson("poisson", "", n_obs, n_exp)

    ds = ROOT.RooDataSet("d", "d", ROOT.RooArgSet(n_obs))
    n_obs.setVal(float(data))

    ds.add(ROOT.RooArgSet(n_obs))
    fres = model.fitTo(ds,ROOT.RooFit.Save(1),ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Verbose(0)) 
    print(r.getErrorHi(), r.getErrorLo())
    return r.getVal(), r.getError()



