from Utils import *
import os
import math

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






parser = input_options()
options = parser.parse_args()

print(options)


outdir = options.outdir
if(not os.path.exists(outdir)): os.system("mkdir %s" % outdir)
jet_str = 'CA'



f_ratio = ROOT.TFile.Open(options.f_ratio)
f_sig = h5py.File(options.fin, "a")

jetR = 1.0
num_excjets = -1



#max_evts = 100
max_evts = None
batch_size = 10000

d = Dataset(f_sig, dtype =1)

nevts = d.get_masked('event_info').shape[0]
if(max_evts is not None): nevts = min(nevts, max_evts)
print("Nevts %i" % nevts)

print("input file", options.fin)

h_ratio = f_ratio.Get("ratio_nom")
f_ratio.cd('pt_extrap')
rdir = ROOT.gDirectory
#rdir = None

keys = [ "lund_weights", "lund_weights_stat_var", "lund_weights_pt_var", "lund_weights_sys_var", "lund_weights_matching_unc",]


#remove old values
for key in keys:
    if(key in list(f_sig.keys())):
        del f_sig[key]


nToys = 100

#Noise used to generated smeared ratio's based on stat unc
rand_noise = np.random.normal(size = (nToys, h_ratio.GetNbinsX(), h_ratio.GetNbinsY(), h_ratio.GetNbinsZ()))
pt_rand_noise = np.random.normal(size = (nToys, h_ratio.GetNbinsY(), h_ratio.GetNbinsZ(), 3))

LP_rw = LundReweighter(jetR = jetR, pt_extrap_dir = rdir, charge_only = options.charge_only)

iters = int(math.ceil(float(nevts)/batch_size))
print("Will run in %i bacthes" % iters)

bad_match = []

for i in range(iters):

    print("batch %i \n" %i)
    start_idx = i*batch_size
    stop_idx = min(nevts, (i+1)*batch_size)

    print("split1")
    subjets1, splittings1, bad_match1 = d.get_matched_splittings(LP_rw, num_excjets = num_excjets, which_j = 1, min_evts = start_idx, max_evts = stop_idx)
    print("split2")
    subjets2, splittings2, bad_match2 = d.get_matched_splittings(LP_rw, num_excjets = num_excjets, which_j = 2, min_evts = start_idx, max_evts = stop_idx)




    print("Rw1")
    weights1, smeared_weights1, pt_smeared_weights1 = d.reweight_LP(LP_rw, h_ratio, num_excjets = num_excjets, 
            min_evts = start_idx, max_evts = stop_idx,  prefix = "", rand_noise = rand_noise, pt_rand_noise = pt_rand_noise, subjets = subjets1, splittings = splittings1)

    print("Rw2")
    weights2, smeared_weights2, pt_smeared_weights2 = d.reweight_LP(LP_rw, h_ratio, num_excjets = num_excjets, 
            min_evts = start_idx, max_evts = stop_idx, prefix = "", rand_noise = rand_noise, pt_rand_noise = pt_rand_noise, subjets = subjets2, splittings = splittings2)

    is_lep = f_sig['event_info'][:,4][:max_evts]
    lep_cut = (is_lep < 0.5)



    weights = weights1 * weights2
    smeared_weights = smeared_weights1 * smeared_weights2
    pt_smeared_weights = pt_smeared_weights1 * pt_smeared_weights2

    #2 jets, so count each bad match as a 50% unc on the evt weight
    bad_match1_val = np.array([0.5*bm for bm in bad_match1])
    bad_match2_val = np.array([0.5*bm for bm in bad_match2])

    bad_match.append(bad_match1_val + bad_match2_val)


    weights /= np.mean(weights)
    smeared_weights = smeared_weights / np.mean(smeared_weights, axis = 1, keepdims=True)
    pt_smeared_weights = smeared_weights / np.mean(pt_smeared_weights, axis = 1, keepdims=True)



    print("sys rw")
    n_sys = 4
    sys_variations = np.ones((len(weights), n_sys))
    if(not options.no_sys):
        #sys_list = list(sys_weights_map.keys())
        sys_list = ["sys_tot_up", "sys_tot_down"]
        for i,sys in enumerate(sys_list):
            sys_ratio = f_ratio.Get("ratio_" + sys)
            sys_str = sys + "_"


            sys_weights1 = d.reweight_LP(LP_rw, sys_ratio, num_excjets = num_excjets, prefix = "", 
                    min_evts = start_idx, max_evts = stop_idx, sys_str = sys_str, subjets = subjets1, splittings = splittings1)

            sys_weights2 = d.reweight_LP(LP_rw, sys_ratio, num_excjets = num_excjets, prefix = "", 
                    min_evts = start_idx, max_evts = stop_idx, sys_str = sys_str, subjets = subjets2, splittings = splittings2)


            sys_weights = sys_weights1 * sys_weights2
            sys_weights /= np.mean(sys_weights)
            sys_variations[:,i] = sys_weights

        #vary weights up/down for b-quark subjets by ratio of b-quark to light quark LP
        b_light_ratio = f_ratio.Get("h_bl_ratio")
        bquark_rw1 = d.reweight_LP(LP_rw, b_light_ratio, num_excjets = num_excjets, prefix = "", 
                min_evts = start_idx, max_evts = stop_idx, sys_str = 'bquark', subjets = subjets1, splittings = splittings1)

        bquark_rw2 = d.reweight_LP(LP_rw, b_light_ratio, num_excjets = num_excjets, prefix = "", 
                min_evts = start_idx, max_evts = stop_idx, sys_str = 'bquark', subjets = subjets2, splittings = splittings2)

        up_bquark_weights = bquark_rw1 * bquark_rw2 * weights
        down_bquark_weights = (1./ (bquark_rw1 * bquark_rw2)) * weights

        up_bquark_weights /= np.mean(up_bquark_weights)
        down_bquark_weights /= np.mean(down_bquark_weights)

        sys_variations[:,2] = up_bquark_weights
        sys_variations[:,3] = down_bquark_weights

    #make sys variations multiplicative factors relative to nom
    sys_variations /= np.expand_dims(weights, -1)

    print(weights.shape)
    print(weights[:10])

    add_dset(f_sig, "lund_weights", data = weights, )
    add_dset(f_sig, "lund_weights_stat_var", data = smeared_weights, )
    add_dset(f_sig, "lund_weights_pt_var", data = pt_smeared_weights, )
    add_dset(f_sig, "lund_weights_sys_var", data = sys_variations, )

bad_match_eff = np.mean(bad_match)
print("Matching eff %.3f" % bad_match_eff)
add_dset(f_sig, "lund_weights_matching_unc", data = [bad_match_eff])


