#from pympler import asizeof
import sys
from Utils import *
import os
import math




def h5_postprocess(f):
    w_min = 0.1
    w_max = 10.
    #fix weight normalizations: first clip outliers, normalize, clip remaining outliers, normalize


    f['lund_weights'][:]  = np.clip(f['lund_weights'][:], 0., w_max)
    f['lund_weights'][:] /= np.mean(f['lund_weights'][:])
    f['lund_weights'][:]  = np.clip(f['lund_weights'][:], w_min, w_max)
    f['lund_weights'][:] /= np.mean(f['lund_weights'][:])
    print('weights_normed', f['lund_weights'][:10])

    f['lund_weights_stat_var'][:]  = np.clip(f['lund_weights_stat_var'][:], 0., w_max)
    f['lund_weights_stat_var'][:]  /= np.mean(f['lund_weights_stat_var'][:], axis = 0, keepdims=True)
    f['lund_weights_stat_var'][:]  = np.clip(f['lund_weights_stat_var'][:], w_min, w_max)
    f['lund_weights_stat_var'][:]  /= np.mean(f['lund_weights_stat_var'][:], axis = 0, keepdims=True)

    f['lund_weights_pt_var'][:]  = np.clip(f['lund_weights_pt_var'][:], 0., w_max)
    f['lund_weights_pt_var'][:]  /= np.mean(f['lund_weights_pt_var'][:], axis = 0, keepdims=True)
    f['lund_weights_pt_var'][:]  = np.clip(f['lund_weights_pt_var'][:], w_min, w_max)
    f['lund_weights_pt_var'][:]  /= np.mean(f['lund_weights_pt_var'][:], axis = 0, keepdims=True)

    f['lund_weights_sys_var'][:]  = np.clip(f['lund_weights_sys_var'][:], 0., w_max)
    f['lund_weights_sys_var'][:]  /= np.mean(f['lund_weights_sys_var'][:], axis = 0, keepdims=True)
    f['lund_weights_sys_var'][:]  = np.clip(f['lund_weights_sys_var'][:], w_min, w_max)
    f['lund_weights_sys_var'][:]  /= np.mean(f['lund_weights_sys_var'][:], axis = 0, keepdims=True)
    #make sys variations multiplicative factors relative to nom
    f['lund_weights_sys_var'][:] /= np.expand_dims(f['lund_weights'][:], -1)
    print('average sys mult fac.' , np.mean(f['lund_weights_sys_var'][:], axis = 0, keepdims=True))
    print( 'average sys weight', np.mean(f['lund_weights_sys_var'][:] * f['lund_weights'][:].reshape(-1,1), axis = 0, keepdims = True))

    #f['lund_weights_sys_var'][:]  = np.clip(f['lund_weights_sys_var'][:], 0.1, 10.)

    bad_match_eff = np.mean(f['lund_weights_matching'][:])


    print("Matching unc %.3f" % bad_match_eff)
    f.create_dataset('lund_weights_matching_unc', data = [bad_match_eff])




def run():

    debug = False
    if(debug):
        import tracemalloc
        tracemalloc.start()


    parser = input_options()
    options = parser.parse_args()

    print(options)


    outdir = options.outdir
    if(not os.path.exists(outdir)): os.system("mkdir %s" % outdir)
    jet_str = 'CA'



    
    f_sig = h5py.File(options.fin, "a")
    if(options.fout != ""):
        f_out = h5py.File(options.fout, "a")
    else:
        f_out = f_sig

    keys = [ "lund_weights", "lund_weights_stat_var", "lund_weights_pt_var", "lund_weights_sys_var", "lund_weights_matching_unc", "lund_weights_matching"]


    #remove old values
    for key in keys:
        if(key in list(f_out.keys())):
            del f_out[key]

    nevts_tot = f_sig['event_info'].shape[0]

    d = Dataset(f_sig, dtype =1)

    jetR = 1.0
    num_excjets = -1



    #max_evts = 100
    max_evts = None
    if(options.max_evts > 0):
        max_evts = options.max_evts

    batch_size = 1000
    #batch_size = 25


    if(max_evts is not None): nevts_tot = min(nevts_tot, max_evts)
    nevts_batch = nevts_tot // options.num_jobs

    global_start_idx = nevts_batch * options.job_idx

    if(options.job_idx == (options.num_jobs -1)):
        nevts_batch += nevts_tot % options.num_jobs


    print("Nevts_tot %i, Nevts_batch %i" % (nevts_tot, nevts_batch))
    print("input file", options.fin)
    print("output file", options.fout)



    nToys = 100

    f_ratio = ROOT.TFile.Open(options.f_ratio)
    h_ratio = f_ratio.Get("ratio_nom")
    f_ratio.cd('pt_extrap')
    rdir = ROOT.gDirectory
    LP_rw = LundReweighter(jetR = jetR, pt_extrap_dir = rdir, charge_only = options.charge_only)

    #Noise used to generated smeared ratio's based on stat unc
    np.random.seed(123)
    rand_noise = np.random.normal(size = (nToys, h_ratio.GetNbinsX(), h_ratio.GetNbinsY(), h_ratio.GetNbinsZ()))
    pt_rand_noise = np.random.normal(size = (nToys, h_ratio.GetNbinsY(), h_ratio.GetNbinsZ(), 3))
    print('rand', rand_noise[0,0,0,:5])

    iters = int(math.ceil(float(nevts_batch)/batch_size))
    print("Will run in %i batches" % iters)

    snapshots = []

    for i in range(iters):
        if(debug): 
            snapshots.append(tracemalloc.take_snapshot())

        sys.stdout.flush()
        print("batch %i \n" %i)





        start_idx = global_start_idx + i*batch_size
        stop_idx = global_start_idx + min(nevts_batch, (i+1)*batch_size)
        print('start, stop', start_idx, stop_idx)

        subjets1, splittings1, bad_match1 = d.get_matched_splittings(LP_rw, num_excjets = num_excjets, which_j = 1, min_evts = start_idx, max_evts = stop_idx)
        subjets2, splittings2, bad_match2 = d.get_matched_splittings(LP_rw, num_excjets = num_excjets, which_j = 2, min_evts = start_idx, max_evts = stop_idx)



        weights1, smeared_weights1, pt_smeared_weights1 = d.reweight_LP(LP_rw, h_ratio, num_excjets = num_excjets, norm = False,
                min_evts = start_idx, max_evts = stop_idx,  prefix = "", rand_noise = rand_noise, pt_rand_noise = pt_rand_noise, subjets = subjets1, splittings = splittings1)

        weights2, smeared_weights2, pt_smeared_weights2 = d.reweight_LP(LP_rw, h_ratio, num_excjets = num_excjets,  norm = False,
                min_evts = start_idx, max_evts = stop_idx, prefix = "", rand_noise = rand_noise, pt_rand_noise = pt_rand_noise, subjets = subjets2, splittings = splittings2)

        print(len(weights1))

        is_lep = f_sig['event_info'][:,4][:max_evts]
        lep_cut = (is_lep < 0.5)



        weights = np.array(weights1) * np.array(weights2)
        smeared_weights = np.array(smeared_weights1) * np.array(smeared_weights2)
        pt_smeared_weights = np.array(pt_smeared_weights1) * np.array(pt_smeared_weights2)

        #2 jets, so count each bad match as a 50% unc on the evt weight
        bad_match1_val = np.array([0.5*bm for bm in bad_match1])
        bad_match2_val = np.array([0.5*bm for bm in bad_match2])

        bad_match = bad_match1_val + bad_match2_val





        n_sys = 4
        sys_variations = np.ones((len(weights), n_sys))
        if(not options.no_sys):
            #sys_list = list(sys_weights_map.keys())
            sys_list = ["sys_tot_up", "sys_tot_down"]
            for i,sys_ in enumerate(sys_list):
                sys_ratio = f_ratio.Get("ratio_" + sys_)
                sys_str = sys_ + "_"


                sys_weights1 = d.reweight_LP(LP_rw, sys_ratio, num_excjets = num_excjets, prefix = "", norm = False,
                        min_evts = start_idx, max_evts = stop_idx, sys_str = sys_str, subjets = subjets1, splittings = splittings1)

                sys_weights2 = d.reweight_LP(LP_rw, sys_ratio, num_excjets = num_excjets, prefix = "",  norm = False,
                        min_evts = start_idx, max_evts = stop_idx, sys_str = sys_str, subjets = subjets2, splittings = splittings2)


                sys_variations[:,i] = np.array(sys_weights1) * np.array(sys_weights2)

            #vary weights up/down for b-quark subjets by ratio of b-quark to light quark LP
            b_light_ratio = f_ratio.Get("h_bl_ratio")
            bquark_rw1 = d.reweight_LP(LP_rw, b_light_ratio, num_excjets = num_excjets, prefix = "", norm = False,
                    min_evts = start_idx, max_evts = stop_idx, sys_str = 'bquark', subjets = subjets1, splittings = splittings1)

            bquark_rw2 = d.reweight_LP(LP_rw, b_light_ratio, num_excjets = num_excjets, prefix = "",  norm = False,
                    min_evts = start_idx, max_evts = stop_idx, sys_str = 'bquark', subjets = subjets2, splittings = splittings2)

            sys_variations[:,2] = np.array(bquark_rw1) * np.array(bquark_rw2) * weights
            sys_variations[:,3] = (1./ (np.array(bquark_rw1) * np.array(bquark_rw2))) * weights



        mjj = d.f['jet_kinematics'][start_idx:stop_idx,0]
        add_dset(f_out, "lund_weights", data = weights, )
        add_dset(f_out, "lund_mjj_check", data = mjj)
        add_dset(f_out, "lund_weights_stat_var", data = smeared_weights, )
        add_dset(f_out, "lund_weights_pt_var", data = pt_smeared_weights, )
        add_dset(f_out, "lund_weights_sys_var", data = sys_variations, )

        add_dset(f_out, "lund_weights_matching", data = bad_match)
        

        del subjets1, splittings1, weights1, smeared_weights1, pt_smeared_weights1
        del subjets2, splittings2, weights2, smeared_weights2, pt_smeared_weights2
        del weights, pt_smeared_weights, smeared_weights, sys_variations

    if(debug): 
        top_stats = snapshots[-1].statistics('lineno')
        top_stats_diff = snapshots[-1].compare_to(snapshots[1], 'lineno')

        print("[ Top 10 overall ]")
        for stat in top_stats[:10]:
            print(stat)


        print("[ Top 10 differences ]")
        for stat in top_stats_diff[:10]:
            print(stat)




    postprocess = (options.num_jobs <= 1)
    if(postprocess):
        h5_postprocess(f_out)

    f_out.close()
    f_sig.close()

if(__name__ == "__main__"):
    run()
