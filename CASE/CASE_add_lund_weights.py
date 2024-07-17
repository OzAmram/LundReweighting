#from pympler import asizeof
import sys, os
sys.path.insert(0, '')
sys.path.append("../")
from utils.Utils import *
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


def run():

    debug = False
    if(debug):
        import tracemalloc
        tracemalloc.start()


    parser = input_options()
    options = parser.parse_args()

    print(options)

    
    f_sig = h5py.File(options.fin, "a")
    if(options.fout != ""):
        os.system("cp " + options.fin + " " + options.fout)
        f_out = h5py.File(options.fout, "a")
    else:
        f_out = f_sig

    keys = [ "lund_weights", "lund_weights_stat_var", "lund_weights_pt_var", "lund_weights_sys_var"]

    #overwrite 
    for key in keys:
        if(key in f_out.keys()):
            del f_out[key]



    max_evts = None
    if(options.max_evts > 0):
        max_evts = options.max_evts


    nevts_tot = f_sig['event_info'].shape[0]

    nums = np.arange(nevts_tot)

    if(max_evts is not None): nevts_tot = min(nevts_tot, max_evts)
    nevts_batch = nevts_tot // options.num_jobs


    if(options.job_idx == (options.num_jobs -1)):
        nevts_batch += nevts_tot % options.num_jobs

    start_idx = nevts_batch * options.job_idx
    stop_idx = start_idx + nevts_batch

    print("Nevts_tot %i, Nevts_batch %i" % (nevts_tot, nevts_batch))
    print("input file", options.fin)
    print("output file", options.fout)

    print('start, stop', start_idx, stop_idx)

    mask_batch = (nums >= start_idx) & ( nums < stop_idx)



    d16 = Dataset(f_sig, dtype =1)
    d17 = Dataset(f_sig, dtype =1)
    d18 = Dataset(f_sig, dtype =1)

    d16.apply_cut(mask_batch)
    d17.apply_cut(mask_batch)
    d18.apply_cut(mask_batch)


    years = d18.get_masked('event_info')[:,-2]

    mask16 = years < 2016.9 # 2016.5 or 2016
    mask17 = np.isclose(years, 2017.)
    mask18 = np.isclose(years, 2018.)

    masks = [mask16, mask17, mask18]

    d16.apply_cut(mask16)
    d17.apply_cut(mask17)
    d18.apply_cut(mask18)

    ds = [d16, d17, d18]
    ns = [ len(d.get_weights()) for d in ds]


    nToys = 100
    nSys = 10

    LP_weights = np.zeros(nevts_batch)
    LP_mjj_check = np.zeros(nevts_batch)
    LP_weights_stat_var = np.zeros((nevts_batch, nToys))
    LP_weights_pt_var = np.zeros((nevts_batch, nToys))
    LP_weights_sys_var = np.zeros((nevts_batch, nSys))


    #f_ratio = ROOT.TFile.Open(options.f_ratio)
    f_ratio18 = ROOT.TFile.Open("data/ratio_2018.root")
    f_ratio17 = ROOT.TFile.Open("data/ratio_2017.root")
    f_ratio16 = ROOT.TFile.Open("data/ratio_2016.root")


    LP_rws = [ LundReweighter(f_ratio = f_ratio16 ), LundReweighter(f_ratio = f_ratio17 ), LundReweighter(f_ratio = f_ratio18 )]

    #Noise used to generated smeared ratio's based on stat unc
    np.random.seed(123)
    rand_noise = np.random.normal(size = (nToys, LP_rws[0].h_ratio.GetNbinsX(), LP_rws[0].h_ratio.GetNbinsY(), LP_rws[0].h_ratio.GetNbinsZ()))
    pt_rand_noise = np.random.normal(size = (nToys, LP_rws[0].h_ratio.GetNbinsY(), LP_rws[0].h_ratio.GetNbinsZ(), 3))
    print('rand', rand_noise[0,0,0,:5])

    snapshots = []
    years = [2016, 2017, 2018]
    sys_keys = ['sys_up', 'sys_down', 'bquark_up', 'bquark_down', 'prongs_up', 'prongs_down', 'unclust_up', 'unclust_down', 'distortion_up', 'distortion_down']

    for i, d in enumerate(ds):
        if(debug): 
            snapshots.append(tracemalloc.take_snapshot())

        sys.stdout.flush()
        print("year %i, %i evts \n" %(years[i], ns[i]))



        LP_weights1 = d.reweight_all(LP_rws[i], which_j = 1, min_evts = start_idx, max_evts = stop_idx, rand_noise = rand_noise, pt_rand_noise = pt_rand_noise)
        LP_weights2 = d.reweight_all(LP_rws[i], which_j = 2, min_evts = start_idx, max_evts = stop_idx, rand_noise = rand_noise, pt_rand_noise = pt_rand_noise)

        LP_weights[masks[i]] = LP_weights1['nom'] * LP_weights2['nom']
        LP_weights_stat_var[masks[i]] = LP_weights1['stat_vars'] * LP_weights2['stat_vars']
        LP_weights_pt_var[masks[i]] = LP_weights1['pt_vars'] * LP_weights2['pt_vars']

        bad_match = 0.5* np.array(LP_weights1['bad_match']) +  0.5 * np.array(LP_weights2['bad_match'])
        print("Bad match %.3f" % np.mean(bad_match))

        #need to be careful with indexing
        mask_idx = np.argwhere(masks[i])
        for j,key in enumerate(sys_keys):
            new = (LP_weights1[key] * LP_weights2[key]).reshape(-1,1)
            LP_weights_sys_var[mask_idx,j]= new

        print(LP_weights_sys_var[masks[i]][:3,0])

        LP_mjj_check[masks[i]] = d.get_masked('jet_kinematics')[:,-0]

        del LP_weights1, LP_weights2

    print(LP_weights[0])
    print(LP_weights_stat_var[0,:3])
    print(LP_weights_pt_var[0,:3])
    print(LP_weights_sys_var[:3,0])
    print(LP_mjj_check[:5])
    print(f_sig['jet_kinematics'][:5,0])
    add_dset(f_out, "lund_weights", data = LP_weights )
    add_dset(f_out, "lund_mjj_check", data = LP_mjj_check)
    add_dset(f_out, "lund_weights_stat_var", data = LP_weights_stat_var )
    add_dset(f_out, "lund_weights_pt_var", data = LP_weights_pt_var )
    add_dset(f_out, "lund_weights_sys_var", data = LP_weights_sys_var)


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
