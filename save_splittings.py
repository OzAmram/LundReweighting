from Utils import *
import math

#d = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_sep29/"
d = "/uscms_data/d3/oamram/CASE_analysis/src/CASE/LundReweighting/Lund_output_files_jan17/"
f_data = d + "SingleMu_2018_merge.h5"
f_ttbar = d + "TT.h5"
f_wjets = d + "QCD_WJets.h5"
f_diboson = d + "diboson.h5"
f_tw = d + "TW.h5"
f_singletop = d + "SingleTop_merge.h5"

#flist = [f_data, f_ttbar, f_wjets, f_diboson, f_tw, f_singletop]
flist = [f_diboson, f_tw, f_singletop]


overwrite = True

num_excjets = 2
jetR = 1.0
tag = "2prong_kt"
batch_size = 10000


split_key = tag + "splittings"
subjet_key = tag + "subjets"
max_splittings = 100


for fname in flist:
    print("running " + fname)


    f = h5py.File(fname, "a")


    all_pf_cands = f['jet1_PFCands'][()].astype(np.float64)

    total_size = all_pf_cands.shape[0]
    iters = int(math.ceil(float(total_size)/batch_size))

    print("Run over %i batches" % iters)


    for it in range(iters):

        print("batch %i \n" %it)
        start_idx = it*batch_size
        end_idx = min(total_size, (it+1)*batch_size)

        pf_cands = all_pf_cands[start_idx:end_idx]


        splittings = np.zeros((pf_cands.shape[0], max_splittings, 3), dtype = np.float32)
        subjets = np.zeros((pf_cands.shape[0], 8), dtype = np.float32)



        for i,pf_cand in enumerate(pf_cands):
            evt_subjets, evt_splittings = get_splittings(pf_cand,  jetR =jetR, num_excjets = num_excjets)

            if(len(evt_splittings) == 0):
                evt_splittings = np.zeros((100,3))
            elif(len(evt_splittings) < max_splittings): #pad with zeros up to 100 splittings
                evt_splittings = np.array(evt_splittings)
                #print("before", evt_splittings.shape)
                evt_splittings = np.pad( evt_splittings, ((0, max_splittings - len(evt_splittings)),(0,0)), 'constant')
                #print("after", evt_splittings.shape)
            elif(len(evt_splittings) > max_splittings):
                print("Warning: splittings greater than max!")
                evt_splittings  = np.array(evt_splittings)[:max_splittings]

            splittings[i] = evt_splittings

            subjets[i] = np.array(evt_subjets).reshape(-1)

                    

        subjets = np.array(subjets)
        splittings = np.array(splittings)
        if(it==0):
            if(split_key in list(f.keys())):
                if(overwrite):
                    print("Overwriting keys %s and %s !" % (split_key, subjet_key))
                    f[split_key].resize((splittings.shape[0]), axis = 0)
                    f[split_key][:] = splittings
                    f[subjet_key].resize((subjets.shape[0]), axis = 0)
                    f[subjet_key][:] = subjets
                else:
                    print("Error: keys exist (%s and %s) but overwrite not specified!" % (split_key, subjet_key))
                    sys.exit(1)

            else:
                print("Creating datasets %s and %s " % (split_key, subjet_key))
                f.create_dataset(split_key, data=splittings, chunks = True, maxshape=(None, splittings.shape[1], splittings.shape[2]), compression = 'gzip')
                f.create_dataset(subjet_key, data=subjets, chunks = True, maxshape=(None, subjets.shape[1]), compression = 'gzip')

        else:
            f[split_key].resize((f[split_key].shape[0] + splittings.shape[0]), axis=0)
            f[split_key][-splittings.shape[0]:] = splittings

            f[subjet_key].resize((f[subjet_key].shape[0] + subjets.shape[0]), axis=0)
            f[subjet_key][-subjets.shape[0]:] = subjets
    f.close()





