from Utils import *
import math

#fname = "Lund_output_files_sep21/TTTo2L2Nu.h5"
fname = "Lund_output_files_sep21/TTToHadronic.h5"
#fname = "Lund_output_files_old/TTbar_semilep_2018_old.h5"
overwrite = True
f = h5py.File(fname, "a")
num_excjets = 2
jetR = -1
tag = "2prong_"
batch_size = 1000


split_key = tag + "splittings"
subjet_key = tag + "subjets"
max_splittings = 100


all_pf_cands = f['jet1_PFCands'][()].astype(np.float64)

total_size = all_pf_cands.shape[0]
iters = int(math.ceil(float(total_size)/batch_size))


for it in range(iters):

    print("batch %i \n" %it)
    start_idx = it*batch_size
    end_idx = min(total_size, (it+1)*batch_size)

    pf_cands = all_pf_cands[start_idx:end_idx]


    splittings = np.zeros((pf_cands.shape[0], max_splittings, 3), dtype = np.float32)
    subjets = np.zeros((pf_cands.shape[0], 8), dtype = np.float32)



    for i,pf_cand in enumerate(pf_cands):
        evt_subjets, evt_splittings = get_splittings(pf_cand,  jetR =jetR, num_excjets = num_excjets)

        if(len(evt_splittings) < max_splittings): #pad with zeros up to 100 splittings
            evt_splittings = np.array(evt_splittings)
            evt_splittings = np.pad( evt_splittings, ((0, max_splittings - len(evt_splittings)),(0,0)), 'constant')
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





