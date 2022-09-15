from Utils import *

fname = "TT_test.h5"
overwrite = True
f = h5py.File(fname, "a")
num_excjets = 2
jetR = -1
tag = "2prong_"


split_key = tag + "splittings"
subjet_key = tag + "subjets"
max_splittings = 150

pf_cands = f['jet1_PFCands'][()].astype(np.float64)[:1000]

all_splittings = np.zeros((pf_cands.shape[0], max_splittings, 3), dtype = np.float32)
all_subjets = np.zeros((pf_cands.shape[0], 8), dtype = np.float32)



for i,pf_cand in enumerate(pf_cands):
    subjets, splittings = get_splittings(pf_cands[i],  jetR =jetR, num_excjets = num_excjets)

    if(len(splittings) < max_splittings): #pad with zeros up to 100 splittings
        splittings = np.array(splittings)
        splittings = np.pad( splittings, ((0, max_splittings - len(splittings)),(0,0)), 'constant')
    elif(len(splittings) > max_splittings):
        print("Warning: splittings greater than max!")
        splittings  = np.array(splittings)[:max_splittings]

    all_subjets[i] = np.array(subjets).reshape(-1)
    all_splittings[i] = splittings

            

if(split_key in list(f.keys())):
    if(overwrite):
        print("Overwriting keys %s and %s !" % (split_key, subjet_key))
        f[split_key].resize((all_splittings.shape[0]), axis = 0)
        f[split_key][:] = all_splittings
        f[subjet_key].resize((all_subjets.shape[0]), axis = 0)
        f[subjet_key][:] = all_subjets
    else:
        print("Error: keys exist (%s and %s) but overwrite not specified!" % (split_key, subjet_key))
        sys.exit(1)

else:
    print("Creating datasets %s and %s " % (split_key, subjet_key))
    f.create_dataset(split_key, data=all_splittings, chunks = True, maxshape=(None, all_splittings.shape[1], all_splittings.shape[2]))
    f.create_dataset(subjet_key, data=all_subjets, chunks = True, maxshape=(None, all_subjets.shape[1]))






