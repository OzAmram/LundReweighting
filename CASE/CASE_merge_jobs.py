from CASE_add_lund_weights import *

# merge multiple h5's together (like hadd)
# syntax is python H5_merge.py output_name.h5 input1.h5 input2.h5 input3.h5 ...

lund_keys = [ "lund_weights", "lund_mjj_check", "lund_weights_stat_var", "lund_weights_pt_var", "lund_weights_sys_var", "lund_weights_matching", "lund_weights_matching_unc"]

def merge(fin_name, fout):
    fin = h5py.File(fin_name, "r")

    for key in lund_keys:
        if('matching_unc' not in key): add_dset(fout, key, fin[key])

    fin.close()
        


def merge_multiple(fout, fs):
    print("Merging H5 files: ", fs)
    print(fs[-1])
    for fin_name in fs:
        print("Merging %s" % fin_name)
        merge(fin_name, fout)



if __name__ == "__main__":

    print("Dest %s" % sys.argv[1])
    fout = h5py.File(sys.argv[1], "r+")

    for key in lund_keys:
        if(key in list(fout.keys())):
            del fout[key]

    merge_multiple(fout, sys.argv[2:])

    h5_postprocess(fout)
    fout.close()

    
    print("Done!")

