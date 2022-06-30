import numpy as np
import h5py
import os

idir = "../CASEUtils/H5_maker/ttbar_output_files_v2/"

ofile = "../CASEUtils/H5_maker/ttbar_output_files_v2/QCD_WJets_merged.h5"

files = [("WJetsToLNu_Pt-100To250.h5", 689.75), 
         ("WJetsToLNu_Pt-250To400.h5", 24.5),
         ("WJetsToLNu_Pt-400To600.h5", 3.11),
         ("WJetsToLNu_Pt-600ToInf.h5", 0.468),
         #("QCD_HT300to500.h5", 347700),
         #("QCD_HT500to700.h5", 3210),
         ("QCD_HT700to1000.h5", 6831),
         ("QCD_HT1000to1500.h5", 1207),
         ("QCD_HT1500to2000.h5", 119.9),
         ("QCD_HT2000toInf.h5", 25.24) ]

#ttbar xsec  364.34

merge_cmd = "python ../CASEUtils/H5_maker/H5_merge.py %s "  % ofile

for fname, xsec in files:
    f = h5py.File(idir + fname, "a")
    gen_weights = f['event_info'][:,3]
    presel_eff = f['preselection_eff'][0]
    rw_factor = xsec * 1000. * presel_eff / np.sum(gen_weights)
    f['event_info'][:,3] *= rw_factor
    print(fname, xsec * 1000. * presel_eff, rw_factor)
    merge_cmd += idir+fname + " " 

print(merge_cmd)
os.system(merge_cmd)

