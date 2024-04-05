from array import array
import numpy as np

#Define constants used in method (ie binning)
jetR = 1.0
n_pt_bins = 6
pt_bins = array('f', [0., 50., 100., 175., 250., 350., 99999.])

pt_max = 1000
dr_bin_min = -1.
dr_bin_max = 8.
kt_bin_min = -5
kt_bin_max = np.log(pt_max)
z_label = "ln(kt/GeV)"
y_label = "ln(0.8/#Delta)"
n_bins_LP = 20
n_bins = 40

kt_bins = array('f', np.linspace(kt_bin_min, kt_bin_max, num = n_bins_LP+1))
dr_bins = array('f', np.linspace(dr_bin_min, dr_bin_max, num = n_bins_LP+1))
