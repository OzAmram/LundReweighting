from array import array
import numpy as np

#Define constants used in method (ie binning)
jetR = 1.0
n_pt_bins = 6
pt_bins_reg = array('f', [10., 65., 110., 175., 240., 300., 99999.])
pt_bins_low = array('f', [0., 10., 20., 30, 50., 70., 100., 130., 175., 250., 350., 99999.])
pt_bins = pt_bins_reg

pt_max = 1000
dr_bin_min = -1.
dr_bin_max = 8.
kt_bin_min = -5
kt_bin_max = np.log(pt_max)
z_label = "ln(kt/GeV)"
y_label = "ln(0.8/#Delta)"
n_bins_LP = 20
n_bins = 40

dR_mu_bjet_cut = 0.4

kt_bins = array('f', np.linspace(kt_bin_min, kt_bin_max, num = n_bins_LP+1))
dr_bins = array('f', np.linspace(dr_bin_min, dr_bin_max, num = n_bins_LP+1))

B_PDG_ID = 5

