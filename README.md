# Code for the Lund Reweighting Method (JME-23-001)

Most of this is personal code used to derive the Lund Plane ratio, plots and results of
JME-23-001.
Files with `CASE` in them are for the application of the reweighting to the
CASE anomaly search (EXO-22-026). 
All of this code is based on running on files using the CASE [h5 file format](https://github.com/case-team/CASEUtils/tree/master/H5_maker)
You can look at `W_SF`, `top_SF` or `CASE_SF` for examples of how a scale factor and uncertainties can be computed. 

For those seeking to just use the Lund Plane ratio for their analysis,
I recomend to use the functions provided in the `LundReweighter` file,
specifically the `reweight_lund_plane` function, which includes
a description of all command line arguments.

The Lund plane weights need to be normalized once they are computed for the
full MC sample (before any substructure cuts).
You can use the `normalize_weights` function to do this.

For a full description of the method and how the different uncertainties should be computed / used,
    please consult the [JME-23-001 documentation](https://cms.cern.ch/iCMS/analysisadmin/cadilines?id=2660&ancode=JME-23-001&tp=an&line=JME-23-001). 
    



