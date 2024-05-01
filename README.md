# Code for the Lund Jet Plane Reweighting Method (JME-23-001)
# Version 0

This is the version of the method corresponding at the time of 'pre-approval' of JME-23-001. It is left here for posterity, but it 
is recommended to switch to the updated version located on the `main` branch


For those seeking to just use the Lund Plane ratio for their analysis,
I recommend to look at the 
at the `example.py` script.
This contains an example of how to apply the method to compute the corrected efficiency
and all components of the uncertainty on a substructure cut.

It makes use of the `LundReweighter` class which does most of the heavy
lifting.
Specifically the `get_splittings_and_matching` and `reweight_lund_plane` functions are the main functions, and they include
a description of the needed inputs.

Keep in mind that Lund plane weights need to be normalized once they are computed for the
full MC sample (before any substructure cuts).
You can use the `normalize_weights` function to do this.

A preliminary version of the data/MC Lund Jet Plane ratio (ie the main ingredient in the correction)
can be found in the `reweighting_files` directory. 


Note that most of the script here are
personal code used to derive the Lund Plane ratio, plots and results of
JME-23-001.
The `CASE` directory contains the script used for the application of the reweighting to the
CASE anomaly search (EXO-22-026). 
All of this code is based on running on files using the CASE [h5 file format](https://github.com/case-team/CASEUtils/tree/master/H5_maker)

For a full description of the method and how the different uncertainties should be computed / used,
    please consult the [JME-23-001 documentation](https://cms.cern.ch/iCMS/analysisadmin/cadilines?id=2660&ancode=JME-23-001&tp=an&line=JME-23-001). 
    



