# Code for the Lund Jet Plane Reweighting Method (JME-23-001)


For those seeking to just use the Lund Plane ratio for their analysis,
I recommend to look at the 
at the `example.py` script.
This contains an example of how to apply the method to compute the corrected efficiency
and all components of the uncertainty on a substructure cut.

It makes use of the `LundReweighter` class which does most of the heavy lifting.
It is recommended to use the `get_all_weights` function to compute weights from the
correction and all the associated systematic variations which encode the uncertainties
This function contains a description of the needed inputs 

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
    
