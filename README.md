# Code for the Lund Jet Plane Reweighting Method (JME-23-001)

For those seeking to just use the Lund Plane ratio for their analysis,
I recommend to look at the at the **`example.py` and `example\_NanoV15.py'** scripts.
Both scripts show step by step how the compute the correcttion, and use it to compute the corrected efficiency
and all components of the uncertainty on a substructure cut.
The first one uses signal inputs based on a custom h5 data format, which has
already been preprocessed to save the needed information. 
The latter directly uses the new [NanoV15](https://gitlab.cern.ch/cms-nanoAOD/nanoaod-doc/-/wikis/Releases/NanoAODv15) 
format for the signal MC, which contains already the needed PFCandidates.

The main part of the correction is implemented in the `LundReweighter` class.
It is recommended to use the `get_all_weights` function to compute weights from the
correction and all the associated systematic variations which encode the uncertainties.
This function contains a description of the needed inputs.

The main ingredients the user needs to provide are: 
- The 4-vectors of the AK8 jets one is reweighting
- The 4-vectors and PDG ID's of the gen-level quarks which are defining the prongs of the jet
- The 4-vectors of PF candidates contained inside the AK8 jet

The main piece of the correction, the data/MC Lund Jet Plane ratio and systematic uncertainties,
can be found as ROOT files in the `data` directory. There is one version for
each data taking year, you should use the version matching the year of the MC
you are applying it to. 

Note that most of the scripts here are
personal code used to derive the Lund Plane ratio, plots and results of
JME-23-001.
The `CASE` directory contains the script used for the application of the reweighting to the
CASE anomaly search (EXO-22-026). 
All of this code is based on running on files using the CASE [h5 file format](https://github.com/case-team/CASEUtils/tree/master/H5_maker)

For a short tutorial on the usage of the correction please see this
[presentation](https://indico.cern.ch/event/1379091/#7-calibrate-jets-with-more-tha) at the CMS DeepDive on boosted jet reconstruction.

For a full description of the method and a description of the different uncertainties,
    please consult the [JME-23-001 documentation](https://cms.cern.ch/iCMS/analysisadmin/cadilines?id=2660&ancode=JME-23-001&tp=an&line=JME-23-001). 
    
