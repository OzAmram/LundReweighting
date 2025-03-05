# Code for the Lund Jet Plane Reweighting Method (JME-23-001)

For those seeking to just use the Lund Plane ratio for their analysis,
I recommend to look at the at the **`example.py` and `example_NanoV15.py`** scripts.
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

The `example_NanoV15.py` script shows how all of these ingredients can be obtained directly from the NanoAOD. 

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
    
## FAQ's

**Can I just compute a scale factor or do I have to propagate all the weights as different shape variations?**

If the final observable you fit is decorrlated with the substructure of the jet you are correcting (eg your search is for a resonance decaying to 2 or more jets and the final observable is the resonance mass), 
then it is fine to use a scale factor which just corrects the efficiency of your substructure cut because changes to the substructure will not change the shape of your final observable.
If your final observable is sensivite to the substructure (eg you are fitting some tagging score) then you need to propagate all the weights to define all the shape variations. 

**Can I apply the method to signals with displaced vertices? Eg. boosted tops or a->bb**

Partially. The method will correct the MC modeling of the substructure. It will not correct the MC modeling of the displaced vertices which likely also impacts the tagging of your jet. 
So if you can find some other way to correct the modeling of these secondary vertices and derive an appropriate uncertainty then you can combine this with the Lund Plane correction to get a full correction.

**Can I apply the method to signals with a lepton inside the jet? Eg. boosted semileptonic H->WW**
Yes this has been done in [HIG-24-008](https://cms.cern.ch/iCMS/analysisadmin/cadilines?line=HIG-24-008&tp=an&id=2807&ancode=HIG-24-008). 
Their procedure was to remove all leptons and photons from the PF candidates in a cone of dR < 0.2 around the gen-level lepton inside the jet.
The correction was then computed based on the remaining PF candidates/prongs. See [this presentation](https://indico.cern.ch/event/1470867/#23-lund-jet-plane-for-boosted) for details

**My systematic uncertainties are very asymmetric, is this typical?** 

Yes many of these systematic uncertainties are generally asymmetric. Particularly the uncertainty on the unclustered subjets is very often strongly asymmetric.
This uncertainty is computed by varying the weights of unclustered subjets up and down by a factor of 5, which often has an asymmetric affect on the efficiency. 

This is best illustrated with a toy example. 
Supposed I have 100 events in my sample which all have a nominal weight of 1. 40 of them pass my substructure cut, 5 of them have an unclustered prong
So my nominal efficiency is 40/100 = 40%. We vary the unclustered weight up and now the total weight of passing the cut is (35 + 5x5)=60 out of a total weight of 95+5x5 =120. So up variation efficiency is 50%
Now for the down variation passing weight is 35+5x0.2 = 36 out of a total weight of 96. So down variation efficiency is 38%. So the uncertainty on the efficiency is +10%/-2%

Now consider the case where the 5 unclustered jets don't pass the cut. Varying the weight up, the effiency is 40/120 = 30%, varying the weight down is 40/96 = 42%. So my uncertainty is +2% / -10% 

So we can see the asymmetric variations are natural to happen, and whether the up/down is larger depends on whether the unclustered jets pass the cut or don't.

**What is a reasonable bad matching fraction?** 

The bad matching fraction depends on the number of prongs in the jet and the pt of the prongs.
The larger number of prongs in the jet the harder it is to cleanly identify them in the reclustering, leading to a higher bad matching fraction.
The relationship between the matching fraction and pt is more subtle. If the jet is too low pt, then it is difficult to identify the prongs inside the jet
and the bad matching fraction is larger. At very high pt the jet becomes too boosted and the prongs start to overlap and it becomes difficult to separate them 
and the bad matching fraction increases. In the middle region there is a sweet spot where the bad matching fraction is lowest.
Due to these competing effects it is difficult to give a precise number on what bad matching fraction is unacceptable, 
but roughly, anything larger than 10% x n_prongs, may be an indication something is going wrong. If the bad matching fraction approaches becomes O(1) the reclustering is clearly no longer working and requires investigation. 
