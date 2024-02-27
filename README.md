# darkQCDplotter

* Takes input hepmc files produced by MadGraph: stored in ```inputs```
* Loops over events in python to fill histograms: written to ```outputs```
* Takes histograms and produces figures: written to ```plots``` 

## Steps to run

```
mkdir inputs; mkdir outputs; mkdir plots
source run.sh # make histos
source merge.sh # merge output to single ROOT file
python plotDists.py # produces figures
``` 
