# darkQCDplotter

* Takes input hepmc files produced by MadGraph: stored in ```inputs```
* Loops over events in python to fill histograms: written to ```outputs```
* Takes histograms and produces figures: written to ```plots``` 

## Steps to run on kdplab01

```
git clone --recurse-submodules git@github.com:kdipetri/darkQCDplotter.git
cd darkQCDplotter
```

to setup do this every time
```
singularity run --nv --bind /cvmfs /local/d1/badea/quick-start/my-image.sif
```

do this first time
```
pip install pyhepmc
pip install fastjet

```
make sure to setup input and output directories correctly

then run
```
source run.sh # make histos with eventlooper.py 
source merge.sh # merge output to single ROOT file
python plotDists.py # produces figures
``` 
