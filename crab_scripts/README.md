# Scripts to submit crab job



## For 2018 ZeroBias
Use lxplus6 (slc6)

### Setup CMSSW & clone the ntuples:
```bash
# Done once to setup environment
cmsrel CMSSW_10_2_0
cd CMSSW_10_2_0/src
git clone git@github.com:cms-lpc-llp/jet_timing_studies.git cms_lpc_llp/jet_timing_studies
scram b
cmsenv
export SCRAM_ARCH=slc6_amd64_gcc700

# setup crab environment
source /cvmfs/cms.cern.ch/crab3/crab.sh

# submit Zerobias for 2018A/B/C
python multi_crab_ntuples.py

```
