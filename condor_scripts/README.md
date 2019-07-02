# Steps to run condor on cmsRun job at Caltech tier2

## Step 1
Make lists on signal samples:

```
python make_list.py
```

Or make lists on qcd samples: 

```
python make_list_qcd.py
```

## Step 2
Make cmsRun scripts from templates.

For signal:
```
python make_cmsrun_script.py
```
with template `signal_template_aod.py`.

For QCD:
```
python make_cmsrun_script_qcd.py
```
with template `qcd_template_aod.py`.

## Step 3
Make condor scripts from templates.

For signal:
```
python make_condor_script.py
```

For QCD:
```
python make_condor_script_qcd.py
```
with template `signal_template_aod_job.jbl` and `signal_template_aod_job.sh`.

