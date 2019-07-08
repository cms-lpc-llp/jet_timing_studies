import os
import subprocess, time, sys, shlex
import fileinput

pt_pairs = [[15,30],[30,50],[50,80],[80,120],[170,300],[470,600],[600,800],[800,1000]]

pwd = os.getcwd()
home_dir = '/mnt/hadoop/store/mc/RunIIFall17DRPremix/'

str1 = 'QCD_Pt_'
str2 = '_TuneCP5_13TeV_pythia8'
for i,pts in enumerate(pt_pairs):
	strpt = str(pts[0])+'to'+str(pts[1])
	name1 = str1 + strpt + str2
	name2 = 'AODSIM/94X_mc2017_realistic_v11_ext1-v1/'
	print(name1)
	list_dir = pwd.replace('condor_scripts','lists/qcd/')
	list_name = name1+'.list'

	script_dir = pwd.replace('condor_scripts','scripts/qcd/')
	os.system('mkdir -p '+script_dir)
	print(script_dir)
	script_name = list_name.replace('.list','_aod.py')
	temp_name = 'qcd_template_aod.py'
	
	exe_dir = pwd+'/qcd_jobs/job_exes/'
	jbl_dir = pwd+'/qcd_jobs/job_jbls/'
	log_dir = pwd+'/qcd_jobs/job_logs/'
	os.system('mkdir -p '+exe_dir)
	os.system('mkdir -p '+jbl_dir)
	os.system('mkdir -p '+log_dir)
	exe_name = script_name.replace('.py','.sh')
	jbl_name = script_name.replace('.py','.jbl')

	label = script_name.replace('.py','')

	exe_temp = 'signal_template_aod_job.sh'
	jbl_temp = 'signal_template_aod_job.jbl'
	
	output_dir = '/store/group/phys_exotica/jmao/jet_timing_studies/samples/qcd/'
	output_name = list_name.replace('.list','.root')
	#script_name = script_dir+script_name
	#print(output_name)
	with open(exe_temp,'r') as fin_exe:
		with open(exe_dir+exe_name,'w') as fout_exe:
			for line in fin_exe:
				if 'condor_scripts' in line:
					new_line = line.replace(line,'cd '+script_dir+'\n')
				elif 'bbbb_ggh_mh125_mx50_pl500_signal_aod_test.py' in line:
					new_line = line.replace('bbbb_ggh_mh125_mx50_pl500_signal_aod_test.py',script_name)
				elif '/store/group/phys_exotica/jmao/jet_timing_studies/samples/signal/' in line:
					new_line = line.replace('/store/group/phys_exotica/jmao/jet_timing_studies/samples/signal/',output_dir)
				elif 'bbbb_ggh_mh125_mx50_pl500.root' in line:
					new_line = line.replace('ntuple_RunIISummer16_bbbb_ggh_mh125_mx50_pl500.root','ntuple_RunIIFall17_'+output_name).replace('/store/user/jmao/output_test.root',output_dir+'ntuple_RunIIFall17_'+output_name)
				else: 
					new_line = line
				fout_exe.write(new_line)

		fout_exe.close()
		fin_exe.close()

	with open(jbl_temp,'r') as fin_jbl:
		with open(jbl_dir+jbl_name,'w') as fout_jbl:
			for line in fin_jbl:
				if 'signal_template_aod_job.sh' in line:
					new_line = line.replace('signal_template_aod_job.sh',exe_dir+exe_name)
				elif 'logs' in line:
					new_line = line.replace('logs/test_job',log_dir+label)
				else: 
					new_line = line
				fout_jbl.write(new_line)

		fout_jbl.close()
		fin_jbl.close()

	cmd = 'condor_submit '+jbl_dir+jbl_name
	print(cmd)
	make_job = subprocess.Popen([cmd],stdout=subprocess.PIPE, shell=True);
        out = make_job.communicate()
        print(out)

	
