import os
import subprocess, time, sys, shlex
import fileinput

model = ['bbbb']
prod = ['vh']
mhmx = [[125,50],[300,125],[500,225],[1000,475],[2000,975]]
ctau = [1000,10000]

pwd = os.getcwd()
home_dir = '/mnt/hadoop/store/user/christiw/RunIISummer16_withISR/'
for i,m in enumerate(model):
#     print(i,m)
    if m=='bbbb':
        str1 = 'SS1Tobb_SS2Tobb'
    else:
        str1 = 'SS1Tobb_SS2Toveve'
        
    for p,h in enumerate(prod):
#         print(p,h)
        
        for j,mm in enumerate(mhmx):
#             print(j,mm)
            str2 = 'mh'+str(mm[0])+'_mx'+str(mm[1])
        
            for k,mmm in enumerate(ctau):
#                 print(k,mmm)
                
#                 print(i,p,j,k,m,h,mm,mmm)

                name1 = 'ppTohToSS1SS2_'+str1+'_'+h+'_withISR_MC_prod'
                name2 = 'ppTohToSS1SS2_'+str1+'_'+h+'_withISR_'+str2+'_pl'+str(mmm)+'_ev100000'
                name3 = 'crab_CMSSW_8_0_21_'+name2.replace('withISR','ISR')+'_AOD_CaltechT2_v1'
                name4 = 'ppTohToSS1SS2_'+str1+'_ggh_withISR_step2_*.root'
		list_dir = pwd.replace('condor_scripts','lists/')+m+'/'+h+'/'
		#list_dir = m+'/'+h+'/'
		script_dir = pwd.replace('condor_scripts','scripts/')+m+'/'+h+'/'
		exe_dir = pwd+'/jobs/job_exes/'+m+'/'+h+'/'
		jbl_dir = pwd+'/jobs/job_jbls/'+m+'/'+h+'/'
		log_dir = pwd+'/jobs/job_logs/'+m+'/'+h+'/'
		os.system('mkdir -p '+list_dir)
		os.system('mkdir -p '+script_dir)
		os.system('mkdir -p '+exe_dir)
		os.system('mkdir -p '+jbl_dir)
		os.system('mkdir -p '+log_dir)
		list_name = m+'_'+h+'_'+str2+'_pl'+str(mmm)+'.list'
		script_name = list_name.replace('.list','_signal_aod.py')
		exe_name = script_name.replace('.py','.sh')
		jbl_name = script_name.replace('.py','.jbl')
	
		label = script_name.replace('.py','')
	
		exe_temp = 'signal_template_aod_job.sh'
		jbl_temp = 'signal_template_aod_job.jbl'
		output_dir = '/store/group/phys_exotica/delayedjets/jet_timing_studies/signal/'+m+'/'+h+'/'
		output_name = list_name.replace('.list','.root')
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
						new_line = line.replace('bbbb_ggh_mh125_mx50_pl500.root',output_name).replace('/store/user/jmao/output_test.root',output_dir+'ntuple_RunIISummer16_'+output_name)
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

	
