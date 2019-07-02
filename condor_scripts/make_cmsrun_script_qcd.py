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

	script_dir = '../scripts/qcd/'
	os.system('mkdir -p '+script_dir)
	script_name = list_name.replace('.list','_aod.py')
	temp_name = 'qcd_template_aod.py'
	script_name = script_dir+script_name
	print(script_name)
	with open(list_dir+list_name,'r') as lin:
		with open(temp_name,'r') as fin:
			with open(script_name,'w') as fout:
				for line in fin:
					output_name = list_name.replace('.list','.root')
					#print(output_name)
					if '.root' in line:
						if 'qcd.root' in line:
							new_line = line.replace('ntuple_RunIIFall17_qcd.root', 'ntuple_RunIIFall17_'+output_name)
							fout.write(new_line)	
						if 'AODSIM' in line:
							for rt in lin:
								rt_name = rt.replace('\n','')
								#new_line = line.replace(line,"\'file:"+rt)
								new_line = line.replace(line,"\'file:"+rt_name+"\',\n")
								fout.write(new_line)	
					else:
						new_line = line
						#print(line.replace('ntuple_RunIISummer16_bbbb_ggH.root', 'ntuple_RunIISummer16_'+output_name))
						fout.write(new_line)	
		fin.close()
		fout.close()
		lin.close()

	#print(cmd)
	#make_list = subprocess.Popen([cmd],stdout=subprocess.PIPE, shell=True);
        #out = make_list.communicate()
        #print(out)

	
