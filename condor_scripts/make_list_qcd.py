import os
import subprocess, time, sys, shlex

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
	os.system('mkdir -p '+list_dir)
	list_name = name1+'.list'
	cmd = 'ls '+home_dir+name1+'/'+name2+'/*/*.root > '+list_dir+list_name

	#print(cmd)
	make_list = subprocess.Popen([cmd],stdout=subprocess.PIPE, shell=True);
        out = make_list.communicate()
        print(out)

		
