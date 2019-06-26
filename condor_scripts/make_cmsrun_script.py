import os
import subprocess, time, sys, shlex
import fileinput

model = ['bbbb','bbmet']
prod = ['ggh','vbfh','vh']
mhmx = [[125,50],[300,125],[500,225],[1000,475],[2000,975]]
ctau = [500,1000,10000]

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
		script_dir = '../scripts/'+m+'/'+h+'/'
		os.system('mkdir -p '+list_dir)
		os.system('mkdir -p '+script_dir)
		list_name = m+'_'+h+'_'+str2+'_pl'+str(mmm)+'.list'
		script_name = list_name.replace('.list','_signal_aod.py')
#                print(name1)
		#print(home_dir+name1+'/'+name2+'/'+name3+'/*/*/'+name4)
		#cmd = 'ls '+home_dir+name1+'/'+name2+'/'+name3+'/*/*/'+name4+' > '+list_dir+list_name
		temp_name = 'signal_template_aod.py'
		#os.system('cp signal_template_aod.py '+script_dir+script_name)
		#list_name = list_dir+list_name
		script_name = script_dir+script_name
		print(script_name)
		with open(list_dir+list_name,'r') as lin:
			with open(temp_name,'r') as fin:
				with open(script_name,'w') as fout:
   					for line in fin:
						output_name = list_name.replace('.list','.root')
						#print(output_name)
						if '.root' in line:
							if 'ggH.root' in line:
        							new_line = line.replace('ntuple_RunIISummer16_bbbb_ggH.root', 'ntuple_RunIISummer16_'+output_name)
								fout.write(new_line)	
        						if 'CMSSW' in line:
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

	
