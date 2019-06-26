import os
import subprocess, time, sys, shlex

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
		os.system('mkdir -p '+list_dir)
		list_name = m+'_'+h+'_'+str2+'_pl'+str(mmm)+'.list'
#                print(name1)
		#print(home_dir+name1+'/'+name2+'/'+name3+'/*/*/'+name4)
		cmd = 'ls '+home_dir+name1+'/'+name2+'/'+name3+'/*/*/'+name4+' > '+list_dir+list_name
		#print(cmd)
		make_list = subprocess.Popen([cmd],stdout=subprocess.PIPE, shell=True);
                out = make_list.communicate()
                print(out)

		
