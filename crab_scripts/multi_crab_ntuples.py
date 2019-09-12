if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    from WMCore.Configuration import Configuration
    config = Configuration()

    config.section_("General")
    config.General.workArea = 'crab'
    config.General.transferOutputs = True
    config.General.transferLogs = True

    config.section_("JobType")
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = '/afs/cern.ch/work/c/christiw/public/LLP/CMSSW_10_2_0/src/cms_lpc_llp/jet_timing_studies/python/jet_timing_studies_data_aod.py'

#    config.JobType.psetName = '/data/christiw/LLP/CMSSW_10_2_0/src/cms_lpc_llp/jet_timing_studies/python/jet_timing_studies_data_aod.py'
    config.JobType.numCores = 1
    config.section_("Data")
    config.Data.inputDBS = 'global'
    config.Data.splitting = 'LumiBased'
    config.Data.unitsPerJob = 10 #when splitting is 'Automatic', this represents jobs target runtime(minimum 180)
    #config.Data.totalUnits = 1
    config.Data.publication = True
    config.Data.ignoreLocality = True
    config.Data.lumiMask = '/afs/cern.ch/work/c/christiw/public/LLP/CMSSW_10_2_0/src/cms_lpc_llp/jet_timing_studies/data/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'

 #   config.Data.lumiMask = '/data/christiw/LLP/CMSSW_10_2_0/src/cms_lpc_llp/jet_timing_studies/data/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'

    config.section_("Site")
    config.Site.storageSite = 'T2_US_Caltech'
    config.Site.whitelist = ['T2_US_Caltech']
    config.Site.ignoreGlobalBlacklist = True
    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################
    dataset = ['/ZeroBias/Run2018A-17Sep2018-v1/AOD', '/ZeroBias/Run2018B-17Sep2018-v1/AOD','/ZeroBias/Run2018C-17Sep2018-v1/AOD']
    name = ['ZeroBias-Run2018A-17Sep2018-v1', 'ZeroBias-Run2018B-17Sep2018-v1','ZeroBias-Run2018C-17Sep2018-v1']


    for i in range(len(dataset)):
        config.General.requestName = 'CMSSW_10_2_0_'+name[i]+'_jettimingstudies_CaltechT2'
        config.Data.inputDataset = dataset[i]
        config.Data.outLFNDirBase = '/store/group/phys_exotica/delayedjets/jet_timing_studies/ZeroBias-17Sep2018-v1/'
        submit(config)
