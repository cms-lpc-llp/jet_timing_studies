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
    config.JobType.psetName = '/data/christiw/LLP/CMSSW_10_2_0/src/cms_lpc_llp/jet_timing_studies/python/jet_timing_studies_data_aod.py'
    config.JobType.numCores = 1
    config.section_("Data")
    config.Data.inputDBS = 'phys03'
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 10 #when splitting is 'Automatic', this represents jobs target runtime(minimum 180)
    config.Data.publication = True
    config.Data.ignoreLocality = True

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

    for i in range(len(dataset)):
        spec = mode+"_mh{}_mx{}_pl{}_ev{}".format(mh,mx,pl,ev)
        config.General.requestName = 'CMSSW_10_2_0_ZeroBias-17Sep2018-v1_jettimingstudies_CaltechT2'
        config.Data.inputDataset = dataset[i]
        config.Data.outLFNDirBase = '/store/group/phys_exotica/delayedjets/jet_timing_studies/ZeroBias-17Sep2018-v1/'
        #        print("gfal-mkdir -p gsiftp://transfer.ultralight.org/"+config.Data.outLFNDirBase)
        submit(config)

                                                                        
