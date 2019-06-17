mport BLT.BLTAnalysis.BatchMaster as bm

import sys

''' Specify parameters '''
cfg        = bm.JobConfig
#path       = '/eos/uscms/store/user/corderom/sync_mc'
path       = 'root://cmseos.fnal.gov//store/user/corderom/sync_mc'
executable = 'execBatch.sh'
#selection  = 'jetjetgmjetjet'
#selection  = 'mumugm'

selection  = 'mumug'
#selection  = 'elelg'

period     = '2016'

#data_samples = ['single_mu']
#data_samples = ['single_el']
if selection == 'mumug':
	data_samples = ['double_mu']
elif selection == 'elelg':	
	data_samples = ['double_el']

#mc_samples   = ['TT']
#mc_samples   = ['DYJets']
#mc_samples   = ['W1Jets','W2Jets','W3Jets','W4Jets' ]
#mc_samples   = ['ZH','WplusH','WminusH']#,'TT','W2Jets','W3Jets' ]
#mc_samples   = ['ZZTo2L2Q','WZTo2L2Q','WZTo1L1Nu2Q']

mc_samples = []
mc_samples += ['Summer16_WW','Summer16_WZ']
mc_samples += ['Summer16_ZZ']
#mc_samples += ['TT']
#mc_samples += ['DYJets']
#mc_samples += ['ZG_ZToLL']

#vbf_samples = ['VBFHToZG_ZToJJ']
#vbf_samples = ['VBFHToZG_ZToJJ_ONE']

''' 
    Set job configurations.  The order of arguments is: (Dataset, path to data,
    number of jobs, arguments to pass to executable, output directory name)
'''
############## DATA ########################
data_dict = {}

##############  DOUBLE MUON DATA ########################
path = '/eos/uscms/store/group/lpcbacon/12d'
data_dict['double_mu'] = [
         cfg(data_name = 'DoubleMuon_2016B_v2',
             path      = '{0}/DoubleMuon_Run2016B-03Feb2017_ver2-v2'.format(path),
             nJobs     = 30,
             suffix    = 'DoubleMuon_2016B'
            ),
         cfg(data_name = 'DoubleMuon_2016C_v1',
             path      = '{0}/DoubleMuon_Run2016C-03Feb2017-v1'.format(path),
             nJobs     = 30,
             suffix    = 'DoubleMuon_2016C'
            ),
         cfg(data_name = 'DoubleMuon_2016D_v1',
             path      = '{0}/DoubleMuon_Run2016D-03Feb2017-v1'.format(path),
             nJobs     = 30,
             suffix    = 'DoubleMuon_2016D'
            ),
         cfg(data_name = 'DoubleMuon_2016E_v1',
            path      = '{0}/DoubleMuon_Run2016E-03Feb2017-v1'.format(path),
            nJobs     = 30,
            suffix    = 'DoubleMuon_2016E'
           ),
         cfg(data_name = 'DoubleMuon_2016F_v1',
             path      = '{0}/DoubleMuon_Run2016F-03Feb2017-v1'.format(path),
             nJobs     = 30,
             suffix    = 'DoubleMuon_2016F'
            ),
         cfg(data_name = 'DoubleMuon_2016G_v1',
             path      = '{0}/DoubleMuon_Run2016G-03Feb2017-v1'.format(path),
             nJobs     = 30,
             suffix    = 'DoubleMuon_2016G'
            ),
         cfg(data_name = 'DoubleMuon_2016H_v2',
             path      = '{0}/DoubleMuon_Run2016H-03Feb2017_ver2-v1'.format(path),
             nJobs     = 30,
             suffix    = 'DoubleMuon_2016H'
            ),
         cfg(data_name = 'DoubleMuon_2016H_v3',
             path      = '{0}/DoubleMuon_Run2016H-03Feb2017_ver3-v1'.format(path),
             nJobs     = 30,
             suffix    = 'DoubleMuon_2016H'
            ),

]

##############  DOUBLE ELECTRON DATA ########################
data_dict['double_el'] = [
        cfg(data_name = 'DoubleEG_2016B_v2',
            path      = '{0}/DoubleEG_Run2016B-03Feb2017_ver2-v2'.format(path),
            nJobs     = 30,
            suffix    = 'DoubleEG_2016B'
           ),
        cfg(data_name = 'DoubleEG_2016C_v1',
            path      = '{0}/DoubleEG_Run2016C-03Feb2017-v1'.format(path),
            nJobs     = 30,
            suffix    = 'DoubleEG_2016C'
           ),
        cfg(data_name = 'DoubleEG_2016D_v2',
            path      = '{0}/DoubleEG_Run2016D-03Feb2017-v1'.format(path),
            nJobs     = 30,
            suffix    = 'DoubleEG_2016D'
           ),
        cfg(data_name = 'DoubleEG_2016E_v2',
            path      = '{0}/DoubleEG_Run2016E-03Feb2017-v1'.format(path),
            nJobs     = 30,
            suffix    = 'DoubleEG_2016E'
           ),
        cfg(data_name = 'DoubleEG_2016F_v2',
            path      = '{0}/DoubleEG_Run2016F-03Feb2017-v1'.format(path),
            nJobs     = 30,
            suffix    = 'DoubleEG_2016F'
           ),
        cfg(data_name = 'DoubleEG_2016G_v2',
            path      = '{0}/DoubleEG_Run2016G-03Feb2017-v1'.format(path),
            nJobs     = 30,
            suffix    = 'DoubleEG_2016G'
           ),
        cfg(data_name = 'DoubleEG_2016H_v2',
            path      = '{0}/DoubleEG_Run2016H-03Feb2017_ver2-v1'.format(path),
            nJobs     = 30,
            suffix    = 'DoubleEG_2016H'
           ),
        cfg(data_name = 'DoubleEG_2016H_v3',
            path      = '{0}/DoubleEG_Run2016H-03Feb2017_ver3-v1'.format(path),
            nJobs     = 30,
            suffix    = 'DoubleEG_2016H'
           ),
]
##############  SINGLE MUON DATA ########################

path = '/eos/uscms/store/group/lpcbacon/12a'
data_dict['single_mu'] = [
        #cfg(data_name = 'muon_2016B_v1',
        #    path      = '{0}/SingleMuon_Run2016B-03Feb2017_ver1-v1'.format(path),
        #    nJobs     = 30,
        #    suffix    = 'muon_2016B'
        #   ),
        #cfg(data_name = 'muon_2016B_v2',
        #    path      = '{0}/SingleMuon_Run2016B-03Feb2017_ver2-v2'.format(path),
        #    nJobs     = 30,
        #    suffix    = 'muon_2016B'
        #   ),
        #cfg(data_name = 'muon_2016C_v1',
        #    path      = '{0}/SingleMuon_Run2016C-03Feb2017-v1'.format(path),
        #    nJobs     = 30,
        #    suffix    = 'muon_2016C'
        #   ),
        #cfg(data_name = 'muon_2016D_v1',
        #    path      = '{0}/SingleMuon_Run2016D-03Feb2017-v1'.format(path),
        #    nJobs     = 30,
        #    suffix    = 'muon_2016D'
        #   ),
        #cfg(data_name = 'muon_2016E_v1',
        #    path      = '{0}/SingleMuon_Run2016E-03Feb2017-v1'.format(path),
        #    nJobs     = 30,
        #    suffix    = 'muon_2016E'
        #   ),
        #cfg(data_name = 'muon_2016F_v1',
        #    path      = '{0}/SingleMuon_Run2016F-03Feb2017-v1'.format(path),
        #    nJobs     = 30,
        #    suffix    = 'muon_2016F'
        #   ),
        cfg(data_name = 'muon_2016G',
            path      = '{0}/SingleMuon_Run2016G-03Feb2017-v1'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016G'
           ),
        cfg(data_name = 'muon_2016G',
            path      = '{0}/SingleMuon_Run2016G-03Feb2017-v1'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016G'
           ),
        cfg(data_name = 'muon_2016H_v2',
            path      = '{0}/SingleMuon_Run2016H-03Feb2017_ver2-v1'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016H'
           ),
        cfg(data_name = 'muon_2016H_v3',
            path      = '{0}/SingleMuon_Run2016H-03Feb2017_ver3-v1'.format(path),
            nJobs     = 30,
            suffix    = 'muon_2016H'
           ),
        ]

##############  SINGLE MUON DATA ########################
data_dict['single_el'] = [
        cfg(data_name = 'electron_2016B',
            nJobs    = 30,
            path     = '{0}/SingleElectron_Run2016B-03Feb2017_ver1-v1'.format(path),
            suffix   = 'electron_2016B'
           ),
        cfg(data_name = 'electron_2016B',
            path     = '{0}/SingleElectron_Run2016B-03Feb2017_ver2-v2'.format(path),
            nJobs    = 30,
            suffix   = 'electron_2016B'
           ),
        cfg(data_name = 'electron_2016C',
            path     = '{0}/SingleElectron_Run2016C-03Feb2017-v1'.format(path),
            nJobs    = 30,
            suffix   = 'electron_2016C'
           ),
        cfg(data_name = 'electron_2016D',
            path     = '{0}/SingleElectron_Run2016D-03Feb2017-v1'.format(path),
            nJobs    = 30,
            suffix   = 'electron_2016D'
           ),
        cfg(data_name = 'electron_2016E',
            path     = '{0}/SingleElectron_Run2016E-03Feb2017-v1'.format(path),
            nJobs    = 30,
            suffix   = 'electron_2016E'
           ),
        cfg(data_name = 'electron_2016F',
            path     = '{0}/SingleElectron_Run2016F-03Feb2017-v1'.format(path),
            nJobs    = 30,
            suffix   = 'electron_2016F'
           ),
        cfg(data_name = 'electron_2016G',
            path     = '{0}/SingleElectron_Run2016G-03Feb2017-v1'.format(path),
            nJobs    = 30,
            suffix   = 'electron_2016G'
           ),
        cfg(data_name = 'electron_2016H',
            path     = '{0}/SingleElectron_Run2016H-03Feb2017_ver2-v1'.format(path),
            nJobs    = 30,
            suffix   = 'electron_2016H'
           ),
        cfg(data_name = 'electron_2016H',
            path     = '{0}/SingleElectron_Run2016H-03Feb2017_ver3-v1'.format(path),
            nJobs    = 30,
            suffix   = 'electron_2016H'
           ),
        ]

#----------------------------------------------------
#---------------------------------------------------
path = '/eos/uscms/store/user/corderom/sync_mc'
mc_dict = {}
mc_dict['ZH'] = [
    cfg(data_name = 'ZH_HToZG_ZToAll',
        path     = '{0}/ZH_HToZG_ZToAll_M-125_13TeV_powheg_pythia8'.format(path),
        nJobs    = 50,
        suffix   = 'ZH_HToZG_ZToAll'
       ),
    ]
mc_dict['WplusH'] = [
    cfg(data_name = 'WplusH_HToZG_WToAll',
        path     = '{0}/WplusH_HToZG_WToAll_M125_13TeV_powheg_pythia8/sync_workarea_1_WplusH_HToZG_WToAll_M125_13TeV_powheg_pythia8/181116_184911/0000'.format(path),
        nJobs    = 50,
        #suffix   = 'WplusH_HToZG_WToAll'
	suffix   = 'wplush' 
       ),
    ]

mc_dict['WminusH'] = [
    cfg(data_name = 'WminusH_HtoZG_WToAll',
        path     = '{0}/WminusH_HToZG_WToAll_M125_13TeV_powheg_pythia8/sync_workarea_1_WminusH_HToZG_WToAll_M125_13TeV_powheg_pythia8/181116_185036/0000'.format(path),
        nJobs    = 40,
        #suffix   = 'WminusH_HtoZG_WToAll'
	suffix   = 'wminush'
       ),
    ]

mc_dict['TT'] = [
    cfg(data_name = 'TT',
        #path     = '{0}/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/'.format(path),
        path     = '{0}/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/sync_workarea_TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/190123_203318/0000/'.format(path),
        nJobs    = 10,
        suffix   = 'tt'
       ),
    ]

mc_dict['DYJets'] = [
    cfg(data_name = 'DYJets',
	path     = '/eos/uscms/store/user/jbueghly/sync_mc/DYJetsToLL_M-50_amcatnlo_all_gen_tmp/',
        #path     = '/eos/uscms/store/group/lpcbacon/12/Summer16_DYJetsToLL_M-50_amcatnlo/',
        #path     = '/eos/uscms/store/group/lpcbacon/12d/DYJetsToLL_M-50_madgraph/',
        nJobs    = 50,
        suffix   = 'dyjets'
	#suffic   = 'DYJets'
       ),
    ]

mc_dict['W1Jets'] = [
    cfg(data_name = 'W1JetsToLNu',
        #path     = '{0}/W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'.format(path),
        path     = '{0}/W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/sync_workarea_W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/190123_200719/0000/'.format(path),
        nJobs    = 50,
        #suffix   = 'W1JetsToLNu'
	suffix   = 'w1jetstolnu'
       ),
    ]

mc_dict['W2Jets'] = [
    cfg(data_name = 'W2JetsToLNu',
        #path     = '{0}/W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'.format(path),
        path     = '{0}/W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/sync_workarea_W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/190104_213056/0000/'.format(path),
        nJobs    = 50,
        #suffix   = 'W2JetsToLNu'
	suffix = 'w2jetstolnu'
       ),
    ]

mc_dict['W3Jets'] = [
    cfg(data_name = 'W3JetsToLNu',
        #path     = '{0}/W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'.format(path),
        path     = '{0}/W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/sync_workarea_W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/190104_213255/0000/'.format(path),
        nJobs    = 10,
        #suffix   = 'W3JetsToLNu'
	suffix   = 'w3jetstolnu'
       ),
    ]

mc_dict['W4Jets'] = [
    cfg(data_name = 'W4JetsToLNu',
        #path     = '{0}/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'.format(path),
        path     = '{0}/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/sync_workarea_W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/190104_214223/0000/'.format(path),
        nJobs    = 10,
        #suffix   = 'W4JetsToLNu'
        suffix   = 'w4jetstolnu'
       ),
    ]

mc_dict['WZTo2L2Q'] = [
    cfg(data_name = 'WZTo2L2Q',
        #path     = '{0}/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'.format(path),
        path     = '/eos/uscms/store/user/jbueghly/sync_mc/WZTo2L2Q_amcatnlo/',
        nJobs    = 10,
        suffix   = 'wzto2l2q'
       ),
    ]


mc_dict['WZTo1L1Nu2Q'] = [
    cfg(data_name = 'WZTo1L1Nu2Q',
        path     = '/eos/uscms/store/user/corderom/sync_mc/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/sync_workarea_WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/190203_065615/0000/',
        nJobs    = 10,
        suffix   = 'wzto1l1nu2q'
       ),
    ]
#############################################################
path = '/eos/uscms/store/group/lpcbacon/12'

######################### ZZ #######################
mc_dict['WJets'] = [
    cfg(data_name = 'WJets',
	path     = '/eos/uscms/store/user/jbueghly/sync_mc/DYJetsToLL_M-50_amcatnlo_all_gen_tmp/',
        nJobs    = 50,
        suffix   = 'wjets'
       ),
    ]
######################### ZZ #######################
mc_dict['Summer16_ZZ'] = [
		    cfg(data_name = 'ZZTo2L2Nu',
			path     = '{0}/Summer16_ZZTo2L2Nu_powheg/'.format(path),
			nJobs    = 10,
			suffix   = 'zzto2l2nu'
		       ),
		    cfg(data_name = 'ZZTo2L2Q',
			path     = '{0}/Summer16_ZZTo2L2Q_amcatnlo/'.format(path),
			nJobs    = 10,
			suffix   = 'zzto2l2q'
		       ),
		    cfg(data_name = 'ZZTo4L',
			path     = '{0}/Summer16_ZZto4L_amcatnlo/'.format(path),
			nJobs    = 10,
			suffix   = 'zzto4l'
		       ),
		    ]
######################### WW #######################mc_dict['Summer16_WW'] = [
mc_dict['Summer16_WW'] = [
		    cfg(data_name = 'WWTo2L2Nu',
			path     = '{0}/Summer16_WWTo2L2Nu_powheg/'.format(path),
			nJobs    = 10,
			suffix   = 'wwto2l2nu'
		       ),
		    ]
########################## WZ  #######################
mc_dict['Summer16_WZ'] = [
		    cfg(data_name = 'WZTo2L2Q',
			path     = '{0}/Summer16_WZTo2L2Q_amcatnlo/'.format(path),
			nJobs    = 10,
			suffix   = 'wzto2l2q'
		       ),
		    cfg(data_name = 'WZTo3LNu',
			#path     = '{0}/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'.format(path),
			path     = '{0}/Summer16_WZTo3LNu_powheg/'.format(path),
			nJobs    = 10,
			suffix   = 'wzto3lnu'
		       ),
		    ]

########################## ZG ######################
mc_dict['ZG_ZToLL']=[	
		 	cfg(data_name = 'ZG_ZToLL',
		 	    path     = '/eos/uscms/store/user/jbueghly/sync_mc/ZGTo2LG_amcatnlo_all_gen/',
			    nJobs    = 10,
			    suffix   = 'zg_ztoll'
			),
		]
######################################################################

VBF_dict = {}

VBF_dict['VBFHToZG_ZToJJ']=[	
			    cfg(data_name = 'VBFHToZG_ZToJJ',
			        path     = '/eos/uscms/store/user/corderom/sync_mc/VBFHToZG_ZToJJ/',
			        nJobs    = 10,
			        suffix   = 'vbfhtozg_ztojj'
			       ),
				]

VBF_dict['VBFHToZG_ZToJJ_ONE']=[	
			    cfg(data_name = 'VBFHToZG_ZToJJ_ONE',
			        path     = '/eos/uscms/store/user/corderom/sync_mc/VBFHToZG_ZToJJ_ONE/',
			        nJobs    = 10,
			        suffix   = 'vbfhtozg_ztojj_one'
			       ),
				]


batch_list = []
batch_list += sum([data_dict[n] for n in data_samples], []) 
#batch_list += sum([mc_dict[n] for n in mc_samples], []) 

#batch_list += sum([VBF_dict[n] for n in vbf_samples], []) 
batch = bm.BatchMaster(config_list = batch_list, 
                       stage_dir   = 'batch',
                       selection   = selection,
                       period      = period,
                       executable  = executable,
                       location    = 'lpc'
                     )
batch.submit_to_batch()

