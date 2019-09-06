import BLT.BLTAnalysis.BatchMaster as bm

import sys

''' Specify parameters '''
cfg        = bm.JobConfig
path       = 'root://cmseos.fnal.gov//store/user/corderom/sync_mc'
executable = 'execBatch.sh'

legacy = True

#selection  = 'jetjetgmjetjet'
selection  = 'mumug'
#selection  = 'elelg'
#selection  = 'ee'

period     = '2016'

#data_samples = ['single_mu']
#data_samples = ['single_el']
if legacy:
	if selection == 'mumug':
		data_samples = ['double_mu_legacy']
	elif selection == 'mumu':
		data_samples = ['double_mu_legacy']
	elif selection == 'ee':	
		data_samples = ['double_el_legacy']
	elif selection == 'elelg':	
		data_samples = ['double_el_legacy']

else:
	if selection == 'mumug':
		data_samples = ['double_mu']
	elif selection == 'mumu':
		data_samples = ['double_mu']
	elif selection == 'ee':	
		data_samples = ['double_el']
	elif selection == 'elelg':	
		data_samples = ['double_el']


mc_samples = []
################################################
#mc_samples += ['JB_TT']
#mc_samples += ['JB_ZG_ZToLL']
#mc_samples += ['JB_WJets']
#mc_samples += ['JB_DYJets']
#mc_samples += ['JB_WZ','JB_ZZ']



################################################
#mc_samples += ['JC_DYJets_Legacy']
#mc_samples += ['JC_WW_Legacy','JC_WZ_Legacy','JC_ZZ_Legacy']
mc_samples += ['JC_TTTo2L2Nu_Legacy']
mc_samples += ['JC_WJets_Legacy']
mc_samples += ['JC_ZGToLLG_Legacy']


''' 
    Set job configurations.  The order of arguments is: (Dataset, path to data,
    number of jobs, arguments to pass to executable, output directory name)
'''
############## DATA ########################
path = '/eos/uscms/store/group/lpcbacon/12d'


data_dict = {}
##############  DOUBLE MUON DATA ########################
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
#######################################################
##################### JC SAMPLES #######################
#######################################################
path = '/eos/uscms/store/user/corderom/sync_mc'

########################## ZH #######################
mc_dict = {}
mc_dict['JC_ZH'] = [
		    cfg(data_name = 'ZH_HToZG_ZToAll',
			path     = '{0}/ZH_HToZG_ZToAll_M-125_13TeV_powheg_pythia8'.format(path),
			nJobs    = 50,
			suffix   = 'ZH_HToZG_ZToAll'
		       ),
		    ]

########################## WH #######################
mc_dict['JC_WH'] = [
		    cfg(data_name = 'WplusH_HToZG_WToAll',
			path     = '{0}/WplusH_HToZG_WToAll_M125_13TeV_powheg_pythia8/sync_workarea_1_WplusH_HToZG_WToAll_M125_13TeV_powheg_pythia8/181116_184911/0000'.format(path),
			nJobs    = 50,
			#suffix   = 'WplusH_HToZG_WToAll'
			suffix   = 'wplush' 
		       ),
		    cfg(data_name = 'JC_WminusH_HtoZG_WToAll',
			path     = '{0}/WminusH_HToZG_WToAll_M125_13TeV_powheg_pythia8/sync_workarea_1_WminusH_HToZG_WToAll_M125_13TeV_powheg_pythia8/181116_185036/0000'.format(path),
			nJobs    = 40,
			#suffix   = 'WminusH_HtoZG_WToAll'
			suffix   = 'wminush'
		       ),
		    ]

########################## TT  #######################
mc_dict['JC_TT'] = [
		    cfg(data_name = 'TT',
			#path     = '{0}/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/'.format(path),
			path     = '{0}/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/sync_workarea_TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/190123_203318/0000/'.format(path),
			nJobs    = 10,
			suffix   = 'tt'
		       ),
		    ]

######################### VBF ZG #################################

VBF_dict = {}

VBF_dict['JC_VBFHToZG_ZToJJ']=[	
			    cfg(data_name = 'VBFHToZG_ZToJJ',
			        path     = '{0}/VBFHToZG_ZToJJ/'.format(path),
			        nJobs    = 10,
			        suffix   = 'vbfhtozg_ztojj'
			       ),
				]

VBF_dict['JC_VBFHToZG_ZToJJ_ONE']=[	
			    cfg(data_name = 'VBFHToZG_ZToJJ_ONE',
			        path     = '{0}/VBFHToZG_ZToJJ_ONE/'.format(path),
			        nJobs    = 10,
			        suffix   = 'vbfhtozg_ztojj_one'
			       ),
				]

########################## ZG ######################
mc_dict['JC_ZGToLLG']=[	
		 	cfg(data_name = 'ZGToLLG',
		 	    path     = '{0}/ZGToLLG_01J_5f_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/'.format(path),
		 	    #path     = '{0}/ZGTo2LG_amcatnlo_all_gen/'.format(path),
			    #path     = '{0}/ZGTo2LG_amcatnlo_calib/'.format(path),
			    nJobs    = 10,
			    suffix   = 'zgtollg'
			),
		]
#######################################################
############### GENERAL SAMPLES #######################
#######################################################

########################## DYJets #######################
path = '/eos/uscms/store/group/lpcbacon/12'
mc_dict['Summer16_DYJets'] = [
    cfg(data_name = 'DYJets',
        path     = '{0}/Summer16_DYJetsToLL_M-50_amcatnlo/'.format(path),
        nJobs    = 50,
        suffix   = 'dyjets'
       ),
    ]

########################## WJets #######################

mc_dict['Summer16_WJets'] = [
			    cfg(data_name = 'W1JetsToLNu',
				path     = '{0}/Summer16_W1JetsToLNu/'.format(path),
				nJobs    = 50,
				suffix   = 'w1jetstolnu'
			       ),
			    cfg(data_name = 'W2JetsToLNu',
				path     = '{0}/Summer16_W2JetsToLNu/'.format(path),
				nJobs    = 50,
				suffix = 'w2jetstolnu'
			       ),
			    cfg(data_name = 'W3JetsToLNu',
				path     = '{0}/Summer16_W3JetsToLNu/'.format(path),
				nJobs    = 10,
				suffix   = 'w3jetstolnu'
			       ),
			    cfg(data_name = 'W4JetsToLNu',
				path     = '{0}/Summer16_W4JetsToLNu/'.format(path),
				nJobs    = 10,
				suffix   = 'w4jetstolnu'
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
######################### WW #######################
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
			path     = '{0}/Summer16_WZTo3LNu_powheg/'.format(path),
			nJobs    = 10,
			suffix   = 'wzto3lnu'
		       ),
		    ]


#######################################################
################## JB SAMPLES #########################
#######################################################
path = '/eos/uscms/store/user/jbueghly/sync_mc/'

#########################  TT  #######################
mc_dict['JB_TT'] = [
		    cfg(data_name = 'TT',
			path     = '{0}/TT_powheg'.format(path),
			nJobs    = 30,
			suffix   = 'tt'
		       ),
			]

########################## DYJets  #######################
mc_dict['JB_DYJets'] = [
		    cfg(data_name = 'DYJets',
			path     = '{0}/DYJetsToLL_M-50_amcatnlo_all_gen_tmp'.format(path),
			nJobs    = 30,
			suffix   = 'dyjets'
		       ),
			]
######################## WJets  #######################
mc_dict['JB_WJets'] = [
		    cfg(data_name = 'W1JetsToLNu',
		        path     = '{0}/W1JetsToLNu_madgraph'.format(path),
		        nJobs    = 30,
		        suffix   = 'w1jetstolnu'
		       ),
		    cfg(data_name = 'W2JetsToLNu',
		        path     = '{0}/W2JetsToLNu_madgraph'.format(path),
		        nJobs    = 30,
		        suffix   = 'w2jetstolnu'
		       ),
		    cfg(data_name = 'W3JetsToLNu',
			path     = '{0}/W3JetsToLNu_madgraph_tmp'.format(path),
			nJobs    = 30,
			suffix   = 'w3jetstolnu'
		       ),
		    cfg(data_name = 'W4JetsToLNu',
			path     = '{0}/W4JetsToLNu_madgraph_tmp'.format(path),
			nJobs    = 30,
			suffix   = 'w4jetstolnu'
		       ),
		    ]
########################## WZ  #######################
mc_dict['JB_WZ'] = [
		    cfg(data_name = 'WZTo2L2Q',
			path     = '{0}/WZTo2L2Q_amcatnlo'.format(path),
			nJobs    = 30,
			suffix   = 'wzto2l2q'
		       ),
		    cfg(data_name = 'WZTo3LNu',
			path     = '{0}/WZTo3LNu_powheg'.format(path),
			nJobs    = 30,
			suffix   = 'wzto3lnu'
		       ),
			]
########################## ZH  #######################
mc_dict['JB_ZH'] = [
                    cfg(data_name = 'ZHToZG',
                        path     = '{0}/ZHtoZG_M-125_powheg'.format(path),
                        nJobs    = 30,
                        suffix   = 'zhtozg'
                       ),
                        ]
########################## ZZ  #######################
mc_dict['JB_ZZ'] = [
                    cfg(data_name = 'ZZTo2L2Q',
                        path     = '{0}/ZZTo2L2Q_amcatnlo'.format(path),
                        nJobs    = 30,
                        suffix   = 'zzto2l2q'
                       ),
                    cfg(data_name = 'ZZTo4L',
                        path     = '{0}/ZZTo4L_amcatnlo'.format(path),
                        nJobs    = 30,
                        suffix   = 'zzto4l'
                       ),
                        ]
########################## ZG ######################
mc_dict['JB_ZG_ZToLL']=[	
		 	cfg(data_name = 'ZG_ZToLL',
		 	    path     = '{0}/ZGTo2LG_amcatnlo_all_gen_tmp/'.format(path),
		 	    #path     = '{0}/ZGTo2LG_amcatnlo_all_gen/'.format(path),
			    #path     = '{0}/ZGTo2LG_amcatnlo_calib/'.format(path),
			    nJobs    = 10,
			    suffix   = 'zg_ztoll'
			),
		]
####################################################################


#----------------------------------------------------
#---------------------------------------------------
#######################################################
##################### JC SAMPLES Legacy #######################
#######################################################
SampleType = 'mc'
path = '/eos/uscms/store/user/corderom/'+SampleType+'_legacy_gen_2016'

########################## TT  #######################
mc_dict['JC_TTTo2L2Nu_Legacy'] = [
		    cfg(data_name = 'TTTo2L2Nu',
			path     = '{0}/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/2016_mc_legacy_gen_TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRu/190824_061004'.format(path),
			nJobs    = 20,
			suffix   = 'ttto2l2nu'
		       ),
		    ]

########################## ZG ######################
mc_dict['JC_ZGToLLG_Legacy']=[	
		 	cfg(data_name = 'ZGToLLG',
		 	    path     = '{0}/ZGToLLG_01J_5f_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/2016_mc_legacy_gen_ZGToLLG_01J_5f_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun/190823_153456'.format(path),
			    nJobs    = 20,
			    suffix   = 'zgtollg'
			),
		]


########################## DYJets ######################
mc_dict['JC_DYJets_Legacy']=[	
		 	cfg(data_name = 'DYJets',
		 	    path     = '{0}/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/2016_mc_legacy_gen_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRu/190824_062451'.format(path),
			    nJobs    = 20,
			    suffix   = 'dyjets'
			),
		]

########################## WW  ######################
mc_dict['JC_WW_Legacy']=[	
		 	cfg(data_name = 'WWToLNuQQ',
		 	    path     = '{0}/WWToLNuQQ_13TeV-powheg/2016_mc_legacy_gen_WWToLNuQQ_13TeV-powheg_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/190824_055930'.format(path),
			    nJobs    = 20,
			    suffix   = 'wwtolnuqq'
			),
		 	cfg(data_name = 'WWTo2L2Nu',
		 	    path     = '{0}/WWTo2L2Nu_13TeV-powheg-herwigpp/2016_mc_legacy_gen_WWTo2L2Nu_13TeV-powheg-herwigpp_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/190824_055749'.format(path),
			    nJobs    = 20,
			    suffix   = 'wwto2l2nu'
			),
		]


########################## WZ  ######################
mc_dict['JC_WZ_Legacy']=[	
		 	cfg(data_name = 'WZTo3LNu',
		 	    path     = '{0}/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/2016_mc_legacy_gen_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic/190824_060249'.format(path),
			    nJobs    = 20,
			    suffix   = 'wzto3lnu'
			),
		]

########################## ZZ  ######################
mc_dict['JC_ZZ_Legacy']=[	
		 	cfg(data_name = 'ZZTo2L2Nu',
		 	    path     = '{0}/ZZTo2L2Nu_13TeV_powheg_pythia8/2016_mc_legacy_gen_ZZTo2L2Nu_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/190824_060811'.format(path),
			    nJobs    = 20,
			    suffix   = 'zzto2l2nu'
			),
		 	cfg(data_name = 'ZZTo2L2Q',
		 	    path     = '{0}/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/2016_mc_legacy_gen_ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptoti/190824_060617'.format(path),
			    nJobs    = 20,
			    suffix   = 'zzto2l2q'
			),
		 	cfg(data_name = 'ZZTo4L',
		 	    path     = '{0}/ZZTo4L_13TeV-amcatnloFXFX-pythia8/2016_mc_legacy_gen_ZZTo4L_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-/190824_060434'.format(path),
			    nJobs    = 20,
			    suffix   = 'zzto4l'
			),
		]

mc_dict['JC_WJets_Legacy'] = [
			    cfg(data_name = 'WJets',
				path     = '{0}/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/2016_mc_legacy_gen_WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_as/190824_062650'.format(path),
				nJobs    = 30,
				suffix   = 'wjets'
			       ),
				]
########################## DoubleMuon  ######################
SampleType = 'data'
path = '/eos/uscms/store/user/corderom/'+SampleType+'_legacy_2016'

########################## DoubleMuon  ######################
data_dict['double_mu_legacy'] = [
			 cfg(data_name = 'DoubleMuon_2016B_v1',
			     path      = '{0}/DoubleMuon/2016_data_legacy_DoubleMuon_Run2016B-17Jul2018_ver2-v1/190812_184322'.format(path),
			     nJobs     = 30,
			     suffix    = 'DoubleMuon_2016B'
			    ),
			 cfg(data_name = 'DoubleMuon_2016C_v1',
			     path      = '{0}/DoubleMuon/2016_data_legacy_DoubleMuon_Run2016C-17Jul2018-v1/190812_184520'.format(path),
			     nJobs     = 30,
			     suffix    = 'DoubleMuon_2016C'
			    ),
			 cfg(data_name = 'DoubleMuon_2016D_v1',
			     path      = '{0}/DoubleMuon/2016_data_legacy_DoubleMuon_Run2016D-17Jul2018-v1/190812_184704'.format(path),
			     nJobs     = 30,
			     suffix    = 'DoubleMuon_2016D'
			    ),
			 cfg(data_name = 'DoubleMuon_2016E_v1',
			    path      = '{0}/DoubleMuon/2016_data_legacy_DoubleMuon_Run2016E-17Jul2018-v1/190812_184916'.format(path),
			    nJobs     = 30,
			    suffix    = 'DoubleMuon_2016E'
			   ),
			 cfg(data_name = 'DoubleMuon_2016F_v1',
			     path      = '{0}/DoubleMuon/2016_data_legacy_DoubleMuon_Run2016F-17Jul2018-v1/190812_185121'.format(path),
			     nJobs     = 30,
			     suffix    = 'DoubleMuon_2016F'
			    ),
			 cfg(data_name = 'DoubleMuon_2016G_v1',
			     path      = '{0}/DoubleMuon/2016_data_legacy_DoubleMuon_Run2016G-17Jul2018-v1/190812_185328'.format(path),
			     nJobs     = 30,
			     suffix    = 'DoubleMuon_2016G'
			    ),
			 cfg(data_name = 'DoubleMuon_2016H_v1',
			     path      = '{0}/DoubleMuon/2016_data_legacy_DoubleMuon_Run2016H-17Jul2018-v1/190812_185517'.format(path),
			     nJobs     = 30,
			     suffix    = 'DoubleMuon_2016H'
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

