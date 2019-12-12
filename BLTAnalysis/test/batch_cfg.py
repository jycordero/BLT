import BLT.BLTAnalysis.BatchMaster as bm

import sys

''' Specify parameters '''
cfg        = bm.JobConfig
path       = 'root://cmseos.fnal.gov//store/user/corderom/sync_mc'
executable = 'execBatch.sh'

legacy = True
#legacy = False

#selection  = 'jetjetgmjetjet'
#selection  = 'mumug'
#selection  = 'mumu'
#selection  = 'elelg'
selection  = 'ee'

#period     = '2016'
period     = '2017'



##########################################################################################
data_samples = []
mc_samples   = []

if period == '2016':
	if legacy:
		############  DATA
		if selection == 'mumug':
			data_samples += ['double_mu_legacy']
		elif selection == 'mumu':
			data_samples += ['double_mu_legacy']
		elif selection == 'ee':	
			data_samples += ['single_el_legacy']
		elif selection == 'elelg':	
			data_samples += ['double_el_legacy']

		############  MC
		if selection == "mumug":
			mc_samples += ['JC_DYJets_Legacy']
			mc_samples += ['JC_WW_Legacy','JC_WZ_Legacy','JC_ZZ_Legacy']
			mc_samples += ['JC_TTTo2L2Nu_Legacy']
			mc_samples += ['JC_WJets_Legacy']
			mc_samples += ['JC_ZGToLLG_Legacy']
		elif selection == "ee":
			mc_samples += ['JC_DYJets_Legacy']
			mc_samples += ['JC_WJets_Legacy']
	else:
		############  DATA
		if selection == 'mumug':
			data_samples += ['double_mu']
		elif selection == 'mumu':
			data_samples += ['double_mu']
		elif selection == 'ee':	
			data_samples += ['single_el']
		elif selection == 'elelg':	
			data_samples += ['double_el']

		############  MC
		if selection == "mumug":
			mc_samples += ['JB_TT']
			mc_samples += ['JB_ZG_ZToLL']
			mc_samples += ['JB_WJets']
			mc_samples += ['JB_DYJets']
			mc_samples += ['JB_WZ','JB_ZZ']
		elif selection == "ee":
			mc_samples += ['JB_DYJets']
			mc_samples += ['JB_WJets']
		elif selection == "mumu":
			mc_samples += ['JB_DYJets']
elif period == '2017':
	############  DATA
	if selection == "mumug":
		data_samples += ['double_mu_2017Reco']
	elif selection == 'ee':
		data_samples += ['single_el_2017Reco']
	elif selection == 'mumu':
		data_samples += ["single_mu_2017Reco"]
	
	############  MC
	if selection == "mumug":
		mc_samples += ['JC_TT_2017Reco']
		mc_samples += ['JC_ZG_ZToLL_2017Reco']
		mc_samples += ['JC_WJets_2017Reco']
		mc_samples += ['JC_DYJets_2017Reco']
		mc_samples += ['JC_WZ_2017Reco','JC_ZZ_2017Reco','JC_WW_2017Reco']

		#mc_samples += ['JC_WZ_2017Reco']
		#mc_samples += ['JC_ZZ_2017Reco']
		#mc_samples += ['JC_WW_2017Reco']

	elif selection == "ee":
		mc_samples += ['JC_DYJets_2017Reco']
		mc_samples += ['JC_WJets_2017Reco']
################################################


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
			cfg(data_name = 'muon_2016B_v1',
			    path      = '{0}/SingleMuon_Run2016B-03Feb2017_ver1-v1'.format(path),
			    nJobs     = 30,
			    suffix    = 'muon_2016B'
			   ),
			cfg(data_name = 'muon_2016B_v2',
			    path      = '{0}/SingleMuon_Run2016B-03Feb2017_ver2-v2'.format(path),
			    nJobs     = 30,
			    suffix    = 'muon_2016B'
			   ),
			cfg(data_name = 'muon_2016C_v1',
			    path      = '{0}/SingleMuon_Run2016C-03Feb2017-v1'.format(path),
			    nJobs     = 30,
			    suffix    = 'muon_2016C'
			   ),
			cfg(data_name = 'muon_2016D_v1',
			    path      = '{0}/SingleMuon_Run2016D-03Feb2017-v1'.format(path),
			    nJobs     = 30,
			    suffix    = 'muon_2016D'
			   ),
			cfg(data_name = 'muon_2016E_v1',
			    path      = '{0}/SingleMuon_Run2016E-03Feb2017-v1'.format(path),
			    nJobs     = 30,
			    suffix    = 'muon_2016E'
			   ),
			cfg(data_name = 'muon_2016F_v1',
			    path      = '{0}/SingleMuon_Run2016F-03Feb2017-v1'.format(path),
			    nJobs     = 30,
			    suffix    = 'muon_2016F'
			   ),
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
######################################################################################################################
############################################## JC SAMPLES 2016 #######################################################
######################################################################################################################
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
mc_dict['JB_DYJets_Madgraph'] = [
		    cfg(data_name = 'DYJets',
			path     = '/eos/uscms/store/group/lpcbacon/15/DYJetsToLL_M_50_TuneCP5_13TeV_10X/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/CRAB3/190213_170924',
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
######################################################################################################################
########################################### JC SAMPLES Legacy ########################################################
######################################################################################################################
SampleType = 'mc'
#path = '/eos/uscms/store/user/corderom/'+SampleType+'_legacy_gen_2016'
path = '/eos/uscms/store/user/lpcbacon/corderom/mc_2016'
########################## TT  #######################
mc_dict['JC_TTTo2L2Nu_Legacy'] = [
		    cfg(data_name = 'TTTo2L2Nu',
			path     = '{0}/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/2016_mc_legacy_gen_TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRu/190824_061004/'.format(path),
			nJobs    = 20,
			suffix   = 'ttto2l2nu'
		       ),
		    ]

########################## ZG ######################
mc_dict['JC_ZGToLLG_Legacy']=[	
		 	cfg(data_name = 'ZGToLLG',
		 	    path     = '{0}/ZGToLLG_01J_5f_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/2016_mc_legacy_gen_ZGToLLG_01J_5f_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun/190823_153456/'.format(path),
			    nJobs    = 20,
			    suffix   = 'zgtollg'
			),
		]


########################## DYJets ######################
mc_dict['JC_DYJets_Legacy']=[	
		 	cfg(data_name = 'DYJets',
		 	    path     = '{0}/Old/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/2016_mc_legacy_gen_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRu/190824_062451/'.format(path),
			    nJobs    = 50,
			    suffix   = 'dyjets'
			),
		]

########################## WW  ######################
mc_dict['JC_WW_Legacy']=[	
		 	cfg(data_name = 'WWToLNuQQ',
		 	    path     = '{0}/WWToLNuQQ_13TeV-powheg/2016_mc_legacy_gen_WWToLNuQQ_13TeV-powheg_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2/190824_055930/'.format(path),
			    nJobs    = 20,
			    suffix   = 'wwtolnuqq'
			),
		 	cfg(data_name = 'WWTo2L2Nu',
		 	    path     = '{0}/WWTo2L2Nu_13TeV-powheg-herwigpp/2016_mc_legacy_gen_WWTo2L2Nu_13TeV-powheg-herwigpp_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/190824_055749/'.format(path),
			    nJobs    = 20,
			    suffix   = 'wwto2l2nu'
			),
		]


########################## WZ  ######################
mc_dict['JC_WZ_Legacy']=[	
		 	cfg(data_name = 'WZTo3LNu',
		 	    path     = '{0}/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/2016_mc_legacy_gen_WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic/190824_060249/'.format(path),
			    nJobs    = 20,
			    suffix   = 'wzto3lnu'
			),
		]

########################## ZZ  ######################
mc_dict['JC_ZZ_Legacy']=[	
		 	cfg(data_name = 'ZZTo2L2Nu',
		 	    path     = '{0}/ZZTo2L2Nu_13TeV_powheg_pythia8/2016_mc_legacy_gen_ZZTo2L2Nu_13TeV_powheg_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2/190824_060811/'.format(path),
			    nJobs    = 20,
			    suffix   = 'zzto2l2nu'
			),
		 	cfg(data_name = 'ZZTo2L2Q',
		 	    path     = '{0}/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/2016_mc_legacy_gen_ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptoti/190824_060617/'.format(path),
			    nJobs    = 20,
			    suffix   = 'zzto2l2q'
			),
		 	cfg(data_name = 'ZZTo4L',
		 	    path     = '{0}/ZZTo4L_13TeV-amcatnloFXFX-pythia8/2016_mc_legacy_gen_ZZTo4L_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-/190824_060434/'.format(path),
			    nJobs    = 20,
			    suffix   = 'zzto4l'
			),
		]

mc_dict['JC_WJets_Legacy'] = [
			    cfg(data_name = 'WJets',
				path     = '{0}/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/2016_mc_legacy_gen_WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_as/190824_062650/'.format(path),
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


data_dict['double_el_legacy'] = [
			 cfg(data_name = 'DoubleEG_2016B_v1',
			     path      = '{0}/DoubleEG/2016_data_legacy_DoubleEG_Run2016B-17Jul2018_ver2-v1/190820_133958'.format(path),
			     nJobs     = 30,
			     suffix    = 'DoubleEG_2016B'
			    ),
			 cfg(data_name = 'DoubleEG_2016C_v1',
			     path      = '{0}/DoubleEG/2016_data_legacy_DoubleEG_Run2016C-17Jul2018-v1/190819_180618'.format(path),
			     nJobs     = 30,
			     suffix    = 'DoubleEG_2016C'
			    ),
			 cfg(data_name = 'DoubleEG_2016D_v1',
			     path      = '{0}/DoubleEG/2016_data_legacy_DoubleEG_Run2016D-17Jul2018-v1/190819_180836'.format(path),
			     nJobs     = 30,
			     suffix    = 'DoubleEG_2016D'
			    ),
			 cfg(data_name = 'DoubleEG_2016E_v1',
			    path      = '{0}/DoubleEG/2016_data_legacy_DoubleEG_Run2016E-17Jul2018-v1/190819_181102'.format(path),
			    nJobs     = 30,
			    suffix    = 'DoubleEG_2016E'
			   ),
			 cfg(data_name = 'DoubleEG_2016F_v1',
			     path      = '{0}/DoubleEG/2016_data_legacy_DoubleEG_Run2016F-17Jul2018-v1/190819_181357'.format(path),
			     nJobs     = 30,
			     suffix    = 'DoubleEG_2016F'
			    ),
			 cfg(data_name = 'DoubleEG_2016G_v1',
			     path      = '{0}/DoubleEG/2016_data_legacy_DoubleEG_Run2016G-17Jul2018-v1/190819_181601'.format(path),
			     nJobs     = 30,
			     suffix    = 'DoubleEG_2016G'
			    ),
			 cfg(data_name = 'DoubleEG_2016H_v1',
			     path      = '{0}/DoubleEG/2016_data_legacy_DoubleEG_Run2016H-17Jul2018-v1/190819_181830'.format(path),
			     nJobs     = 30,
			     suffix    = 'DoubleEG_2016H'
			    ),
			]
path = "/eos/uscms/store/user/corderom/data_legacy_gen_2016"
data_dict['single_el_legacy'] = [
			 cfg(data_name = 'Electron_2016B_v1',
			     path      = '{0}/SingleElectron/2016_data_legacy_Electron_Run2016B-17Jul2018_ver2-v1/190820_133958'.format(path),
			     nJobs     = 30,
			     suffix    = 'Electron_2016B'
			    ),
			 cfg(data_name = 'Electron_2016C_v1',
			     path      = '{0}/SingleElectron/2016_data_legacy_SingleElectron_Run2016C-17Jul2018-v1/191002_063637'.format(path),
			     nJobs     = 30,
			     suffix    = 'Electron_2016C'
			    ),
			 cfg(data_name = 'Electron_2016D_v1',
			     path      = '{0}/SingleElectron/2016_data_legacy_SingleElectron_Run2016D-17Jul2018-v1/191002_063848'.format(path),
			     nJobs     = 30,
			     suffix    = 'Electron_2016D'
			    ),
			 cfg(data_name = 'Electron_2016E_v1',
			    path      = '{0}/SingleElectron/2016_data_legacy_SingleElectron_Run2016E-17Jul2018-v1/191002_064057'.format(path),
			    nJobs     = 30,
			    suffix    = 'Electron_2016E'
			   ),
			 cfg(data_name = 'Electron_2016F_v1',
			     path      = '{0}/SingleElectron/2016_data_legacy_Electron_Run2016F-17Jul2018-v1/190819_181357'.format(path),
			     nJobs     = 30,
			     suffix    = 'Electron_2016F'
			    ),
			 cfg(data_name = 'Electron_2016G_v1',
			     path      = '{0}/SingleElectron/2016_data_legacy_SingleElectron_Run2016G-17Jul2018-v1/191002_064504'.format(path),
			     nJobs     = 30,
			     suffix    = 'Electron_2016G'
			    ),
			 cfg(data_name = 'Electron_2016H_v1',
			     path      = '{0}/SingleElectron/2016_data_legacy_SingleElectron_Run2016H-17Jul2018-v1/191002_064704'.format(path),
			     nJobs     = 30,
			     suffix    = 'Electron_2016H'
			    ),
			]


######################################################################################################################
############################################## SAMPLES 2017 ##########################################################
######################################################################################################################
path = '/eos/uscms/store/user/lpcbacon/corderom/data_2017/DoubleMuon/'


##############  DOUBLE MUON DATA ########################
data_dict['double_mu_2017Reco'] = [
			 cfg(data_name = 'DoubleMuon_2017B_v1',
			     path      = '{0}/2017_data_legacy_DoubleMuon_Run2017B-31Mar2018-v1/190923_203833/'.format(path),
			     nJobs     = 30,
			     suffix    = 'DoubleMuon_2017B'
			    ),
			 cfg(data_name = 'DoubleMuon_2017C_v1',
			     path      = '{0}/2017_data_legacy_DoubleMuon_Run2017C-31Mar2018-v1/190923_204104/'.format(path),
			     nJobs     = 30,
			     suffix    = 'DoubleMuon_2017C'
			    ),
			 cfg(data_name = 'DoubleMuon_2017D_v1',
			     path      = '{0}/2017_data_legacy_DoubleMuon_Run2017D-31Mar2018-v1/190923_204317/'.format(path),
			     nJobs     = 30,
			     suffix    = 'DoubleMuon_2017D'
			    ),
			 cfg(data_name = 'DoubleMuon_2017E_v1',
			    path      = '{0}/2017_data_legacy_DoubleMuon_Run2017E-31Mar2018-v1/190923_204527/'.format(path),
			    nJobs     = 30,
			    suffix    = 'DoubleMuon_2017E'
			   ),
			 cfg(data_name = 'DoubleMuon_2017F_v1',
			     path      = '{0}/2017_data_legacy_DoubleMuon_Run2017F-31Mar2018-v1/190923_204752/'.format(path),
			     nJobs     = 30,
			     suffix    = 'DoubleMuon_2017F'
			    ),
			]
path = "/eos/uscms/store/user/lpcbacon/corderom/data_2017/"
data_dict['single_el_2017Reco'] = [
			 #cfg(data_name = 'Electron_2017B_v1',
			 #    path      = '{0}/SingleElectron/2017_data_legacy_SingleElectron_Run2017B-31Mar2018-v1/191127_043138/'.format(path),
			 #    nJobs     = 30,
			 #    suffix    = 'Electron_2017B'
			 #   ),
			 cfg(data_name = 'Electron_2017C_v1',
			     path      = '{0}/SingleElectron/2017_data_legacy_SingleElectron_Run2017C-31Mar2018-v1/191127_043440/'.format(path),
			     nJobs     = 30,
			     suffix    = 'Electron_2017C'
			    ),
			 #cfg(data_name = 'Electron_2017D_v1',
			 #    path      = '{0}/SingleElectron/2017_data_legacy_SingleElectron_Run2017D-31Mar2018-v1/191127_043843/'.format(path),
			 #    nJobs     = 30,
			 #    suffix    = 'Electron_2017D'
			 #   ),
			 #cfg(data_name = 'Electron_2017E_v1',
			 #   path      = '{0}/SingleElectron/2017_data_legacy_SingleElectron_Run2017E-31Mar2018-v1/191127_044217/'.format(path),
			 #   nJobs     = 30,
			 #   suffix    = 'Electron_2017E'
			 #  ),
			 #cfg(data_name = 'Electron_2017F_v1',
			 #   path      = '{0}/SingleElectron/2017_data_legacy_SingleElectron_Run2017F-31Mar2018-v1/191127_044441/'.format(path),
			 #   nJobs     = 30,
			 #   suffix    = 'Electron_2017F'
			 #  ),
			]



path = "/eos/uscms/store/user/lpcbacon/15/"
data_dict['single_mu_2017Reco'] = [
			 cfg(data_name = 'Muon_2017B_v1',
			     path      = '{0}/SingleMuonRun2017B_17Nov2017_v1'.format(path),
			     nJobs     = 30,
			     suffix    = 'Muon_2017B'
			    ),
			 cfg(data_name = 'Muon_2017C_v1',
			     path      = '{0}/SingleMuonRun2017C_17Nov2017_v1'.format(path),
			     nJobs     = 30,
			     suffix    = 'Muon_2017C'
			    ),
			 cfg(data_name = 'Muon_2017D_v1',
			     path      = '{0}/SingleMuonRun2017D_17Nov2017_v1'.format(path),
			     nJobs     = 30,
			     suffix    = 'Muon_2017D'
			    ),
			 cfg(data_name = 'Muon_2017E_v1',
			    path      = '{0}/SingleMuonRun2017E_17Nov2017_v1'.format(path),
			    nJobs     = 30,
			    suffix    = 'Muon_2017E'
			   ),
			 cfg(data_name = 'Muon_2017F_v1',
			    path      = '{0}/SingleMuonRun2017F_17Nov2017_v1'.format(path),
			    nJobs     = 30,
			    suffix    = 'Muon_2017F'
			   ),
			]



path = "/eos/uscms/store/user/corderom/mc_legacy_gen_2017"
mc_dict['JC_DYJets_2017Reco']=[	
		 	cfg(data_name = 'DYJets',
		 	    path      = '{0}/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_mc_legacy_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017/190925_185738'.format(path),
			    nJobs     = 50,
			    suffix    = 'dyjets'
			),
		]

mc_dict['JC_ZG_ZToLL_2017Reco']=[	
		 	cfg(data_name = 'ZGToLLG',
		 	    #path      = '{0}/ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_mc_legacy_ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_/191018_203100/'.format(path),
			    path      = '/eos/uscms/store/user/lpchzg/corderom/mc_2017/ZGToLLG_01J_LoosePtlPtg_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_mc_legacy_ZGToLLG_01J_LoosePtlPtg_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018/191205_135218/',
			    nJobs     = 50,
			    suffix    = 'zgtollg'
			),
		]

mc_dict['JC_TT_2017Reco']=[	
		 	cfg(data_name = 'TTTo2L2Nu',
		 	    path      = '{0}/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/2017_mc_legacy_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v/191018_204042/0000/'.format(path),
			    nJobs     = 30,
			    suffix    = 'ttto2l2nu'
			),
			]

mc_dict['JC_WJets_2017Reco']=[	
				cfg(data_name = 'WJets',
				    path      = '{0}/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/2017_mc_legacy_WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_reali/191011_040655/'.format(path),
				    nJobs     = 30,
				    suffix    = 'wjets'
				),
				]

	
mc_dict['JC_WZ_2017Reco']=[	
				cfg(data_name = 'WZTo2L2Q',
				    path      = '{0}/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/2017_mc_legacy_WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realis/191011_035512/'.format(path),
				    nJobs     = 30,
				    suffix    = 'wzto2l2q'
				),
				cfg(data_name = 'WZTo3LNu',
				    path      = '{0}/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_mc_legacy_WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realis/191018_203257/'.format(path),
				    nJobs     = 30,
				    suffix    = 'wzto3lnu'
				),
				]


mc_dict['JC_WW_2017Reco']=[	
				cfg(data_name = 'WWTo2L2Nu',
				    path      = '{0}/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/2017_mc_legacy_WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_rea/191011_035302/0000/'.format(path),
				    nJobs     = 10,
				    suffix    = 'wwto2l2nu'
				),
				]

mc_dict['JC_ZZ_2017Reco']=[	
				cfg(data_name = 'ZZTo4L',
				    path      = '{0}/ZZTo4L_TuneCP5_13TeV-amcatnloFXFX-pythia8/2017_mc_legacy_ZZTo4L_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realisti/191018_203812/'.format(path),
				    nJobs     = 10,
				    suffix    = 'zzto4l'
				),
				cfg(data_name = 'ZZTo2L2Nu',
				    path      = '{0}/ZZTo2L2Nu_13TeV_powheg_pythia8/2017_mc_legacy_ZZTo2L2Nu_13TeV_powheg_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/191011_040301/0000/'.format(path),
				    nJobs     = 10,
				    suffix    = 'zzto2l2nu'
				),
				cfg(data_name = 'ZZTo2L2Q',
				    path      = '{0}/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/2017_mc_legacy_ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realis/191018_203518/'.format(path),
				    nJobs     = 10,
				    suffix    = 'zzto2l2q'
				),
				]
######################################################################################################


batch_list = []
batch_list += sum([data_dict[n] for n in data_samples], []) 
#batch_list += sum([mc_dict[n]   for n in   mc_samples], []) 

batch = bm.BatchMaster(config_list = batch_list, 
                       stage_dir   = 'batch',
                       selection   = selection,
                       period      = period,
                       executable  = executable,
                       location    = 'lpc'
                     )
batch.submit_to_batch()

