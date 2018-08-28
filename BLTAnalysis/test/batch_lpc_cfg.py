#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg        = bm.JobConfig
selection  = 'emu'
period     = '2016'
executable = 'execBatch.sh'
location   = 'lpc'

#data_samples = ['single_mu', 'single_el']
#mc_samples   = ['selected', 'zjets']

#data_samples = ['double_mu', 'double_eg']
#mc_samples   = ['zjets_12d', 'ttbar_12d', 'diboson_12d']

data_samples = ['double_eg']


''' 
Set job configurations.  
'''

# DATA #
data_dict = {}

path = '/eos/uscms/store/group/lpcbacon/12d'
data_dict['double_mu'] = \
[
    cfg(data_name = 'muon_2016B_v2',
        path     = '{0}/DoubleMuon_Run2016B-03Feb2017_ver2-v2'.format(path),
        nJobs    = 51,
        suffix   = 'muon_2016B'
       ),
    cfg(data_name = 'muon_2016C_v1',
        path     = '{0}/DoubleMuon_Run2016C-03Feb2017-v1'.format(path),
        nJobs    = 67,
        suffix   = 'muon_2016C'
       ),
    cfg(data_name = 'muon_2016D_v1',
        path     = '{0}/DoubleMuon_Run2016D-03Feb2017-v1'.format(path),
        nJobs    = 51,
        suffix   = 'muon_2016D'
       ),
    cfg(data_name = 'muon_2016E_v1',
        path     = '{0}/DoubleMuon_Run2016E-03Feb2017-v1'.format(path),
        nJobs    = 51,
        suffix   = 'muon_2016E'
       ),
    cfg(data_name = 'muon_2016F_v1',
        path     = '{0}/DoubleMuon_Run2016F-03Feb2017-v1'.format(path),
        nJobs    = 51,
        suffix   = 'muon_2016F'
       ),
    cfg(data_name = 'muon_2016G_v1',
        path     = '{0}/DoubleMuon_Run2016G-03Feb2017-v1'.format(path),
        nJobs    = 51,
        suffix   = 'muon_2016G'
       ),
    cfg(data_name = 'muon_2016H_v2',
        path     = '{0}/DoubleMuon_Run2016H-03Feb2017_ver2-v1'.format(path),
        nJobs    = 51,
        suffix   = 'muon_2016H'
       ),
    cfg(data_name = 'muon_2016H_v3',
        path     = '{0}/DoubleMuon_Run2016H-03Feb2017_ver3-v1'.format(path),
        nJobs    = 28,
        suffix   = 'muon_2016H'
       ),
]

data_dict['double_eg'] = \
[
#   cfg(data_name = 'electron_2016B_v2',
#       path     = '{0}/DoubleEG_Run2016B-03Feb2017_ver2-v2'.format(path),
#       nJobs    = 51,
#       suffix   = 'electron_2016B'
#       ),
#   cfg(data_name = 'electron_2016C_v1',
#       path     = '{0}/DoubleEG_Run2016C-03Feb2017-v1'.format(path),
#       nJobs    = 36,
#       suffix   = 'electron_2016C'
#       ),
#   cfg(data_name = 'electron_2016D_v1',
#       path     = '{0}/DoubleEG_Run2016D-03Feb2017-v1'.format(path),
#       nJobs    = 51,
#       suffix   = 'electron_2016D'
#       ),
#   cfg(data_name = 'electron_2016E_v1',
#       path     = '{0}/DoubleEG_Run2016E-03Feb2017-v1'.format(path),
#       nJobs    = 51,
#       suffix   = 'electron_2016E'
#       ),
#   cfg(data_name = 'electron_2016F_v1',
#       path     = '{0}/DoubleEG_Run2016F-03Feb2017-v1'.format(path),
#       nJobs    = 51,
#       suffix   = 'electron_2016F'
#       ),
#   cfg(data_name = 'electron_2016G_v1',
#       path     = '{0}/DoubleEG_Run2016G-03Feb2017-v1'.format(path),
#       nJobs    = 51,
#       suffix   = 'electron_2016G'
#       ),
    cfg(data_name = 'electron_2016H_v2',
        path     = '{0}/DoubleEG_Run2016H-03Feb2017_ver2-v1'.format(path),
        nJobs    = 51,
        suffix   = 'electron_2016H'
        ),
    cfg(data_name = 'electron_2016H_v3',
        path     = '{0}/DoubleEG_Run2016H-03Feb2017_ver3-v1'.format(path),
        nJobs    = 29,
        suffix   = 'electron_2016H'
        ),
]

path = '/eos/uscms/store/group/lpcbacon/12a'
data_dict['single_mu'] = \
[
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
   cfg(data_name = 'muon_2016C',
       path      = '{0}/SingleMuon_Run2016C-03Feb2017-v1'.format(path),
       nJobs     = 30,
       suffix    = 'muon_2016C'
       ),
   cfg(data_name = 'muon_2016D',
       path      = '{0}/SingleMuon_Run2016D-03Feb2017-v1'.format(path),
       nJobs     = 30,
       suffix    = 'muon_2016D'
       ),
   cfg(data_name = 'muon_2016E',
       path      = '{0}/SingleMuon_Run2016E-03Feb2017-v1'.format(path),
       nJobs     = 30,
       suffix    = 'muon_2016E'
       ),
   cfg(data_name = 'muon_2016F',
       path      = '{0}/SingleMuon_Run2016F-03Feb2017-v1'.format(path),
       nJobs     = 30,
       suffix    = 'muon_2016F'
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

data_dict['single_el'] = \
[
   cfg(data_name = 'electron_2016B_v1',
       nJobs    = 30,
       path     = '{0}/SingleElectron_Run2016B-03Feb2017_ver1-v1'.format(path),
       suffix   = 'electron_2016B'
       ),
   cfg(data_name = 'electron_2016B_v2',
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
   cfg(data_name = 'electron_2016H_v2',
       path     = '{0}/SingleElectron_Run2016H-03Feb2017_ver2-v1'.format(path),
       nJobs    = 30,
       suffix   = 'electron_2016H'
       ),
   cfg(data_name = 'electron_2016H_v3',
       path     = '{0}/SingleElectron_Run2016H-03Feb2017_ver3-v1'.format(path),
       nJobs    = 30,
       suffix   = 'electron_2016H'
       ),
]

data_dict['mueg'] = \
[
    cfg(data_name = 'mueg_2016B_v1',
        path     = '{0}/MuonEG_Run2016B-03Feb2017_ver1-v1'.format(path),
        nJobs    = 30,
        suffix   = 'mueg_2016B'
        ),
    cfg(data_name = 'mueg_2016B_v2',
        path     = '{0}/MuonEG_Run2016B-03Feb2017_ver2-v2'.format(path),
        nJobs    = 30,
        suffix   = 'mueg_2016B'
        ),
    cfg(data_name = 'mueg_2016C',
        path     = '{0}/MuonEG_Run2016C-03Feb2017-v1'.format(path),
        nJobs    = 30,
        suffix   = 'mueg_2016C'
        ),
    cfg(data_name = 'mueg_2016D',
        path     = '{0}/MuonEG_Run2016D-03Feb2017-v1'.format(path),
        nJobs    = 30,
        suffix   = 'mueg_2016D'
        ),
    cfg(data_name = 'mueg_2016E',
        path     = '{0}/MuonEG_Run2016E-03Feb2017-v1'.format(path),
        nJobs    = 30,
        suffix   = 'mueg_2016E'
        ),
    cfg(data_name = 'mueg_2016F',
        path     = '{0}/MuonEG_Run2016F-03Feb2017-v1'.format(path),
        nJobs    = 30,
        suffix   = 'mueg_2016F'
        ),
    cfg(data_name = 'mueg_2016G',
        path     = '{0}/MuonEG_Run2016G-03Feb2017-v1'.format(path),
        nJobs    = 30,
        suffix   = 'mueg_2016G'
        ),
    cfg(data_name = 'mueg_2016H_v2',
        path     = '{0}/MuonEG_Run2016H-03Feb2017_ver2-v1'.format(path),
        nJobs    = 30,
        suffix   = 'mueg_2016H'
        ),
    cfg(data_name = 'mueg_2016H_v3',
        path     = '{0}/MuonEG_Run2016H-03Feb2017_ver3-v1'.format(path),
        nJobs    = 30,
        suffix   = 'mueg_2016H'
        ),
]


# MONTE CARLO #
mc_dict = {}
path = '/eos/uscms/store/group/lpcbacon/12d'
mc_dict['zjets_12d'] = \
[
    cfg(data_name = 'DYJetsToLL_M-50',
        path     = '{0}/DYJetsToLL_M-50_madgraph'.format(path),
        nJobs    = 50,
        suffix   = 'zjets_50'
       ),
    cfg(data_name = 'DYJetsToLL_M-10to50',
        path     = '{0}/DYJetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'zjets_10to50'
       ),
    cfg(data_name = 'DY1JetsToLL_M-50',
        path     = '{0}/DY1JetsToLL_M-50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z1jets_50'
       ),
    cfg(data_name = 'DY1JetsToLL_M-10to50',
        path     = '{0}/DY1JetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z1jets_10to50'
       ),
    cfg(data_name = 'DY2JetsToLL_M-50',
        path     = '{0}/DY2JetsToLL_M-50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z2jets_50'
       ),
    cfg(data_name = 'DY2JetsToLL_M-10to50',
        path     = '{0}/DY2JetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z2jets_10to50'
       ),
    cfg(data_name = 'DY3JetsToLL_M-50',
        path     = '{0}/DY3JetsToLL_M-50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z3jets_50'
       ),
    cfg(data_name = 'DY3JetsToLL_M-10to50',
        path     = '{0}/DY3JetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z3jets_10to50'
       ),
    cfg(data_name = 'DY4JetsToLL_M-50',
        path     = '{0}/DY4JetsToLL_M-50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z4jets_50'
       ),
    cfg(data_name = 'DY4JetsToLL_M-10to50',
        path     = '{0}/DY4JetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z4jets_10to50'
       ),
]

mc_dict['ttbar_12d'] = \
[
        cfg(data_name = 'ttbar_inclusive',
            path     = '{0}/TT_powheg'.format(path),
            nJobs    = 51,
            suffix   = 'ttbar'
            ),
]

mc_dict['diboson_12d'] = \
[
    cfg(data_name = 'WW',
        path     = '{0}/WWTo2L2Nu_powheg'.format(path),
        nJobs    = 10,
        suffix   = 'ww'
       ),
    cfg(data_name = 'WZJetsTo2L2Q',
        path     = '{0}/WZTo2L2Q_amcatnlo'.format(path),
        nJobs    = 51,
        suffix   = 'wz_2l2q'
        ),
    cfg(data_name = 'WZJetsTo3LNu',
        path     = '{0}/WZTo3LNu_powheg'.format(path),
        nJobs    = 27,
        suffix   = 'wz_3lnu'
        ),
    cfg(data_name = 'ZZJetsTo2L2Nu',
        path     = '{0}/ZZTo2L2Nu_powheg'.format(path),
        nJobs    = 10,
        suffix   = 'zz_2l2nu'
       ),
    cfg(data_name = 'ZZJetsTo2L2Q',
        path     = '{0}/ZZTo2L2Q_amcatnlo'.format(path),
        nJobs    = 10,
        suffix   = 'zz_2l2q'
        ),
    cfg(data_name = 'ZZJetsTo4L',
        path     = '{0}/ZZTo4L_amcatnlo'.format(path),
        nJobs    = 10,
        suffix   = 'zz_4l'
        ),
]

mc_dict['higgs_12d'] = \
[
    cfg(data_name = 'GluGluHToZZTo4L',
        path     = '{0}/GluGlu_HToZZTo4L'.format(path),
        nJobs    = 50,
        suffix   = 'ggH_zz_4l'
       ),
    cfg(data_name = 'VBFHToZZTo4L',
        path     = '{0}/VBF_HToZZTo4L'.format(path),
        nJobs    = 50,
        suffix   = 'H_zz_4l'
       ),
]


path = '/eos/uscms/store/group/lpcbacon/12'
mc_dict['selected'] = \
[
    cfg(data_name = 'DYJetsToLL_M-10to50',
        path     = '{0}/Summer16_DYJetsToLL_M-10to50_amcatnlo'.format(path),
        nJobs    = 50,
        suffix   = 'dy_m-10to50'
        ),
    cfg(data_name = 'DYJetsToLL_M-50',
        path     = '{0}/Summer16_DYJetsToLL_M-50_amcatnlo'.format(path),
        nJobs    = 50,
        suffix   = 'dy_m-50'
        ),
    cfg(data_name = 'GluGluHToZZTo4L',
        path     = '{0}/Summer16_GluGluHToZZTo4L_powheg'.format(path),
        nJobs    = 10,
        suffix   = 'ggH_zz_4l'
        ),
    cfg(data_name = 'HToZZTo4L',
        path     = '{0}/Summer16_VBF_HToZZTo4L_powheg'.format(path),
        nJobs    = 10,
        suffix   = 'H_zz_4l'
        ),
    cfg(data_name = 'ttbar_inclusive',
        path     = '{0}/Summer16_TT_powheg'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar'
        ),
    cfg(data_name = 'TTZToLLNuNu_M-10',
        path     = '{0}/Summer16_TTZToLLNuNu_M-10'.format(path),
        nJobs    = 10,
        suffix   = 'ttz_2l2nu'
        ),
    cfg(data_name = 'WW',
        path     = '{0}/Summer16_WWTo2L2Nu_powheg'.format(path),
        nJobs    = 10,
        suffix   = 'ww_2l2nu'
        ),
    cfg(data_name = 'WZJetsTo3LNu',
        path     = '{0}/Summer16_WZTo3LNu_powheg'.format(path),
        nJobs    = 10,
        suffix   = 'wz_3lnu'
        ),
    cfg(data_name = 'ZZJetsTo4L',
        path     = '{0}/Summer16_ZZto4L_amcatnlo'.format(path),
        nJobs    = 10,
        suffix   = 'zz_4l'
        ),
]

path = '/eos/uscms/store/group/lpcbacon/12d'
mc_dict['zjets'] = \
[
#   cfg(data_name = 'DYJetsToLL_M-50',
#       path     = '{0}/Summer16_DYJetsToLL_M-50_amcatnlo'.format(path),
#       nJobs    = 50,
#       suffix   = 'zjets_m-50'
#       ),
#   cfg(data_name = 'DYJetsToLL_M-10to50',
#       path     = '{0}/Summer16_DYJetsToLL_M-10to50_amcatnlo'.format(path),
#       nJobs    = 10,
#       suffix   = 'zjets_m-10to50'
#       ),
#   cfg(data_name = 'DYJetsToLL_M-50',
#       path     = '{0}/DYJetsToLL_M-50_madgraph'.format(path),
#       nJobs    = 50,
#       suffix   = 'zjets_m-50'
#      ),
#   cfg(data_name = 'DYJetsToLL_M-10to50',
#       path     = '{0}/DYJetsToLL_M-10to50_madgraph'.format(path),
#       nJobs    = 10,
#       suffix   = 'zjets_m-10to50'
#      ),
    cfg(data_name = 'DY1JetsToLL_M-50',
        path     = '{0}/DY1JetsToLL_M-50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z1jets_m-50'
       ),
    cfg(data_name = 'DY1JetsToLL_M-10to50',
        path     = '{0}/DY1JetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z1jets_m-10to50'
       ),
    cfg(data_name = 'DY2JetsToLL_M-50',
        path     = '{0}/DY2JetsToLL_M-50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z2jets_m-50'
       ),
    cfg(data_name = 'DY2JetsToLL_M-10to50',
        path     = '{0}/DY2JetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z2jets_m-10to50'
       ),
    cfg(data_name = 'DY3JetsToLL_M-50',
        path     = '{0}/DY3JetsToLL_M-50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z3jets_m-50'
       ),
    cfg(data_name = 'DY3JetsToLL_M-10to50',
        path     = '{0}/DY3JetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z3jets_m-10to50'
       ),
    cfg(data_name = 'DY4JetsToLL_M-50',
        path     = '{0}/DY4JetsToLL_M-50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z4jets_m-50'
       ),
    cfg(data_name = 'DY4JetsToLL_M-10to50',
        path     = '{0}/DY4JetsToLL_M-10to50_madgraph'.format(path),
        nJobs    = 10,
        suffix   = 'z4jets_m-10to50'
       ),
]

mc_dict['wjets'] = \
[
    cfg(data_name = 'W1JetsToLNu',
        path     = '{0}/Summer16_W1JetsToLNu'.format(path),
        nJobs    = 40,
        suffix   = 'w1jets'
        ),
    cfg(data_name = 'W2JetsToLNu',
        path     = '{0}/Summer16_W2JetsToLNu'.format(path),
        nJobs    = 40,
        suffix   = 'w2jets'
        ),
    cfg(data_name = 'W3JetsToLNu',
        path     = '{0}/Summer16_W3JetsToLNu'.format(path),
        nJobs    = 40,
        suffix   = 'w3jets'
        ),
    cfg(data_name = 'W4JetsToLNu',
        path     = '{0}/Summer16_W4JetsToLNu'.format(path),
        nJobs    = 40,
        suffix   = 'w4jets'
        ),
]

mc_dict['qcd'] = \
[
    cfg(data_name = 'QCD_HT50to100',
        path     = '{0}/Summer16_QCD_HT50to100'.format(path),
        nJobs    = 10,
        suffix   = 'qcd_ht50to100'
        ),
    cfg(data_name = 'QCD_HT100to200',
        path     = '{0}/Summer16_QCD_HT100to200'.format(path),
        nJobs    = 10,
        suffix   = 'qcd_ht100to200'
        ),
    cfg(data_name = 'QCD_HT200to300',
        path     = '{0}/Summer16_QCD_HT200to300'.format(path),
        nJobs    = 10,
        suffix   = 'qcd_ht200to300'
        ),
    cfg(data_name = 'QCD_HT300to500',
        path     = '{0}/Summer16_QCD_HT300to500'.format(path),
        nJobs    = 10,
        suffix   = 'qcd_ht300to500'
        ),
    cfg(data_name = 'QCD_HT500to700',
        path     = '{0}/Summer16_QCD_HT500to700'.format(path),
        nJobs    = 10,
        suffix   = 'qcd_ht500to700'
        ),
    cfg(data_name = 'QCD_HT700to1000',
        path     = '{0}/Summer16_QCD_HT700to1000'.format(path),
        nJobs    = 10,
        suffix   = 'qcd_ht700to1000'
        ),
    cfg(data_name = 'QCD_HT1000to1500',
        path     = '{0}/Summer16_QCD_HT1000to1500'.format(path),
        nJobs    = 10,
        suffix   = 'qcd_ht1000to1500'
        ),
    cfg(data_name = 'QCD_HT1500to2000',
        path     = '{0}/Summer16_QCD_HT1500to2000'.format(path),
        nJobs    = 10,
        suffix   = 'qcd_ht1500to2000'
        ),
    cfg(data_name = 'QCD_HT2000toInf',
        path     = '{0}/Summer16_QCD_HT2000toInf'.format(path),
        nJobs    = 10,
        suffix   = 'qcd_ht2000'
        ),
]

mc_dict['ttbar'] = \
[
    cfg(data_name = 'ttbar_inclusive',
        path     = '{0}/Summer16_TT_powheg'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_inclusive'
        ),
    cfg(data_name = 'ttbar_inclusive_tunedown',
        path     = '{0}/Summer16_TT_powheg_TuneCUETP8M2T4down'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_inclusive_down'
       ),
    cfg(data_name = 'ttbar_inclusive_tuneup',
        path     = '{0}/Summer16_TT_powheg_TuneCUETP8M2T4up'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_inclusive_up'
       ),
    cfg(data_name = 'ttbar_inclusive_isrdown',
        path     = '{0}/Summer16_TT_powheg_isrdown'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_inclusive_isrdown'
       ),
    cfg(data_name = 'ttbar_inclusive_isrup',
        path     = '{0}/Summer16_TT_powheg_isrup'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_inclusive_isrup'
       ),
    cfg(data_name = 'ttbar_inclusive_fsrdown',
        path     = '{0}/Summer16_TT_powheg_fsrdown'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_inclusive_fsrdown'
       ),
    cfg(data_name = 'ttbar_inclusive_fsrup',
        path     = '{0}/Summer16_TT_powheg_fsrup'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_inclusive_fsrup'
       ),
    cfg(data_name = 'ttbar_inclusive_hdampdown',
        path     = '{0}/Summer16_TT_powheg_hdampDOWN'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_inclusive_hdampdown'
       ),
    cfg(data_name = 'ttbar_inclusive_hdampup',
        path     = '{0}/Summer16_TT_powheg_hdampUP'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_inclusive_hdampup'
        ),
    cfg(data_name = 'ttbar_leptonic',
        path     = '{0}/Summer16_TTTo2L2Nu_powheg'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_lep'
       ),
    cfg(data_name = 'ttbar_semileptonic',
        path     = '{0}/Summer16_TTToSemilepton_powheg'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_semilep'
       ),
    cfg(data_name = 'ttbar_leptonic',
        path     = '{0}/Summer16_TTJets_DiLept_madgraph'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar_lep'
       ),
]

mc_dict['t'] = \
[
    cfg(data_name = 'T_s-channel',
        path     = '{0}/Summer16_ST_s-channel_4f_leptonDecays_amcatnlo'.format(path),
        nJobs    = 10,
        suffix   = 't_s'
       ),
    cfg(data_name = 'Tbar_s-channel',
        path     = '{0}/'.format(path),
        nJobs    = 10,
        suffix   = 'tbar_s'
       ),
    cfg(data_name = 'T_t-channel',
        path     = '{0}/Summer16_ST_t_channel_top_4f_inclusiveDecays_TuneCUETP8M2T4'.format(path),
        nJobs    = 10,
        suffix   = 't_t'
       ),
    cfg(data_name = 'Tbar_t-channel',
        path     = '{0}/Summer16_ST_t_channel_antitop_4f_inclusiveDecays_TuneCUETP8M2T4'.format(path),
        nJobs    = 10,
        suffix   = 'tbar_t'
       ),
    cfg(data_name = 'T_tW-channel',
        path     = '{0}/Summer16_ST_tW_top_5f_inclusiveDecays_TuneCUETP8M2T4'.format(path),
        nJobs    = 10,
        suffix   = 't_tw'
        ),
    cfg(data_name = 'Tbar_tW-channel',
        path     = '{0}/Summer16_ST_tW_antitop_5f_inclusiveDecays_TuneCUETP8M2T4'.format(path),
        nJobs    = 10,
        suffix   = 'tbar_tw'
        ),
]

mc_dict['diboson'] = \
[
    cfg(data_name = 'WW',
        path     = '{0}/Summer16_WWTo2L2Nu_powheg'.format(path),
        nJobs    = 10,
        suffix   = 'ww'
        ),
    cfg(data_name = 'WZJetsTo2L2Q',
        path     = '{0}/Summer16_WZTo2L2Q_amcatnlo'.format(path),
        nJobs    = 10,
        suffix   = 'wz_2l2q'
        ),
    cfg(data_name = 'WZJetsTo3LNu',
        path     = '{0}/Summer16_WZTo3LNu_powheg'.format(path),
        nJobs    = 10,
        suffix   = 'wz_3lnu'
        ),
    cfg(data_name = 'ZZJetsTo2L2Nu',
        path     = '{0}/Summer16_ZZTo2L2Nu_powheg'.format(path),
        nJobs    = 10,
        suffix   = 'zz_2l2nu'
        ),
    cfg(data_name = 'ZZJetsTo2L2Q',
        path     = '{0}/Summer16_ZZTo2L2Q_amcatnlo'.format(path),
        nJobs    = 10,
        suffix   = 'zz_2l2q'
        ),
    cfg(data_name = 'ZZJetsTo4L',
        path     = '{0}/Summer16_ZZto4L_amcatnlo'.format(path),
        nJobs    = 10,
        suffix   = 'zz_4l'
        ),
    cfg(data_name = 'GluGluHToZZTo4L',
        path     = '{0}/Summer16_GluGluHToZZTo4L_powheg'.format(path),
        nJobs    = 10,
        suffix   = 'ggH_zz_4l'
        ),
]



batch_list = []
batch_list += sum([data_dict[n] for n in data_samples], []) 
#batch_list += sum([mc_dict[n] for n in mc_samples], []) 

batch = bm.BatchMaster(config_list = batch_list, 
                      stage_dir   = 'batch',
                      selection  = selection,
                      period     = period,
                      executable = executable,
                      location   = 'lpc'
                     )
batch.submit_to_batch()

