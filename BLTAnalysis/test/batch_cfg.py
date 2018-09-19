#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg        = bm.JobConfig
executable = 'execBatch.sh'
selection  = 'double'
period     = '2012'

#data_samples = ['single_mu', 'single_el']
data_samples = ['double_mu', 'double_el']
mc_samples   = ['zjets', 'ttbar', 'diboson']



''' 
Set job configurations.  
'''

# DATA #
data_dict = {}

path       = '/tthome/share/bacon/production/04'
data_dict['single_mu'] = \
[
    cfg(data_name = 'muon_2012A',
        path     = '{0}/SingleMu_2012A-22Jan2013'.format(path),
        nJobs    = 1,
        suffix   = 'muon_2012A'
        ),
    cfg(data_name = 'muon_2012B',
        path     = '{0}/SingleMu_2012B-22Jan2013'.format(path),
        nJobs    = 7,
        suffix   = 'muon_2012B'
        ),
    cfg(data_name = 'muon_2012C',
        path     = '{0}/SingleMu_2012C-22Jan2013'.format(path),
        nJobs    = 11,
        suffix   = 'muon_2012C'
        ),
    cfg(data_name = 'muon_2012D',
        path     = '{0}/SingleMu_2012D-22Jan2013'.format(path),
        nJobs    = 13,
        suffix   = 'muon_2012D'
        ),
]

data_dict['single_el'] = \
[
    cfg(data_name = 'electron_2012A',
        path     = '{0}/SingleElectron_2012A-22Jan2013'.format(path),
        nJobs    = 3,
        suffix   = 'electron_2012A'
        ),
    cfg(data_name = 'electron_2012B',
        path     = '{0}/SingleElectron_2012B-22Jan2013'.format(path),
        nJobs    = 14,
        suffix   = 'electron_2012B'
        ),
    cfg(data_name = 'electron_2012C',
        path     = '{0}/SingleElectron_2012C-22Jan2013'.format(path),
        nJobs    = 20,
        suffix   = 'electron_2012C'
        ),
    cfg(data_name = 'electron_2012D',
        path     = '{0}/SingleElectron_2012D-22Jan2013'.format(path),
        nJobs    = 20,
        suffix   = 'electron_2012D'
        ),
]

data_dict['double_mu'] = \
[
    cfg(data_name = 'muon_2012A',
        path     = '{0}/DoubleMuParked_2012A-22Jan2013'.format(path),
        nJobs    = 1,
        suffix   = 'muon_2012A'
        ),
    cfg(data_name = 'muon_2012B',
        path     = '{0}/DoubleMuParked_2012B-22Jan2013'.format(path),
        nJobs    = 7,
        suffix   = 'muon_2012B'
        ),
    cfg(data_name = 'muon_2012C',
        path     = '{0}/DoubleMuParked_2012C-22Jan2013'.format(path),
        nJobs    = 11,
        suffix   = 'muon_2012C'
        ),
    cfg(data_name = 'muon_2012D',
        path     = '{0}/DoubleMuParked_2012D-22Jan2013'.format(path),
        nJobs    = 13,
        suffix   = 'muon_2012D'
        ),
]

data_dict['double_el'] = \
[
    cfg(data_name = 'electron_2012A',
        path     = '{0}/2012_Data_Multicrab_DoubleElectron_Run2012A-22Jan2013-v1/161006_183851/0000'.format(path),
        nJobs    = 3,
        suffix   = 'electron_2012A'
        ),
    cfg(data_name = 'electron_2012B',
        path     = '{0}/2012_Data_Multicrab_DoubleElectron_Run2012B-22Jan2013-v1/161006_184140/0000'.format(path),
        nJobs    = 14,
        suffix   = 'electron_2012B'
        ),
    cfg(data_name = 'electron_2012C',
        path     = '{0}/2012_Data_Multicrab_DoubleElectron_Run2012C-22Jan2013-v1/161006_184501/0000'.format(path),
        nJobs    = 20,
        suffix   = 'electron_2012C'
        ),
    cfg(data_name = 'electron_2012D',
        path     = '{0}/2012_Data_Multicrab_DoubleElectron_Run2012D-22Jan2013-v1/161006_184923/0000'.format(path),
        nJobs    = 20,
        suffix   = 'electron_2012D'
        ),
]



# MONTE CARLO #
mc_dict = {}

mc_dict['zjets'] = \
[
    cfg(data_name = 'DYJetsToLL_M-50',
        path     = '{0}/Summer12_DYJetsToLL_M-50_TuneZ2Star'.format(path),
        nJobs    = 50,
        suffix   = 'zjets_m-50'
       ),
    cfg(data_name = 'DYJetsToLL_M-10to50',
        path     = '{0}/Summer12_DYJetsToLL_M-10To50_TuneZ2Star'.format(path),
        nJobs    = 10,
        suffix   = 'zjets_m-10to50'
       ),
]

mc_dict['ttbar'] = \
[
    cfg(data_name = 'TTJets_FullLept',
        path     = '{0}/Summer12_TTJets_FullLeptMGDecays'.format(path),
        nJobs    = 6,
        suffix   = 'ttbar_full'
       ),
    cfg(data_name = 'TTJets_SemiLept',
        path     = '{0}/Summer12_TTJets_SemiLeptMGDecays'.format(path),
        nJobs    = 13,
        suffix   = 'ttbar_semi'
       ),
    cfg(data_name = 'TTZJets',
        path     = '{0}/Summer12_TTZJets'.format(path),
        nJobs    = 10,
        suffix   = 'ttz'
       ),
]

mc_dict['diboson'] = \
[
    cfg(data_name = 'WW',
        path     = '{0}/Summer12_WW_TuneZ2star'.format(path),
        nJobs    = 10,
        suffix   = 'ww'
       ),
    cfg(data_name = 'WZJetsTo2L2Q',
        path     = '{0}/Summer12_WZJetsTo2L2Q_TuneZ2star'.format(path),
        nJobs    = 10,
        suffix   = 'wz_2l2q'
       ),
    cfg(data_name = 'WZJetsTo3LNu',
        path     = '{0}/Summer12_WZJetsTo3LNu_TuneZ2'.format(path),
        nJobs    = 10,
        suffix   = 'wz_3lnu'
       ),
    cfg(data_name = 'ZZJetsTo2L2Nu',
        path     = '{0}/Summer12_ZZJetsTo2L2Nu_TuneZ2star'.format(path),
        nJobs    = 10,
        suffix   = 'zz_2l2nu'
       ),
    cfg(data_name = 'ZZJetsTo2L2Q',
        path     = '{0}/Summer12_ZZJetsTo2L2Q_TuneZ2star'.format(path),
        nJobs    = 10,
        suffix   = 'zz_2l2q'
       ),
    cfg(data_name = 'ZZ',
        path     = '{0}/Summer12_ZZ_TuneZ2star'.format(path),
        nJobs    = 10,
        suffix   = 'zz'
       ),
]



batch_list = []
batch_list += sum([data_dict[n] for n in data_samples], []) 
batch_list += sum([mc_dict[n] for n in mc_samples], []) 

batch = bm.BatchMaster(config_list = batch_list, 
                       stage_dir   = 'batch',
                       selection   = selection,
                       period      = period,
                       executable  = executable,
                       location    = 'nut3'
                     )
batch.submit_to_batch()
