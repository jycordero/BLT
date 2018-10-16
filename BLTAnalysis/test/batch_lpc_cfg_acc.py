#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg        = bm.JobConfig
selection  = 'gen'
period     = '2017'
executable = 'execBatchAcc.sh'
location   = 'lpc'

mc_samples   = ['signal']



''' 
Set job configurations.  
'''

# MONTE CARLO #
mc_dict = {}

path = '/eos/uscms/store/user/jrainbol/Bacon/'
mc_dict['signal'] = \
[
    cfg(data_name = 'ZZTo4L',
        path     = '{0}/Fall17_ZZTo4L_powheg'.format(path),
        nJobs    = 20,
        suffix   = 'zz_4l'
        ),
]



batch_list = []
batch_list += sum([mc_dict[n] for n in mc_samples], []) 

batch = bm.BatchMaster(config_list = batch_list, 
                      stage_dir   = 'batch',
                      selection  = selection,
                      period     = period,
                      executable = executable,
                      location   = 'lpc'
                     )
batch.submit_to_batch()

