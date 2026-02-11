#!/usr/bin/python3

import os
import subprocess


################################################################################
################################################################################
################################################################################
class SCEE_Gaussian(object):

    def __init__(self):

        self.param = {}
        self.nproc = 8
        self.basis = 'aug-cc-pvtz'
        self.method = 'b3lyp'

        self.rundir = './opt_b3lypaug-cc-pvtz_spb3lyp/'
        self.g09root='/home/mjksill/Research/PhD/MacPherson/SCEE/Gaussian/g09_pgi/'
        self.GAUSS_SCRDIR = '/home/mjksill/idrive/ProcessModelling/mjksill/scratch/'


################################################################################
################################################################################
    def __setitem__(self, key, value):
        self.param[key] = value

        
    def __getitem__(self, key):
        return self.param[key]

    
################################################################################
################################################################################
################################################################################

