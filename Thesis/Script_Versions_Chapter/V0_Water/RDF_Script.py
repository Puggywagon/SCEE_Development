import gromacs
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
################################################################################
class RDF_Script(object):
    def __init__(self):
        pass
################################################################################
    def make_ndx(self,system_title):
        input_1='del 2 \n del 1 \n del 0 \n a O* \n name 0 Oxygen \n a H* \n name 1 Hydrogen \n q'
        tmp = gromacs.tools.Make_ndx(f=f'{system_title}_QMMM_md3',
                                o='index.ndx',
                                input=input_1)
        chk, makendx_output, stderr = tmp.run()
################################################################################
    def rdf_OO(self,system_title):
        input_str='0\n0\n'
        tmp = gromacs.tools.Rdf(f=f'{system_title}_QMMM_md3', 
                                s=f'{system_title}_QMMM_md3',
                                n='index.ndx',
                                excl=[],
                                bin=0.004,
                                o='OO_RDF.xvg',
                                input=input_str)
        chk, rdf_output, stderr = tmp.run()
################################################################################
    def rdf_OH(self,system_title):
        input_str='0\n 1\n'
        tmp = gromacs.tools.Rdf(f=f'{system_title}_QMMM_md3', 
                                s=f'{system_title}_QMMM_md3',
                                n='index.ndx',
                                bin=0.004,
                                excl=[],
                                o='OH_RDF.xvg',
                                input=input_str)
        chk, rdf_output, stderr = tmp.run()
################################################################################
    def rdf_HH(self,system_title):
        input_str='1\n 1\n'
        tmp = gromacs.tools.Rdf(f=f'{system_title}_QMMM_md3', 
                                s=f'{system_title}_QMMM_md3',
                                n='index.ndx',
                                bin=0.004,
                                excl=[],
                                o='HH_RDF.xvg',
                                input=input_str)
        chk, rdf_output, stderr = tmp.run()
##############################################################################################

                                
