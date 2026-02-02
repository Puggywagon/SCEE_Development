#!/usr/bin/python3

import numpy as np

import Particle
import SimBox


################################################################################
################################################################################
################################################################################
class Molecule(object):
    
    def __init__(self, ID):
        self.ID = ID
    
        angle = 150.0 * np.pi/180.0
        self.HBcutoff = 0.24
        self.HBdot = np.cos(angle)
        self.particle_list = []
        self.acceptor_list = []
        self.donor_list = []
        self.SimBox = None

        
################################################################################
################################################################################
    def is_acceptor(self, mj):

        for D in self.donor_list:
            Di, D0 = D
            for Aj in mj.acceptor_list:
                dpos = self.box.get_dposition(Di.position, Aj.position)
                dpos0 = self.box.get_dposition(Di.position, D0.position)
        dot=0.0; d=0.0; d0=0.0;
        for i in range(len(dpos)):
            dot += dpos[i]*dpos0[i]
            d += dpos[i]*dpos[i]
            d0 += dpos0[i]*dpos0[i]
        d = np.sqrt(d)
        d0 = np.sqrt(d0)
        dot /= d*d0
        if (d < self.HBcutoff and dot < self.HBdot):
            return True
        
        return False

    
    def is_donor(self, mj):
        for Ai in self.acceptor_list:
          for D in mj.donor_list:
              Dj, D0 = D
              dpos = self.box.get_dposition(Dj.position, Ai.position)
              dpos0 = self.box.get_dposition(Dj.position, D0.position)
              dot=0.0; d=0.0; d0=0.0;
          for i in range(len(dpos)):
              dot += dpos[i]*dpos0[i]
              d += dpos[i]*dpos[i]
              d0 += dpos0[i]*dpos0[i]
          d = np.sqrt(d); d0 = np.sqrt(d0); dot /= d*d0;
          if (d < self.HBcutoff and dot < self.HBdot):
              return True;        
        return False

    
    def add_particle(self, p):
        self.particle_list.append(p)

        
    def add_acceptor(self, p):
        self.acceptor_list.append(p)
            
    def add_donor(self, p, p0):
        self.donor_list.append( (p, p0) )

            
################################################################################
################################################################################
################################################################################
