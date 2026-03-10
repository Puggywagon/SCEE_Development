#!/usr/bin/python3

import numpy as np


################################################################################
################################################################################
################################################################################
class SimBox(object):
    def __init__(self, x, y, z):
        self.box = np.array([x, y, z])
        
        angle = 150.0 * np.pi/180.0
        self.HBcutoff = 0.24
        self.HBdot = np.cos(angle)


    def get_dposition(self, vec1, vec2):
        tmp = []
        for k in range(len(vec1)):
            dx = vec1[k] - vec2[k]
            if (dx >   0.5*self.box[k]): dx -= self.box[k]
            if (dx <= -0.5*self.box[k]): dx += self.box[k]
            tmp.append(dx) 
        dpos = np.array(tmp)
        return dpos
        

################################################################################
################################################################################
    def is_acceptor(self, m0, mj):

        for D in m0.donor_list:
            Di, D0 = D
            for Aj in mj.acceptor_list:
                dpos = m0.box.get_dposition(Di.position, Aj.position)
                dpos0 = m0.box.get_dposition(Di.position, D0.position)
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

    
    def is_donor(self, m0, mj):
        for Ai in m0.acceptor_list:
          for D in mj.donor_list:
              Dj, D0 = D
              dpos = m0.box.get_dposition(Dj.position, Ai.position)
              dpos0 = m0.box.get_dposition(Dj.position, D0.position)
              dot=0.0; d=0.0; d0=0.0;
          for i in range(len(dpos)):
              dot += dpos[i]*dpos0[i]
              d += dpos[i]*dpos[i]
              d0 += dpos0[i]*dpos0[i]
          d = np.sqrt(d); d0 = np.sqrt(d0); dot /= d*d0;
          if (d < self.HBcutoff and dot < self.HBdot):
              return True;        
        return False

