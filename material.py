import numpy as np
import warnings

'''
doping in cm-3
'''

E_MASS = 5.6933648125e-16 # eV cm-2 s2
kT = 0.026   # eV @ 300 K
PLANCK_H = 4.135667e-15  # eV s

class layer:
    def __init__(self, thickness, material):
        self.thickness = thickness
        self.material = material

class material:
    def __init__(self, epsilon, Ec, Ev, work):
        self.epsilon = epsilon
        self.levels = {
            'Ec': Ec,
            'Ev': Ev,
            'E0': Ec + work
        }

class metal(material):
    def __init__(self, work_function):
        super().__init__(1e10, 0, 0, work_function) # epsilon = 1e10 = np.inf

class semiconductor(material):
    def __init__(self, Eg, electron_affinity, epsilon=1., doping_type=None, doping=None, effective_m_e=1, effective_m_h=None):
        if effective_m_h is None:
            effective_m_h = effective_m_e
        
        Nc = 2*(2*np.pi*E_MASS*effective_m_e*kT/PLANCK_H**2)**1.5
        Nv = 2*(2*np.pi*E_MASS*effective_m_h*kT/PLANCK_H**2)**1.5
        ni = np.sqrt(Nc*Nv*np.exp(-Eg/kT))

        if doping is not None and doping < ni:
            warnings.warn("You have chosen a too low doping (< intrinsic doping), units are cm-3. Assuming no doping")

        if doping_type == 'n':
            self.n = doping
            self.p = ni**2/doping
        elif doping_type == 'p':
            self.p = doping
            self.n = ni**2/doping
        else:
            if doping_type is not None:
                warnings.warn("You have selected a doping type different than 'n' or 'p'. Assuming no doping")
                doping_type = None
            self.n = ni
            self.p = ni
            
        Ec = -kT*np.log(self.n/Nc)
        Ev = kT*np.log(self.p/Nv)
        super().__init__(epsilon, Ec, Ev, electron_affinity)
    
    def is_p_doped(self):
        return self.p > self.n            
    
    def is_n_doped(self):
        return self.n > self.p