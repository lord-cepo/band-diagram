import numpy as np
import matplotlib.pyplot as plt
import material as mat

EPSILON_0 = 8.8541878128e-14  # C V-1 cm-1
E_CHARGE = 1.602176634e-19    # C


class mesh:
    def __init__(self, layers, n_points=1000):
        self.thickness = sum([layer.thickness for layer in layers]) 
        self.applied_voltage = 0.
        self.grid = np.linspace(0, self.thickness, n_points)
        self.levels = {
            'Ec': np.zeros(n_points),
            'Ev': np.zeros(n_points),
            'E0': np.zeros(n_points)
        }
        partial_thickness = 0.
        self.interfaces = []
        for layer in layers:
            start = round(partial_thickness/self.thickness*n_points)
            self.interfaces.append(start)
            end = round((partial_thickness+layer.thickness)/self.thickness*n_points)
            for k in self.levels.keys():
                self.levels[k] += np.concatenate((np.zeros(start), 
                                                 layer.material.levels[k]*np.ones(end-start), 
                                                 np.zeros(n_points-end)
                ))
            partial_thickness += layer.thickness
        self.Ef = np.zeros(n_points) # Fermi level is common at 0 applied voltage
        self.layers = layers
        del self.interfaces[0] # remove -1:0-th interface
        
    
    def bend(self):
        densities = []
        n_points = self.Ef.size
        for material in [layer.material for layer in self.layers]:
            n = 0
            if type(material) is mat.metal:
                n = 1e23 # e densities of noble metals (cm-3)
            if type(material) is mat.semiconductor:
                if material.is_p_doped():
                    n = material.p
                else:
                    n = material.n
            densities.append(n)
        
        for i, position in enumerate(self.interfaces):
            left_overflow = False
            right_overflow = False
            
            V_bi = self.levels['E0'][position]-self.levels['E0'][position-1]
            a = np.sqrt(2*EPSILON_0*self.layers[i].material.epsilon*self.layers[i+1].material.epsilon /
                        (   densities[i]*self.layers[i].material.epsilon +
                            densities[i+1]*self.layers[i+1].material.epsilon) *
                        np.abs((V_bi-self.applied_voltage))/E_CHARGE ) * 1000 # cm to um
            depletion_left = a*np.sqrt(densities[i+1]/densities[i])
            depletion_right = a*np.sqrt(densities[i]/densities[i+1])
            
            if depletion_left > self.layers[i].thickness:
                left_overflow = True
                start = round(position - self.layers[i].thickness/self.thickness*n_points)
            else:
                start = round(position - depletion_left/self.thickness*n_points)
                
            if depletion_right > self.layers[i+1].thickness:
                right_overflow = True
                end = round(position + self.layers[i+1].thickness/self.thickness*n_points)
            else:
                end = round(position + depletion_right/self.thickness*n_points)
            
            for k in self.levels.keys():
                self.levels[k] += np.concatenate((np.zeros(start),
                    E_CHARGE*densities[i]/(2*EPSILON_0*self.layers[i].material.epsilon)*(self.grid[start:position]-position + depletion_left)**2,
                    E_CHARGE*densities[i+1]/(2*EPSILON_0*self.layers[i+1].material.epsilon)*(self.grid[end:position]-position - depletion_right)**2,
                    np.zeros(n_points-end)
                ))
                
                
    
    def plot(self, display_E0=False):
        for key, level in self.levels.items():
            if key != 'E0' or display_E0:
                plt.plot(self.grid, level, label=key)
        plt.plot(self.grid, self.Ef)
        plt.legend()
        plt.tight_layout()
        plt.show()
        
        
    def apply_voltage(volt, right=True):
        pass
        