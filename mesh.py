import numpy as np
import matplotlib.pyplot as plt
import material as mat

EPSILON_0 = 8.8541878128e-14  # C V-1 cm-1
E_CHARGE = 1.602176634e-19    # C


class mesh:
    def __init__(self, layers, n_points=1000):
        self.thickness = sum([layer.thickness for layer in layers]) 
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
        self.applied_voltage = np.zeros(len(self.interfaces))
        
    
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
            epsilon_left = self.layers[i].material.epsilon
            epsilon_right = self.layers[i+1].material.epsilon
            thickness_left = self.layers[i].thickness
            thickness_right = self.layers[i+1].thickness
            delta_V = self.levels['E0'][position]-self.levels['E0'][position-1]
            
            a = np.sqrt(2*EPSILON_0*epsilon_left*epsilon_right /
                        (   densities[i]*epsilon_left +
                            densities[i+1]*epsilon_right) *
                        np.abs(delta_V)/E_CHARGE )
            depletion_left = a*np.sqrt(densities[i+1]/densities[i])*1e4 # um
            depletion_right = a*np.sqrt(densities[i]/densities[i+1])*1e4 # um
            
            if depletion_left > thickness_left:
                left_overflow = True
                start = round(position - thickness_left/self.thickness*n_points)
            else:
                start = round(position - depletion_left/self.thickness*n_points)
                
            if depletion_right > thickness_right:
                right_overflow = True
                end = round(position + thickness_right/self.thickness*n_points)
            else:
                end = round(position + depletion_right/self.thickness*n_points)
            
            for k in self.levels.keys():
                self.levels[k] += np.concatenate((
                    np.zeros(start),
                    np.sign(delta_V)*densities[i]*1e-8*E_CHARGE/(2*EPSILON_0*epsilon_left)*((self.grid[start:position] - position*self.thickness/n_points + depletion_left))**2,
                    -np.sign(delta_V)*densities[i+1]*1e-8*E_CHARGE/(2*EPSILON_0*epsilon_right)*((self.grid[position:end] - position*self.thickness/n_points - depletion_right))**2,
                    np.zeros(n_points-end)
                ))
        
    
    
    def plot(self, display_E0=False):
        for key, level in self.levels.items():
            if key != 'E0' or display_E0:
                plt.plot(self.grid, level, label=key)
        plt.plot(self.grid, self.Ef)
        plt.legend()
        # plt.tight_layout()
        # plt.show()
        
        
    def apply_voltage(self, volt, ith_layer=None):
        n_of_interfaces = len(self.interfaces)
        n_points = self.Ef.size
        if ith_layer is None:
            ith_layer = n_of_interfaces
        
        positions = [0] + self.interfaces + [n_points]
        arr = np.concatenate((
            np.zeros(positions[ith_layer]),
            volt*np.ones(positions[ith_layer+1]-positions[ith_layer]),
            np.zeros(n_points-positions[ith_layer+1])
        ))
        
        self.Ef -= arr
        for k in self.levels.keys():
            self.levels[k] -= arr
