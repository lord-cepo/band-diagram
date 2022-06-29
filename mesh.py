import numpy as np
import matplotlib.pyplot as plt
import material

class mesh:
    def __init__(self, layers, n_points=1000) -> None:
        thickness = sum([layer.thickness for layer in layers]) 
        self.grid = np.linspace(0, thickness, n_points)
        self.levels = {
            'Ec': np.zeros(n_points),
            'Ev': np.zeros(n_points),
            'E0': np.zeros(n_points)
        }
        partial_thickness = 0.
        for layer in layers:
            start = round(partial_thickness/thickness*n_points)
            end = round((partial_thickness+layer.thickness)/thickness*n_points)
            for k in self.levels.keys():
                self.levels[k] += np.concatenate(np.zeros(start), 
                                                 layer.material.levels[k]*np.ones(end-start), 
                                                 np.zeros(n_points-end)
                                                 )
            partial_thickness += layer.thickness
        self.levels['Ef'] = np.zeros()
        
    
    def plot(self, display_E0=False):
        for key, level in self.levels.items():
            plt.plot(self.grid, level, label=key)
        plt.show()
        
        
    def apply_voltage(volt, right=True):
        pass
        