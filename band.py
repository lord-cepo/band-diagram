import numpy as np
import matplotlib.pyplot as plt
import material as mat
import warnings

EPSILON_0 = 8.8541878128e-14  # C V-1 cm-1
E_CHARGE = 1.602176634e-19    # C

class band_diagram:
    """
    Description:
    ----------------------------------------------------------------------------
    Describes the entire band diagram (energy levels as a function of the
    spatial extension of the device). Bands are bent according to the Anderson's
    rule: vacuum level is continuous; Fermi level is unique at equilibrium.
    Band bending follows the full depletion region approximation, valid (with a
    correction of the order of the thermal energy of 26 meV - not implemented -)
    for bendings far greater than thermal energy. 
    
    Methods:
    ----------------------------------------------------------------------------
    bend()
        does the dirty job of bending bands in both sides of each interface,
        allowing E0 to be continuous.
    plot(display_E0=False, show=True)
        plots Ef, Ec, Ev and E0 if display_E0 is True. show should be set to
        false if we want to insert plt.show() manually (multiple plots in same
        pane)
    
    Attributes:
    ----------------------------------------------------------------------------
    thickness : float
        sum of thicknesses of various layers
    grid : np.array()
        grid of points, from zero to thickness
    levels : dict
        contains the three levels, Ec (conduction), Ev (valence) and E0 (vacuum)
        . Fermi is defined as zero, levels bend equally, mantaining the same
        distance between each others
    layers : list
        copy of the input array of N layers
    interfaces : list[int]
        list of N-1 integers, in which each interface point is given as an index
        of the grid array. 
    Ef : np.array()
        Fermi level, separate from the other levels. It's initialized at zero
        and changed only by applying voltage
    """
    def __init__(self, layers, n_points=1000) -> None:
        """
        Description:
        ------------------------------------------------------------------------
        From the layers list, this constructor creates builds flat Ec, Ev, E0
        levels, keeping the Fermi level constant.

        Arguments:
        ------------------------------------------------------------------------
        layers : Iterable[layer]
            list of layer instances. Example: [layer(thickness1, material1),
            layer(thickness2, material2)]. Layer is inside material module.
        n_points : int (optional)
            Number of points in the actual plot. Defaults to 1000.
        """
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
        
    
    def bend(self) -> None:
        """
        Generates the right amount of bending to cancel out discontinuities on
        E0. Formulas are found in docs. Once device is initializated, first 
        apply a voltage (optional) and then type device.bend()
        """
        # electronic density is obtained for each layer, giving huge density to
        # metals (1e23 cm-3)
        densities = []
        n_points = self.Ef.size
        for material in [layer.material for layer in self.layers]:
            density = 0
            if type(material) is mat.metal:
                density = 1e23 # el. densities of noble metals (cm-3)
            elif type(material) is mat.semiconductor:
                if material.is_p_doped():
                    density = material.p
                else:
                    density = material.n
            else:
                warnings.warn("Material is not metal nor semic, behaviour\
                    not implemented yet")
            densities.append(density)
        
        # in every interface, bending depends on "left" and "right" material
        # properties. These are just full depletion calculations.
        epsilons = [layer.material.epsilon for layer in self.layers]
        thicknesses = [layer.thickness for layer in self.layers]
        for i, position in enumerate(self.interfaces):
            overflow = []
            delta_V = self.levels['E0'][position]-self.levels['E0'][position-1]
            
            a = np.sqrt(2*EPSILON_0*epsilons[i]*epsilons[i+1] /
                        (   densities[i]*epsilons[i] +
                            densities[i+1]*epsilons[i+1]) *
                        np.abs(delta_V)/E_CHARGE )
            depletion_left = a*np.sqrt(densities[i+1]/densities[i])*1e4 # um
            depletion_right = a*np.sqrt(densities[i]/densities[i+1])*1e4 # um
            
            # overflow check (warns if depletion > thickness)
            if depletion_left > thicknesses[i]:
                overflow.append(i)
                start = round(position - thicknesses[i]/self.thickness*n_points)
            else:
                start = round(position - depletion_left/self.thickness*n_points)
                
            if depletion_right > thicknesses[i+1]:
                overflow.append(i+1)
                end = round(position + thicknesses[i+1]/self.thickness*n_points)
            else:
                end = round(position + depletion_right/self.thickness*n_points)          
            if len(overflow) != 0:
                warnings.warn(f"Depletion region larger than thickness of \
                    materials: {overflow}")
            
            # left bending start in position - left_dep and ends in position +
            # right_dep. See calculations for details
            c = np.sign(delta_V)*1e-8*E_CHARGE/(2*EPSILON_0)
            bending = np.concatenate((
                np.zeros(start),
                c*densities[i]/epsilons[i]*(self.grid[start:position] - position*self.thickness/n_points + depletion_left)**2,
                -c*densities[i+1]/epsilons[i+1]*(self.grid[position:end] - position*self.thickness/n_points - depletion_right)**2,
                np.zeros(n_points-end)
            ))
            
            # Ec, Ev, E0 are mantained parallel
            for k in self.levels.keys():
                self.levels[k] += bending
        
    
    
    def plot(self, display_E0=False, show=True):
        """
        Plotter of band_diagram. Plots conduction band Ec, valence band Ev,
        Fermi level Ef. 

        Arguments:
        ------------------------------------------------------------------------
        display_E0 : bool (optional)
            Displays vacuum level. Defaults to False.
        show : bool (optional)
            If False, plt.show() has to be performed manually. Useful if we want
            multiple plots on one figure. Defaults to True.
        """
        for key, level in self.levels.items():
            if key != 'E0' or display_E0:
                plt.plot(self.grid, level, label=key)
        plt.plot(self.grid, self.Ef)
        plt.legend()
        if show:
            plt.tight_layout()
            plt.show()
        
        
    def apply_voltage(self, volt, ith_layer=None) -> None:
        """
        Translate rigidly all the levels of the specified layer by minus voltage

        Arguments:
        ------------------------------------------------------------------------
        volt : float
            Voltage applied (in Volts). If positive, the energy shift will be
            negative: energy = (-e) * voltage
        ith_layer : int (optional) 
            Specifies which layer will receive voltage. Every other layer should
            be considered as grounded. Naturally, we are under the (not so 
            realistic) hypothesis that voltage drops, and hence jumps in Ef,
            can exist only at the interfaces. If not specified, the rightmost
            layer will receive the voltage. Defaults to None
        """
        # equals to N_layers - 1
        n_of_interfaces = len(self.interfaces)
        # equals to n_points inside __init__
        n_points = self.Ef.size
        if ith_layer is None:
            # if not specified, voltage is applied to last layer (N_layers-1)th
            ith_layer = n_of_interfaces
        
        # first and last elements are added to avoid overflow
        positions = [0] + self.interfaces + [n_points]
        
        # creation of array of volts that will be added to each level
        volts = np.concatenate((
            np.zeros(positions[ith_layer]),
            volt*np.ones(positions[ith_layer+1]-positions[ith_layer]),
            np.zeros(n_points-positions[ith_layer+1])
        ))
        
        # volts SUBTRACTED from Ef, because energy = (-e) * voltage
        self.Ef -= volts
        # volts subtracted from Ec, Ev, E0 
        for k in self.levels.keys():
            self.levels[k] -= volts
