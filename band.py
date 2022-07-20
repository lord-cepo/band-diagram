import numpy as np
import matplotlib.pyplot as plt # type: ignore
import material
import warnings

from collections.abc import Iterable

EPSILON_0 = 8.8541878128e-14  # C V-1 cm-1
E_CHARGE = 1.602176634e-19    # C

class band_diagram:
    """
    Description
    -----------
    Describes the entire band diagram (energy levels as a function of the
    spatial extension of the device). Bands are bent according to the Anderson's
    rule: vacuum level is continuous; Fermi level is unique at equilibrium.
    Band bending follows the full depletion region approximation, valid (with a
    correction of the order of the thermal energy of 26 meV - not implemented -)
    for bendings far greater than thermal energy. 
    
    Methods
    -------
    reset()
        just resets the device
    bend(fermi=False)
        does the dirty job of bending bands in both sides of each interface,
        allowing E0 to be continuous.
    plot(title=None, display_E0=False, show=True, display_eh=False)
        plots Ef, Ec, Ev and E0 if display_E0 is True. show should be set to
        false if we want to insert plt.show() manually (multiple plots in same
        pane)
    apply_voltage(volt, index=None)
        moves all four levels of -volt of specified layer. Layers are numbered
        from 0 and from left to right, by default voltage is applied to 
        rightmost layer 
    
    Attributes
    ----------
    thickness : float
        sum of thicknesses of various layers (um)
    grid : np.ndarray
        grid of points, from zero to thickness (um)
    levels : dict[str, numpy.ndarray]
        contains the three levels, Ec (conduction), Ev (valence) and E0 (vacuum)
        . Fermi is defined as zero, levels bend equally, mantaining the same
        distance between each others. Levels are in eV.
    layers : List
        copy of the input array of N layers
    interfaces : list[int]
        list of N-1 integers, in which each interface point is given as an index
        of the grid array. 
    Ef : np.ndarray
        Fermi level, separate from the other levels. It's initialized at zero
        and changed only by applying voltage
    """
    def __init__(self, 
        layers: Iterable[material.layer],
        n_points: int = 1000
        ) -> None:
        """
        From the layers list, this constructor creates builds flat Ec, Ev, E0
        levels, keeping the Fermi level constant, passing from reset() method

        Parameters
        ----------
        layers : Iterable[material.layer]
            list of layer instances. Example: [layer(thickness1, material1),
            layer(thickness2, material2)]. Layer is inside material module.
        n_points : int, optional
            Number of points in the actual plot, by default 1000.
        """
        if layers != []:
            self.reset(layers, n_points)
        else:
            print('Device cannot be initialized if layers is []')
            
    
    def reset(self,
        layers: Iterable[material.layer] = [],
        n_points: int = 1000
        ) -> None:
        """
        User can call it without parameters (after having initialized the
        device) to reset to initial conditions, otherwise is called by __init__

        Parameters
        ----------
        layers : Iterable[material.layer], optional
            list of layer instances. Example: [layer(thickness1, material1),
            layer(thickness2, material2)]. Layer is inside material module,
            by default []
        n_points : int, optional
            Number of points in the actual plot, by default 1000.
        """
        
        # executed if reset is called by __init__
        if layers != []:
            self.layers = layers
        
        self.thickness = sum([layer.thickness for layer in self.layers]) 
        self.grid = np.linspace(0, self.thickness, n_points)
        self.levels = {
            'Ec': np.zeros(n_points),
            'Ev': np.zeros(n_points),
            'E0': np.zeros(n_points)
        }
        partial_thickness = 0.
        self.interfaces = []
        for layer in self.layers:
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
        del self.interfaces[0] # remove -1:0-th interface
    
    # TODO: design with bias between ends
    def bend(self,
        fermi: bool = False
        ) -> None:
        """
        Generates the right amount of bending to cancel out discontinuities on
        E0. Formulas are found in docs. Once device is initializated, first 
        apply a voltage (optional) and then type device.bend()
        
        Parameters
        ----------
        fermi : bool, optional
            If True, fermi level is also bended to remove discontinuities, 
            occurring in presence of applied voltage.
            This procedure is maybe a bit unphysical, use it only to cancel 
            discontinuities in populations during plotting, by default False
        """
        self.__bend(fermi=False)
        if fermi:
            self.__bend(fermi=True)
             
    
    def __bend(self,
        fermi: bool
        ) -> None:
        """
        Private method, does the hard work on bending
        """
        # electronic density is obtained for each layer, giving huge density to
        # metals (1e23 cm-3)
        densities = []
        n_points = self.Ef.size
        for mat in [layer.material for layer in self.layers]:
            density = 0.0
            if type(mat) is material.metal:
                density = 1e23 # el. densities of noble metals (cm-3)
            elif type(mat) is material.semiconductor:
                if mat.is_p_doped():
                    density = mat.p
                else:
                    density = mat.n
            else:
                warnings.warn("Material is not metal nor semic, behaviour not implemented yet")
            densities.append(density)
        
        # in every interface, bending depends on "left" and "right" material
        # properties. These are just full depletion calculations.
        epsilons = [layer.material.epsilon for layer in self.layers]
        thicknesses = [layer.thickness for layer in self.layers]
        for i, position in enumerate(self.interfaces):
            overflow = []
            if fermi:
                delta_V = self.Ef[position+1]-self.Ef[position-1]
            else:
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
                warnings.warn(f"Depletion region larger than thickness of materials: {overflow}")
            
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
            if fermi:
                self.Ef += bending
            else:
                for k in self.levels.keys():
                    self.levels[k] += bending
    
    
    def plot(self,
        title: str = None, 
        display_E0: bool = False, 
        show: bool = True,
        display_eh: bool = False,
        smearing: float = 0.5
        ) -> None:
        """
        Plotter of band_diagram. Plots conduction band Ec, valence band Ev,
        Fermi level Ef. 

        Parameters
        ----------
        title : str, optional
            title of the plot, by default no title
        display_E0 : bool, optional
            Displays vacuum level. Defaults to False.
        show : bool, optional
            If False, plt.show() has to be performed manually. Useful if we want
            multiple plots on one figure, by default True
        display_eh: bool, optional
            If True, displays density of electrons (red) and holes (blue),
            according to fictious smearing, by default False
        smearing: float, optional
            At 300 K should be 0.026 eV, but can be set to higher values for
            representative purposes of the two populations, by default 0.5
        """
        if title is not None:
            plt.title(title)
        if display_eh:
            n_points = self.Ef.size
            list_of_levels = [self.levels['Ec'], self.levels['Ev'], self.Ef]
            E_min = min([min(E) for E in list_of_levels])
            E_max = max([max(E) for E in list_of_levels])
            x, y = np.meshgrid(self.grid, np.linspace(E_min-0.1, E_max+0.1, 1000))
            
            z_el = np.exp(-(y-self.Ef[(x/self.thickness*(n_points-1)).round().astype(int)])/smearing)
            z_el = np.where(z_el>1, 0, z_el)
            for i in range(n_points):
                for j, _y in enumerate(y):
                    if _y[0] < self.levels['Ec'][i]:
                        z_el[j][i] = 0.
                        
            z_min, z_max = z_el.min(), z_el.max()
            plt.imshow(z_el, cmap='Reds', vmin=z_min, vmax=z_max,
                extent=[x.min(), x.max(), y.min(), y.max()],
                interpolation='bilinear', origin='lower', alpha=1.0)
            
            z_ho = np.exp((y-self.Ef[(x/self.thickness*(n_points-1)).round().astype(int)])/smearing)
            z_ho = np.where(z_ho>1, 0, z_ho)
            for i in range(n_points):
                for j, _y in enumerate(y):
                    if _y[0] > self.levels['Ev'][i]:
                        z_ho[j][i] = 0.

            plt.imshow(z_ho, cmap='Blues', vmin=z_min, vmax=z_max,
                extent=[x.min(), x.max(), y.min(), y.max()],
                interpolation='bilinear', origin='lower', alpha=0.70)
        
        plt.plot(self.grid, self.levels['Ec'], label='Ec', color='maroon')
        plt.plot(self.grid, self.levels['Ev'], label='Ev', color='darkslategrey')
        if display_E0:
            plt.plot(self.grid, self.levels['E0'], label='E0', color='darkmagenta')
        plt.plot(self.grid, self.Ef, label='Ef', linestyle='--', color='black')
        plt.ylabel('energy ($eV$)')
        plt.xlabel('layer thickness ($\mu m$)')
        plt.legend()
        if show:
            plt.tight_layout()
            # saves to plots folder (executing some_device it saves everything)
            # plt.savefig('plots/'+ title + '.png', bbox_inches='tight', dpi=200)
            plt.show()
        
        
    def apply_voltage(self, 
        volt: float, 
        index: int = None
        ) -> None:
        """
        Translate rigidly all the levels of the specified layer by minus voltage
        
        Parameters
        ----------
        volt : float
            Voltage applied (in Volts). If positive, the energy shift will be 
            negative: energy = (-e) * voltage
        index : int, optional
            Specifies which layer will receive voltage. Every other layer should
            be considered as grounded. Naturally, we are under the (not so 
            realistic) hypothesis that voltage drops, and hence jumps in Ef,
            can exist only at the interfaces. If not specified, the rightmost
            layer will receive the voltage, by default None
        """
        # equals to N_layers - 1
        n_of_interfaces = len(self.interfaces)
        # equals to n_points inside __init__
        n_points = self.Ef.size
        if index is None:
            # if not specified, voltage is applied to last layer (N_layers-1)th
            index = n_of_interfaces
        
        # first and last elements are added to avoid overflow
        positions = [0] + self.interfaces + [n_points]
        
        # creation of array of volts that will be added to each level
        volts = np.concatenate((
            np.zeros(positions[index]),
            volt*np.ones(positions[index+1]-positions[index]),
            np.zeros(n_points-positions[index+1])
        ))
        
        # volts SUBTRACTED from Ef, because energy = (-e) * voltage
        self.Ef -= volts
        # volts subtracted from Ec, Ev, E0 
        for k in self.levels.keys():
            self.levels[k] -= volts
