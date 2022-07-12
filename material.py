import numpy as np
import warnings

E_MASS = 5.6933648125e-16 # eV cm-2 s2
kT = 0.026   # eV @ 300 K
PLANCK_H = 4.135667e-15  # eV s

class layer:
    """
    Description:
    ----------------------------------------------------------------------------
    Once materials are defined, layer is only used to couple material type
    with thickness. An ordered list of layers (left to right) is used to 
    initialize band_diagram object (top level class)
    
    Attributes:
    ----------------------------------------------------------------------------
    thickness : float
        thickness of layer in micrometer (um)
    material : material
        material type. Some general use materials are already defined in
        tutorials
    """
    def __init__(self, thickness, material):
        """
        Description:
        ------------------------------------------------------------------------
        Once materials are defined, layer is only used to couple material type
        with thickness. An ordered list of layers (left to right) is used to 
        initialize mesh object (top level class)
        
        Arguments:
        ------------------------------------------------------------------------
        thickness : float
            thickness of layer in micrometer (um)
        material : material
            material type. Some general use materials are already defined in
            tutorials
        """
        self.thickness = thickness
        self.material = material

class material:
    """
    Description:
    ----------------------------------------------------------------------------
    General use class, use it if you want to create an artificial material,
    choosing directly the energy levels. Otherwise, use "semiconductor" or
    "metal" constructor for correct initialization.
    
    material Attributes:
    ----------------------------------------------------------------------------
    epsilon : float
        DC dielectric constant, used inside poisson equations. Large dielectric
        means low band bending.
    levels : dict
        Contains conduction and valence levels ('Ec', 'Ev') and vacuum level E0.
        Fermi Energy is implicitely equal to zero.
    """
    def __init__(self, epsilon, Ec, Ev, work):
        """
        The use of this constructor is not indicated. Use semiconductor() or
        metal() instead.

        Arguments:
        ------------------------------------------------------------------------
        epsilon : float
            Dc dielectric constant
        Ec : float
            Conduction band
        Ev : float
            Valence band
        work : float
            Generalization of work_function and electron_affinity
        """
        self.epsilon = epsilon
        self.levels = {
            'Ec': Ec,
            'Ev': Ev,
            'E0': Ec + work
        }

class metal(material):
    """
    Description:
    ----------------------------------------------------------------------------
    Child of material class. Conduction and valence bands coincide with
    Fermi level, DC dielectric constant goes to infinity (a large number is
    chosen because (np.inf+ float)/np.inf = nan instead of one)
    """
    def __init__(self, work_function):
        """
        Conduction and valence bands coincide with Fermi level, DC dielectric 
        constant is a very large number.

        Arguments:
        ------------------------------------------------------------------------
        work_function : float
            Difference between vacuum level (E0) and Fermi level (Ef). It is
            generally expected to be positive.
        """
        super().__init__(1e30, 0, 0, work_function) # epsilon = 1e30 = np.inf

class semiconductor(material):
    """
    Description:
    ----------------------------------------------------------------------------
    Child of material class. Calculates effective densities in conduction and 
    valence bands (Nc and Nv) and intrinsic density (ni). Sets n and p according 
    to doping type, (strong doping hypothesis, if doping < intrinsic density it 
    throws a warning), calculates conduction and valence levels and calls the 
    parent constructor.
    
    semiconductor Attributes:
    ----------------------------------------------------------------------------
    n : float
        equals to doping if n-type, this parameter is set so that mass action
        law is respected (np = ni^2). Doping is always in cm-3
    p : float
        equals to doping if p-type, this parameter is set so that mass action
        law is respected (np = ni^2). Doping is always in cm-3
    
    material Attributes:
    ----------------------------------------------------------------------------
    epsilon : float
        DC dielectric constant, used inside poisson equations. Large dielectric
        means low band bending.
    levels : dict
        Contains conduction and valence levels ('Ec', 'Ev') and vacuum level E0.
        Fermi Energy is implicitely equal to zero.
    """
    def __init__(self, Eg, electron_affinity, epsilon=1., doping_type=None, 
                 doping=None, effective_m_e=1, effective_m_h=None):
        """
        Description:
        ------------------------------------------------------------------------
        Calculates effective densities in conduction and valence bands (Nc and 
        Nv) and intrinsic density (ni). Sets n and p according to doping type, 
        (strong doping hypothesis, if doping < intrinsic density it throws a 
        warning), calculates conduction and valence levels and calls the parent 
        constructor.
        
        Arguments:
        ------------------------------------------------------------------------
        Eg : float
            Energy gap of the semiconductor. Should be > 0
        
        electron_affinity : float
            Electron affinity, defined as the difference between vacuum level E0 
            and conduction band Ec. Generally, bulk affinities aredifferent from 
            the isolated atom ones, found in periodic table. Should be > 0. 
        
        epsilon : float (optional)
            Relative DC dielectric constant (permittivity). Defaults to 1.
        
        doping_type : str (optional)
            Can be 'n' or 'p'. Defaults to None.
        
        doping : float (optional)
            Doping density of the type specified in doping_type, expressed in 
            cm-3. Only one (strong) doping type supported. Defaults to None.
            
        effective_m_e : int (optional)
            effective electron mass. It's used to calculate density of states. 
            If n_v is the number of conduction band minima valleys in whole band
            structure, effective_m_e = (n_v^2 mx my mz)^(1/3), where the three 
            effective masses around principal axes. This software is just for 
            plotting, rigour in the choice of this parameter can be easily 
            avoided without any consequences. Defaults to 1.
        
        effective_m_h : float (optional)
            effective hole mass. Same as electron case. If None, it's set equal 
            to effective_m_e. Defaults to None.
        """
        if effective_m_h is None:
            effective_m_h = effective_m_e
        
        Nc = 2*(2*np.pi*E_MASS*effective_m_e*kT/PLANCK_H**2)**1.5
        Nv = 2*(2*np.pi*E_MASS*effective_m_h*kT/PLANCK_H**2)**1.5
        ni = np.sqrt(Nc*Nv*np.exp(-Eg/kT))

        if doping is not None and doping < ni:
            warnings.warn("You have chosen a too low doping (< intrinsic \
                doping), units are cm-3. Assuming no doping")
            doping_type = None
        
        if doping_type == 'n':
            self.n = doping
            self.p = ni**2/doping
        elif doping_type == 'p':
            self.p = doping
            self.n = ni**2/doping
        else:
            if doping_type is not None:
                warnings.warn("You have selected a doping type other than 'n' \
                    or 'p'. Assuming no doping")
                doping_type = None
            self.n = ni
            self.p = ni
            
        Ec = kT*np.log(Nc/self.n)
        Ev = -kT*np.log(Nv/self.p)
        super().__init__(epsilon, Ec, Ev, electron_affinity)
    
    def is_p_doped(self):
        return self.p > self.n            
    
    def is_n_doped(self):
        return self.n > self.p