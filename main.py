import material
import mesh
import matplotlib.pyplot as plt

Si_n = material.semiconductor(Eg=1.12, electron_affinity=4.05, epsilon=11.7, doping_type='n', doping=1e15)
Si_p = material.semiconductor(Eg=1.12, electron_affinity=4.05, epsilon=11.7, doping_type='p', doping=1e15)
GaAs = material.semiconductor(Eg=1.424, 
                              electron_affinity=4.07, 
                              epsilon=12.9, 
                              effective_m_e=0.063, 
                              effective_m_h=0.51,
                              doping_type='n', doping=1e15)
x = 0.8
AlGaAs = material.semiconductor(Eg=1.424+1.247*x if x<0.45 else 1.9+0.125*x+0.143*x**2,
                                electron_affinity=4.07-1.1*x if x<0.45 else 3.64-0.14*x,
                                epsilon=12.9-2.84*x,
                                effective_m_e=(0.063+0.083*x)**1.5 if x<0.41 else (0.85-0.14*x)**1.5,
                                effective_m_h=(0.51+0.25*x)**1.5 if x<0.41 else (0.85-0.14*x)**1.5,
                                doping_type='n',
                                doping=GaAs.n
)
Al = material.metal(4.1)

def example_diode():
    list_layers = [(10, Si_n), (10, Si_p)]
    layers = [ material.layer(*el) for el in list_layers]
    diode = mesh.mesh(layers=layers)
    # diode.plot()
    # diode.apply_voltage(volt=-5.)
    diode.bend()
    diode.plot(display_E0=True)

def example_npn():
    list_layers = [(5, Si_n), (5, Si_p), (5, Si_n)]
    layers = [ material.layer(*el) for el in list_layers]
    npn = mesh.mesh(layers=layers)
    
    # npn.bend()
    # npn.plot(display_E0=False) # display_E0=True
    npn.apply_voltage(volt=1, ith_layer=1)
    npn.bend()
    npn.apply_voltage(volt=1, ith_layer=1)
    npn.bend()
    npn.plot(display_E0=True) # display_E0=True
    
def example_schottky_contact():
    list_layers = [(1, Al), (5, Si_n)]
    layers = [ material.layer(*el) for el in list_layers]
    schottky_contact = mesh.mesh(layers=layers)
    
    schottky_contact.bend()
    schottky_contact.plot()

def example_heterojunction():
    list_layers = [(1, GaAs), (1, AlGaAs)]
    layers = [material.layer(*el) for el in list_layers]
    htj = mesh.mesh(layers)
    htj.bend()
    htj.plot(display_E0=True)

# example_diode()
# example_npn()
# example_schottky_contact()
example_heterojunction()