import material
import mesh

Si_n = material.semiconductor(Eg=1.12, electron_affinity=4.05, epsilon=11.7, doping_type='n', doping=1e16)
Si_p = material.semiconductor(Eg=1.12, electron_affinity=4.05, epsilon=11.7, doping_type='p', doping=1e16)
Al = material.metal(4.1)

def example_diode():
    list_layers = [(1, Si_p), (1, Si_n)]
    layers = [ material.layer(*el) for el in list_layers]
    diode = mesh.mesh(layers=layers)
    diode.bend()
    diode.plot(display_E0=False) # display_E0=True

def example_schottky_contact():
    list_layers = [(1, Al), (1, Si_n)]
    layers = [ material.layer(*el) for el in list_layers]
    schottky_contact = mesh.mesh(layers=layers)
    schottky_contact.bend()
    schottky_contact.plot()

example_diode()
example_schottky_contact()
