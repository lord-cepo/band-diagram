import material
import mesh

Si_n = material.semiconductor(Eg=1.5, electron_affinity=4.05, epsilon=11.7)
Al = material.metal(4.1)
list_layers = [(50, Al), (100, Si_n)]
layers = [ material.layer(*el) for el in list_layers]

diode = mesh.mesh(layers=layers)
diode.plot(display_E0=True) # display_E0=True
