# band-diagram

1D Plotter of band diagrams of various devices. Energy of conduction, valence, Fermi and vacuum levels are displayed as a function of the dimension of the device.
Voltage can be applied to single layers of material. The script is intended as a visual helper and results should not be taken quantitatively. I've made the following assumptions:
- applied voltage drops only at the interface between two materials, hence Fermi level is constant in bulk material;
- the existence of surface states in neglected;
- band bending occurs under full depletion region approximation, boltzmann smearing correction on the population at the edges of the depletion region
has not been implemented yet;
- metals are designed as semiconductor with zero band gap, high fixed electronic density (![equation](https://latex.codecogs.com/svg.image?\inline&space;10^{23}) ![equation](https://latex.codecogs.com/svg.image?\inline&space;\text{cm}^{-3}))
, infinite DC dielectric constant;

# equations
At the interface between two materials (1 and 2), depletion regions $x_1$ and $x_2$ are:
$$x_1 = \left(\frac{2}{q}\epsilon_1\epsilon_2 \frac{n_2}{n_1}\frac{|V_{bi}-V|}{\epsilon_1 n_1 + \epsilon_2 n_2} \right)^{1/2}$$ $$x_2 = \left(\frac{2}{q}\epsilon_1\epsilon_2 \frac{n_1}{n_2}\frac{|V_{bi}-V|}{\epsilon_1 n_1 + \epsilon_2 n_2}\right)^{1/2}$$
where $q$ is the elementary charge, $\epsilon$ the dielectric constant, $n$ the carrier density, $V_{bi}$ the built-in potential, $V$ the applied voltage.
Potential $\phi$ grows and decays quadratically across the interface, bending in opposite manner all the energy levels ($E=-q \phi$):
$$\phi(x) =\phi_{-\infty} \pm \frac{q n_1 (x + x_1)^2}{2 \epsilon_1} \quad \text{if $-x_1 < x < 0$}$$ 
$$\phi(x) =\phi_{+\infty} \mp \frac{q n_2 (x - x_2)^2}{2 \epsilon_2} \quad \text{if $0 < x < x_2$}$$
Sign depends on the sign of the interface potenital ($\phi_{+\infty} - \phi_{-\infty} = V - V_{bi})$, such potential is automatically continuous at the interface thanks to the definition of $x_1$ and $x_2$.
# tutorial


# plots


