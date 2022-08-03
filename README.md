# band-diagram
1D Plotter of band diagrams of various devices. Energy of conduction, valence, Fermi and vacuum levels are displayed as a function of the dimension of the device.
Voltage can be applied to single layers of material. The script is intended as a visual helper and results should not be taken quantitatively. I've made the following assumptions:
- applied voltage drops only at the interface between two materials, hence Fermi level is constant in bulk material;
- the existence of surface states in neglected;
- band bending occurs under full depletion region approximation, boltzmann smearing correction on the population at the edges of the depletion region
has not been implemented yet;
- metals are designed as semiconductor with zero band gap, high fixed electronic density $(10^{23}\text{ cm}^{-3})$, infinite DC dielectric constant;

# tutorials
**Installation**: clone this repository. **Requirements**: `numpy` and `matplotlib`. There are three core modules: 
- [`material.py`](https://github.com/lord-cepo/band-diagram/blob/master/material.py), containing the definition of metal and semiconductor classes; 
- [`mat_data.py`](https://github.com/lord-cepo/band-diagram/blob/master/mat_data.py), a tiny database of semiconductors and metals that can be used as a playground;
- [`band.py`](https://github.com/lord-cepo/band-diagram/blob/master/band.py), calculates band bending and plots the results.

A [`quick_start`](https://github.com/lord-cepo/band-diagram/blob/master/notebooks/quick_start.ipynb) guide can be found inside `nootebooks` folder. On [`some_device`](https://github.com/lord-cepo/band-diagram/blob/master/notebooks/some_device.ipynb) notebook you can find some ideas on device design

# equations
At the interface between two materials (1 and 2), depletion regions $x_1$ and $x_2$ are:
$$x_1 = \left(\frac{2}{q}\epsilon_1\epsilon_2 \frac{n_2}{n_1}\frac{|V_{bi}-V|}{\epsilon_1 n_1 + \epsilon_2 n_2} \right)^{1/2}$$ $$x_2 = \left(\frac{2}{q}\epsilon_1\epsilon_2 \frac{n_1}{n_2}\frac{|V_{bi}-V|}{\epsilon_1 n_1 + \epsilon_2 n_2}\right)^{1/2}$$
where $q$ is the elementary charge, $\epsilon$ the dielectric constant, $n$ the carrier density, $V_{bi}$ the built-in potential, $V$ the applied voltage.
Potential $\phi$ grows and decays quadratically across the interface, bending in opposite manner all the energy levels $(E=-q \phi)$:
$$\phi(x) =\phi_{-\infty} \pm \frac{q n_1 (x + x_1)^2}{2 \epsilon_1} \quad \text{if $-x_1 < x < 0$}$$ 
$$\phi(x) =\phi_{+\infty} \mp \frac{q n_2 (x - x_2)^2}{2 \epsilon_2} \quad \text{if $0 < x < x_2$}$$
Sign depends on the sign of the interface potenital $(\phi_{+\infty} - \phi_{-\infty} = V - V_{bi})$, such potential is automatically continuous at the interface thanks to the definition of $x_1$ and $x_2$.
If bended, Fermi level follows same equations, given $V_{bi}=0$.

# plots
Inside [plots folder](https://github.com/lord-cepo/band-diagram/tree/master/plots) you will find all the outputs of `some_device.ipynb`. Here are some of them:
<img src="./plots/Reverse%20bias%20(with%20more%20smearing).png" alt="pn_junction" width="380"/>
<img src="./plots/Si_Ge%20heterojunction.png" alt="heterojunction" align="top" width="620"/>

