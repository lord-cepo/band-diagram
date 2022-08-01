import material
import band
import mat_data
from hypothesis import given, strategies as st
import random
random.seed(4169)
import numpy as np

########################################
# material.py
########################################

@given(work_function=st.floats(min_value=1e-2, max_value=1e2))
def test_metal(work_function):
    """
    Given: 
        material is metal, with various work functions
    Tests: 
        if energy gap is zero and vacuum level corresponds with work function 
    """
    obj = material.metal(work_function=work_function)
    assert obj.levels['Ec'] == 0
    assert obj.levels['Ev'] == 0
    assert obj.levels['E0'] == work_function


@given(Eg=st.floats(min_value=1e-4, max_value=1e2),
    effective_m_e=st.floats(min_value=1e-10, max_value=1e10),
    effective_m_h=st.floats(min_value=1e-10, max_value=1e10),
    doping_type=st.sampled_from(['n', 'p', 'others']),
    doping=st.floats(min_value=1e0, max_value=1e20))
def test_semiconductor(Eg, doping_type, doping, effective_m_e, effective_m_h):
    """
    Given: 
        semiconductor with a broad range of Eg and doping 
        - 1e-4 < Eg < 1e2
        - 1e-10 < effective_m < 1e10
        - doping type: n, p, undoped
        - doping < 1e20
    Tests: 
        if energy gap is positive
    """
    obj = material.semiconductor(
        Eg, electron_affinity=Eg+1., 
        doping_type=doping_type, doping=doping,
        effective_m_e=effective_m_e, effective_m_h=effective_m_h
    )
    assert obj.levels['Ec'] > obj.levels['Ev']


@given(Eg=st.floats(min_value=0.01, max_value=10.0),
    doping_type=st.sampled_from(['n', 'p', 'others']),
    doping=st.floats(min_value=1e0, max_value=1e17))
def test_nondegenerate(Eg, doping_type, doping):
    """
    Given: 
        semiconductor with narrower range of Eg and mild doping
        - 0.01 < Eg < 10
        - doping type: n, p, undoped
        - doping < 1e17

    Tests: 
        if semiconductor is nondegenerate, i.e. if Fermi level is within the gap
    """
    obj = material.semiconductor(Eg, electron_affinity=Eg+1., doping_type=doping_type, doping=doping)
    assert obj.levels['Ec'] >= 0
    assert obj.levels['Ev'] <= 0

########################################
# band.py and mat_data.py
########################################

def E_is_continuous(level, device):
    """
    Support function, tests if vacuum level is continuous. This should be always
    true immediately after calling bend(). Discontinuity occurs if the energy
    distance between two consecutive points is more than 10 meV.
    """
    device.bend()
    continuous = True
    if level == 'Ef':
        for i, el in enumerate(device.Ef):
            if i>0 and abs(el-device.Ef[i-1]) > 1e-2:
                continuous = False
                break
        return continuous
    else:
        for i, el in enumerate(device.levels[level]):
            if i>0 and abs(el-device.levels[level][i-1]) > 1e-2:
                continuous = False
                break
        return continuous
             

@given(layer=st.floats(min_value=1e-3, max_value=1e1))
def test_layer_thick(layer):
    """
    Given: 
        band diagram with n/p-Si with various thicknesses
        - 1 nm < thickness < 10 um
        - points on band diagram grid = 10_000
    Tests: 
        if E0 is continuous, to find safe thickness range over which the image is 
        too grainy. Increase grid resolution if you want a clearer image
    """
    Si_n = mat_data.Si('n', 1e16)
    Si_p = mat_data.Si('p', 1e16)
    device = band.band_diagram(
        [
            material.layer(layer, Si_p),
            material.layer(layer, Si_n)
        ],
        n_points=10000
    )
    assert E_is_continuous('E0', device)


@given(semi1=st.sampled_from(mat_data.list_of_semiconductors),
    semi2=st.sampled_from(mat_data.list_of_semiconductors))
def test_semiconductor_types(semi1, semi2):
    """
    Given:
        a random semiconductor/semiconductor junction, 2 um thick
    Tests:
        if E0 is continuous, after bending
    """
    device = band.band_diagram(
        [
            material.layer(1, semi1()),
            material.layer(1, semi2())
        ]
    )
    assert E_is_continuous('E0', device)


@given(metal=st.sampled_from(mat_data.list_of_metals),
    semi=st.sampled_from(mat_data.list_of_semiconductors))
def test_semi_metal_types(metal, semi):
    """
    Given:
        a random metal/semiconductor junction, 2 um thick
    Tests:
        if E0 is continuous, after bending. Grid must be of 10_000 points, metals
        introduce a more discontinuous behaviour
    """
    device = band.band_diagram(
        [
            material.layer(1, metal),
            material.layer(1, semi())
        ],
        n_points=10000
    )
    assert E_is_continuous('E0', device)


def test_densities():
    """
    Given:
        ms contact
    Tests:
        if density values are all positive
    """
    device = band.band_diagram(
        [
            material.layer(1, mat_data.Si('p', 1e15)),
            material.layer(1, mat_data.Al)
        ]
    )
    assert all([el > 0 for el in device._densities()])    


def test_epsilons():
    """
    Given:
        ms contact
    Tests:
        if epsilon values are all positive
    """
    device = band.band_diagram(
        [
            material.layer(1, mat_data.Si('p', 1e15)),
            material.layer(1, mat_data.Al)
        ]
    )
    epsilons = [layer.material.epsilon for layer in device.layers]
    assert all([el > 0 for el in epsilons])    


def test_apply_voltage():
    """
    Given:
        fixed junction on which random voltages are applied
    Tests:
        if the discontinuity of Fermi level across the interface is minus the
        applied voltage
    """
    device = band.band_diagram(
        [
            material.layer(1, mat_data.Si('p', 1e15)),
            material.layer(1, mat_data.Ge('n', 1e16))
        ]
    )
    s = 0
    n_points = device.Ef.size
    for i in range(100):
        volt = i*random.random()
        s += volt
        device.apply_voltage(volt)
        assert abs(device.Ef[n_points//2]-device.Ef[n_points//2-1] + s) < 1e-10


def test_multi_bend_Ev_Ec():
    """
    Given:
        device on which we call bend() a random amount of times
    Tests:
        if the device has not changed drastically, i.e. if Ec and Ev has not moved
        with respect to the first bend() call. bend() effectively make a change once
        called multiple times, but only slightly and without being graphically
        noticed (bending after first shift each point by no more than 1 meV).
    """
    device = band.band_diagram(
        [
            material.layer(1, mat_data.Si('p', 1e15)),
            material.layer(1, mat_data.Si('n', 1e15))
        ]
    )
    device.bend()
    initial_Ev = np.copy(device.levels['Ev'])
    initial_Ec = np.copy(device.levels['Ec'])
    for _ in range(random.randint(1,10)):
        device.bend()
        assert np.sum(np.abs(initial_Ev-device.levels['Ev'])) < 1
        assert np.sum(np.abs(initial_Ec-device.levels['Ec'])) < 1


def test_multi_bend_E0():
    """
    Given:
        device on which we call bend() a random amount of times
    Tests:
        if E0 is continuous
    """
    device = band.band_diagram(
        [
            material.layer(1, mat_data.Si('p', 1e15)),
            material.layer(1, mat_data.Si('n', 1e15))
        ]
    )
    for _ in range(random.randint(1,10)):
        device.bend()
    assert E_is_continuous('E0', device)


def test_random_calls():
    """
    Given:
        device in which random methods of band_diagram are called a random amount of
        times.
    Tests:
        if E0 is continuous
    """
    device = band.band_diagram(
        [
            material.layer(1, mat_data.Si('p', 1e15)),
            material.layer(1, mat_data.Si('n', 1e15))
        ]
    )
    for _ in range(random.randint(1, 10)):
        r = random.randint(1, 3)
        match r:
            case 1:
                device.apply_voltage(r)
            case 2:
                device.bend()
            case 3:
                device.reset()
    E_is_continuous('E0', device)
    

def test_reset():
    """
    Given:
        device in a random state
    Tests:
        if reset() is able to bring back the device to original state.
    """
    device = band.band_diagram(
        [
            material.layer(1, mat_data.Si('p', 1e15)),
            material.layer(1, mat_data.Si('n', 1e15))
        ]
    )
    initial_Ev = np.copy(device.levels['Ev'])
    initial_Ec = np.copy(device.levels['Ec'])
    for _ in range(random.randint(1, 10)):
        r = random.randint(1, 2)
        match r:
            case 1:
                device.apply_voltage(r)
            case 2:
                device.bend()
    device.reset()
    assert np.sum(np.abs(initial_Ev - device.levels['Ev'])) < 1e-10
    assert np.sum(np.abs(initial_Ec - device.levels['Ec'])) < 1e-10


def test_bend_fermi():
    """
    Given:
        device with random voltage applied and with Fermi level bended
    Tests:
        if Fermi level is continuous
    """
    device = band.band_diagram(
        [
            material.layer(1, mat_data.Si('p', 1e15)),
            material.layer(1, mat_data.Si('n', 1e15))
        ]
    )
    device.apply_voltage(random.randint(-10, 10))
    device.bend(fermi=True)
    assert E_is_continuous('Ef', device)


def test_Eg_constant():
    """
    Given:
        bended homojucntion
    Tests:
        if Eg is constant across the device
    """
    device = band.band_diagram(
        [
            material.layer(1, mat_data.Si('p', 1e15)),
            material.layer(1, mat_data.Si('n', 1e15))
        ]
    )
    device.bend()
    assert np.std(device.levels['Ec']-device.levels['Ev']) < 1e-10


def test_interface_continuity():
    """
    Given:
        bended homojunction Si/Si
    Tests:
        if Ev, Ec, E0 are continuous only at interface, with and without an applied
        voltage
    """
    device = band.band_diagram(
        [
            material.layer(1, mat_data.Si('p', 1e15)),
            material.layer(1, mat_data.Si('n', 1e16))
        ]
    )
    device.bend()
    i = device.interfaces[0]
    for E in device.levels.values():
        assert np.abs(E[i]-E[i-1]) < 1e-2
    device.apply_voltage(1)
    device.bend()
    for E in device.levels.values():
        assert np.abs(E[i]-E[i-1]) < 1e-2


def test_interface_discontinuity():
    """
    Given:
        bended heterojunction Si/Ge
    Tests:
        if Ev, Ec are discontinuous only at interface, with and without an 
        applied voltage
    """
    device = band.band_diagram(
        [
            material.layer(1, mat_data.Si('p', 1e15)),
            material.layer(1, mat_data.Ge('n', 1e16))
        ]
    )
    device.bend()
    i = device.interfaces[0]
    assert np.abs(device.levels['Ev'][i]-device.levels['Ev'][i-1]) > 1e-2
    assert np.abs(device.levels['Ec'][i]-device.levels['Ec'][i-1]) > 1e-2
    device.apply_voltage(1)
    device.bend()
    assert np.abs(device.levels['Ev'][i]-device.levels['Ev'][i-1]) > 1e-2
    assert np.abs(device.levels['Ec'][i]-device.levels['Ec'][i-1]) > 1e-2


def test_apply_voltage_symmetry():
    """
    Given:
        junction in which voltage is applied to last layer and then to first
        layer changing sign
    Tests:
        if the difference in Ef at the interface is equal
    """
    device = band.band_diagram(
        [
            material.layer(1, mat_data.Si('p', 1e15)),
            material.layer(1, mat_data.Si('n', 1e15))
        ]
    )
    device.apply_voltage(volt=10, index=1)
    i = device.interfaces[0]
    delta_Ef_1 = device.Ef[i]- device.Ef[i-1]
    device.reset()
    
    device.apply_voltage(volt=-10, index=0)
    delta_Ef_2 = device.Ef[i]- device.Ef[i-1]
    assert np.abs(delta_Ef_1 - delta_Ef_2) < 1e-10


def test_apply_voltage_translation():
    """
    Given:
        junction in which same voltage is applied to both ends
    Tests:
        if Ef is the same across the device (if Ef is simply translated)
    """
    device = band.band_diagram(
        [
            material.layer(1, mat_data.Si('p', 1e15)),
            material.layer(1, mat_data.Si('n', 1e15))
        ]
    )
    device.apply_voltage(volt=10, index=0)
    device.apply_voltage(volt=10, index=1)
    assert np.std(device.Ef) < 1e-10

