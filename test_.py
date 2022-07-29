import material
import band
import mat_data
from hypothesis import given, strategies as st

########################################
# material.py
########################################

"""
Given: 
    material is metal, with various work functions
Tests: 
    if energy gap is zero and vacuum level corresponds with work function 
"""
@given(work_function=st.floats(min_value=1e-2, max_value=1e2))
def test_metal(work_function):
    obj = material.metal(work_function=work_function)
    assert obj.levels['Ec'] == 0
    assert obj.levels['Ev'] == 0
    assert obj.levels['E0'] == work_function

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
@given(Eg=st.floats(min_value=1e-4, max_value=1e2),
    effective_m_e=st.floats(min_value=1e-10, max_value=1e10),
    effective_m_h=st.floats(min_value=1e-10, max_value=1e10),
    doping_type=st.sampled_from(['n', 'p', 'others']),
    doping=st.floats(min_value=1e0, max_value=1e20))
def test_semiconductor(Eg, doping_type, doping, effective_m_e, effective_m_h):
    obj = material.semiconductor(
        Eg, electron_affinity=Eg+1., 
        doping_type=doping_type, doping=doping,
        effective_m_e=effective_m_e, effective_m_h=effective_m_h
    )
    assert obj.levels['Ec'] > obj.levels['Ev']

"""
Given: 
    semiconductor with narrower range of Eg and mild doping
    - 0.01 < Eg < 10
    - doping type: n, p, undoped
    - doping < 1e17

Tests: 
    if semiconductor is nondegenerate, i.e. if Fermi level is within the gap
"""
@given(Eg=st.floats(min_value=0.01, max_value=10.0),
    doping_type=st.sampled_from(['n', 'p', 'others']),
    doping=st.floats(min_value=1e0, max_value=1e17))
def test_nondegenerate(Eg, doping_type, doping):
    obj = material.semiconductor(Eg, electron_affinity=Eg+1., doping_type=doping_type, doping=doping)
    assert obj.levels['Ec'] >= 0
    assert obj.levels['Ev'] <= 0

########################################
# band.py and mat_data.py
########################################

"""
Support function, tests if vacuum level is continuous. This should be always
true immediately after calling bend(). Discontinuity occurs if the energy
distance between two consecutive points is more than 10 meV.
"""
def E0_is_continuous(device):
    continuous = True
    for i, el in enumerate(device.levels['E0']):
        if i>0 and abs(el-device.levels['E0'][i-1]) > 1e-2:
            continuous = False
            break
    return continuous
             
"""
Given: 
    band diagram with n/p-Si with various thicknesses
    - 1 nm < thickness < 10 um
    - points on band diagram grid = 10_000
Tests: 
    if E0 is continuous, to find safe thickness range over which the image is 
    too grainy. Increase grid resolution if you want a clearer image
"""
@given(layer=st.floats(min_value=1e-3, max_value=1e1))
def test_layer_thick(layer):
    Si_n = mat_data.Si('n', 1e16)
    Si_p = mat_data.Si('p', 1e16)
    device = band.band_diagram(
        [
            material.layer(layer, Si_p),
            material.layer(layer, Si_n)
        ],
        n_points=10000
    )
    device.bend()
    assert E0_is_continuous(device)

"""
Given:
    a random semiconductor/semiconductor junction, 2 um thick
Tests:
    if E0 is continuous, after bending
"""
@given(semi1=st.sampled_from(mat_data.list_of_semiconductors),
    semi2=st.sampled_from(mat_data.list_of_semiconductors))
def test_semiconductor_types(semi1, semi2):
    device = band.band_diagram(
        [
            material.layer(1, semi1()),
            material.layer(1, semi2())
        ]
    )
    device.bend()
    assert E0_is_continuous(device)

"""
Given:
    a random metal/semiconductor junction, 2 um thick
Tests:
    if E0 is continuous, after bending. Grid must be of 10_000 points, metals
    introduce a more discontinuous behaviour
"""
@given(metal=st.sampled_from(mat_data.list_of_metals),
    semi=st.sampled_from(mat_data.list_of_semiconductors))
def test_semi_metal_types(metal, semi):
    device = band.band_diagram(
        [
            material.layer(1, metal),
            material.layer(1, semi())
        ],
        n_points=10000
    )
    device.bend()
    assert E0_is_continuous(device)
