import material
import mesh
from hypothesis import given, strategies as st

@given(work_function=st.floats(min_value=1e-2, max_value=1e2))
def test_metal(work_function):
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
       obj = material.semiconductor(Eg, electron_affinity=Eg+1., doping_type=doping_type, doping=doping,
                                    effective_m_e=effective_m_e, effective_m_h=effective_m_h)
       assert obj.levels['Ec'] > obj.levels['Ev']

@given(Eg=st.floats(min_value=0.01, max_value=10.0),
       doping_type=st.sampled_from(['n', 'p', 'others']),
       doping=st.floats(min_value=1e0, max_value=1e17))
def test_nondegenerate(Eg, doping_type, doping):
       obj = material.semiconductor(Eg, electron_affinity=Eg+1., doping_type=doping_type, doping=doping)
       assert obj.levels['Ec'] >= 0
       assert obj.levels['Ev'] <= 0
