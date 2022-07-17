import material

Si_n = material.semiconductor(Eg=1.12, electron_affinity=4.05, epsilon=11.7, doping_type='n', doping=1e15)
Si_p = material.semiconductor(Eg=1.12, electron_affinity=4.05, epsilon=11.7, doping_type='p', doping=1e15)
GaAs = material.semiconductor(Eg=1.424, 
                              electron_affinity=4.07, 
                              epsilon=12.9, 
                              effective_m_e=0.063, 
                              effective_m_h=0.51,
)

# Intended as Al(x) Ga(1-x) As
def AlGaAs(x):
    return material.semiconductor(Eg=1.424+1.247*x if x<0.45 else 1.9+0.125*x+0.143*x**2,
                                electron_affinity=4.07-1.1*x if x<0.45 else 3.64-0.14*x,
                                epsilon=12.9-2.84*x,
                                effective_m_e=(0.063+0.083*x)**1.5 if x<0.41 else (0.85-0.14*x)**1.5,
                                effective_m_h=(0.51+0.25*x)**1.5 if x<0.41 else (0.85-0.14*x)**1.5,
                                doping_type='n',
                                doping=GaAs.n
    )

Al = material.metal(work_function=4.1)