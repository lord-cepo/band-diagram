import numpy as np

class layer:
    def __init__(self, thickness, material) -> None:
        self.thickness = thickness
        self.material = material

class material:
    def __init__(self, epsilon, Ec, Ev, work) -> None:
        self.epsilon = epsilon
        self.levels = {
            'Ec': Ec,
            'Ev': Ev,
            'E0': Ec + work
        }
    
    def is_metal(self):
        pass

    def is_semiconductor(self):
        pass

class metal(material):
    def __init__(self, work_function) -> None:
        super().__init__(np.inf, 0, 0, work_function)

class semiconductor(material):
    def __init__(self, epsilon, electron_affinity, doping_type=None, doping=0) -> None:
        if doping_type == 'n':
            pass
        elif doping_type == 'p':
            pass
        else:
            pass
        Ec = 0
        Ev = 0        
        super().__init__(epsilon, Ec, Ev, electron_affinity)