import numpy as np

def get_params(ranges=None):
    # Default values from file (in exact order used in model.rhs)
    defaults = np.array([
        1e-15, 1e-8, 1e-7, 1e-8, 1e-9, 1e-10,  # beta* (6)
        0.1282, 0.1282, 0.7, 0.01, 0.15, 1.0, 2.0, 2.0, 0.3,  # gamma* (9)
        0.8055, 19.757, 6.1212, 8.664339, 5.37e-5, 1.02, 0.05, 2e-7, 5.2939, 2.0, 2.0, 1.0,  # delta* (12)
        1e8, 1e6, 1e5, 5e5, 5e10, 1e5, 1e5, 1e5, 1e7,  # lambda* (9)
        0.75, 0.9, 0.17, 0.18, 1e-10, 1.5e-5, 1e-9, 0.1245, 1e-7,  # mu* (9)
        1e10, 1e10,  # max* (2)
        10.0, 0.001, 10.0, 2.0531, 3.02, 6.7979, 6.9937, 0.01, 0.001,  # k*, ktc* (9)
        1e9, 1e8, 1e9, 1e9, 0.01, 4e-7, 0.2, 0.05, 0.0001, 1e-5, 0.1, 0.01, 0.01, 0.2, 0.01 # some other params (15)
    ])
    
    params = defaults.copy()
    if ranges is not None:
        for i, (low_mult, high_mult) in ranges.items():
            params[i] *= np.random.uniform(low_mult, high_mult)
    
    return params

def get_default_ic():
    # Placeholder ICs (adjust as needed; not in file)
    return np.array([1, 0, 0, 0, 85000, 15000, 71000, 12000, 56000, 8000, 0.0085, 0.12, 0.0094])
