import json

def 

def generate_json(
    N: int,
    sigma,
    sides,
    bumps_per_side,
    k = 1,
    dt = 0.001,
    max_energy_iter = 1000000,
    E_thresh = 0.0000000001,
    P_t = 0.00001,
    P_tol_log = 0.001,
    r_scale = 0.999
):
    return json.dumps(
        {"N": N,
        "sigma": sigma,
        "sides": sides,
        "bumps_per_side": bumps_per_side,
        "k": k,
        "dt": dt,
        "max_energy_iter": max_energy_iter,
        "E_thresh": E_thresh,
        "P_t": P_t,
        "P_tol_log": P_tol_log,
        "r_scale": r_scale
        })
def main():
