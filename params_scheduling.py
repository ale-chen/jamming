import json
import os.path
from datetime import datetime
import argparse

def param_ranges(min_sides, max_sides, min_bumps_per_side, max_bumps_per_side, N, sigma):
    jsons = []
    filenames = []
    for sides in range(min_sides,max_sides+1):
        for bumps in range(min_bumps_per_side, max_bumps_per_side+1):
            jsons.append(generate_json(N,sigma,sides,bumps))
            filenames.append(str(N)+'polygons'+str(sides)+'sides'+str(bumps)+'perSide')
    return jsons,filenames

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
def main(json_directory, min_sides = 3, max_sides = 12, min_bumps_per_side = 1, max_bumps_per_side = 20, N = 512, sigma = 1):
    command_file_string = ''
    command_string = """module load MATLAB/2023a && matlab -nodisplay -r "generate_initial_state('_output','./_param_files','{filename}.json','false'); exit;
    """""

    jsons, filenames = param_ranges(min_sides,max_sides,min_bumps_per_side, max_bumps_per_side, N, sigma);
    # 10 copies per param
    for i in range(0,len(jsons)):
        command_file_string+=command_string.format(filename=filenames[i])
        with open(os.path.join(json_directory,filenames[i]+'.json'), 'w', encoding='utf-8') as f:
            f.write(jsons[i])
            now = datetime.now()
    dt_string = now.strftime("%d_%m_%Ytime%H_%M_%S")
    with open(dt_string+'script.txt', 'w') as txt_file:
        txt_file.write(command_file_string)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate parameter files and MATLAB commands')
    parser.add_argument('json_directory', type=str, help='Directory to output JSON files')
    parser.add_argument('--min_sides', type=int, default=3, help='Minimum number of sides (default: 3)')
    parser.add_argument('--max_sides', type=int, default=12, help='Maximum number of sides (default: 12)')
    parser.add_argument('--min_bumps_per_side', type=int, default=1, help='Minimum bumps per side (default: 1)')
    parser.add_argument('--max_bumps_per_side', type=int, default=20, help='Maximum bumps per side (default: 20)')
    parser.add_argument('--N', type=int, default=512, help='Number of polygons (default: 512)')
    parser.add_argument('--sigma', type=float, default=1, help='Sigma value (default: 1)')

    args = parser.parse_args()

    main(
        args.json_directory,
        args.min_sides,
        args.max_sides,
        args.min_bumps_per_side,
        args.max_bumps_per_side,
        args.N,
        args.sigma
    )
