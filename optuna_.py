
from typing import Dict
import optuna
import subprocess
import re
from multiprocessing import Pool
import argparse
from functools import partial
from pathlib import Path


def create_model(trial):
    # n_layers = trial.suggest_int('n_layers', 1, 2)
    n_layers = 2
    envs = {"n_layers": str(n_layers),
            "time": str(trial.suggest_float('time', 0.0, 1.0))
            }

    # two_opt_is_annealing = trial.suggest_categorical(f"{i}_two_opt_is_annealing", [0, 1])
    # replace_is_annealing = trial.suggest_categorical(f"{i}_replace_is_annealing", [0, 1])

    two_opt_is_annealing = 0
    replace_is_annealing = 0

    # two_opt_temp_alpha = 0.1
    # two_opt_time_offset = 1.0

    two_opt_temp_alpha = 0.1
    two_opt_time_offset = 1.0
    replace_temp_alpha = 0.1
    replace_time_offset = 1.0

    envs.update({
        f"n0_two_opt_is_annealing": str(two_opt_is_annealing),
        f"n0_replace_is_annealing": str(replace_is_annealing),
        f"n0_two_opt_temp_alpha": str(two_opt_temp_alpha),
        f"n0_two_opt_time_offset": str(two_opt_time_offset),
        f"n0_replace_temp_alpha": str(replace_temp_alpha),
        f"n0_replace_time_offset": str(replace_time_offset),
    })

    # two_opt_is_annealing = trial.suggest_categorical(f"{i}_two_opt_is_annealing", [0, 1])
    # replace_is_annealing = trial.suggest_categorical(f"{i}_replace_is_annealing", [0, 1])

    two_opt_is_annealing = 1
    replace_is_annealing = 1

    if two_opt_is_annealing:
        two_opt_temp_alpha = trial.suggest_float(f"1_two_opt_temp_alpha", 1e-6, 1e-1, log=True)
        two_opt_time_offset = trial.suggest_float(f"1_two_opt_time_offset", 0, 10)
    else:
        two_opt_temp_alpha = 0.1
        two_opt_time_offset = 1.0

    if replace_is_annealing:
        replace_temp_alpha = trial.suggest_float(f"1_replace_temp_alpha", 1e-6, 1e-1, log=True)
        replace_time_offset = trial.suggest_float(f"1_replace_time_offset", 0, 10)
    else:
        replace_temp_alpha = 0.1
        replace_time_offset = 1.0

    envs.update({
        f"n1_two_opt_is_annealing": str(two_opt_is_annealing),
        f"n1_replace_is_annealing": str(replace_is_annealing),
        f"n1_two_opt_temp_alpha": str(two_opt_temp_alpha),
        f"n1_two_opt_time_offset": str(two_opt_time_offset),
        f"n1_replace_temp_alpha": str(replace_temp_alpha),
        f"n1_replace_time_offset": str(replace_time_offset),
    })
    return envs


# d18512.tsp
# fnl4461.tsp
# nrw1379.tsp
# d1291.tsp
# d657.tsp
# rd400.tsp
# ch150.tsp
# eil51.tsp
# test.tsp

def run(env: Dict[str, str], file: Path) -> int:
    ret = subprocess.check_output(f"cat  {file} | ./build/seminar-competition", env=env, shell=True)
    # ret = b'recomputed tour length = 1897\ntime for the search:     12.00 seconds\ntime to read the instance:    0.00 seconds\n'
    length = int(re.match(r"recomputed tour length = (\d+)", ret.decode()).group(1))
    return length


def objective(trial, n: int, file: Path):
    model = create_model(trial)
    p = Pool(n)
    result = p.map(partial(run, file=file), [model] * n)
    return sum(result) / n


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str)
    parser.add_argument("-n", type=int, default=4)

    args = parser.parse_args()

    study = optuna.create_study(
        direction="minimize",
        pruner=optuna.pruners.MedianPruner()
    )

    file = Path(args.file)
    if not file.exists():
        print(f"{file} is not found")
        return

    study.optimize(lambda x: objective(x, args.n,  file), n_trials=100)

    pruned_trials = [t for t in study.trials if t.state == optuna.structs.TrialState.PRUNED]
    complete_trials = [t for t in study.trials if t.state == optuna.structs.TrialState.COMPLETE]
    print("Study statistics: ")
    print("  Number of finished trials: ", len(study.trials))
    print("  Number of pruned trials: ", len(pruned_trials))
    print("  Number of complete trials: ", len(complete_trials))

    # 結果の表示
    print("Best trial:")
    trial = study.best_trial
    print("  Value: ", trial.value)
    print("  Params: ")
    for key, value in trial.params.items():
        print(f"    {key}: {value}")


if __name__ == "__main__":
    main()
