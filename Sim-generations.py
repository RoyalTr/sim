import numpy as np
import os
import sys
import time
import multiprocessing

# Variables user can change
Repetitions = 20  # No. of times to rerun using the same parameter values
generations = 100000000000  # Prevent a possible runaway situation

# Prevent system sleep on Windows (when running overnight)
if os.name == 'nt':
    import ctypes
    ES_CONTINUOUS = 0x80000000
    ES_SYSTEM_REQUIRED = 0x00000001
    ctypes.windll.kernel32.SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED)

start_time = time.time()

def column_headings():
    return "Population size (N);selection coeff (s);attempts"

# Modified: Removed 'Max gen' from the headings.
def results_headings():
    return "Simul Nr;Rep;Pop size;Select coeff;Attempts;Prob allele loss;Prob allele fix;St. dev. allele fix;Aver gen to fix;St. dev. gen to fix"

output_headings_avg = "Reps;SimulationNr;Pop size (N);Select coeff;Attempts;Prob allele loss;Prob of allele fix;St. dev. allele fix;Aver gen to fix;St. dev. gen to fix"
results_filename_avg = "results_data_avg.txt"

example_rows = [
    "1000;0.001;5000",
    "10000;0.005;20000"
]

input_filename = "input_data.txt"
headings = column_headings()

if not os.path.exists(input_filename):
    with open(input_filename, "w") as f:
        f.write(headings + "\n")
        for row in example_rows:
            f.write(row + "\n")
    print("Please enter the parameters to run in file input_data.txt (see example data), then rerun the program.")
    sys.exit(0)
    
with open(input_filename, "r") as f:
    lines = [line.strip() for line in f if line.strip()]

if not lines:
    with open(input_filename, "w") as f:
        f.write(headings + "\n")
        for row in example_rows:
            f.write(row + "\n")
    print("Please enter the parameters to run in file input_data.txt (see example data), then rerun the program.")
    sys.exit(0)

if lines[0] != headings:
    os.remove(input_filename)
    with open(input_filename, "w") as f:
        f.write(headings + "\n")
        for row in example_rows:
            f.write(row + "\n")
    print("Please enter the parameters you want in file input_data.txt (see example data), then rerun the program.")
    sys.exit(0)

if len(lines) == 1:
    with open(input_filename, "a") as f:
        for row in example_rows:
            f.write(row + "\n")
    print("Please enter the parameters you want in file input_data.txt (see example data), then rerun the program.")
    sys.exit(0)

valid_data = []
error_found = False
for line_num, line in enumerate(lines[1:], start=2):
    parts = line.split(";")
    if len(parts) != 3:
        print(f"The data in file input_data.txt in line {line_num} is wrong. Please correct and rerun.")
        error_found = True
        continue
    try:
        N_val = int(parts[0])
        if not (1 <= N_val <= 1000000000000):
            print(f"The data in line {line_num} is wrong. Please correct")
            error_found = True
            continue
    except:
        print(f"The data in line {line_num} is wrong. Please correct")
        error_found = True
        continue
    try:
        s_val = float(parts[1])
        if not (-2 <= s_val <= 2):
            print(f"The data in line {line_num} is wrong. Please correct")
            error_found = True
            continue
    except:
        print(f"The data in line {line_num} is wrong. Please correct")
        error_found = True
        continue
    try:
        attempts_val = int(parts[2])
        if not (1 <= attempts_val <= 1000000000000):
            print(f"The data in line {line_num} is wrong. Please correct")
            error_found = True
            continue
    except:
        print(f"The data in line {line_num} is wrong. Please correct")
        error_found = True
        continue

    valid_data.append((N_val, s_val, attempts_val))

if error_found:
    print("Please correct the data and rerun the program")
    sys.exit(0)

def simulate_population(N, s, p0, generations, attempts):
    # Initial counters for the beginning of each simulation run (affects each CPU running in parallel)
    losses = 0
    fixation = 0
    sum_fixation_gens = 0.0
    sum_fixation_gens_sq = 0.0
    p_fix  = 1 - (1 / (2 * N))   # When the prob. of the new allele p_t > p_fix then only that allele is present in the pop.
    # p_fix = 0.999
    # The simulation will see if an allele fixes many times, namely 'attempts' times
    for i in range(attempts):
        p_t = p0
        for gen in range(generations):
            if p_t == 0.0:  # Random segregation removed all the new alleles in the pop.
                losses += 1
                break
            # elif p_t == 1.0: # Bad idea, a round-off error could prevent this from ever occurring
            elif p_t > p_fix:
                fixation += 1  # Increase the counter of allele fixed for 'attempts' times tried
                # If the current attempt to fix an allele succeeded add gen (the number of generations required) to sum_fixation_gens
                sum_fixation_gens += gen  # In Python 3 there is no upper limit for integers so sum_fixation_gens can be used safely
                sum_fixation_gens_sq += gen * gen
                break
            homozygous = p_t * p_t
            heterozygous = 2.0 * p_t * (1.0 - p_t)
            inheritance_prob = homozygous + 0.5 * heterozygous  # Prob. of inheriting allele A before natural selection
            # fit = ((1.0 + s) * inheritance_prob) / (((1.0 + s) * inheritance_prob) + (1.0 - inheritance_prob))
            fit = ((1.0 + s) * inheritance_prob) / (1.0 + inheritance_prob * s)
            # Sample 2N allele copies from a binomial distribution having a median of fit, which adjusts for allele relative fitness
            C = np.random.binomial(2 * N, float(fit))
            p_t = C / (2 * N)
    loss_probability = losses / attempts
    fixation_probability = fixation / attempts
    fixation_std = np.sqrt(fixation_probability * (1.0 - fixation_probability) / attempts)
    # If some alleles fixed, calculate the average of the generations required over all those that fixed
    if fixation > 0:
        avg_fixation_gen = sum_fixation_gens / fixation
        variance = (sum_fixation_gens_sq / fixation) - (avg_fixation_gen * avg_fixation_gen)
        std_fixation_gen = np.sqrt(variance) if variance > 0 else 0.0
    else:
        # Document that none of the alleles fixed during this simulation run
        avg_fixation_gen = np.nan
        std_fixation_gen = np.nan
    return loss_probability, fixation_probability, fixation_std, avg_fixation_gen, std_fixation_gen

# Define a worker function that wraps a single simulation run.
def worker(job):
    idx, rep, N, s, p0, generations, attempts = job
    result = simulate_population(N, s, p0, generations, attempts)
    return idx, rep, N, s, attempts, result[0], result[1], result[2], result[3], result[4]

# Multiprocessing section.
if __name__ == '__main__':
    max_processes = multiprocessing.cpu_count()
    print(f"Maximum number of processes supported: {max_processes}")
    
    # Build a list of simulation jobs.
    jobs = []
    for idx, (N, s, attempts) in enumerate(valid_data, start=1):
        for rep in range(1, Repetitions + 1):
            # Diploid organism p0 = 1/(2N) possible alleles
            # p0 = 1.0 / (2 * N)
            p0 = 100.0 / (2 * N)
            jobs.append((idx, rep, N, s, p0, generations, attempts))
    
    with multiprocessing.Pool(processes=max_processes) as pool:
        results = pool.map(worker, jobs)
    
    # Sort individual simulation results by simulation number and repetition.
    individual_results_sorted = sorted(results, key=lambda x: (x[0], x[1]))
    
    # Group results by parameter row for averaging.
    from collections import defaultdict
    grouped_results = defaultdict(list)
    for res in individual_results_sorted:
        grouped_results[res[0]].append(res)
    
    results_by_param = []
    for idx, group in grouped_results.items():
        loss_probs = [r[5] for r in group]
        fixation_probs = [r[6] for r in group]
        fixation_std_devs = [r[7] for r in group]
        fixation_gens = [r[8] for r in group]
        fixation_std_gens = [r[9] for r in group]
        avg_loss_prob = sum(loss_probs) / len(loss_probs)
        avg_fixation_prob = sum(fixation_probs) / len(fixation_probs)
        avg_fixation_std_dev = sum(fixation_std_devs) / len(fixation_std_devs)
        avg_fixation_gens_value = sum(fixation_gens) / len(fixation_gens)
        avg_std_fixation_gens_value = sum(fixation_std_gens) / len(fixation_std_gens)
        N_val = group[0][2]
        s_val = group[0][3]
        attempts_val = group[0][4]
        results_by_param.append((idx, N_val, s_val, attempts_val,
                                  avg_loss_prob, avg_fixation_prob, avg_fixation_std_dev,
                                  avg_fixation_gens_value, avg_std_fixation_gens_value))
    
    # Write individual simulation results to file "results_data.txt"
    results_filename = "results_data.txt"
    output_headings = results_headings()
    lines_to_write = []
    # Format: Simul Nr;Rep;Pop size;Select coeff;Attempts;Prob allele loss;Prob allele fix;St. dev. allele fix;Aver gen to fix;St. dev. gen to fix
    for rec in individual_results_sorted:
        line = f"{rec[0]};{rec[1]};{rec[2]};{rec[3]};{rec[4]};{rec[5]:.8f};{rec[6]:.8f};{rec[7]:.10f};{rec[8]:.2f};{rec[9]:.2f}"
        lines_to_write.append(line)
    
    if os.path.exists(results_filename):
        with open(results_filename, "a") as f:
            for line in lines_to_write:
                f.write(line + "\n")
    else:
        with open(results_filename, "w") as f:
            f.write(output_headings + "\n")
            for line in lines_to_write:
                f.write(line + "\n")
    
    print(f"Results stored in file {results_filename}.")
    
    # If Repetitions > 1, write averaged results to file "results_data_avg.txt"
    if Repetitions > 1:
        skip_creation = False
        try:
            if os.path.exists(results_filename_avg):
                os.remove(results_filename_avg)
        except Exception as e:
            print(f"Warning: Could not delete {results_filename_avg} because: {e}")
            skip_creation = True
        if not skip_creation:
            with open(results_filename_avg, "w") as f:
                f.write(output_headings_avg + "\n")
        avg_lines = []
        # Format: Reps;SimulationNr;Pop size (N);Select coeff;Attempts;Prob allele loss;Prob of allele fix;St. dev. allele fix;Aver gen to fix;St. dev. gen to fix
        for rec in results_by_param:
            line = f"{Repetitions};{rec[0]};{rec[1]};{rec[2]};{rec[3]};{rec[4]:.8f};{rec[5]:.8f};{rec[6]:.10f};{rec[7]:.2f};{rec[8]:.2f}"
            avg_lines.append(line)
        if not skip_creation:
            with open(results_filename_avg, "a") as f:
                for line in avg_lines:
                    f.write(line + "\n")
            print(f"The average values based on the repetitions were stored in file: {results_filename_avg}")
        else:
            print(f"Skipped creating average results file {results_filename_avg} due to file access issues.")
    
    end_time = time.time()
    execution_time = end_time - start_time
    print(f"\nExecution time required: {execution_time:.2f} seconds")
