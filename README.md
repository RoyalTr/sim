# sim.py

## Overview

`sim.py` simulates the stochastic fate (loss or fixation) of a new allele in a diploid population under selection and genetic drift. It runs multiple replicates in parallel, tracks allele-frequency trajectories, and outputs both per‐replicate and averaged statistics.

## Features

* Diploid model: samples 2N gene copies each generation
* Hardy–Weinberg genotype weighting (homozygote AA vs. heterozygote Aa) with heterozygote fitness = ½·s
* Configurable population size (N), selection coefficient (s), and replicate count
* Tracks allele‐frequency trajectories for the first 10 generations of the first parameter set
* Parallel execution via multiprocessing
* Input file validation with helpful error messages
* Output files:

  * `results_data.txt`: line per replicate
  * `results_data_avg.txt`: averaged statistics when `Repetitions > 1`

## Requirements

* Python 3.7+
* NumPy
* Matplotlib

## Usage

1. Populate `input_data.txt` (see **Input Format** below).
2. Configure parameters at top of `sim.py` (e.g., `Repetitions`, `max_generations`).
3. Run:

   ```bash
   python sim.py
   ```
4. A trajectory plot will appear for the first data row.
5. Two output files are generated:

   * `results_data.txt`
   * `results_data_avg.txt` (if `Repetitions > 1`)

## Input Format (`input_data.txt`)

* **Header (must match exactly):**

  ```
  Population size (N);selection coeff (s);attempts
  ```
* **Data rows:** semicolon‐delimited

  ```
  1000;0.001;5000
  10000;0.005;20000
  ```

  | Column       | Type    | Range        | Description                        |
  | ------------ | ------- | ------------ | ---------------------------------- |
  | **N**        | Integer | 1 to 1e12    | Diploid population size            |
  | **s**        | Float   | −2.0 to +2.0 | Selection coefficient (homozygote) |
  | **attempts** | Integer | 1 to 1e12    | Number of replicate simulations    |

## Output Files

* **`results_data.txt`**: columns:

  ```
  SimulNr;Rep;N;s;attempts;P_loss;P_fix;SD_fix;Mean_gen_fix;SD_gen_fix
  ```
* **`results_data_avg.txt`** (if `Repetitions>1`):

  ```
  Reps;SimulNr;N;s;attempts;Avg_P_loss;Avg_P_fix;Avg_SD_fix;Avg_mean_gen_fix;Avg_SD_gen_fix
  ```

## License

MIT License. Feel free to adapt and extend.
