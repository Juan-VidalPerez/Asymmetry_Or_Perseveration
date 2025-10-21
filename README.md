# Code and Data for: Learning Asymmetry Or Perseveration? A Critical Re-Evaluation And Solution To A Pervasive Confound

This repository contains the MATLAB code and data necessary to reproduce the computational modeling analyses and main figures presented in the paper.

**Paper:** [https://osf.io/preprints/psyarxiv/xdse5_v1](https://osf.io/preprints/psyarxiv/xdse5_v1)

---
## 📂 Repository Structure

* **`/` (Root Directory):** Contains the main MATLAB functions (`.m` files) for model fitting, simulation, and plotting, including the `slanCM.m` colormap function.
* **`/data/`:** Contains all necessary MATLAB data files (`.mat`).

---
## 🛠️ Requirements & Dependencies

* **MATLAB:** Tested on **R2024a**. Earlier versions might work but are untested.
* **MATLAB Toolboxes:**
    * Statistics and Machine Learning Toolbox (for `fitglme`, `ttest`, `corr`, etc.)
    * Optimization Toolbox (for `fit_models.m`, `fit_models_Chambon.m`)
    * Parallel Computing Toolbox (recommended for feasible fitting times using `parfor` loops in fitting functions)
* **`slanCM.m` Colormap Function:** Included in the root directory. Used for specific color schemes in Figures 7 & 8.

---
## 📊 Datasets (`/data/` directory)

All data is stored in MATLAB `.mat` files within the `/Data/` subdirectory.

### 1. Behavioral Task Data

Raw behavioral data from various experiments, typically structured as cell arrays (each cell = one participant).

* `S2021.mat`: Sugawara & Katahira (2021) - Used for S1a & S1b.
* `P2017a.mat`: Palminteri et al. (2017) - Used for P1.
* `P2017b.mat`: Palminteri et al. (2017) - Used for P2.
* `L2017a.mat`: Lefebvre et al. (2017) - Used for L1.
* `L2017b.mat`: Lefebvre et al. (2017) - Used for L2.
* `C2020a.mat`: Chambon et al. (2020) - Used for C1.
* `C2020b.mat`: Chambon et al. (2020) - Used for C2.
* `C2020c.mat`: Chambon et al. (2020) - Used for C3.
* `C2020d.mat`: Chambon et al. (2020) - Used for C4.

**Structure within task data files:**

* `sta`: {1 x nSubs} cell array of state vectors (bandit pair id).
* `cho`: {1 x nSubs} cell array of choice vectors (1 or 2).
* `out`: {1 x nSubs} cell array of outcome vectors (-1 or 1).
* `cou`: {1 x nSubs} cell array of counterfactual outcome vectors (0 if partial feedback).
* `obs`: {1 x nSubs} cell array (only in C1-C3) indicating trial type (1=free, 0=forced).

### 2. Pre-computed Model Fits & Simulation Results

Contain results allowing direct figure generation without re-running fits/simulations.

* `MAP_fits.mat` & `MLE_fits.mat`: Results from model fitting (MAP/MLE).
    * `experiments`: {'L1',...,'S1b'} - Order of experiments.
    * `parameters_[MODEL]`: {1 x 10} cell array of fitted parameters.
    * `LPP_[MODEL]`: {1 x 10} cell array of log posterior probabilities / log likelihoods.
    * `parameters_[SIM_MODEL]sim_[FIT_MODEL]fit`: Parameters from fitting `FIT_MODEL` to data simulated using `SIM_MODEL`. Data is a {1 x 10} cell array, where each cell contains a `[subject x parameter x simulation]` matrix (1001 simulations per subject).
    * `LPP_[SIM_MODEL]sim_[FIT_MODEL]fit`: Corresponding LPP/LL values.
* `CBPERS_recoveries.mat`: Results for CBPERS model recovery (Fig 3a/b).
    * `parameters_CBPERSgener_MAP`/`MLE`: Generative parameters used (cell array per experiment).
    * `parameters_CBPERSfitted_MAP`/`MLE`: Recovered parameters after fitting (cell array per experiment).
* `figure4_fits.mat`: Data for recovery vs. session length (Fig 4, S5, S11).
    * `parameters_PSL_MAP`/`MLE`: Generative PSL parameters.
    * `parameters_PSLsim_CBPERSfit_MAP`/`MLE`: Fitted CBPERS parameters. Data is a {1 x 10} cell array where each cell contains a `[session_length x simulation x parameter]` matrix (200 simulations for MAP, 4000 for MLE).
    * `session_length_MAP`/`MLE`: Vector indicating the number of sessions simulated.
* `figure5_data.mat`: Data for negative perseveration analysis (Fig 5c,d,e).
    * `phi_MAP`/`MLE`: Generative negative phi values used.
    * `parameters_negphi_MAP`/`MLE`: Recovered CBPERS parameters.
* `parameters_sweep.mat`: Data from parameter sweeps (Fig S4, plotted by `plot_parametersweep`).
    * `parameters_sweep_MAP`/`MLE`: Cell array {swept param} of fitted params `[simulation x parameter x swept_value]`.
    * `swept_MAP`/`MLE`: Cell array {swept param} of swept values.
    * `generative_MAP`/`MLE`: Base generative parameters.

---
## ⚙️ MATLAB Functions (Root Directory)

### Core Fitting & Simulation

* **`fit_models.m`**: Fits models (RW, CB, PSL, CBPERS) to standard task data (all tasks except C1, C2, C3). *(Formerly `fit_pers_MAP.m`)*
* **`fit_models_Chambon.m`**: Fits models to observational task data (C1, C2, C3). *(Formerly `fit_pers_Chambon.m`)*
* `simulate_signatures.m`: Simulates data for behavioral signature paradigms (Fig 7).
* `simulate_newtask.m`: Simulates data for the 4-condition task (Fig 8).

### Plotting Functions

These functions load pre-computed data from `/Data/` and generate figures.

* `plot_fig2.m`: Figure 2 & S10 (Parameter estimates).
* `plot_fig3.m`: Figure 3 (Recovery & MAP vs MLE).
* `plot_fig4.m`: Figure 4, S5 & S11 (Recovery vs session length).
* `plot_fig5_cde.m`: Figure 5c,d,e (Negative phi effects).
* `plot_parametersweep.m`: Figure 5a-b, S6 & S13 (Parameter sweep recovery).
* `plot_fig6.m`: Figure 6 (Bootstrap analysis).
* `plot_fig7.m`: Figure 7 (Behavioral signatures).
* `plot_fig8.m`: Figure 8 (New task signatures).

### Parameter Order in Fitting Functions

**`fit_models.m`**
* Model 1 (RW): `[beta, lr1]`
* Model 2 (CB): `[beta, lr1, lr2]`
* Model 3 (PSL): `[beta, lr1, tau, phi]`
* Model 4 (CBPERS): `[beta, lr1, lr2, tau, phi]`

**`fit_models_Chambon.m`**
* Model 1: `[beta, lr1, lr3]`
* Model 2: `[beta, lr1, lr2, lr3]`
* Model 3: `[beta, lr1, lr3, tau, phi]`
* Model 4: `[beta, lr1, lr2, lr3, tau, phi]`

*(Note: `lr1` corresponds to $\alpha_c$, `lr2` to $\alpha_d$)*

---
## 🚀 Usage

### Generating Figures

1.  **Setup:** Ensure MATLAB (R2024a recommended) is installed with the required toolboxes (Statistics, Optimization, Parallel Computing). Place all `.mat` data files in the `/Data/` subdirectory.
2.  **Run:** Open MATLAB, navigate to the root directory containing the `.m` files, and execute the desired plotting function from the command window.

    * **Figure 2 (MAP):** `plot_fig2()`
    * **Figure S10 (MLE):** `plot_fig2('MLE')`
    * **Figure 3:** `plot_fig3()`
    * **Figure 4 & S5 (MAP):** `plot_fig4()`
    * **Figure S11 (MLE):** `plot_fig4('MLE')`
    * **Figure 5a-b, S13 (MLE):** `plot_parametersweep()`
    * **Figure 5c,d,e (MLE):** `plot_fig5_cde()`
    * **Figure 6:** `plot_fig6()`
    * **Figure 7:** `plot_fig7(50000, 1000)`
    * **Figure 8:** `plot_fig8(10000)`

    *(Note: Figures 7 & 8 require `simulate_signatures.m` and `simulate_newtask.m`)*

### Running Model Fits (Example)

You can re-run the model fitting procedures using the provided functions. This requires the behavioral data files in `/Data/`.

```matlab
% Example: Fit the CBPERS model (model 4) to the P2 dataset using MAP
[parameters, LPP] = fit_models('P2', 4, 'MAP');

% Example: Fit model 2 to the C1 dataset using MLE
[parameters_C1, LPP_C1] = fit_models_Chambon('C1', 2, 'MLE');

