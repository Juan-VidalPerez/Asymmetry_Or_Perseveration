# Code and Data for: Learning asymmetry or perseveration? Disentangling overlapping mechanisms in behavioral modeling

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
    * Statistics and Machine Learning Toolbox
    * Optimization Toolbox
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
*Note: the LA model is sometimes refered in code as the CB model, and the hybrid model as the CBPERS model.*

* `MAP_fits.mat`, `MAP_fits_wide.mat` & `MLE_fits.mat`: Results from model fitting (MAP with N(0,1) phi prior, MAP with N(0,5) phi prior and MLE).
    * `experiments`: {'L1',...,'S1b'} - Order of experiments.
    * `parameters_[MODEL]`: {1 x 10} cell array of fitted parameters.
    * `LPP_[MODEL]`: {1 x 10} cell array of log posterior probabpailities / log likelihoods.
    * `parameters_[SIM_MODEL]sim_[FIT_MODEL]fit`: Parameters from fitting `FIT_MODEL` to data simulated using `SIM_MODEL`. Data is a {1 x 10} cell array, where each cell contains a `[subject x parameter x simulation]` matrix (1001 simulations per subject).
    * `LPP_[SIM_MODEL]sim_[FIT_MODEL]fit`: Corresponding LPP/LL values.
*`figure3a_data.mat`: Data for the PSL simulation analysis (Fig 3a).
    * `parameters_P2`: Average parameters used for the PSL simulation.
    * `sd`: Vector containing the standard deviation of the $\phi$ prior for each fit.
    * `parameters_PSL_sd`: Resulting fits of simulations `[sd x simulation x parameters]`.
* `figure3bc_data.mat`: Results for hybrid model recovery (Fig 3b,c).
    * `parameters_CBPERSgener_MAP`/`MLE`: Generative parameters used (cell array per experiment).
    * `parameters_CBPERSfitted_MAP`/`MLE`: Recovered parameters after fitting (cell array per experiment).
* `figure4a_fits.mat`: Data for recovery vs. session length (Fig 4a, S9, S17).
    * `parameters_PSL_MAP`/`MLE`: Generative PSL parameters.
    * `parameters_PSLsim_CBPERSfit_MAP`/`MLE`: Fitted CBPERS parameters. Data is a {1 x 10} cell array where each cell contains a `[session_length x simulation x parameter]` matrix (200 simulations for MAP, 4000 for MLE).
    * `session_length_MAP`/`MLE`: Vector indicating the number of sessions simulated.
*`figure5_data.mat`: Data for negative perseveration analysis and APPS (Fig 5, S13).
    * `parameters_CBPERS`: Hybrid fits of the experimental data `[participant x parameter]`.
    * `parameters_CBPERSsim_CBPERSfit`: Fits used to generate the APPS null distribution `[participant x parameter x anti-perseveration x simulation]`.
* `parameters_sweep.mat`: Data from parameter sweeps (Fig 4b, S10, S20).
    * `parameters_sweep_MAP`/`MLE`: Cell array {swept param} of fitted params `[simulation x parameter x swept_value]`.
    * `swept_MAP`/`MLE`: Cell array {swept param} of swept values.
    * `generative_MAP`/`MLE`: Base generative parameters.

---
## ⚙️ MATLAB Functions (Root Directory)

### Core Fitting & Simulation

* `fit_models.m`: Fits reinforcement learning models (RW, LA, PSL, HYBRID) to behavioral data using Maximum a Posteriori (MAP) or Maximum Likelihood Estimation (MLE). This function is adapted from Palminteri et al. (2023), *Nat Rev Neurosci*.
* `simulate_models.m`: Simulates behavioral data for reinforcement learning models (RW, LA, PSL, HYBRID). This function routes simulations to specific sub-functions based on the experimental design (handling standard vs. C1-3 designs).
* `simulate_signatures.m`: Simulates data for behavioral signature paradigms (Fig 7).
* `simulate_newtask.m`: Simulates data for the 4-condition task (Fig 8).

### Plotting Functions

These functions load pre-computed data from `/Data/` and generate figures.

* `plot_fig2.m`: Parameter estimates.
* `plot_fig3a.m`, `plot_fig3bcde.m`: Recovery analysis (Figure 3).
* `plot_fig4a.m`: Recovery vs. session length (Figure 4a).
* `plot_fig4b.m`: Negative phi effects (Figure 4b).
* `plot_fig5.m`: APPS analysis on the dataset (Figure 5).
* `plot_fig6.m`: Behavioral signatures (Figure 6).
* `plot_figS17.m`: Parameter recovery sweep (Supplementary Figure 17).
* `plot_figS22.m`: Supplementary analysis (Supplementary Figure 22).
* `plot_figS6_S16.m`: Combined plots for Supplementary Figures 6 and 16.
* `plot_CBsweep.m`: Independent confirmation bias sweep analysis.

### Parameter Order in Fitting Functions

**`fit_models.m`** 
(experiments L1,L2,P1,P2,C4,S1a and S1b)
* Model 1 (RW): `[beta, lr1]`
* Model 2 (CB): `[beta, lr1, lr2]`
* Model 3 (PSL): `[beta, lr1, tau, phi]`
* Model 4 (CBPERS): `[beta, lr1, lr2, tau, phi]`

(experiments C1,C2,C3)
* Model 1: `[beta, lr1, lr3]`
* Model 2: `[beta, lr1, lr2, lr3]`
* Model 3: `[beta, lr1, lr3, tau, phi]`
* Model 4: `[beta, lr1, lr2, lr3, tau, phi]`

*(Note: `lr1` corresponds to* $\alpha_{c}$*, `lr2` to* $\alpha_{d}$ *and* `lr3` *to* $\alpha$ *for forced choices in (C1-3)).*

---
## 🚀 Usage

### Generating Figures

1.  **Setup:** Ensure MATLAB (R2024a recommended) is installed with the required toolboxes (Statistics, Optimization, Parallel Computing). Place all `.mat` data files in the `/Data/` subdirectory.
2.  **Run:** Open MATLAB, navigate to the root directory containing the `.m` files, and execute the desired plotting function from the command window.

    * **Figure 2 (MAP N(0,1) phi prior):** `plot_fig2()`
    * **Figure 2 (MAP N(0,5) phi prior):** `plot_fig2('MAP_wide')`
    * **Figure 3a:** `plot_fig3a()`
    * **Figure 3bcde:** `plot_fig3bcde()`
    * **Figure 4a and S17 (MLE):** `plot_fig4a()`
    * **Figure 4b and S19 (MLE):** `plot_fig4b()`
    * **Figure 5 and S18 (APPS):** `plot_fig5()`
    * **Figure 6:** `plot_fig6()`
    * **Figure 7:** `plot_fig7(50000, 1000)`
    * **Figure 8:** `plot_fig8(10000)`
    *(Note: Figures 7 & 8 require `simulate_signatures.m` and `simulate_newtask.m`)*

    Other SI figures:
    * **Figure S6:** `plot_figS6_S16()`
    * **Figure S10 (MAP):** `plot_fig4a('MAP')`
    * **Figure S12 (MLE):** `plot_fig4b('MLE')`
    * **Figure S14 (MLE):** `plot_fig2('MLE')`
    * **Figure S17:** `plot_figS6_S17('MLE')`
    * **Figure S22:** `plot_figS22()`


### Running Model Fits (Example)

You can re-run the model fitting procedures using the provided functions. This requires the behavioral data files in `/Data/`.

```matlab
% Example: Fit the CBPERS(hybrid) model (model 4) to the P2 dataset using MAP
[parameters, LPP] = fit_models('P2', 4, 'MAP');

% Example: Fit model 2(LA) model to the C1 dataset using MLE
[parameters_C1, LPP_C1] = fit_models_Chambon('C1', 2, 'MLE');

