# Causal Discovery in Exchangeable Data

This repository contains the code and resources associated with the thesis "Causal Discovery in Exchangeable Data," authored by Tiago Brogueira with the support of Prof. Mário Figueiredo, conducted at Instituto de Telecomunicações, Técnico Lisboa, Portugal. This work has led to the submitted paper, "Rethinking Causal Discovery through the lens of Exchangeability," currently under review.

***

## 📄 Overview

This thesis introduces a novel **synthetic dataset** specifically designed for bivariate causal discovery. The dataset, along with its associated weights, is packaged in the `exchangeable_synthetic_dataset.zip` file. This archive includes:
* Individual examples as `.txt` files, organized with respect to the controllable design choices.
* A `data` folder containing a `.csv` file with dataset weights.
* An auxiliary Python file with functions: `run_synthetic()`, `analyse_predictions()`, and `analyse_exchangeable_dataset()`, which are identical to those found in `test_synthetic.py`.

***

## 📁 Repository Structure

This repository has the following structure:

```
Thesis-Causal-Discovery-In-Exchangeable-Data/
├── code/
│   ├── allfigs/
│   ├── data/
│   ├── extra/
│   ├── generate_scenarios/
│   ├── generate_txt/
│   ├── interesting/
│   ├── predictions/
│   ├── anm.ipynb
│   ├── cepairsimplementation.ipynb
│   ├── cgnn.ipynb
│   ├── emd.ipynb
│   ├── exchangeable_synthetic_dataset.zip
│   ├── igci.ipynb
│   ├── lingam.ipynb
│   ├── onehitwonders.ipynb
│   ├── pnl.ipynb
│   ├── reci.ipynb
│   ├── synthetic_nn_keras.ipynb
│   ├── test_synthetic.ipynb
│   ├── testing_functions.ipynb
│   └── testR.R
├── other implementations/
│   ├── bqcd/
│   │   └── bqcd.R
│   └── slope-0181208/
│       └── Slope.R
├── Pairs/
├── previous_results.xlsx
├── README.md
├── requirements.txt
└── TuebingenAnalysis.xlsx
```

It is organized into the following main directories:

* **`code/`**: Contains all original code, along with generated data (images and raw results).
* **`other_implementations/`**:
    * Includes implementations of **bCQD** and **Slope** methods, adapted from their original R repositories ([tagas/bQCD](https://github.com/tagas/bQCD) and [eda.rg.cispa.io/prj/slope/](https://eda.rg.cispa.io/prj/slope/)).
    * Houses four synthetic causal discovery datasets obtained from [Harvard Dataverse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/3757KX).
* **`Pairs/`**: Contains information from the **Tuebingen Cause-Effect Pairs Collection**, sourced from [webdav.tuebingen.mpg.de/cause-effect/](https://webdav.tuebingen.mpg.de/cause-effect/).

In addition to the folders, two key Excel files provide valuable metadata and results:

* **`Tuebingen_analysis.xlsx`**: This Excel file offers a detailed analysis for each Tuebingen pair, including:
    * `source_new`, `destination_new`, `Datapoints`, and `nature?` (continuous, discrete, or categorical) as described in the original Tuebingen dataset.
    * `Function`, `Aspect?`, and `Is Invertible?`: These columns reflect an analysis of the relationship between variables. `Function` indicates if the relationship follows a known function (e.g., linear, exponential), `Aspect?` assesses if the points visually adhere to these functions, and `Is Invertible?` is a visual "eye test" for invertibility.
    * `Dimensions?`: Indicates the dimensionality of the cause and effect variables. A value of `1` means both are one-dimensional.
    * `Cause Dist?`: An "eye test" analysis of the cause's distribution, categorized as "uniform," "normal," "bipolar," or "skewed."
    * `weights`: Original weights from the dataset metadata file.
    * `exchangeability_finally`: Our theoretical analysis of whether each pair has a reasonable latent variable.
    * `latent examples:`: Provides an example of such a latent variable.
    
* **`previous_results.xlsx`**: Serves as a logbook for tested methods, detailing the method name, specifications (if any), specific dataset used, and the obtained AUROC and accuracy.

***

## 📦 Requirements and Installation

This repository includes a `requirements.txt` file that specifies all the necessary Python dependencies for running the experiments and analyses.

To install the required packages, first ensure you are using **Python 3.11.9** (recommended), then run:

```bash
# Create and activate a virtual environment (recommended)
python -m venv venv
source venv/bin/activate   # On macOS/Linux
venv\Scripts\activate      # On Windows

# Install all dependencies
pip install -r requirements.txt
```


***

## 💻 Implementations of Causal Discovery Methods

Within the `code/` folder, various causal discovery method implementations are included:

* **RECI** and **IGCI**: Taken verbatim from [CausalDiscoveryToolbox](https://fentechsolutions.github.io/CausalDiscoveryToolbox/html/index.html).
* **CGNN** and **ANM**: Based on the CausalDiscoveryToolbox repository, with minor hyperparameter adjustments (for example, ANM uses Hoeffding's D test for computational efficiency).
* **EMD**: Sourced from [tagas/bQCD](https://github.com/tagas/bQCD).
* **LiNGAM** and **PNL**: Obtained from [ssamot/causality](https://github.com/ssamot/causality/tree/master).

All method implementations adhere to a homogenized structure: they accept an array `[x, y]` representing the variables. Most methods currently support only 1D variables and return `NaN` otherwise. The return value ranges from `-inf` to `inf`; a higher value (further from 0) indicates greater confidence that X causes Y, and vice versa. Many functions also support hyperparameter tuning via keyword arguments.

***

## ⚙️ Key Scripts and Functionality

* **`cepairsimplementation.ipynb`**: An auxiliary Jupyter Notebook containing essential functions:
    * `data, weights = getTuebingen()`: Loads the Tuebingen dataset.
    * `auroc, accuracy = testTuebingen(func, kwargs)`: Tests any causal discovery method on the Tuebingen dataset.
    * `score = test_independence(x, y, method_name, kwargs)`: Performs one of 8 different independence tests between two 1D variables. For homogenity, a higher score indicates greater dependence.

* **`exchangeabledatasetcreation.ipynb`**: This notebook is responsible for generating the synthetic dataset presented in the paper.
    * The core function is `create_dataset(func, disttheta, distpsi, eps, N, dist_samples, **kwargs)`.
        * `func`: Specifies the causal function for dataset generation.
        * `disttheta` and `distpsi`: Define the prior shapes for the latent variables theta and psi.
        * `eps`: Noise power.
        * `N`: Number of examples to generate.
        * `dist_samples`: Distribution from which the number of points in each pair will be sampled.
    * After generation, each `(x, y)` example (where `x` is the cause and `y` is the effect) is saved to a `.txt` file, and its visual representation is saved as an image.

* **`onehitwonders.py`**: Contains various one-time experiments, including the different normalized prior distributions and the code for generating Figure 2 in the paper.

* **`testR.R`**: Serves as a wrapper for running **Slope** and **bCQD** (named `qcd_function` and `Slope_` respectively; original Slope returns negative scores when X causes Y) on different datasets.
    * `run_tuebingen(method)`: Runs a method on the Tuebingen dataset.
    * `run_old(dataset_name, method)`: Runs a method on the old existing synthetic datasets from the literature.
    * `test_all(method)`: Runs a method on the newly proposed datasets.
    * **Note**: These functions save results to the `predictions` folder for later analysis and do not directly return evaluation metrics.

* **`test_synthetic.py`**: The primary file for analyzing different methods on both existing literature synthetic datasets and the new dataset proposed in this paper.
    * `run_synthetic(function, dataset_name)`: Runs a function on a specified dataset and saves results to the `predictions` folder.
    * `analyse_prediction(function_name, dataset_name)`: Computes AUROC and accuracy for a given function and dataset, saving the results to the `previous_results` Excel file.

***

## 📊 New Dataset Construction and Analysis

The second part of `test_synthetic.ipynb` focuses on constructing the full new dataset (Subsection 4.3 in the paper) by computing the optimal weight for each specific dataset based on its design choices.

* `compute_matrices()`: Outputs the **A** and **b** matrices as described in Equation 9 of the paper ($A_{paper} = \text{diag(weights)} \times A_{code}$). Keyword arguments allow choosing which methods to consider and whether to include noisy samples.
* `compute_coefficients_general()`: Receives the computed matrices and outputs the coefficients (weights, `w` in the paper) for all analyzed metrics (AUROC and accuracy). It also returns the minimized average error for both metrics. The type of regression (least squares, ridge, lasso) can be selected via the `regression_type` keyword argument.
* `compute_full_metric()`: Computes the precise result of any metric in the constructed dataset.
* `leave_one_out_validation()`: Computes the average cross-validated error as explained in Subsection 4.4 of the paper, with a selectable penalty (`"l1"` or `"l2"`).
* `analyse_exchangeable_dataset()`: Returns both the average AUROC and accuracy, considering their respective dataset weights.
* `test_exchangeable_synthetic()`: Receives a causal discovery function, runs it on the designed synthetic dataset, and outputs the AUROC and accuracy.

***

## 🧠 Neural Network Implementation

The **SynthNN** model is implemented in `synthetic_nn_keras.ipynb`.

* `build_binary_classifier()`: Allows the creation of a neural network with specified configurations. Default settings follow those used in the paper.
* `test_network()`: Performs N simulations of the neural network. Hyperparameters controlling network architecture, image preprocessing, outlier removal, and Gaussian filtering can be passed as keyword arguments.
