# CHIRON EFFECT GASTRULOIDS BULK RNAseq ANALYSIS

The structure of the project is:

  - `setup.sh`: Bash script to setup the folder.
  - `requirements.txt`: File containing the python packages used for the analysis and to setup a reproducible environment.
  - Analysis files:
    - `bulk_analysis.ipynb`: Analysis of the bulk dataset.
    - `chir_conditions_analysis.ipynb`: Analysis of the chiron conditions and comparison with mouse data.
    - `custom_functions.py`: An auxiliary file with custom functions used to produce some of the plots in the analysis and that are reused between analysis.

# Setup

To setup the analysis, follow the following steps:

 1. **Clone** the repository:

    ```
    git clone https://github.com/stembryo-lab/chiron_gastruloids_bulk_rna_seq
    ```

 2. **Conditions** This analysis should not be affected by specific versions of the packages as the package dependency is quite simple and does not depend on stochstic algorithms. Nonetheless you can:
    
    - Check you have the installed packages from `requirements.txt` in your python environment.
    
    - Create a nother environment using conda executing in a terminal and from the repository folder:

        ```
        chiron_gastruloids_bulk_rna_seq>>> conda create env_chiron_gastruloids
        chiron_gastruloids_bulk_rna_seq>>> conda activate env_chiron_gastruloids
        chiron_gastruloids_bulk_rna_seq>>> conda install --file requirements.txt
        ```

    - Use venv instead of conda. Here it is a link with [documentation](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/).

 3. Inside the repository, in a terminal (with the virtual environment activated if you have created it), execute the setup script:

    ```
    chiron_gastruloids_bulk_rna_seq>>> conda activate env_chiron_gastruloids
    chiron_gastruloids_bulk_rna_seq>>> ./setup.sh
    ```

    This step will download all the data used in the analysis and will prepare it in the folder `data`. A second folder `results` is created, where the figures from the analysis will be saved.

# Analysis

For executing the analysis, the different notebooks can be executed in indistinct order. To execute them just activate the environment and open jupyter lab.