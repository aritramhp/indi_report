# indi_report
This program generates reports (in PDF) for individual patient from the microbiome signatures

# Installation and Setup:
## List of packages:
1) Install VS Code
2) Install Anaconda
3) Open Anaconda Prompt
4) Navigate to the `indi_report` folder
5) Run the following commands:
  - `conda create --name indi_report --file requirements.txt`
  - `conda activate indi_report`
  - `python -m ipykernel install --user --name=indi_report`

# Setup Reference DB
The tool comes with a default reference DB by including KORA healthy, KORA obese and Student cohort of 2016, 2017, 2018, 2019, 2023.
The data is kept in `ref_db` folder.
To update the DB,
1) Execute IMNGS2 / NGSToolkit on the reference samples: KORA-Healthy, KORA-Obese, and Student cohort -> Generate zOTU table
2) Execute Rhea:
  - Alpha diverisity: required file `alpha-diversity.tab`
  - Beta diversity on KORA-Healthy and Students cohort: required file `samples-Tree.nwk`
  - Taxonomy binning on KORA-Healthy, KORA-Obese, and Students cohort: required files `1.Phyla.all.tab` and `5.Genera.all.tab`
3) Keep the files in `ref_db` folder
4) Modify the `clade.tab` file to match with the clades in the tree

# Execution:
1) Open `patientReport_v1.ipynb` on VS Code
2) Select `indi_report` environment on VS Code
3) You are ready to generate reports

# Files
- `patientReport_v1`: 
- `patientReport_v2`: It is able to run Rhea from python. However, not all the computers support this feature.
- `functions.py`: All the necessary functions are present in the file. 
