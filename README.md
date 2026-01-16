# SensinSilico
**In Silico Simulation of Sensor Response Time Mismatch in Dual-Sensor Systems**

SensinSilico is a graphical user interface (GUI)–based simulation software published together with the peer-reviewed article  
**“Timing Matters: The Overlooked Issue of Response Time Mismatch in pH-Dependent Analyte Sensing using Multiple Sensors”**,  
published in the journal *Analyst*.  
DOI: https://doi.org/10.1039/D3AN01207G

The software enables a systematic investigation of **error propagation caused by mismatched sensor response times** when combining two (electrochemical) sensors for total-parameter monitoring.

---

## Overview

In environmental and analytical sensing, it is often required to monitor **multiple analytes simultaneously** and ideally at the same location. In practice, however, sensors frequently differ in key characteristics such as:

- Sensor response time  
- Acquisition rate  
- Signal dynamics  

A common approach is to adapt the measurement procedure to the **faster sensor**, implicitly accepting a small inaccuracy introduced by the slower one. While this limitation is well known, its quantitative impact on calculated total parameters has rarely been assessed.

**SensinSilico** addresses this gap by providing a simulation framework that allows users to:

- Examine dual-sensor systems with different response times  
- Quantify systematic bias and error propagation  
- Optimize measurement protocols prior to experimental deployment  

---

## Scientific Background

### Response Time Mismatch in Dual-Sensor Measurements

When two sensors respond at different rates to changing analyte concentrations, their combined signal can introduce a **systematic bias** into calculated total parameters (e.g. total acidity or other composite chemical quantities).

SensinSilico simulates:

- Time-dependent sensor response behavior  
- Dynamic changes in analyte concentration  
- Mathematical reconstruction of total parameters  

This enables users to assess whether commonly accepted inaccuracies remain tolerable or fundamentally compromise measurement validity.

---

## Key Features

### Software Capabilities

- **Dual-Sensor Simulation**  
  Simulates two sensors with independently defined response times and kinetics

- **Gompertz-Based Response Modeling**  
  Graphical implementation of Gompertz response curves for realistic sensor dynamics

- **Error Propagation Analysis**  
  Quantifies bias introduced by delayed sensor response

- **Measurement Protocol Optimization**  
  Explore acquisition timing strategies to minimize total-parameter error

- **Interactive Graphical User Interface**  
  No programming knowledge required for standard operation

### Exploratory & Advanced Use

- **Custom Response Curves**  
  Extend beyond Gompertz kinetics using the included Jupyter Notebook

- **Transparent and Reproducible Simulations**  
  Deterministic workflow suitable for method development, teaching, and reproducibility studies

---

## Prerequisites

### System Requirements

- **Operating System:** Windows, macOS, or Linux  
- **Python:** Version 3.8 or higher  
- **Memory:** 4 GB RAM minimum (8 GB recommended)  
- **Storage:** ~500 MB available disk space  

### Knowledge Requirements

- Basic understanding of chemical or electrochemical sensing  
- Familiarity with sensor response behavior (recommended)  
- No programming knowledge required for GUI use  

---

## Installation

### 1. Clone or Download the Repository

```bash
git clone https://github.com/yourusername/SensorResponse-inSilico.git
cd SensorResponse-inSilico
```
Alternatively, download the repository as a ZIP archive and extract it locally.

```bash
python3 -m venv sensinsilico-env
```
Then, activate the environment.

### 2. Install Dependencies
```bash
pip install -r requirements.txt
```
The current release has been tested with the following package versions.
The listed versions correspond to the tested release environment.
Older versions may cause errors, while newer versions should generally remain compatible.

---

## Directory Structure
Please ensure that the directory structure remains unchanged:
```bash
SensorResponse-inSilico/
│
├── SensinSilico.py            # Main GUI application
├── requirements.txt           # Python dependencies
├── Instruction.pdf            # User manual
│
├── picture/                   # GUI graphics and figures
├── scripts/                   # Supplementary scripts
│   ├── Code_snippes.ipynb
│   └── dbs_CodeSnippes.py
```
Relocating or renaming files may break internal file references.

---

## Running the Application
To start the GUI, execute:
```bash
python SensinSilico.py
```
Make sure, that you are using python3. 

---

## Exploring Alternative Response Curves
In addition to the GUI (which supports Gompertz response curves), `SensinSilico` includes tools for exploratory modeling:
```bash
Jupyter Notebook: scripts/Code_snippes.ipynb
Supporting Script: scripts/dbs_CodeSnippes.py
```

To execute the notebook, ensure both files are located in the same directory. 
Launch Jupyter Notebook or JupyterLab and execute the notebook cells step by step.
These tools allow users to implement and evaluate alternative sensor response models.

---

## Manual & Documentation
A detailed step-by-step user guide is provided in: [Instruction.pdf](Instruction.pdf)
The manual includes:
  * Theoretical background
  * Description of all GUI elements
  * Example simulation workflows
  * Guidance on result interpretation
  * Validation & Compatibility
    * Tested on macOS (Apple M1, version 13.4)
    * Tested on Windows 10 (version 10.0.19045)

---

## Citation
If you use SensinSilico in your work, please cite:
```bash
@article{Zieger2023Timing,
  title={Timing Matters: The Overlooked Issue of Response Time Mismatch in pH-Dependent Analyte Sensing using Multiple Sensors},
  journal={Analyst},
  year={2023},
  doi={10.1039/D3AN01207G}
}
```

---

## License
Copyright © 2023–2026
SilviaE. Zieger
This software is released under the MIT License.
See the LICENSE file for full license text.
