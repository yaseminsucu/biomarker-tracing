# Biomarker Tracing Pipeline — VM-Compatible Version (AD Workbench)

This repository contains a **VM-compatible adaptation** of **Phuc Nguyen's biomarker tracing pipeline**, originally designed to run on an HPC environment.  
The goal of this fork is to enable running the pipeline on the **AD Workbench virtual machine** with minimal changes to the original analysis logic.

- **Original pipeline author:** Phuc Nguyen 
- **Original repository / reference:**  https://github.com/Mustardburger/biomarker-tracing
- This repo is maintained by **Yasemin Sucu** and focuses on VM portability (paths, execution wrappers, and environment setup).

## What’s Different in This VM Version
Compared to the HPC version, this repository includes modifications to:
- Replace HPC job submission (e.g., LSF/Slurm) with local/VM execution scripts
- Update filesystem paths and I/O assumptions for AD Workbench VM
- Add reproducible environment setup for the VM (conda/pip)


> **Note:** The analysis intent of the pipeline is unchanged; modifications are primarily to support the VM runtime environment.

## Requirements
- OS: Linux (tested on AD Workbench VM)
- Dependencies: see `packages.yml` in the original repository 



