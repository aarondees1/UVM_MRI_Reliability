# UVM_MRI_Reliability

## UVM MRI Reliability Study: Cartilage Thickness Processing and Comparison

This repository contains a suite of MATLAB scripts designed to process, analyze, and compare subchondral bone and cartilage data from MRI scans, specifically focusing on the tibia and femur. The scripts establish anatomical coordinate systems, transform digitized points, calculate cartilage thicknesses, and perform comparisons across different MRI scan types (T1FFE, T1rho, and T2S).

The workflow is sequential, with each script building upon the outputs of the preceding ones.

---

## Workflow and Order of M-Files

The scripts in this repository must be executed in the following order to ensure proper data flow and dependencies are met:

```
tibias08b_AD.m
tibias08c_AD.m
tcart08_sag_AD.m
tcthk_cmp_AD.m
femurs08b_AD.m
femurs08c_AD.m
fcart08_sag_AD.m
fcthk_cmp_AD.m
```

Here is a flowchart illustrating the complete workflow for the UVM MRI Reliability Study, from initial MRI scans and digitization to the final comparison of cartilage thicknesses using your MATLAB scripts, updated with the new .mat file names for the sd structure.

```mermaid
graph TD
    A[MRI Scans] --> B{OsiriX Digitization}
    B --> C[Raw Data CSVs]

    C --> D(Create Tibia Bone Meshes - tibias08b_AD.m)
    D --> E(Create Tibia Cartilage Meshes - tibias08c_AD.m)
    E --> F(Calculate Tibia Cartilage Thicknesses - tcart08_sag_AD.m)
    F --> G(Compare Tibia Cartilage Thicknesses and Export Results in Excel - tcthk_cmp_AD.m)
    
    F --> I(Create Femur Bone Meshes - femurs08b_AD.m)
    I --> J(Create Femur Cartilage Meshes - femurs08c_AD.m)
    J --> K(Calculate Femur Cartilage Thicknesses - fcart08_sag_AD.m)
    K --> L(Compare Femur Cartilage Thicknesses and Export Results in Excel- fcthk_cmp_AD.m)
    L --> M[Femur Analysis Results]

    style A fill:#f9f,stroke:#333,stroke-width:2px
    style B fill:#f9f,stroke:#333,stroke-width:2px
    style C fill:#f9f,stroke:#333,stroke-width:2px
    style H fill:#bbf,stroke:#333,stroke-width:2px
    style M fill:#bbf,stroke:#333,stroke-width:2px
``` 

The flowchart visually represents the entire data processing pipeline, highlighting the sequential nature of the MATLAB scripts and their key inputs and outputs.

---

## Script Descriptions

Each script performs a specific set of tasks, contributing to the overall data processing pipeline:

### 1. `tibias08b_AD.m`
**Purpose**: Processes digitized tibia bone data to establish the tibial coordinate system.

**Functionality**: Reads raw bone digitization points and transforms the data into the tibia's anatomical coordinate system. It plots the coordinate system and the transformed bone data.

**Outputs**:
- PDF Plots: Plots of the tibia coordinate system and subchondral bone, following the format `***_tibias08b.pdf`.
- MAT File: `***_tibiaCS.mat`
- sd Structure: Updates the sd structure (`All Subjects Tibia Bone Data.mat`).

### 2. `tibias08c_AD.m`
**Purpose**: Processes digitized tibia cartilage data.

**Dependencies**: Relies on the `***_tibiaCS.mat` files.

**Outputs**:
- PDF Plots: `***_tibias08c.pdf`
- MAT File: `***_tibiaCart.mat`
- sd Structure: Updates `All Subjects Tibia Bone And Cartilage Data.mat`

...

*(Note: For brevity in this code cell, remaining script descriptions are omitted but will be retained in the actual file.)*

---

## Common Dependencies

Several custom MATLAB M-files are required in the current directory or MATLAB path for these scripts to run correctly. These include (but are not limited to):

```
fix_pts_AD.m
tibia_cs8.m
f_cs_14.m
rd_roi6.m
mk_tri6.m
mk_tri4f.m
tri_fix2.m
car_thk8.m
gridproj.m
comb_dat.m
nod2tri.m
quadconn.m
freg_axpf2_AD.m
```

---

## Data Structure (sd)

The `sd` (subject data) structure is a critical component of this workflow. It is initialized in the first script and progressively updated by subsequent scripts.

---

## Usage

To use these scripts:

1. Organize your raw data: Ensure all necessary raw data files (e.g., CSV digitizations) are placed within a directory named `Data`.
2. Ensure all necessary custom M-files are in your MATLAB path.
3. Run the scripts sequentially as listed.
4. Each script will prompt you to select the base data directory.

---

## Authors

- Mack Gardner-Morse
- Aaron Dees

**Last Edited**: July 7, 2025