# Molecular Modelling and Quantum Chemistry — Master Course

## Complete Module Overview

**Course:** Molecular Modelling and Quantum Chemistry (Master)  
**Instructors:** Conrad Hübler  
**Duration:** 3 × 90 minutes (Sessions) + Pre-Course  
**Format:** Online with VirtualBox / OpenSUSE KDE  
**Language:** English  
**Last Updated:** April 29, 2026

---

## 🎯 Module Objectives

By the end of this module, students will be able to:

1. **Understand molecular modelling foundations:** Born-Oppenheimer approximation, force fields, ensembles
2. **Perform geometry optimization** on small molecules using Curcuma
3. **Run molecular dynamics simulations** using both Curcuma and Gromacs
4. **Analyze MD results:** energy, temperature, RMSD, trajectories
5. **Use AlphaFold** to predict protein structures
6. **Compare predictions with simulations:** when each method works, when it fails
7. **Work confidently with the Linux terminal** and scientific software
8. **Create production-ready MD input files** for cluster submission
9. **Visualize and interpret molecular structures** using PyMOL and Avogadro

---

## 📁 Course Structure

### **Pre-Course: Linux Console Basics** (Self-paced)
- **File:** `BashRPG.md` (LiaScript)
- **Topics:** File navigation, text editing, pipes, wildcards, grep, find
- **Deliverable:** Students should be comfortable with terminal before Session 1

### **Session 1: Geometry Optimization** (Self-paced)
- **File:** `session1_geometry_optimization.md` (LiaScript)
- **Topics:**
  - Born-Oppenheimer approximation
  - Force fields: UFF vs GFN-FF
  - Convergence criteria
  - Using Curcuma for optimization
- **Practical Exercises:**
  1. Glucose optimization (Monosaccharide, 2 methods)
  2. Sucrose optimization (Disaccharide)
  3. Peptide optimization (AAAAA and WRKLQ)
  4. (Optional) Fructose as alternative isomer
- **Analysis:** RMSD calculation, trajectory visualization (Avogadro)

### **Session 2A: Molecular Dynamics — Thermostats** (Self-paced)
- **File:** `session2a_md_curcuma.md` (LiaScript)
- **Topics:**
  - Molecular dynamics theory (NVE, NVT, NPT ensembles)
  - Thermostat concepts: Berendsen vs CSVR
  - MD with Curcuma
  - Data extraction and analysis
  - Gnuplot plotting
- **Practical Exercises:**
  1. Glucose MD: Berendsen vs CSVR comparison
  2. RMSD analysis during MD
  3. Peptide comparison: AAAAA vs WRKLQ
- **Analysis Methods:**
  - Extract temperature, energy, RMSD with AWK
  - Create publication-quality plots with Gnuplot
  - Interpret thermostat effects

### **Session 2B: Introduction to Gromacs** (Self-paced)
- **File:** `session2b_gromacs_intro.md` (LiaScript)
- **Topics:**
  - Gromacs file formats (.mdp, .top, .gro, .tpr, .edr, .xtc)
  - Workflow compiled by Lemkul J. Phys. Chem. B 2024, 128, 39, 9418–9435 (https://pubs.acs.org/doi/full/10.1021/acs.jpcb.4c04901)
  - Step-by-step MD protocol
  - MDP file parameters
  - Student role: prepare inputs, instructor runs on cluster
- **Workflow Steps:**
  1. Energy minimization
  2. NVT equilibration
  3. NPT equilibration
  4. Production MD
  5. Analysis (energy, RMSD, flexibility)

### **Session 3: AlphaFold vs Molecular Dynamics** (Live + analysis)
- **File:** `session3_alphafold_and_md.md` (LiaScript)
- **Topics:**
  - AlphaFold2 predictions
  - Comparison with MD simulations
  - When each method succeeds/fails
  - Protein dynamics interpretation
- **Instructor Actions:**
  - Collect .tpr files from Session 2B
  - Run production MD on cluster (~2-24 hours)
  - Generate .xtc, .edr, .log outputs
- **Student Analysis:**
  - Compare AlphaFold structure vs MD
  - Calculate RMSD (trajectory vs prediction)
  - Analyze per-residue flexibility (RMSF)
  - Interpret results
- **Deliverable:** Online presentation (10-15 min)
  - Protein choice and motivation
  - AlphaFold prediction quality metrics
  - MD simulation results (energy, RMSD, stability)
  - Structural comparison (visualizations)
  - Biological interpretation

---

## 🛠️ Tools Required

### **On your system**

| Tool | Version | Purpose | Installation |
|------|---------|---------|--------------|
| **Curcuma** | Latest | Geometry optimization, MD, RMSD | https://github.com/conradhuebler/curcuma |
| **Gromacs** | 2024.1 (or compatible) | Production MD simulations | https://www.gromacs.org/Downloads |
| **Avogadro** | 1.2+ | Structure visualization | https://two.avogadro.cc/install/index.html#install |
| **PyMOL** | 2.x | Advanced visualization | https://www.pymol.org/ |
| **Open Babel** | 3.x | SMILES ↔ 3D conversion | https://openbabel.org/index.html |
| **Gnuplot** | 5.x | Data plotting | http://www.gnuplot.info/ |
| **nano** | (included) | Text editing | Already on system |
| **bash/grep/sed/awk** | (included) | Command-line utilities | Already on system |

### **Online Resources**

| Resource | Purpose | Link |
|----------|---------|------|
| **ChemSpider** | SMILES/structure lookup | https://www.chemspider.com/ |
| **PDB** | Protein structures | https://www.rcsb.org/ |
| **AlphaFold Server** | Protein structure prediction | https://alphafoldserver.com/ |
| **Lemkul Tutorials** | Gromacs protocols | https://pubs.acs.org/doi/full/10.1021/acs.jpcb.4c04901 |

---

## 📚 References & Further Reading

### **Primary Literature**

1. **Lemkul, J.A.** "Introductory Tutorials for Simulating Protein Dynamics with GROMACS" *J. Phys. Chem. B* **2024**, 128, 9418−9435. https://doi.org/10.1021/acs.jpcb.4c04901

2. **Bussi, G.; Donadio, D.; Parrinello, M.** "Canonical sampling through velocity rescaling" *J. Chem. Phys.* **2007**, 126, 014101. https://doi.org/10.1063/1.2408420

3. **Rappe, A.K.; Casewit, C.J.; et al.** "UFF, a Full Periodic Table Force Field for Molecular Mechanics and Molecular Dynamics Simulations" *J. Am. Chem. Soc.* **1992**, 114(25), 10024–10035.

### **Software Documentation**

- **Curcuma GitHub:** https://github.com/conradhuebler/curcuma
- **GROMACS Manual:** https://manual.gromacs.org/
- **Open Babel:** https://openbabel.org/
- **Avogadro:** https://avogadro.cc/
- **PyMOL:** https://pymol.org/
- **Gnuplot:** http://www.gnuplot.info/

### **Educational Resources**

- **Frenkel, D.; Smit, B.** *Understanding Molecular Simulation*, 2nd ed. Academic Press.

---

## Funding

The development of this course material is funded by:

- TUBAFdigital
- European Union (Erasmus+ National Agency for Higher Education)

Funded by the European Union. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or Erasmus+ National Agency for Higher Education (German Academic Exchange Service). Neither the European Union nor the granting authority can be held responsible for them.

![Co-funded by the European Union](https://github.com/conradhuebler/TeacherTwinMolecularModelling/raw/main/EN%20Co-funded%20by%20the%20EU_POS.jpg)

## 📞 Contact & Support

**Course Instructor:** Conrad Hübler
