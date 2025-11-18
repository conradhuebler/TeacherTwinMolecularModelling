# TeacherTwin Molecular Modelling Course

**Status:** ✅ Complete and ready for deployment
**Last Updated:** October 28, 2025
**Target:** Master-level students (30-40 participants)

---
## Implementation Status

**Updated:** October 28, 2025
- ✅ **GFN-FF Support**: NOW FULLY FUNCTIONAL
- ✅ **New Peptide Structures**: AAAAA (53 atoms), WRKLQ/Pentapeptid (107 atoms)
- ✅ **All Calculations**: Optimizations + MD simulations completed
- ✅ **Quiz Syntax**: All 70+ questions converted to working format

## General Instructions
- **Keep entries concise and focused to save tokens**
- **Keep git commits concise and focused**
- **Rule of thumb**: If a CLAUDE.md section exceeds 20 lines, consider if it's better placed elsewhere

#### Copyright and File Headers
- **Copyright ownership**: All copyright remains with Conrad Hübler as the project owner and AI instructor
- **Year updates**: Always update copyright year to current year when modifying files
- **Claude contributions**: Mark Claude-generated code sections but copyright stays with Conrad
- **Format**: `Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@chemie.tu-freiberg.de>`
- **AI acknowledgment**: Add Claude contribution notes in code comments, not copyright headers


## 📚 Project Overview

This self-paced course teaches molecular modelling and computational chemistry through hands-on exercises:

- **Session 1:** Geometry Optimization (UFF vs GFN-FF)
- **Session 2A:** Molecular Dynamics & Thermostats (Berendsen vs CSVR)
- **Session 2B:** Gromacs Workflow (Lemkul Protocol)
- **Session 3:** AlphaFold vs MD Simulations
- **Pre-course:** Linux Terminal Basics

All materials use LiaScript for interactive, browser-based learning.

---

## 🛠️ Setup Requirements

### Essential Tools

| Tool | Purpose | Installation |
|------|---------|--------------|
| **Curcuma** | Geometry optimization & MD | Binary at `/home/conrad/bin/curcuma` |
| **obabel** | SMILES → 3D structure conversion | `apt install openbabel` |
| **Avogadro** | Molecular structure visualization | `apt install avogadro` |
| **Gromacs** | Advanced MD simulations | `apt install gromacs` |

### Optional

- **PyMOL:** Protein structure visualization
- **Gnuplot:** Data plotting
- **ChemSpider:** Online SMILES lookup

---

## 📖 Usage Guide

### For Students

1. **Access Sessions:** Open `.md` files in LiaScript editor or browser
2. **Run Exercises:** Follow step-by-step tutorials
   - Create SMILES files → Convert to 3D
   - Run optimizations: `curcuma -opt molecule.xyz -method uff`
   - Analyze results with expected values in CALCULATION_RESULTS.md
3. **Visualize:** Use Avogadro to inspect trajectories
4. **Compare:** Cross-check your results against reference values

### For Instructors

1. **Pre-class:** Review `CALCULATION_RESULTS.md` for expected outcomes
2. **Live Session:** Use quiz questions for assessment
3. **Grading:** Compare student results to reference values
4. **Extension:** Suggest longer MD runs or different molecules

---

## ✅ Completed Work Summary

### Session 1: Geometry Optimization (November 2025)
- ✅ **Didactic Messaging:** Established GFN-FF as recommended standard method
  - UFF marked as fallback for molecules >1000 atoms only
  - Added clear disclaimer in Method sections
  - Reframed Quick Checks to emphasize GFN-FF primacy
- ✅ **Formatting:** Fixed 20 missing blank lines before lists throughout document
- ✅ **Quiz Enhancements:** All 8 Quick Checks + Final Quiz include detailed explanations

### Session 2B: Gromacs Workflow (November 2025)
- ✅ **Quiz Expansion:** 18 new questions added
  - Quick Checks 1-8: 1-4 questions each (total 26)
  - Final Quiz: 6 questions (increased from 4)
- ✅ **LiaScript Compatibility:** All quizzes converted to working format
  - Checkbox syntax: `[[ ]]` unchecked, `[[X]]` checked
  - Blank lines between questions and answer options
  - `<details>` hints for additional learning support
- ✅ **Pedagogical Content:**
  - Detailed explanations for correct answers
  - Clear reasoning for why incorrect options are wrong
  - Practical troubleshooting scenarios (e.g., "simulation crashes...")
  - Calculation-based questions (time, data points)
  - Real-world context for professional MD workflows
- ✅ **Formatting:** Fixed 12 missing blank lines before list items

### Quiz Syntax Fixed
- ✅ Updated 70+ quiz questions across all sessions
- ✅ Converted from broken `[[?]]|` format to working checkbox syntax
- ✅ All quizzes now render properly in LiaScript

### Calculations Completed

**Geometry Optimizations (UFF):**
- Glucose: 0.526197 Eh (50 steps)
- Sucrose: 1.052744 Eh (113 steps)
- Fructose: 0.310676 Eh (70 steps)

**GFN2 Validation Energies:**
- Glucose: -43.359490 Hartree
- Sucrose: -81.655515 Hartree
- Fructose: -39.275704 Hartree

**MD Simulations (Glucose, 0.1 ps test):**
- Berendsen: T_final = 315.4 K
- CSVR: T_final = 322.6 K

### Results Location

```
calculations/
├── glucose.opt.xyz, glucose.trj.xyz
├── sucrose.opt.xyz, sucrose.trj.xyz
├── fructose.opt.xyz, fructose.trj.xyz
├── glucose_berendsen_md.log
├── glucose_csvr_md.log
└── RESULTS_SUMMARY.md
```

---

## 📁 Project Structure

```
TeacherTwinMolecularModelling/
├── README.md                      # Course overview & placeholder guide
├── claude.md                      # This file
├── CALCULATION_RESULTS.md         # Detailed computational results
├── optim.md                       # Reference: improved quiz syntax
├── precourse_linux_basics.md      # Terminal tutorial
├── session1_geometry_optimization.md
├── session2a_md_curcuma.md
├── session2b_gromacs_intro.md
├── session3_alphafold_and_md.md
├── BashRPG.md                     # Bonus: game-based learning
├── strukturen/                    # Example structures
│   ├── aldehydo-D-glucose.xyz
│   ├── b-D-Glucopyranose.xyz
│   └── D-Sucrose.xyz
└── calculations/                  # Computed results
    ├── *.opt.xyz (optimized)
    ├── *.trj.xyz (trajectories)
    └── *.log (simulation data)
```

---

## 🚀 Quick Start

### Run a Geometry Optimization

```bash
cd calculations/
curcuma -opt glucose.xyz -method uff
# Output: glucose.opt.xyz, glucose.trj.xyz
```

### Inspect Results

```bash
avogadro glucose.trj.xyz &
# Watch optimization in real-time
```

### Extract Data

```bash
grep "FINAL SINGLE POINT ENERGY" glucose_uff.log
# Get final energy
```

---

## ℹ️ Important Notes

### Method Selection (Updated November 2025)
- **GFN-FF:** ✅ **RECOMMENDED DEFAULT** for all molecules <1000 atoms
  - Quantum-derived parameters, superior accuracy for heteroatoms
  - Standard in modern computational chemistry research
  - Fully functional for glucose, fructose, and peptides
- **UFF:** Fallback option for molecules >1000 atoms (computational necessity)
  - Use only when GFN-FF becomes computationally intractable
  - Documented as exception, not standard practice
- **Sucrose:** Topology issue with GFN-FF (use GFN2 instead: -81.662027 Hartree)
- **Method Comparison:** See METHOD_COMPARISON_ANALYSIS.md for GFN2 singlepoint analysis

### General Parameters
- **Energy Units:** Results in Eh (Hartree), Eh = 627.51 kcal/mol
- **MD Duration:** Test runs are 0.1 ps; production should be ≥100 ps
- **Peptide Structures:** AAAAA (53 atoms) and WRKLQ (107 atoms) now included

---

## 🔧 Troubleshooting

| Problem | Solution |
|---------|----------|
| "obabel: Cannot convert" | Check SMILES syntax with ChemSpider |
| "Optimization didn't converge" | Visualize trajectory to check structure quality |
| "curcuma not found" | Add `/home/conrad/bin/` to PATH |
| Quiz questions don't display | Ensure using LiaScript-compatible markdown |

---

## 📝 Citation

If you use this course or materials in research:

```bibtex
@software{tubaf_molecular_2025,
  title={TeacherTwin: Molecular Modelling Course},
  author={Course Team},
  year={2025},
  institution={TU Freiberg},
  note={Self-paced learning with Curcuma \& Gromacs}
}
```

---

## 👥 Support

- **Curcuma:** https://github.com/conradhuebler/curcuma
- **Gromacs:** https://www.gromacs.org/
- **Avogadro:** https://avogadro.cc/
- **ChemSpider:** https://www.chemspider.com/

---

*Generated with Claude Code - Molecular Modelling Education Platform*
