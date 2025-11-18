<!--
author: Molecular Modelling Course Team
language: en
narrator: US English Female
version: 1.0

Session 2B: Introduction to Gromacs
Part of: Molecular Modelling and Quantum Chemistry (Master)
-->

# Session 2B: Introduction to Gromacs — Professional Molecular Dynamics

## Welcome to the Industry Standard

> **Session 2A** taught you MD concepts with **Curcuma** (simple, educational).
>
> **Session 2B** introduces **Gromacs** — the most widely used MD software in research.
>
> We follow the **Lemkul workflow** (Lemkul 2024, *J. Phys. Chem. B*):
> a proven, step-by-step approach to running protein simulations.

---

## 🎯 Learning Objectives

By the end of Session 2B, you will:

- ✅ Understand **Gromacs file formats** (.mdp, .top, .gro, .tpr, .edr, .xtc)
- ✅ Know the **Lemkul workflow:** prep → energy min → equil → production
- ✅ Learn to prepare input files for **production-scale** simulations
- ✅ Submit jobs to a cluster (you send inputs, instructor runs them)
- ✅ Analyze Gromacs results (energy, temperature, trajectory)

---

## Part 1️⃣: Gromacs Overview

### What is Gromacs?

**Gromacs** = Groningen Machine for Chemical Simulations

- Powerful open-source MD package
- Used in > 50,000 publications
- Highly optimized for performance
- Supports multiple force fields (CHARMM36, Amber, OPLS, etc.)
- Standard in pharma, academia, biochemistry

**Key advantage:** Gromacs does everything (prep, equil, production, analysis) in one ecosystem.

---

### When to Use Gromacs vs. Curcuma

| Task | Curcuma | Gromacs |
|------|---------|---------|
| **Teaching MD concepts** | ✓✓✓ | ✓✓ |
| **Quick prototyping** | ✓✓✓ | ✓✓ |
| **Production protein MD** | ✓ | ✓✓✓ |
| **Large-scale simulations** | ✗ | ✓✓✓ |
| **Advanced sampling** | ✓ (partial) | ✓✓✓ |
| **Ease of use** | ✓✓✓ | ✓✓ |

**In this course:** You'll prepare Gromacs inputs (which you understand from Curcuma), then the instructor runs production simulations on a cluster.

---

### ✅ Quick Check 1: Gromacs Use Cases

**Question 1:** Which task is Gromacs best suited for?

- [[ ]] Teaching MD to beginners
- [[X]] Running 100 ns simulations of proteins on a cluster
- [[ ]] Quick energy minimization tests
- [[ ]] Simple molecules in vacuum
********************
**Explanation:**
Gromacs is **production-scale MD software**, designed for large, long simulations. **Why the answers:**

- ✅ **Correct:** 100 ns protein simulations on clusters require specialized software with parallelization, energy efficiency, and advanced analysis tools—Gromacs excels here.
- ❌ Teaching beginners: Use **Curcuma** or **xtb** instead (simpler, more intuitive).
- ❌ Quick tests: **Curcuma** or **xtb** are faster for prototyping; Gromacs has setup overhead.
- ❌ Small Molecules in vacuum: For single molecules without solvation, **QM/SQM programs (ORCA, Gaussian, xtb, GAMESS)** or **Curcuma** are better choices. Gromacs expects explicit solvent for realistic dynamics.

**Key insight:** Gromacs is overkill for small systems but essential for production biochemistry simulations.

********************

**Question 2:** Which is a disadvantage of Gromacs compared to Curcuma?

- [[X]] Steeper learning curve with more setup complexity
- [[ ]] Much less accurate results
- [[ ]] Slower for small molecules
- [[ ]] Cannot handle proteins
********************
**Explanation:**

- ✅ **Correct:** Gromacs requires understanding of MDP files, topology generation, solvation, force fields, and ensemble control. Curcuma: single command with reasonable defaults.
- ❌ Less accurate: Gromacs is more accurate as it has well tested force fields.
- ❌ Slower for small molecules: **Curcuma** or **xtb** *are* faster for prototyping, but that's a limitation of starting overhead, not speed.
- ❌ Cannot handle proteins: Gromacs is designed specifically for proteins and other biomolecules.

**Practical tip:** Curcuma for teaching + exploration; Gromacs for publication-ready results.

********************

**Question 3:** Why does this course use BOTH Curcuma and Gromacs?

- [[X]] Learn MD concepts with Curcuma, then apply them professionally with Gromacs
- [[ ]] They use different force fields, and are complementary
- [[ ]] Gromacs cannot do geometry optimization
- [[ ]] Curcuma is required before Gromacs
********************

**Explanation:**
- ✅ **Correct:** Session 2A (Curcuma) builds conceptual understanding; Session 2B (Gromacs) shows how real researchers run simulations. This progression ensures you understand *why* before using complex *how*.

- ❌ Different force fields: While this is true, that is not the key of using both. Curcuma is for small molecules, with ability to use force fields and semiemprical methods. But curcuma is used as easy way to get to know molecular dynamics
- ❌ Gromacs cannot do optimization: Gromacs can minimize, NVT, NPT, and production. Limitation is that this course focuses you on prep → submit workflow.
- ❌ Curcuma required: Technically not required, but pedagogically essential for understanding.

********************

## Part 2️⃣: Gromacs File Formats

Gromacs works with specific file types. Understanding them is **critical**.

### Input Files

| File | Purpose | Format | Example |
|------|---------|--------|---------|
| **.pdb** | Protein structure (from PDB or AlphaFold) | Text | `protein.pdb` |
| **.top** | System topology (atoms, bonds, interactions) | Text | `topol.top` |
| **.itp** | Included topology (force field parameters) | Text | `amber99sb.itp` |
| **.mdp** | MD parameters (all simulation settings) | Text | `minim.mdp` |

### Processing Files

| File | Purpose | Format | Created by |
|------|---------|--------|------------|
| **.gro** | Coordinate file (XYZ with velocities, box) | Text (fixed format) | `gmx pdb2gmx` |
| **.tpr** | Binary run input (combines gro + top + mdp) | Binary | `gmx grompp` |

### Output Files

| File | Purpose | Format | Size |
|------|---------|--------|------|
| **.xtc** | Trajectory (compressed, no velocities) | Binary | Large (~1-100 GB) |
| **.trr** | Trajectory (full precision, with velocities) | Binary | Huge (~10-1000 GB) |
| **.edr** | Energy/temperature/pressure data | Binary | Analyzed with `gmx energy` |
| **.log** | Simulation log (human-readable) | Text | Small (~1-10 MB) |

---

### ✅ Quick Check 2: File Formats

**Question 1:** Which file contains the force field parameters?

- [[ ]] .pdb
- [[ ]] .gro
- [[X]] .top / .itp
- [[ ]] .xtc
********************

**Explanation:**

- ✅ **Correct:** Topology files (.top = main topology, .itp = included topology modules) contain force field parameters: bond strengths, angles, dihedrals, nonbonded interactions. .itp files are often imported from force field databases (AMBER, CHARMM).
- ❌ .pdb: Just coordinates + element types (from PDB database or AlphaFold).
- ❌ .gro: Coordinates + velocities, created *after* topology processing.
- ❌ .xtc: Compressed trajectory (coordinates only, no parameters).

**Key insight:** .top is the "recipe" for your molecule; .gro is just the "photo" of it.

********************

**Question 2:** What does .tpr stand for?

- [[X]] Topology, Positions, and Run parameters
- [[ ]] Trajectory print result
- [[ ]] Total potential representation
- [[ ]] Time-point recording
********************

**Explanation:**

- ✅ **Correct:** .tpr is a **binary preprocessed file** created by `gmx grompp`. It combines:

  1. **Topology** (bonds, parameters from .top)
  2. **Positions** (coordinates from .gro)
  3. **Run parameters** (everything from .mdp file)
  This single file is fed to `gmx mdrun`.
- ❌ Trajectory result: That would be .trr or .xtc (created *after* running mdrun).
- ❌ Total potential: Confusing terminology; .tpr is binary, not a potential-specific file.
- ❌ Time-point recording: Misunderstands file purpose entirely.

**Practical:** Think of .tpr as the "ready-to-run" package. It's binary because Gromacs optimizes it for fast reading during simulation.

********************

**Question 3:** When would you use .trr instead of .xtc?

- [[X]] When you need velocities and forces (full precision MD data)
- [[ ]] Always, because .trr is smaller
- [[ ]] Never, .xtc is always better
- [[ ]] Only for proteins, not other molecules
********************

**Explanation:**

- ✅ **Correct:**

  - **.trr** = full precision, includes velocities & forces, HUGE file size (10-1000 GB for long sims). Use when:
  
    - Restarting simulation from specific frame (need velocities)
    - Advanced analysis requiring forces
    - Publication-quality archival
  - **.xtc** = compressed, coordinates only, manageable size (1-100 GB). Standard for distribution & analysis.
- ❌ .trr always smaller: Opposite! .trr is 10× larger than .xtc.
- ❌ .xtc never good: .xtc is the standard for analysis; .trr is backup/restart.
- ❌ Only for proteins: Both work for any system.

**Real practice:** Most papers use .xtc for analysis, .trr as archive.

********************

## Part 3️⃣: The Lemkul Workflow

The **Lemkul workflow** (from the 2024 paper) is a proven, step-by-step process:

```
1. PREPARATION
   Input: protein.pdb
   ↓
   gmx pdb2gmx  →  topology (.top), gro file
   
2. ENERGY MINIMIZATION
   Input: minim.mdp, topol.top, conf.gro
   ↓
   gmx grompp  →  em.tpr
   gmx mdrun   →  em.trr, em.edr, em.log
   
3. EQUILIBRATION (NVT)
   Input: nvt.mdp, topol.top, em.gro
   ↓
   gmx grompp  →  nvt.tpr
   gmx mdrun   →  nvt.trr, nvt.edr, nvt.log
   
4. EQUILIBRATION (NPT)
   Input: npt.mdp, topol.top, nvt.gro
   ↓
   gmx grompp  →  npt.tpr
   gmx mdrun   →  npt.trr, npt.edr, npt.log
   
5. PRODUCTION MD
   Input: prod.mdp, topol.top, npt.gro
   ↓
   gmx grompp  →  prod.tpr
   gmx mdrun   →  prod.trr, prod.edr, prod.log
   
6. ANALYSIS
   Input: prod.edr, prod.xtc
   ↓
   gmx energy   →  energy.xvg (plot)
   gmx rms      →  rmsd.xvg (plot)
   gmx rmsf     →  rmsf.xvg (per-residue flexibility)
```

Each step generates input for the next step. **Very systematic.**

---

### Key Workflow Principles

1. **Energy minimization first** — Remove bad contacts, relax starting structure
2. **NVT equilibration** — Heat system to target T, keep volume constant
3. **NPT equilibration** — Allow volume to equilibrate, adjust box size
4. **Production run** — Long, unbiased simulation at correct conditions
5. **Analysis** — Extract energies, RMSDs, fluctuations from outputs

---

### ✅ Quick Check 3: Workflow Order

**Question 1:** The correct order of Lemkul workflow is:

- [[X]] Minimize → NVT equil → NPT equil → Production
- [[ ]] NVT → Minimize → NPT → Production
- [[ ]] Production → Minimize → Equilibration
- [[ ]] Minimize → Production directly
********************

**Explanation:**

- ✅ **Correct:** This is the thermodynamic sequence:

  1. **Minimize:** Remove bad contacts, relax structure
  2. **NVT:** Heat system to 300K at constant volume (temperature equilibration)
  3. **NPT:** Allow box to adjust to 1 bar at constant temperature (pressure equilibration)
  4. **Production:** Long, unbiased simulation at correct T & P
- ❌ NVT before minimize: Impossible—system still has strain; NVT would crash. Why?
- ❌ Production first: Meaningless—system isn't equilibrated.
- ❌ Skip NPT: Box size would be wrong; pressure artifacts in results.

**Key:** Each step feeds the next. Skip one = simulation fails or unphysical results.

********************

**Question 2:** What happens if you skip energy minimization?

- [[X]] Bad contacts cause simulation to crash or become unstable
- [[ ]] Nothing—it's optional, just slower
- [[ ]] Simulation runs faster
- [[ ]] You cannot create .tpr file
********************

**Explanation:**

- ✅ **Correct:** Unminimized structures have atoms overlapping or too close:

  - Repulsive forces become **enormous** (van der Waals)
  - Simulation explodes (atoms fly apart)
  - NVT/NPT cannot equilibrate
  - Common error message: "Steepest descent not converged," "simulation blew up"
- ❌ It's optional: Wrong—minimization is mandatory in every published workflow.
- ❌ Runs faster: Missing minimization adds 10-20 failed attempts; total time increases.
- ❌ Cannot create .tpr: grompp will still create .tpr even with bad geometry; problem appears during mdrun.

**Lesson:** Minimization is cheap (minutes); failed simulations are expensive (hours wasted).

********************

**Question 3:** Why must NVT equilibration run BEFORE NPT?

- [[X]] System must reach target temperature first, then adjust volume (pressure)
- [[ ]] NPT is more expensive computationally
- [[ ]] The order doesn't matter; they're independent
- [[ ]] NVT creates files that NPT depends on
********************

**Explanation:**

- ✅ **Correct:** Thermodynamic hierarchy:

  - **NVT (constant V):** System can't expand/contract; all heat input increases temperature
  - **NPT (constant P & T):** If you skip NVT, system tries to equilibrate T *and* P simultaneously = chaotic, often unstable
  - **Proper sequence:** NVT first → system reaches 300K → *then* NPT adjusts box → both T & P stable
  - Analogy: Fill balloon with air (heat, NVT), then let it reach pressure equilibrium (NPT).
- ❌ NPT more expensive: Actually the other way—NPT barostat adds some cost, not the reason for order.
- ❌ Order doesn't matter: Would immediately see pressure instability, temperature overshoot.
- ❌ NVT creates files for NPT: Both generate .gro files; order is physics-driven, not file-driven.

**Pro tip:** If your NPT crashes, 90% of the time it's because initial structure wasn't properly NVT-equilibrated.

********************

## Part 4️⃣: MDP Files — The Heart of Gromacs

The **.mdp file** controls everything about your simulation. It's **crucial**.

### Example: Energy Minimization (minim.mdp)

```
; minim.mdp - Steepest Descent Energy Minimization
; Run parameters
integrator              = steep      ; steepest descent
emtol                   = 100.0      ; stop when F < 100 kJ/mol/nm
emstep                  = 0.01       ; initial step size (nm)
nsteps                  = 50000      ; max iterations

; Output control
nstenergy               = 500        ; save energy every 500 steps
nstlog                  = 500        ; save log every 500 steps
nstxout-compressed      = 500        ; save trajectory every 500 steps
```

**Key parameters:**

- `integrator`: Algorithm (steep = steepest descent, CG = conjugate gradient)
- `emtol`: Convergence criterion (energy gradient threshold)
- `nsteps`: Maximum iterations
- `nstenergy`: Output frequency

---

### Example: NVT Equilibration (nvt.mdp)

```
; nvt.mdp - NVT Equilibration (300 K, 100 ps)
; Run parameters
integrator              = md         ; normal MD
dt                      = 0.002      ; time step (2 fs)
nsteps                  = 50000      ; 2 fs × 50000 = 100 ps
nstcomm                 = 100        ; remove COM motion every 100 steps

; Output control
nstenergy               = 1000       ; save energy every 1000 steps (2 ps)
nstlog                  = 1000       
nstxout-compressed      = 1000       ; trajectory every 2 ps

; Temperature coupling
tcoupl                  = V-rescale  ; thermostat (like CSVR)
tau-t                   = 0.1        ; coupling time (ps)
ref-t                   = 300        ; target temperature (K)
tc-grps                 = System     ; couple entire system

; Pressure coupling (constant V)
pcoupl                  = no         ; NO pressure coupling
```

**Key differences from minimization:**

- `integrator = md` (instead of steep)
- `dt = 0.002` (time step in ps, usually 2 fs)
- `nsteps = 50000` (run duration)
- `tcoupl = V-rescale` (CSVR-like thermostat!)
- `pcoupl = no` (constant volume)

---

### Example: NPT Equilibration (npt.mdp)

```
; npt.mdp - NPT Equilibration (300 K, 1 bar, 100 ps)
; (Same as nvt.mdp, but with pressure coupling)

integrator              = md
dt                      = 0.002
nsteps                  = 50000

nstenergy               = 1000
nstlog                  = 1000
nstxout-compressed      = 1000

; Temperature coupling (same as NVT)
tcoupl                  = V-rescale
tau-t                   = 0.1
ref-t                   = 300
tc-grps                 = System

; Pressure coupling (allow volume to adjust)
pcoupl                  = Parrinello-Rahman  ; barostat
pcoupltype              = isotropic          ; uniform scaling
tau-p                   = 2.0                ; pressure coupling time
ref-p                   = 1.0                ; target pressure (bar)
compressibility         = 4.5e-5             ; water compressibility
```

**New parameters:**

- `pcoupl = Parrinello-Rahman` — Barostat (pressure control)
- `ref-p = 1.0` — Target pressure (1 bar = atmospheric)
- `tau-p = 2.0` — Pressure coupling time

---

### Example: Production Run (prod.mdp)

```
; prod.mdp - Production MD (300 K, 1 bar, 500 ps)
; Most parameters same as NPT, but longer

integrator              = md
dt                      = 0.002
nsteps                  = 250000     ; 2 fs × 250000 = 500 ps (longer!)

nstenergy               = 1000
nstlog                  = 1000
nstxout-compressed      = 500        ; trajectory every 1 ps

; Temperature coupling
tcoupl                  = V-rescale
tau-t                   = 0.1
ref-t                   = 300
tc-grps                 = System

; Pressure coupling
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau-p                   = 2.0
ref-p                   = 1.0
compressibility         = 4.5e-5
```

**Key difference:** Longer `nsteps` (500 ps instead of 100 ps).

---

### ✅ Quick Check 4: MDP Parameters

**Question 1:** What does `dt = 0.002` mean?

- [[ ]] 0.002 seconds per step
- [[X]] 0.002 picoseconds per step (2 femtoseconds)
- [[ ]] 0.002 nanometers
- [[ ]] Not important
********************

**Explanation:**

- ✅ **Correct:** dt is **time per integration step**. 0.002 ps = 2 fs. Physics requires small timesteps to accurately integrate forces without errors.
- ❌ Seconds: That would be absurdly long (2 million femtoseconds per step!).
- ❌ Nanometers: That's distance, not time. Confusing units.
- ❌ Not important: Wrong—timestep is critical. Too large → simulation explodes; too small → simulation crawls.

**Key:** dt typically ranges 0.5–2 fs for proteins. 2 fs is standard (fast bond vibrations ~0.01 fs, so 2 fs captures them safely).

********************

**Question 2:** In NVT ensemble, which should be "no"?

- [[X]] pcoupl (pressure coupling)
- [[ ]] tcoupl (temperature coupling)
- [[ ]] nstenergy (output frequency)
- [[ ]] integrator
********************

**Explanation:**

- ✅ **Correct:** NVT = **constant Volume**. So:

  - **pcoupl = no** (no pressure control allowed)
  - **tcoupl = MUST be on** (need thermostat to maintain 300K)
  - Volume is fixed → box size doesn't change
- ❌ tcoupl: Wrong—you NEED thermostat in NVT to control temperature.
- ❌ nstenergy: That's output frequency, unrelated to NVT constraint.
- ❌ integrator: That's the algorithm choice (md, steep), doesn't define NVT.

**Tip:** Remember: **NVT = No Pressure coupling. NPT = Pressure coupling ON.**

********************

**Question 3:** If `dt=0.002` ps and you want a 200 ps simulation, what is `nsteps`?

- [[ ]] 200
- [[ ]] 2,000
- [[ ]] 50,000
- [[X]] 100,000
********************

**Explanation:**

- ✅ **Correct:** nsteps = duration / dt = 200 ps / 0.002 ps = **100,000 steps**

  - Verification: 100,000 steps × 0.002 ps/step = 200 ps ✓
- ❌ 200: Confusing nsteps with total time (way too short!).
- ❌ 2,000: Off by factor of 50 (would give only 4 ps).
- ❌ 50,000: Off by factor of 2 (would give 100 ps, not 200 ps).

**Pro skill:** Always verify: **Actual time = nsteps × dt**. This calculation appears in every workflow.

********************

**Question 4:** What does `nstenergy = 1000` control?

- [[X]] Energy data is saved every 1000 steps (every 2 ps)
- [[ ]] Maximum energy allowed in simulation
- [[ ]] Timestep size
- [[ ]] Force calculation frequency
********************

**Explanation:**

- ✅ **Correct:** **nstenergy** = output frequency for energy file (.edr). With dt=0.002:

  - 1000 steps × 0.002 ps/step = 2 ps per energy output
  - You get ~100 energy points for a 200 ps run
  - Too frequent (nstenergy=10) → huge .edr file; too rare (nstenergy=10000) → can't see behavior detail
- ❌ Maximum energy: That would be `emtol` (minimization threshold), not nstenergy.
- ❌ Timestep size: That's `dt`, not nstenergy.
- ❌ Force calculation: Forces recalculated every step regardless of nstenergy; this just controls output.

**Practical:** Use nstenergy=1000 (default). For analysis, you need ~100-200 data points; this gives the right resolution.

********************

## Part 5️⃣: Gromacs Commands — The Workflow in Practice

### Step 1: Prepare Structure (pdb2gmx)

**Input:** Protein structure from PDB or AlphaFold

```bash
gmx pdb2gmx -f protein.pdb -o conf.gro -p topol.top -ff amber99sb
```

**Parameters:**

- `-f protein.pdb` — Input protein file
- `-o conf.gro` — Output coordinate file
- `-p topol.top` — Output topology
- `-ff amber99sb` — Force field (amber99sb, charmm36, etc.)

**Output files:**

- `conf.gro` — Protein coordinates (gromacs format)
- `topol.top` — System topology with bonding info
- `posre.itp` — Position restraints file

---

### Step 2: Add Water Box (editconf + solvate)

```bash
# Define box around protein
gmx editconf -f conf.gro -o conf_newbox.gro -c -d 1.0 -bt cubic

# Add water molecules
gmx solvate -cp conf_newbox.gro -cs spc216.gro -o conf_solv.gro -p topol.top
```

**Result:** `conf_solv.gro` — Protein + water molecules

---

### Step 3: Add Ions (genion)

```bash
# Create tpr for ion addition
gmx grompp -f ions.mdp -c conf_solv.gro -p topol.top -o ions.tpr

# Add counterions
gmx genion -s ions.tpr -o conf_ions.gro -p topol.top -pname NA -nname CL -neutral
```

**Result:** Neutralized system ready for minimization

---

### Step 4: Energy Minimization

```bash
# Prepare minimization tpr file
gmx grompp -f minim.mdp -c conf_ions.gro -p topol.top -o em.tpr

# Run minimization
gmx mdrun -v -deffnm em
```

**Input:**

- `minim.mdp` — Minimization parameters
- `topol.top` — Topology
- `conf_ions.gro` — Starting structure

**Output:**

- `em.tpr` — Binary input
- `em.gro` — Minimized structure
- `em.edr` — Energy data
- `em.trr` — Trajectory

---

### Step 5: NVT Equilibration

```bash
gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt
```

**Input:**

- `nvt.mdp` — NVT parameters
- `em.gro` — Minimized structure from previous step

**Output:** `nvt.gro`, `nvt.edr`, `nvt.trr`

---

### Step 6: NPT Equilibration

```bash
gmx grompp -f npt.mdp -c nvt.gro -p topol.top -t nvt.cpt -o npt.tpr
gmx mdrun -v -deffnm npt
```

**Note:** `-t nvt.cpt` continues from checkpoint (preserves velocities)

**Output:** `npt.gro`, `npt.edr`, `npt.trr`

---

### Step 7: Production MD

```bash
gmx grompp -f prod.mdp -c npt.gro -p topol.top -t npt.cpt -o prod.tpr
gmx mdrun -v -deffnm prod
```

**This is the LONG run** — takes hours or days depending on system size.

**Output:** `prod.gro`, `prod.edr`, `prod.xtc` (compressed trajectory)

---

### ✅ Quick Check 5: Workflow Commands

**Question 1:** The command `gmx grompp` does what?

- [[ ]] Runs MD simulation
- [[X]] Preprocesses inputs (.mdp, .top, .gro) and creates binary .tpr file
- [[ ]] Analyzes trajectory results
- [[ ]] Adds water to the system
********************

**Explanation:**

- ✅ **Correct:** `grompp` = **GROMacs PreProcessor**. It:

  1. Reads .mdp file (simulation settings)
  2. Reads .top file (topology)
  3. Reads .gro file (coordinates)
  4. Creates .tpr file (combined binary ready for mdrun)
  - It validates inputs, checks for errors, pre-processes topology
  - Cannot run without .tpr file from grompp
- ❌ Runs MD: That's `gmx mdrun`, not grompp.
- ❌ Analyzes: That's `gmx energy`, `gmx rms`, etc.
- ❌ Adds water: That's `gmx solvate`, not grompp.

**Key:** `grompp` = prep, `mdrun` = execute.

********************

**Question 2:** What does `gmx mdrun` do?

- [[X]] Runs the actual MD simulation (reads .tpr, creates trajectory)
- [[ ]] Creates topology files
- [[ ]] Analyzes results
- [[ ]] Preprocesses inputs
********************

**Explanation:**

- ✅ **Correct:** `mdrun` = **MD RUN**. It:

  1. Reads .tpr file (from grompp)
  2. Performs integration loop: calculate forces → move atoms → repeat
  3. Produces: .trr/.xtc (trajectory), .edr (energy), .log (output), .cpt (checkpoint)
  - mdrun is the CPU-intensive bottleneck
  - On clusters, mdrun is parallelized (GPU/CPU)
- ❌ Creates topology: Topology created during grompp, not mdrun.
- ❌ Analyzes: That's separate gmx commands (energy, rms, etc.).
- ❌ Preprocesses: That's grompp.

**Pro:** mdrun takes hours/days; grompp takes seconds. If mdrun crashes, problem is usually in structure or .tpr, not mdrun itself.

********************

**Question 3:** If `gmx grompp` fails, which file is most likely missing?

- [[X]] .mdp or .top file
- [[ ]] .xtc trajectory file
- [[ ]] .edr energy file
- [[ ]] .log output file
********************

**Explanation:**

- ✅ **Correct:** grompp **needs inputs**: .mdp, .top, .gro. Missing any one → error.

  - Missing .mdp: "cannot find minim.mdp"
  - Missing .top: "cannot find topol.top"
  - Missing .gro: "cannot find conf.gro"
- ❌ .xtc, .edr, .log: These are **outputs** created *after* mdrun, not needed by grompp.

**Troubleshooting tip:** "grompp failed" → check inputs exist first (`ls -la minim.mdp topol.top conf.gro`).

********************

## Part 6️⃣: Your Role in This Course

### What YOU Will Do

1. **Prepare inputs** (MDP + topology using templates)
2. **Create `.tpr` files** (using `gmx grompp`)
3. **Send `.tpr` + `.gro` files** to the instructor
4. **Instructor runs on cluster** (Gromacs mdrun)
5. **Get back results** (.xtc, .edr, .log)

### MDP Templates (Copy-Paste Ready)

**You'll get templates like these** — fill in blanks, use for your protein:

```bash
# Template: energy minimization
gmx grompp -f minim.mdp -c conf_ions.gro -p topol.top -o em.tpr

# Template: NVT (100 ps, 300 K)
gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr

# Template: NPT (100 ps, 300 K, 1 bar)
gmx grompp -f npt.mdp -c nvt.gro -p topol.top -t nvt.cpt -o npt.tpr

# Template: Production (500 ps, 300 K, 1 bar)
gmx grompp -f prod.mdp -c npt.gro -p topol.top -t npt.cpt -o prod.tpr
```

**Your submission:**

```
session3_inputs/
├── protein_name.pdb          (from AlphaFold or PDB)
├── minim.mdp
├── nvt.mdp
├── npt.mdp
├── prod.mdp
└── em.tpr, nvt.tpr, npt.tpr, prod.tpr  (created by gmx grompp)
```

---

### Instructor's Role

1. Takes your `.tpr` files
2. Runs production MD on cluster
3. Returns `.xtc`, `.edr`, `.log`
4. You analyze the results in Session 3

**Why this division?** MD runs are **computationally expensive**. You learn the workflow; cluster runs the heavy computation.

---

### ✅ Quick Check 6: Your Workflow

**Question 1:** In this course, what do you do with Gromacs?

- [[X]] Prepare inputs and create .tpr files
- [[ ]] Run full production simulations
- [[ ]] Only analyze results
- [[ ]] Nothing (only Curcuma)
********************

**Explanation:**

- ✅ **Correct:** Your role is **input preparation**:

  1. Write/modify MDP files
  2. Run `gmx grompp` to create .tpr files
  3. Submit .tpr + .gro files to instructor
  4. Instructor runs `gmx mdrun` on cluster
  5. You analyze results (Session 3)
- ❌ Run full simulations: That's on the cluster (computationally expensive, requires HPC resources).
- ❌ Only analyze: You need to understand prep—analysis is just the final step.
- ❌ Nothing: Gromacs input preparation is 30% of the learning.

**Why this split?** Production MD for 500+ ps takes 8-48 CPU hours. Teaching lab computers can't spare that. You learn the workflow; cluster runs the heavy compute.

********************

**Question 2:** What files do you submit to the instructor?

- [[X]] .tpr files (and MDP templates for documentation)
- [[ ]] Only .pdb files
- [[ ]] .xtc trajectory files
- [[ ]] Analyzed results
********************

**Explanation:**

- ✅ **Correct:** You submit:

  - **.tpr files** (em.tpr, nvt.tpr, npt.tpr, prod.tpr) — ready to run with mdrun
  - **.mdp files** (optional but good for documentation)
  - **.gro file** (em.gro or nvt.gro, starting coords for next step)
  - **.top file** (topology, so instructor can verify system is correct)
- ❌ .pdb: Already used to create .gro; not needed for submission.
- ❌ .xtc: You don't have it yet—it's created *by* mdrun on the cluster.
- ❌ Analyzed results: That's for *after* you get results back.

**Submission checklist:** `ls -la *.tpr *.mdp *.gro topol.top` — should see 4 .tpr files + topology.

********************

**Question 3:** Why doesn't the course have you run mdrun?

- [[X]] Computationally expensive; requires HPC cluster resources
- [[ ]] mdrun is too difficult to use
- [[ ]] Gromacs is unavailable on your computer
- [[ ]] It's not important for learning
********************

**Explanation:**

- ✅ **Correct:** Production MD is computationally expensive:

  - 500 ps simulation = 250,000 timesteps
  - Each timestep = calculate forces for ~100,000 atoms
  - Total: ~25 billion force calculations
  - On single CPU: 24-48 hours
  - On cluster GPU: 2-8 hours (still significant)
  - Teaching lab: finite compute budget; can't waste on individual runs
- ❌ Too difficult: mdrun is actually simple—`gmx mdrun -v -deffnm prod`—one command!
- ❌ Unavailable: Gromacs can be on your computer; limitation is CPU time, not access.
- ❌ Not important: This is where the real computation happens; you're missing nothing conceptually, just not waiting 8 hours.

**Industry practice:** You prep inputs locally; submit to HPC cluster; get results back. You're learning the real workflow.

********************

## Part 7️⃣: Parameter Guide — Choosing Values

### Time Step (dt)

**Rule of thumb:** Smaller dt = more stable, but slower

- **2 fs (0.002 ps)** — Default for proteins, usually stable
- **1 fs (0.001 ps)** — If system unstable or RATTLE constraints
- **0.5 fs** — Only if serious problems

**For this course:** Use **2 fs** (standard).

---

### Coupling Times

**Thermostat coupling time (tau-t):**

- Small (0.05 ps) — Strong coupling, artificial but stable
- Medium (0.1 ps) — Good balance
- Large (0.5 ps) — Weak coupling, more realistic but slower

**For this course:** Use **0.1 ps** (recommended in literature).

---

### Run Duration

**How long should you simulate?**

| System | Typical Duration |
|--------|------------------|
| Small peptide (5 AA) | 10-50 ps |
| Medium protein (100 AA) | 100-500 ps |
| Large protein (500 AA) | 500 ps - 1 ns |
| Protein complex | 1-10 ns |
| Real production | 10 ns - 1 µs |

**For this course:**

- Equilibration: **100 ps each** (NVT + NPT)
- Production: **[EXPECTED_PRODUCTION_LENGTH_PLACEHOLDER]** ps

---

### Temperature & Pressure

**Temperature:** Always **300 K** (room temperature, ~27°C)

**Pressure:** Always **1.0 bar** (1 atmosphere)

---

### ✅ Quick Check 7: Parameters

**Question 1:** Which dt is most stable for protein MD?

- [[ ]] 0.1 ps
- [[X]] 0.002 ps (2 fs)
- [[ ]] 0.5 ps
- [[ ]] 0.0001 ps (overkill)
********************

**Explanation:**

- ✅ **Correct:** 0.002 ps (2 fs) is the **standard** for protein MD:

  - Bond vibrations ~0.01 fs frequency → need dt << 0.01 fs to capture
  - 2 fs = 200× smaller than frequency → safe margin
  - Faster than needed (gives room for error), but not wasteful
- ❌ 0.1 ps: Way too large! Forces change dramatically in 100 fs. Simulation becomes unstable, diverges.
- ❌ 0.5 ps: Even worse (500× larger than needed).
- ❌ 0.0001 ps: Technically stable but 20× overkill. Simulation runs 20× slower for no added accuracy.

**Trade-off:** dt=2fs = **Goldilocks zone** → stable, reasonable speed.

********************

**Question 2:** If `dt=0.002` ps and `nsteps=50000`, how long is the simulation?

- [[ ]] 10 ps
- [[ ]] 50 ps
- [[X]] 100 ps
- [[ ]] 250 ps
********************

**Explanation:**

- ✅ **Correct:** Duration = nsteps × dt = 50,000 × 0.002 ps = **100 ps**
  - Verification: 50,000 × 0.002 = 100 ✓
  - This is exactly NVT equilibration time in Lemkul workflow
- ❌ 10 ps: Off by factor of 10 (nsteps would be 5,000).
- ❌ 50 ps: Off by factor of 2 (nsteps would be 25,000).
- ❌ 250 ps: That would be nsteps=125,000 (production run length).

**Pro skill:** **Time = nsteps × dt**. Always verify. Mistakes here waste huge CPU time.

********************

**Question 3:** What does `tau-t = 0.1` ps control?

- [[X]] How strongly temperature is controlled (thermostat coupling strength)
- [[ ]] Duration of simulation
- [[ ]] Pressure coupling strength
- [[ ]] Output frequency
********************

**Explanation:**

- ✅ **Correct:** **tau-t** = thermostat **coupling time**:

  - Small tau-t (0.05 ps): **Strong coupling**—temperature tightly controlled, artificial but stable. Good for equilibration.
  - Medium tau-t (0.1 ps): **Balance**—realistic dynamics + stability. Recommended for this course.
  - Large tau-t (0.5 ps): **Weak coupling**—more realistic but takes longer to reach target T.
  - Analogy: Like a thermostat in your house: small response time = tight control, large = sluggish response.
- ❌ Duration: That's nsteps × dt, not tau-t.
- ❌ Pressure: That's tau-p, not tau-t.
- ❌ Output: That's nstenergy, nstlog, etc.

**Practice:** Use tau-t = 0.1 ps (default). Changing it affects how quickly system equilibrates, not length of run.

********************

**Question 4:** Decreasing tau-t (e.g., 0.1 → 0.05) means:

- [[X]] Stronger temperature control; system reaches 300K faster
- [[ ]] Weaker temperature control; system heats slower
- [[ ]] Simulation runs shorter
- [[ ]] No effect on equilibration
********************

**Explanation:**

- ✅ **Correct:** Smaller tau-t = faster response to temperature:

  - tau-t=0.05 ps: Thermostat acts every 0.05 ps → tightly constrained
  - tau-t=0.1 ps: Thermostat acts every 0.1 ps → looser control
  - Result: tau-t=0.05 reaches 300K in ~1-2 ps; tau-t=0.1 takes ~5-10 ps
  - Trade-off: Small tau-t is faster but less realistic; large tau-t is slower but more physical
- ❌ Weaker control: Opposite! Smaller tau-t = stronger, tighter coupling.
- ❌ Shorter simulation: tau-t doesn't change nsteps; duration = nsteps × dt, unchanged.
- ❌ No effect: tau-t has major effect on equilibration speed and quality.

**Decision:** For teaching (fast equilibration) → tau-t=0.1 ps (balance). For production (realistic) → tau-t=0.1-0.5 ps (literature values).

********************

## Part 8️⃣: MDP Template for Your Use

Here's a **complete set of templates** you can copy and modify:

### minim.mdp

```
; Energy minimization
integrator              = steep
emtol                   = 100.0
emstep                  = 0.01
nsteps                  = 50000

nstenergy               = 500
nstlog                  = 500
nstxout-compressed      = 500

; Non-bonded
cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
```

### nvt.mdp

```
; NVT equilibration (100 ps)
integrator              = md
dt                      = 0.002
nsteps                  = 50000

nstenergy               = 1000
nstlog                  = 1000
nstxout-compressed      = 1000

; Temperature
tcoupl                  = V-rescale
tau-t                   = 0.1
ref-t                   = 300
tc-grps                 = System

; Pressure
pcoupl                  = no

; Non-bonded
cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
```

### npt.mdp

```
; NPT equilibration (100 ps)
integrator              = md
dt                      = 0.002
nsteps                  = 50000

nstenergy               = 1000
nstlog                  = 1000
nstxout-compressed      = 1000

; Temperature
tcoupl                  = V-rescale
tau-t                   = 0.1
ref-t                   = 300
tc-grps                 = System

; Pressure
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau-p                   = 2.0
ref-p                   = 1.0
compressibility         = 4.5e-5

; Non-bonded
cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
```

### prod.mdp

```
; Production MD (500 ps)
integrator              = md
dt                      = 0.002
nsteps                  = 250000

nstenergy               = 1000
nstlog                  = 1000
nstxout-compressed      = 500

; Temperature
tcoupl                  = V-rescale
tau-t                   = 0.1
ref-t                   = 300
tc-grps                 = System

; Pressure
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau-p                   = 2.0
ref-p                   = 1.0
compressibility         = 4.5e-5

; Non-bonded
cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0
```

**All templates use:**

- Time step: **2 fs**
- Temperature: **300 K**
- Pressure: **1 bar**
- Thermostat: **V-rescale** (like CSVR from Session 2A)
- Barostat: **Parrinello-Rahman**

---

## Part 9️⃣: Key Differences from Lemkul Paper

**This course simplifies the Lemkul workflow slightly:**

| Step | Lemkul | This Course |
|------|--------|------------|
| 1. Prep | pdb2gmx, editconf, solvate, genion | [Will be provided] |
| 2. Min | Yes | Yes |
| 3. NVT | Yes (100 ps) | Yes (100 ps) |
| 4. NPT | Yes (100 ps) | Yes (100 ps) |
| 5. Prod | 1 ns | [EXPECTED_PRODUCTION_LENGTH_PLACEHOLDER] |
| 6. Analysis | gmx energy, gmx rms, gmx rmsf | You do this |

**Simplifications:**

- We skip some of Lemkul's advanced options (e.g., LINCS constraints)
- Use default Gromacs settings (more robust)
- Shorter runs for teaching (1-2 hours CPU instead of 8+ hours)

---

### ✅ Quick Check 8: Course vs Lemkul

**Question 1:** This course is based on the Lemkul workflow but:

- [[ ]] Skips energy minimization entirely
- [[ ]] Uses Berendsen thermostat only
- [[ ]] Doesn't follow any published protocol
- [[X]] Uses shorter equilibration times for teaching
********************

**Explanation:**

- ✅ **Correct:** Course follows Lemkul's proven 5-step workflow but with teaching modifications:

  - Lemkul: NVT = 100 ps, NPT = 100 ps, Production = 1 ns (appropriate for publication)
  - Course: NVT = 100 ps, NPT = 100 ps, Production = 500 ps (shorter to fit lab time)
  - Rationale: Even 500 ps gives meaningful data for learning; 1 ns is research-grade overkill
- ❌ Skips minimization: Wrong—minimization is **essential** first step.
- ❌ Uses Berendsen: We use **V-rescale** (CSVR), not Berendsen. (Berendsen is older, less accurate.)
- ❌ Doesn't follow protocol: We **strictly** follow Lemkul 2024—one of the most cited MD tutorial papers.

**Key:** Shorter runs ≠ different method. Same workflow, same physics, just adjusted for pedagogical time.

********************

**Question 2:** What does LINCS (mentioned in Lemkul) do?

- [[ ]] Analyzes trajectory results
- [[X]] Constrains bond lengths, allows larger timesteps
- [[ ]] Couples temperature and pressure
- [[ ]] Compresses trajectory files
********************

**Explanation:**
- ✅ **Correct:** **LINCS** = Linear Constraint Solver. It:

  - Fixes bond lengths at their equilibrium values (e.g., C-H always 0.109 nm)
  - Allows dt → 4-5 fs (instead of 2 fs) because fast bond vibrations don't need calculating
  - Reduces computation by ~20%
  - Trade-off: Slightly more artificial dynamics, but standard in production MD
- ❌ Couples T/P: That's thermostats and barostats, not LINCS.
- ❌ Compresses trajectory: That's .xtc format (grompp -o option), not LINCS.
- ❌ Analyzes: That's gmx energy, gmx rms, etc.

**Why we skip LINCS:** It's advanced; default approach (no constraints) is simpler for teaching. Production workflows use LINCS.

********************

**Question 3:** Why does this course use shorter production runs than Lemkul?

- [[ ]] Gromacs can't handle long simulations
- [[ ]] Shorter production is scientifically equivalent to Lemkul
- [[X]] Teaching time constraints; 500 ps still shows meaningful dynamics and equilibration
- [[ ]] Shorter runs are more accurate than longer ones
********************

**Explanation:**

- ✅ **Correct:** Practical trade-off:

  - **500 ps (course):** Shows equilibration of temperature, pressure, energy. Good for understanding dynamics. CPU: ~2 hours.
  - **1 ns (Lemkul):** Shows longer-timescale fluctuations, better statistics. CPU: ~4+ hours.
  - Both are **scientifically valid**; the difference is statistical convergence (more data = lower noise).
  - For learning the workflow: 500 ps is sufficient. For publication: 1 ns is minimum.
- ❌ Shorter = more accurate: Opposite! Longer = more data = lower noise = more reliable averages.
- ❌ Gromacs limitation: Gromacs easily handles 10 ns, 1 µs, 1 ms sims. No technical limit.
- ❌ Equivalent: They're not equivalent—500 ps has ~2× statistical noise vs. 1 ns. But both are educationally valid.

**Insight:** In research, you run **as long as you can afford (CPU time)**. In this course, you learn the **workflow**, not achieve publication-grade convergence.

********************

## Part 🔟: Summary & Next Steps

### What You Learned

- ✅ Gromacs file formats and workflow  
- ✅ MDP parameter files (minimization, NVT, NPT, production)  
- ✅ The Lemkul step-by-step protocol  
- ✅ How to prepare inputs for cluster submission  
- ✅ Your role: prepare inputs, instructor runs on cluster, you analyze results  

### What You'll Do Next (Session 3)

1. **Get a protein** (from AlphaFold or PDB)
2. **Prepare Gromacs inputs** (MDP files + grompp commands)
3. **Send to instructor** (session3_inputs/ folder)
4. **Instructor runs MD** on cluster (~2-24 hours depending on protein size)
5. **Analyze results** (energy, RMSD, trajectories)
6. **Compare:** AlphaFold prediction vs. MD simulation

---

### ✅ Final Quiz: Session 2B

**Question 1:** What does .tpr file contain?

- [[ ]] Just the protein coordinates
- [[X]] Binary combination of topology, coordinates, and all MD settings
- [[ ]] Energy data from simulation
- [[ ]] Trajectory information
********************

**Explanation:**

- ✅ **Correct:** .tpr (Topology, Positions, Run parameters) is the preprocessed binary package:

  1. **Topology** from .top file (bond parameters, force field)
  2. **Positions** from .gro file (atomic coordinates)
  3. **Run parameters** from .mdp file (dt, nsteps, thermostat, etc.)
  - Created by `gmx grompp`, consumed by `gmx mdrun`
  - Single self-contained file; everything mdrun needs is inside
- ❌ Just coordinates: That's .gro or .pdb, not .tpr.
- ❌ Energy data: That's .edr (created *after* mdrun), not input.
- ❌ Trajectory: That's .trr or .xtc (output), not .tpr.

**Key:** Think of .tpr as the "recipe ready for execution."

********************

**Question 2:** In the Lemkul workflow, what comes immediately AFTER energy minimization?

- [[ ]] NPT equilibration (constant pressure phase)
- [[X]] NVT equilibration (constant temperature heating)
- [[ ]] Adding the protein to water
- [[ ]] Production MD
********************


**Explanation:**
- ✅ **Correct:** The Lemkul sequence is:

  1. **Minimize** (em.gro output)
  2. **NVT** (heat minimized structure to 300K, fixed volume) ← comes right after minimization
  3. **NPT** (adjust box to 1 bar pressure)
  4. **Production** (long unbiased simulation)
  - Thermodynamic reason: Must reach target T *before* allowing volume changes
- ❌ NPT: NPT comes *after* NVT (step 4, not step 3).
- ❌ Add water: Water added at the beginning (solvate step), long before minimization.
- ❌ Production: Production comes last (step 5), not immediately after minimize.

**Sequence matters:** Minimize → NVT → NPT → Production. This order is physics-driven, not flexible.

********************

**Question 3:** V-rescale thermostat is equivalent to:

- [[X]] CSVR (from Session 2A)
- [[ ]] Berendsen thermostat
- [[ ]] NVE (no thermostat)
- [[ ]] Langevin dynamics
********************

**Explanation:**

- ✅ **Correct:** V-rescale = **Velocity rescaling** = same as CSVR (Canonical Sampling via Velocity Rescaling) from Session 2A:

  - Both rescale atomic velocities to match target temperature
  - Both are fast, stable, good for equilibration
  - Used in every modern MD protocol
  - Gromacs calls it "V-rescale"; literature calls it "CSVR"—same thing
- ❌ Berendsen: Older thermostat (still used but less accurate). Different algorithm.
- ❌ NVE: That's *no* thermostat (microcanonical, energy conserving). Opposite of what we want.
- ❌ Langevin: Different thermostat (stochastic, adds friction). Different algorithm and dynamics.

**Terminology note:** Name variation between codes is common. V-rescale = CSVR—identical.

********************

**Question 4:** NPT ensemble means constant:

- [[X]] Number of particles, Pressure, and Temperature
- [[ ]] Number, Pressure, and Time
- [[ ]] Number, Volume, and Temperature
- [[ ]] Nothing is constant
********************

**Explanation:**

- ✅ **Correct:** **NPT** = Canonical Isobaric-Isothermal Ensemble:

  - **N** = Number of atoms (fixed; no atoms added/removed)
  - **P** = Pressure (constant at 1 bar; barostat controls)
  - **T** = Temperature (constant at 300K; thermostat controls)
  - **Volume is NOT constant** (allowed to fluctuate to maintain pressure)
  - Contrast: **NVT** = Number, **Volume** (not Pressure), Temperature
- ❌ Number, Pressure, Time: "Time" isn't a thermodynamic variable in ensembles.
- ❌ Number, Volume, Temp: That's **NVT**, not NPT.
- ❌ Nothing constant: Wrong—N, P, T are all controlled.

**Memory aid:** **N** = atomic number (fixed), **P** = pressure, **T** = temperature. Volume adjusts to maintain P.

********************

**Question 5:** Your simulation crashes after 500 MD steps with "LINCS warning: distance constraint." First thing to try?

- [[X]] Reduce dt (timestep), e.g., 0.002 → 0.001 ps
- [[ ]] Increase nsteps
- [[ ]] Change force field (e.g., amber99sb → charmm36)
- [[ ]] Skip NVT equilibration
********************

**Explanation:**

- ✅ **Correct:** LINCS warning = **bond constraint violation** caused by timestep too large:

  - dt=0.002 ps might be too big for your system (forces very large, atoms moving too far)
  - Reduce dt → smaller steps → smoother integration → fixes constraint violations
  - Usually crashes after 500 steps = system instability building up from day 1
  - Solution: Restart with dt=0.001 ps, try again
- ❌ Increase nsteps: Won't help—problem is per-step integration, not number of steps.
- ❌ Change FF: Force fields don't cause LINCS violations (unless structure is wrong); this wastes time.
- ❌ Skip NVT: Wrong—NVT is necessary for proper equilibration. Skipping makes it worse.

**Real troubleshooting:** If it crashes, always check: (1) timestep, (2) starting geometry (was it minimized?), (3) simulation conditions.

********************

**Question 6:** Calculate: If `dt=0.002` ps, `nsteps=50000`, and `nstenergy=500`, how many energy data points do you get?

- [[ ]] 50 data points
- [[X]] 100 data points (simulation is 100 ps, energy saved every 1 ps)
- [[ ]] 100 data points
- [[ ]] 500 data points
********************

**Explanation:**

- ✅ **Correct:** Multiple calculations:

  - Duration: nsteps × dt = 50,000 × 0.002 ps = **100 ps**
  - Energy frequency: nstenergy × dt = 500 × 0.002 ps = **1 ps per save**
  - Data points: Duration / frequency = 100 ps / 1 ps = **100 points**
  - This gives ~100 independent energy values for analysis
- ❌ 50 points: Would need nstenergy=1000 (1000 × 0.002 = 2 ps per save).
- ❌ 500 points: That's nstenergy=100 (100 × 0.002 = 0.2 ps per save)—overkill, huge file.

**Practical:** With ~100 data points, you can see equilibration trends, temperature/pressure stability. For statistical averaging, 100 pts is reasonable (not amazing, but educational).

**Calculation skill:** Know **Duration = nsteps × dt** and **Frequency = nstenergy × dt**. These appear in every MDP file.

********************

## 📚 Resources

- **Gromacs Documentation:** https://manual.gromacs.org/
- **Lemkul GROMACS Tutorials:** https://www.bevanlab.org/static/papers/
- **GROMACS Manual (Nonbonded):** https://manual.gromacs.org/current/user-guide/mdrun-features/nonbonded-interactions.html
- **Force Field References:** CHARMM36, Amber99SB papers

---

## Questions for Live Session

1. Why do we need both NVT and NPT equilibration?
2. What would happen if you skipped energy minimization?
3. How do you choose the coupling times for thermostat/barostat?
4. What's the difference between .xtc and .trr trajectories?
5. If your simulation crashes, what are common reasons?

---

*Session 2B — Introduction to Gromacs (Lemkul Workflow)*  
*Last updated: October 27, 2025*  
*Course: Molecular Modelling and Quantum Chemistry (Master)*
