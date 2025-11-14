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

**Gromacs** = **GRO**tesque **M**olecular **A**nd **C**hemical **S**imulation

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

Which task is Gromacs best suited for?
- [ ] Teaching MD to beginners
- [x] Running 100 ns simulations of proteins on a cluster
- [ ] Quick energy minimization tests
- [ ] Simple molecules in vacuum

---

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

Which file contains the force field parameters?
- [ ] .pdb
- [ ] .gro
- [x] .top / .itp
- [ ] .xtc

What does .tpr stand for?
- [x] Topology, Positions, and Run input
- [ ] Trajectory print result
- [ ] Total potential representation
- [ ] Time-point recording

---

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

The correct order of Lemkul workflow is:
- [x] Minimize → NVT equil → NPT equil → Production
- [ ] NVT → Minimize → NPT → Production
- [ ] Production → Minimize → Equilibration
- [ ] Minimize → Production directly

---

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

What does `dt = 0.002` mean?
- [ ] 0.002 seconds per step
- [x] 0.002 picoseconds per step (2 femtoseconds)
- [ ] 0.002 nanometers
- [ ] Not important

In NVT ensemble, which should be "no"?
- [x] pcoupl (pressure coupling)
- [ ] tcoupl (temperature coupling)
- [ ] nstenergy (output frequency)
- [ ] integrator

---

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

The command `gmx grompp` does what?
- [ ] Runs MD simulation
- [x] Preprocesses inputs (.mdp, .top, .gro) and creates binary .tpr file
- [ ] Analyzes trajectory results
- [ ] Adds water to the system

---

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

In this course, what do you do with Gromacs?
- [x] Prepare inputs and create .tpr files
- [ ] Run full production simulations
- [ ] Only analyze results
- [ ] Nothing (only Curcuma)

---

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

Which dt is most stable for protein MD?
- [ ] 0.1 ps
- [x] 0.002 ps (2 fs)
- [ ] 0.5 ps
- [ ] 0.0001 ps (overkill)

---

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

This course is based on the Lemkul workflow but:
- [x] Uses shorter equilibration times for teaching
- [ ] Skips energy minimization entirely
- [ ] Uses Berendsen thermostat only
- [ ] Doesn't follow any published protocol

---

## Part 🔟: Summary & Next Steps

### What You Learned

✅ Gromacs file formats and workflow  
✅ MDP parameter files (minimization, NVT, NPT, production)  
✅ The Lemkul step-by-step protocol  
✅ How to prepare inputs for cluster submission  
✅ Your role: prepare inputs, instructor runs on cluster, you analyze results  

### What You'll Do Next (Session 3)

1. **Get a protein** (from AlphaFold or PDB)
2. **Prepare Gromacs inputs** (MDP files + grompp commands)
3. **Send to instructor** (session3_inputs/ folder)
4. **Instructor runs MD** on cluster (~2-24 hours depending on protein size)
5. **Analyze results** (energy, RMSD, trajectories)
6. **Compare:** AlphaFold prediction vs. MD simulation

---

### ✅ Final Quiz: Session 2B

What does .tpr file contain?
- [x] Binary combination of topology, coordinates, and all MD settings
- [ ] Just the protein coordinates
- [ ] Energy data from simulation
- [ ] Trajectory information

In the Lemkul workflow, minimization runs BEFORE:
- [x] NVT and NPT equilibration
- [ ] Adding the protein to water
- [ ] Creating the topology
- [ ] All the above

V-rescale thermostat is equivalent to:
- [x] CSVR (from Session 2A)
- [ ] Berendsen thermostat
- [ ] NVE (no thermostat)
- [ ] Langevin dynamics

NPT ensemble means constant:
- [x] Number of particles, Pressure, and Temperature
- [ ] Number, Pressure, and Time
- [ ] Number, Volume, and Temperature
- [ ] Nothing is constant

---

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
