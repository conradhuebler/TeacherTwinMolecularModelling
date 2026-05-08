<!--
author: Molecular Modelling Course Team
language: en
narrator: US English Female
version: 1.0

Session 2A: Molecular Dynamics with Curcuma
Part of: Molecular Modelling and Quantum Chemistry (Master)
-->

[Open this course in LiaScript](https://liascript.github.io/course/?https://raw.githubusercontent.com/conradhuebler/TeacherTwinMolecularModelling/main/session2a_md_curcuma.md)

# Session 2A: Molecular Dynamics — Understanding Thermostats with Curcuma

## Part A: From Optimization to Dynamics

> **Welcome to Session 2A!**
>
> You've optimized molecular structures. Now they're going to **move**.
>
> In this session, you'll run **Molecular Dynamics (MD) simulations** using Curcuma,
> compare different **thermostats**, and analyze how temperature and energy 
> evolve during a simulation.

---

## 🎯 Learning Objectives

By the end of Part A, you will:

- ✅ Understand **molecular dynamics** and **ensembles** (NVE, NVT, etc.)
- ✅ Know what **thermostats** do and why they matter
- ✅ Compare **Berendsen** vs. **CSVR** thermostats practically
- ✅ Extract and analyze MD data (temperature, energy, RMSD)
- ✅ Plot MD results with **Gnuplot**
- ✅ Interpret trajectories and energetics

---

## Part 1️⃣: Concept — Molecular Dynamics

### What is Molecular Dynamics?

**Molecular Dynamics (MD)** simulates how molecules move and evolve over time by integrating Newton's equations of motion:

$$\vec{F}_i = m_i \vec{a}_i = -\nabla E(\vec{r})$$

**In plain English:**

1. Calculate forces on each atom (from force field gradient)
2. Update velocities and positions (Newton's laws)
3. Repeat many times (thousands to millions of steps)

**Result:** A **trajectory** showing how atomic positions change over nanoseconds or microseconds.

**Why it matters:**

- Proteins don't sit still—they move, flex, breathe
- MD captures dynamics that optimization cannot
- We can observe transitions between conformations
- Extract thermodynamic properties (temperature, energy, density)

---

### Ensembles: NVE, NVT, NPT

Different simulation conditions define different **ensembles**:

| Ensemble | Constant | Variable | Thermostat | Real system |
|----------|----------|----------|-----------|------------|
| **NVE** | N, V, E | T, P | ✗ | Isolated system (vacuum) |
| **NVT** | N, V, T | E, P | ✓ | System in thermal bath |
| **NPT** | N, P, T | V, E | ✓ | System in box, room conditions |

**Most common:** **NVT** (constant number of particles, volume, temperature) — like a test tube in a water bath.

In **NVE**, energy is conserved but temperature fluctuates (unrealistic for most biochemistry).

### ✅ Quick Check 1: Ensembles

**Question 1:** Which ensemble is most realistic for a protein in aqueous solution?

- [[ ]] NVE (isolated system)
- [[X]] NVT (constant temperature, like in a water bath)
- [[ ]] NPT (constant pressure)
- [[ ]] None of the above
********************

**Explanation:**

- ✅ **Correct:** NVT (constant N, V, T) mimics a real protein in a water bath at lab conditions (298 K, 1 atm ambient). The aqueous environment acts as a heat reservoir, keeping temperature constant. The protein is confined in its hydration shell (constant volume at MD timescales).
- ❌ **NVE:** Energy-conserving but unrealistic—proteins in solution exchange heat with surroundings. NVE is only appropriate for vacuum simulations or as a theoretical reference.
- ❌ **NPT:** This is used for production MD of pure liquids/bulk simulations, not proteins in solution (which have a constant volume container).
- ❌ **None of the above:** All options listed are valid ensembles; NVT is simply the best match.

**Key insight:** NVT = "protein in a water bath at constant temperature." This is what happens in your lab!

********************

**Question 2:** In NVE ensemble, what is conserved?

- [[ ]] Temperature
- [[X]] Total energy
- [[ ]] Volume
- [[ ]] Pressure
********************

**Explanation:**

- ✅ **Correct:** NVE conserves total energy (E = KE + PE = constant). No heat exchange with environment, so energy cannot be gained or lost. This is why it's called the **microcanonical ensemble** in statistical mechanics.
- ❌ **Temperature:** Not conserved in NVE. As atoms move, kinetic energy oscillates while potential energy changes, so T fluctuates significantly.
- ❌ **Volume:** NVE keeps volume fixed, but this is assumed, not conserved as an energy consequence.
- ❌ **Pressure:** Not conserved—follows from temperature and density fluctuations.

**Key insight:** "NVE" literally means N (particle number), V (volume), and E (energy) are constant. Temperature and pressure are *derived properties* that fluctuate.

********************

---

### Thermostats: Maintaining Temperature

A **thermostat** keeps temperature constant by rescaling or reassigning velocities. It couples your system to a "heat bath" at temperature T.

Without a thermostat (NVE), temperature fluctuates wildly and is unrealistic for biochemistry.

**Why we need one:** A real molecule in solution is at constant T. MD simulations should reflect this.

---

## Part 2️⃣: Thermostat Comparison — Berendsen vs. CSVR

### Berendsen Thermostat

**What it does:** Weakly couples the system to a heat bath. If current temperature differs from target, velocities are **rescaled**.

**Equation:**
$$T_{\text{new}} = T_{\text{current}} \left(1 + \frac{\Delta t}{\tau_T}\left(\frac{T_0}{T_{\text{current}}} - 1\right)\right)$$

**Characteristics:**

- ✅ Simple, fast, easy to use
- ✅ Converges to target temperature quickly
- ❌ Does NOT generate proper NVT ensemble
- ❌ Underestimates fluctuations
- ❌ "Artificial" temperature control

**When to use:** Initial equilibration, not for production runs.

**Coupling time:** How strongly coupled is the system?

- Small τ = strong coupling (fast temperature control, but artificial)
- Large τ = weak coupling (more natural, but slower equilibration)

### ✅ Quick Check 2: Berendsen Thermostat

**Question:** Berendsen thermostat rescales velocities. What does this accomplish?

- [[X]] Adjusts kinetic energy to match target temperature
- [[ ]] Changes potential energy
- [[ ]] Moves atoms to new positions
- [[ ]] Adds random forces
********************

**Explanation:**

- ✅ **Correct:** Velocity rescaling directly changes kinetic energy (KE = ½m·v²). By scaling all velocities by a factor λ, kinetic energy (and thus temperature, since T ∝ KE) is adjusted toward the target. If T_current > T_target, velocities are reduced (λ < 1). If T_current < T_target, velocities are increased (λ > 1).
- ❌ **Changes potential energy:** No—Berendsen doesn't modify atomic positions, only velocities. PE depends only on atomic coordinates, not velocities.
- ❌ **Moves atoms to new positions:** Berendsen is purely a velocity thermostat. Positions are updated via normal MD integration in the next step.
- ❌ **Adds random forces:** That's CSVR or Langevin thermostats. Berendsen is **deterministic** (no randomness).

**Practical tip:** Velocity rescaling is like "speeding up" or "slowing down" the atomic motion globally to change the temperature. Simple but artificial!

**Hint for study:** Remember KE = ½m·v². If you want to change T without moving atoms, you must change their velocities.

********************

---

### CSVR Thermostat (Canonical Sampling via Velocity Rescaling)

**What it does:** More sophisticated approach. Applies velocity rescaling **stochastically** (with randomness) to generate a proper NVT ensemble.

**Reference:** Bussi, Donadio, Parrinello, *J. Chem. Phys.* **2007**, 126, 014101

**Characteristics:**

- ✅ Generates **correct NVT ensemble** (proper statistics)
- ✅ Velocities include random fluctuations (more realistic)
- ✅ Better for production runs
- ❌ Slightly more complex
- ❌ Older computers: longer to run (negligible nowadays)

**When to use:** Production MD where you care about correct thermodynamic properties.

**Coupling time:** Similar meaning as Berendsen, but applied stochastically.

---

### Quick Comparison

| Property | Berendsen | CSVR |
|----------|-----------|------|
| **Theory** | Weak coupling | Stochastic rescaling |
| **Ensemble** | ✗ (approximate NVT) | ✓ (true NVT) |
| **Speed** | Slightly faster | Slightly slower |
| **Temperature fluctuations** | Underestimated | Realistic |
| **Use case** | Equilibration | Production |

### ✅ Quick Check 3: When to Use Which

**Question:** You're running a 10 ns production MD. Best thermostat?

- [[ ]] Berendsen (faster equilibration)
- [[X]] CSVR (proper NVT ensemble for correct statistics)
- [[ ]] Both simultaneously
- [[ ]] NVE (no thermostat needed)
********************

**Explanation:**

- ✅ **Correct:** CSVR is **the industry standard** for production MD because it generates a proper canonical (NVT) ensemble with correct Boltzmann statistics. This is essential if you want reliable thermodynamic properties (free energy, entropy, equilibrium populations). A 10 ns production run deserves rigorous statistics!
- ❌ **Berendsen:** Fast but only for equilibration (initial ~1 ns). It underestimates fluctuations and distorts statistical distributions. Using it for production would bias your results.
- ❌ **Both simultaneously:** Mixing thermostats is unphysical and generates artifacts.
- ❌ **NVE:** Microcanonical ensemble doesn't match lab conditions. Better only for benchmark tests or if your system is already equilibrated and you want to monitor energy conservation (but even then, use it sparingly).

**Practical tip:** **Equilibration = Berendsen. Production = CSVR.** This is the workflow used in every modern MD paper.

**Hint:** "Production" data is what you publish. You need the best thermostat (CSVR) to make your data reliable!

********************

---

## Part 3️⃣: Practical MD with Curcuma

### Basic Syntax

```bash
curcuma -md input.xyz -method methodname -thermostat thermostat_name [options]
```

**Example:**

```bash
# MD with Berendsen thermostat for 100 fs
curcuma -md structure.xyz -method gfnff -thermostat berendsen -maxtime 100

# MD with CSVR thermostat for 100 ps
curcuma -md structure.xyz -method gfnff -thermostat csvr -maxtime 1e5
```

---

### Output Format

Curcuma prints results to **stdout** in **tabular format**:

```
Time(ps)  e_tot  <e_tot>  e_kin  <e_kin>  e_pot  <e_pot>  T(K)  <T(K)>  ...
0.000000  1.254  1.254    0.166  0.166    1.420  1.420    298.7  298.7
0.010000  1.248  1.255    0.265  0.267    1.514  1.420    477.8  480.0
0.020000  1.235  1.249    0.249  0.269    1.484  1.465    448.0  484.6
...
```

**Columns:**

- `Time` = Simulation time in picoseconds (ps)
- `e_tot` = Total energy (instant)
- `<e_tot>` = Average total energy
- `e_kin` = Kinetic energy (instant)
- `<e_kin>` = Average kinetic energy
- `e_pot` = Potential energy (instant)
- `<e_pot>` = Average potential energy
- `T` = Temperature (instant, in Kelvin)
- `<T>` = Average temperature

**Files created:**

- `structure.trj.xyz` — Full trajectory (all atoms, all steps)
- Energy/temp data on stdout (can be piped to file)

---

### ✅ Quick Check 4: MD Output

**Question:** What does the `<e_tot>` column represent?

- [[ ]] Total energy at a single snapshot
- [[X]] Running average of total energy (time-averaged)
- [[ ]] The final energy after MD
- [[ ]] The initial energy before MD
********************

**Explanation:**

- ✅ **Correct:** `<e_tot>` is a **time-averaged (moving average)** of total energy up to that point. This smooths out instantaneous fluctuations and reveals the true energy trend. If the simulation is stable, `<e_tot>` should plateau after equilibration (no drift).
- ❌ **Single snapshot:** That's `e_tot` (without angle brackets). The angle brackets `< >` always denote averaging in physics/chemistry.
- ❌ **Final energy after MD:** `<e_tot>` evolves throughout the simulation. The *last* value is approximately the final energy.
- ❌ **Initial energy before MD:** No—`<e_tot>` starts small and builds up as the average includes more and more frames.

**Practical tip:** When analyzing output, look at both `e_tot` (noisy) and `<e_tot>` (smooth). The trend in `<e_tot>` tells you if your system is stable!

**Hint:** If you see a zigzag in `e_tot` but a flat line in `<e_tot>`, your simulation is working correctly—high-frequency oscillations are expected, but the average should be stable.

**Hint:** In Gnuplot, use column 3 (`<e_tot>`) for cleaner plots, not column 2 (`e_tot`).

********************

---

## Part 4️⃣: Exercise 1 — Compare Thermostats on Glucose

### Objective

Run the **optimized glucose** from Session 1 through two MD simulations:

1. **Berendsen thermostat** 
2. **CSVR thermostat** 

Compare temperature stability and energy behavior.

---

### Setup

```bash
mkdir -p session2a_md
cd session2a_md

# Use optimized glucose from Session 1
cp ../session1_glucose/glucose.opt.xyz glucose_start.xyz
```

---

### Step 1: MD with Berendsen (100 ps)

```bash
# Run MD with Berendsen
curcuma -md glucose_start.xyz -method gfnff -thermostat berendsen \
    -maxtime 100 -print 10 -dump 1 > md_berendsen.dat

# This creates:
# - glucose_start.trj.xyz (trajectory)
# - md_berendsen.dat (energy/temp data)
cp glucose_start.trj.xyz glucose_berendson.trj.xyz # store trajectory for later comparison
```

**Important flags:**

- `-maxtime 1e5` = Run for 100 picoseconds
- `-print 10` = Print output every 10 fs (0.01 ps)
- `-dump 1` = Write trajectory every 1 fs

**Expected output file:**

```
Time(ps)  e_tot  <e_tot>  e_kin  <e_kin>  e_pot  <e_pot>  T(K)  <T(K)>
0.000000  ...
0.010000  ...
...
```

**Save this file!** 💾

---

### Step 2: MD with CSVR (100 ps)

**Start fresh from the optimized geometry** (don't use the Berendsen trajectory):

```bash
# Run MD with CSVR
curcuma -md glucose_start.xyz -method gfnff -thermostat csvr \
    -maxtime 1e5 -print 10 -dump 1 > md_csvr.dat

# Creates:
# - glucose_start.trj.xyz (overwrites previous!)
# - md_csvr.dat (energy/temp data)
cp glucose_start.trj.xyz glucose_csvr.trj.xyz
```

**Rename the trajectory so we don't lose it:**

```bash
mv glucose_start.trj.xyz glucose_csvr.trj.xyz
```

---

### Step 3: Extract Data with AWK

Now we have two data files: `md_berendsen.dat` and `md_csvr.dat`

**Extract specific columns using `awk`:**

```bash
# Extract: Time (col 1) and Temperature (col 8)
awk '{print $1, $8}' md_berendsen.dat > temp_berendsen.txt

# Extract: Time (col 1) and Total Energy (col 2)
awk '{print $1, $2}' md_berendsen.dat > energy_berendsen.txt

# Same for CSVR
awk '{print $1, $8}' md_csvr.dat > temp_csvr.txt
awk '{print $1, $2}' md_csvr.dat > energy_csvr.txt
```

**Verify the files:**

```bash
head temp_berendsen.txt
# Should show: Time  Temperature
```

---

### Step 4: Plot with Gnuplot

**Basic Gnuplot script** (`plot_md.gnu`):

```gnuplot
set terminal png size 1200,800
set output 'md_comparison.png'
set multiplot layout 2,2

# Plot 1: Temperature (Berendsen)
set title 'Temperature: Berendsen Thermostat'
set xlabel 'Time (ps)'
set ylabel 'Temperature (K)'
plot 'temp_berendsen.txt' u 1:2 w l title 'Berendsen'

# Plot 2: Temperature (CSVR)
set title 'Temperature: CSVR Thermostat'
set xlabel 'Time (ps)'
set ylabel 'Temperature (K)'
plot 'temp_csvr.txt' u 1:2 w l title 'CSVR'

# Plot 3: Energy (Berendsen)
set title 'Total Energy: Berendsen'
set xlabel 'Time (ps)'
set ylabel 'Energy (kcal/mol)'
plot 'energy_berendsen.txt' u 1:2 w l title 'Berendsen'

# Plot 4: Energy (CSVR)
set title 'Total Energy: CSVR'
set xlabel 'Time (ps)'
set ylabel 'Energy (kcal/mol)'
plot 'energy_csvr.txt' u 1:2 w l title 'CSVR'

unset multiplot
```

**Run Gnuplot:**

```bash
gnuplot plot_md.gnu
# Creates: md_comparison.png
```

**View the result:**

```bash
display md_comparison.png
# or: feh md_comparison.png
# or open in image viewer
```

---

### Step 5: Analyze Results

**Look at your plots. Answer these questions:**

1. **Temperature stability:**

   - Which thermostat keeps T closer to 298 K?
   - Which has bigger fluctuations?
   - Expected: CSVR should be more stable (that's the point)

2. **Energy behavior:**

   - Is total energy constant or drifting?
   - In NVT, total energy should fluctuate around a mean (but not drift)
   - Which thermostat looks "cleaner"?

3. **Differences between thermostats:**

   - Can you see artifacts in Berendsen (e.g., jerky T changes)?
   - Does CSVR look more "natural"?
   
4. **Check the trajectories:**

   - Use Avogadro to visualise both trajectories. Do you observe differences?
   - Does CSVR or Berendson look more "resonable"?
   
**⚠️ Important Observation: Flying Ice Cube Artefact**

---

### ✅ Quick Check 5: Interpreting MD Results

**Question 1:** If total energy drifts downward throughout the MD, this suggests:

- [[ ]] The thermostat is working perfectly
- [[X]] There might be numerical instabilities or poor parameters
- [[ ]] It's normal and expected
- [[ ]] The force field is wrong
********************

**Explanation:**

- ✅ **Correct:** Continuous energy drift (monotonic decrease or increase) indicates **problems**: (a) time step too large (numerical integration errors accumulate), (b) bad starting geometry (atoms close to each other, repulsive forces dominate), (c) cutoff distance too short (missing long-range interactions), or (d) integration algorithm mismatch. In well-behaved simulations, `<e_tot>` plateaus after equilibration.
- ❌ **Thermostat working perfectly:** A good thermostat maintains T but doesn't prevent energy drift caused by integration errors.
- ❌ **Normal and expected:** **Drift is a red flag.** A healthy simulation shows `<e_tot>` stability (within ±5% after equilibration).
- ❌ **Force field is wrong:** Force field issues cause *high* total energy or bad dynamics, but not necessarily monotonic drift. Drift is usually numerical.

**Practical tip:** Always plot `<e_tot>` during equilibration. If it's still drifting after 1 ns, reduce your time step or check your starting structure!

**Hint:** Numerical instability shows up as drift. A bad force field shows up as wrong *values*, not trends.
********************

**Question 2:** Temperature fluctuations in CSVR should be:

- [[ ]] Zero (perfectly constant)
- [[X]] Non-zero but realistic (Boltzmann distribution)
- [[ ]] Larger than Berendsen
- [[ ]] Unpredictable
********************

**Explanation:**

- ✅ **Correct:** CSVR generates stochastic velocity rescaling, so T oscillates around the target value. These fluctuations follow the **Boltzmann distribution** (canonical ensemble). In a 300-atom system at 298 K, you might see ±20 K swings—this is physically realistic, not a bug!
- ❌ **Zero (perfectly constant):** Impossible in any real system. Even NVE shows T fluctuations (from KE ↔ PE exchange). A perfectly flat T indicates artificial thermostat damping (like Berendsen) or too small a system.
- ❌ **Larger than Berendsen:** Opposite! Berendsen artificially suppresses fluctuations. CSVR allows larger, more realistic T swings.
- ❌ **Unpredictable:** CSVR is stochastic but **reproducible** if you set a random seed. It's not random in an uncontrolled way.

**Practical tip:** A healthy CSVR simulation shows T oscillating ±10-20% around the target. If T is stone-flat, check if Berendsen is still active or if your system is too small.

**Hint:** In real life, your lab temperature isn't exactly 298.000000 K—it fluctuates slightly. CSVR captures that!

**Hint:** Plot both `T` (noisy) and `<T>` (smooth). The average `<T>` should match your target; scatter in `T` is fine and expected.

********************

---

## Part 5️⃣: Exercise 2 — RMSD Analysis During MD

### Objective

Track how much the **structure changes** from its starting point during MD.

**RMSD (Root-Mean-Square Displacement)** measures average atomic displacement:
$$\text{RMSD} = \sqrt{\frac{1}{N} \sum_{i=1}^{N} (\vec{r}_i - \vec{r}_i^{\text{ref}})^2}$$

Small RMSD = structure stays stable
Large RMSD = structure changes a lot (unfolds, denatures, etc.)

---

### Step 1: Extract Snapshots from Trajectory

The `.trj.xyz` file contains many snapshots. We need to extract individual frames.

**Using `grep` and `sed` (or Curcuma's trajectory tools):**

```bash
# Count how many frames in trajectory
grep -c "^[0-9]" glucose_csvr.trj.xyz
# Let's say it's 10000 frames (100 ps at 1 fs/step)

# Extract frame 0 (beginning)
# Frames are separated by blank lines; each frame has atom count on first line
# This is a bit complex; Curcuma might have a tool for this
```

**Simpler approach: Use Avogadro to export key frames**

```bash
avogadro glucose_csvr.trj.xyz &
# In Avogadro: Export first frame, last frame as separate XYZ files
# glucose_frame0.xyz (start)
# glucose_frame_last.xyz (end)
```

---

### Step 2: Calculate RMSD vs. Starting Structure

```bash
curcuma -rmsd glucose_start.xyz glucose_frame_last.xyz
# Output: RMSD = 0.656 Angstrom (Berendsen) or 0.911 Angstrom (CSVR)
```

**Interpretation:**

- RMSD < 0.5 Å: Structure stable (good!)
- RMSD 0.5 - 2 Å: Significant change (flexible)
- RMSD > 2 Å: Major refolding or unfolding


---

### ✅ Quick Check 6: RMSD Interpretation

**Question:** High RMSD at the end of MD means:

- [[ ]] The simulation failed
- [[X]] The structure moved/changed significantly
- [[ ]] Better thermostat control
- [[ ]] Lower temperature
********************

**Explanation:**

- ✅ **Correct:** RMSD measures atomic displacement from the starting structure. High RMSD (>2 Å) means atoms have moved far from their initial positions. This could indicate: (a) the structure is flexible/floppy, (b) it's exploring multiple conformations, or (c) it's unfolding/denaturing. None of these are failures—they tell you how the molecule behaves dynamically!
- ❌ **Simulation failed:** High RMSD doesn't mean failure. A protein that flexes has naturally high RMSD. Only if RMSD explodes to 10+ Å *and* energy explodes should you suspect a crash.
- ❌ **Better thermostat control:** Thermostat doesn't directly control RMSD. RMSD is a structural property that depends on the force field and temperature.
- ❌ **Lower temperature:** Temperature and RMSD are related but opposite: lower T → less atomic motion → lower RMSD. Higher T → more motion → higher RMSD.

**Practical tip:** Expected RMSD ranges:
- **< 0.5 Å:** Very stable structure (like a rigid molecule or protein in a crystal).
- **0.5–2 Å:** Flexible but stable (typical for proteins in solution).
- **> 2 Å:** Major structural rearrangement (unfolding, domain motion, or system instability).

**Hint:** RMSD is a *symptom*, not a diagnosis. Plot it to understand *how* the structure changes, then investigate *why*.

**Hint for comparison:** If Berendsen gives RMSD = 1.0 Å and CSVR gives RMSD = 1.5 Å, that's normal—different thermostats sample phase space differently, and CSVR explores more conformational space!

**Reference:** See the RMSD calculation formula in Part 5️⃣. It's the root-mean-square of individual atomic displacements from the starting positions.

********************

---

## Part 6️⃣: Exercise 3 — Compare Both Peptides

### Objective

Run MD on both peptides from Session 1:

- **AAAAA** (homogeneous)
- **WRKLQ** (heterogeneous)

Use **only CSVR thermostat** (the good one).

**Predict:** Which peptide will be more stable (lower RMSD)?

---

### Run MD on AAAAA

```bash
curcuma -md aaaaa.opt.xyz -method gfnff -thermostat csvr -t 300 -maxtime 10000

# Extract frames and calculate RMSD
head -55 aaaaa.opt.trj.xyz > aaaaa_start.xyz
tail -55 aaaaa.opt.trj.xyz > aaaaa_end.xyz
curcuma -rmsd aaaaa_start.xyz aaaaa_end.xyz
```

💾 **Results:**

- Final energy: **-8.884406 Eh**
- RMSD: **2.77 Å**

---

### Run MD on WRKLQ (pentapeptid)

```bash
curcuma -md pentapeptid.opt.xyz -method gfnff -thermostat csvr -t 300 -maxtime 10000

# Extract frames and calculate RMSD
head -109 pentapeptid.opt.trj.xyz > wrklq_start.xyz
tail -109 pentapeptid.opt.trj.xyz > wrklq_end.xyz
curcuma -rmsd wrklq_start.xyz wrklq_end.xyz
```

💾 **Results:**

- Final energy: **-18.628899 Eh**
- RMSD: **3.78 Å**

---

### Compare

| Property | AAAAA | WRKLQ |
|----------|-------|-------|
| **Final Energy** | -8.884406 Eh | -18.628899 Eh |
| **Energy per Atom** | -0.1677 Eh/atom | -0.1742 Eh/atom |
| **RMSD (start to end)** | 2.77 Å | 3.78 Å |
| **Flexibility** | Low (homogeneous) | High (heterogeneous) |

---

### Analysis Questions

1. Which peptide is more stable structurally (lower RMSD)?
2. Is AAAAA more stable because it's simpler?
3. Does heterogeneity (WRKLQ) lead to more flexibility?
4. Any relationship between sequence and dynamics?

### ✅ Quick Check 7: Peptide Dynamics

**Question:** Why might AAAAA be more stable than WRKLQ?

- [[X]] Homogeneous sequences have fewer degrees of freedom
- [[ ]] CSVR doesn't work on diverse sequences
- [[ ]] Small peptides are always stable
- [[ ]] GFN-FF favors uniform structures
********************

**Explanation:**

- ✅ **Correct:** Alanine (A) is a small, non-polar amino acid with minimal side-chain interactions. A homopolypeptide like AAAAA can only form simple secondary structures (α-helix or extended). WRKLQ, by contrast, has bulky (W, R, K) and charged (R, K, Q) residues that create multiple electrostatic interactions, hydrogen bonding networks, and conformational preferences. More interaction diversity = more degrees of freedom = higher RMSD. Think of it as "a simple structure is easier to lock in place than a complex one."
- ❌ **CSVR doesn't work on diverse sequences:** CSVR is a thermostat, not a method that "favors" certain sequences. It works equally well on all peptides.
- ❌ **Small peptides are always stable:** Size alone doesn't determine stability. WRKLQ is also small (5 residues) but more flexible than AAAAA.
- ❌ **GFN-FF favors uniform structures:** The force field doesn't discriminate based on sequence uniformity. It just calculates forces from the atoms present.

**Practical tip:** **Sequence complexity = conformational flexibility.** Uniform sequences (poly-A, poly-G) are more rigid; diverse sequences explore more space.

**Biological insight:** In real proteins, diverse residues allow *folding* (structure formation through interactions). A protein made entirely of alanine would be boring and wouldn't fold properly!

**Hint for your analysis:** When comparing AAAAA and WRKLQ RMSD values, WRKLQ will likely show higher RMSD. That's not a failure—it's expected biochemistry!

**Hint:** Plot RMSD(t) for both peptides. AAAAA's curve should plateau faster (reaches equilibrium sooner); WRKLQ's curve may rise more gradually (still exploring conformations even at 100 ps).

**Connection:** This is why studying peptide dynamics is interesting—sequence determines flexibility!

********************

---

## Part 7️⃣: Gnuplot Deep Dive

### Why Gnuplot?

Gnuplot is a **command-line plotting tool**:

- ✅ Scriptable (reproducible)
- ✅ Powerful (publication-quality plots)
- ✅ Standard in computational chemistry
- ✅ Free and open-source

### Script Structure

Basic Gnuplot script (`.gnu` file):

```gnuplot
# 1. Set output format and file
set terminal png size 1000,600
set output 'myplot.png'

# 2. Set labels and title
set title 'My Simulation Results'
set xlabel 'Time (ps)'
set ylabel 'Temperature (K)'

# 3. Set ranges (optional)
set xrange [0:100]
set yrange [250:350]

# 4. Plot data
plot 'data.txt' u 1:2 w l title 'Simulation Data'

# 5. Save
# (automatic with set output)
```

**Alternative: Interactive plot (X11)**

```gnuplot
set terminal x11
set title 'Temperature'
plot 'data.txt' u 1:2 w l
# Plot appears in a window; close to exit
```

---

### Useful Gnuplot Options

| Option | Meaning |
|--------|---------|
| `u 1:2` | Use columns 1 (x) and 2 (y) |
| `u 1:3:4` | Columns with error bars |
| `w l` | Plot with lines |
| `w p` | Plot with points |
| `w lp` | Plot with lines and points |
| `title 'Label'` | Legend label |
| `set xrange [a:b]` | Limit x-axis |
| `set grid` | Add grid |
| `set logscale y` | Logarithmic y-axis |

---

### Advanced: Multiple Datasets

```gnuplot
set terminal png size 1200,800
set output 'comparison.png'

set title 'Temperature Comparison'
set xlabel 'Time (ps)'
set ylabel 'T (K)'

plot 'temp_berendsen.txt' u 1:2 w l title 'Berendsen', \
     'temp_csvr.txt' u 1:2 w l title 'CSVR'
```

---

### Create Custom Plots for YOUR Data

**Template for temperature comparison:**

```gnuplot
set terminal png size 1200,600
set output 'temp_all.png'

set title 'Temperature Evolution in MD'
set xlabel 'Time (ps)'
set ylabel 'Temperature (K)'
set grid

plot 'temp_berendsen.txt' u 1:2 w l lw 2 title 'Berendsen', \
     'temp_csvr.txt' u 1:2 w l lw 2 title 'CSVR'
```

**Template for energy comparison:**

```gnuplot
set terminal png size 1200,600
set output 'energy_all.png'

set title 'Total Energy During MD'
set xlabel 'Time (ps)'
set ylabel 'Energy (kcal/mol)'
set grid

plot 'energy_berendsen.txt' u 1:2 w l lw 2 title 'Berendsen', \
     'energy_csvr.txt' u 1:2 w l lw 2 title 'CSVR'
```

---

### ✅ Quick Check 8: Gnuplot Syntax

**Question:** What does `plot 'data.txt' u 1:2 w l` mean?

- [[ ]] Use only column 1 for plotting
- [[X]] Plot column 1 vs column 2 with lines
- [[ ]] Plot data with points
- [[ ]] Use logarithmic scale
********************

**Explanation:**

- ✅ **Correct:** This Gnuplot command:

  - `plot`: Start a plot
  - `'data.txt'`: Read from this file
  - `u 1:2`: **U**se columns 1 (x-axis) and 2 (y-axis)
  - `w l`: **W**ith **l**ines (not points)
  - Result: A line plot of column 2 vs column 1
- ❌ **Only column 1:** `u 1` alone would be meaningless. You always need two columns for 2D plotting (one for x, one for y).
- ❌ **Plot with points:** That would be `w p`. `w l` specifically means lines.
- ❌ **Logarithmic scale:** That requires `set logscale y` or similar. This command doesn't set any scales.

**Gnuplot vocabulary breakdown:**

- `u` = using (which columns to use)
- `w` = with (line style)
- `l` = lines
- `p` = points
- `lp` = lines and points

**Practical example:**

```gnuplot
plot 'temp_csvr.txt' u 1:8 w l title 'Temperature'
```
This plots column 1 (time) vs column 8 (temperature) as lines, labeled "Temperature."

**Hint:** If you want both lines and points:
```gnuplot
plot 'data.txt' u 1:2 w lp
```

**Hint:** To use columns 1 and 3 (skip column 2):
```gnuplot
plot 'data.txt' u 1:3 w l
```

**Common mistake:** Forgetting the file name or the column specification. Both are required!

**Connection:** In this course, `u 1:8` (time vs temperature) and `u 1:2` (time vs energy) are your main plots.

**Reference:** See the Gnuplot section (Part 7️⃣) for more examples and templates you can copy-paste.

********************

---

## Part 8️⃣: Synthesis & Summary

### What We Did

- ✅ **Concept:** NVT ensembles, thermostats, MD workflow  
- ✅ **Practice:** Run Curcuma MD with Berendsen and CSVR  
- ✅ **Analysis:** Extract data, analyze with AWK, plot with Gnuplot  
- ✅ **Comparison:** Glucose (Berendsen vs CSVR) and peptides (AAAAA vs WRKLQ)  

### Key Insights

1. **Berendsen is artificial** — good for first steps and understanding
2. **CSVR is realistic** — use for production
  
  - Alternatives could be Anderson Thermostat, Nose-Hover, Nose-Hover-Chain 
  
3. **Temperature should fluctuate around target** — not constant!
4. **RMSD shows structural stability** — small = stable, large = flexible
5. **Visualization is essential** — plots reveal patterns that numbers hide

### ✅ Quick Check 9: Session Summary

**Question 1:** A thermostat maintains temperature by:

- [[ ]] Directly moving atoms to cooler regions
- [[X]] Rescaling or reassigning atomic velocities (kinetic energy)
- [[ ]] Changing the force field parameters
- [[ ]] Adding random forces only to aromatic residues
********************

**Explanation:**

- ✅ **Correct:** All thermostats work by modifying **velocities**, which directly control kinetic energy and thus temperature. Berendsen rescales all velocities uniformly; CSVR rescales stochastically. Neither moves atoms or changes the force field—they only change *how fast* atoms are moving.
- ❌ **Moving atoms to cooler regions:** Thermostats don't have spatial preferences. They work globally on all atoms.
- ❌ **Changing force field parameters:** That would change the dynamics entirely. Thermostats leave the force field untouched.
- ❌ **Random forces on aromatic residues:** Thermostats treat all atoms equally, regardless of chemistry. Also, CSVR's randomness is applied to all atoms, not just aromatics.

**Practical tip:** Think of a thermostat like a **speed limiter** on a highway. It doesn't move the cars; it just prevents them from going too fast or too slow!
********************

**Question 2:** Which thermostat should you use for production MD where accurate thermodynamic properties matter?

- [[ ]] Berendsen (it's older, so it's more standard)
- [[ ]] NVE (no thermostat, purest physics)
- [[X]] CSVR (generates proper canonical ensemble)
- [[ ]] All are equally good for production
********************

**Explanation:**

- ✅ **Correct:** CSVR is the gold standard for production MD because it generates a **proper NVT ensemble** with correct Boltzmann statistics. This is essential if you want publishable thermodynamic data (free energies, fluctuations, equilibrium populations).
- ❌ **Berendsen:** Older code often used it for legacy reasons, but modern standards prefer CSVR. Berendsen is fine for equilibration only.
- ❌ **NVE:** Microcanonical ensemble, not representative of lab conditions (proteins in solution). Use only for special purposes (e.g., vacuum dynamics benchmark).
- ❌ **All equally good:** They serve different purposes. Thermostat choice significantly affects statistics.

**Industry standard:** **Phase 1 (Equilibration): Berendsen. Phase 2 (Production): CSVR.** This workflow is used in nearly every published MD paper.
********************

**Question 3:** To plot temperature vs time from MD output using Gnuplot, what command structure would you use?

- [[ ]] `plot 'data.txt' u 8 w l`
- [[ ]] `plot 'data.txt' u 2:8 w l`
- [[X]] `plot 'data.txt' u 1:8 w l title 'Temperature'`
- [[ ]] `plot 'data.txt' u 8:1 w l`
********************

**Explanation:**

- ✅ **Correct:** Gnuplot plots x-axis first, then y-axis: `u 1:8` means use column 1 (time) for x, column 8 (T) for y. The command `plot 'data.txt' u 1:8 w l title 'Temperature'` creates a line plot of temperature over time with a legend label.
- ❌ `u 8 w l`: Column 8 alone (no x-axis specified). Gnuplot won't know what x-values to use.
- ❌ `u 2:8 w l`: Column 2 is `<e_tot>` (energy), not time. This plots energy vs temperature, which is the wrong axes.
- ❌ `u 8:1 w l`: Reversed axes—time on y-axis, temperature on x-axis. Nonsensical for time-series data.

**Gnuplot column reference (from Curcuma output):**
- Column 1: Time (ps)
- Column 2: e_tot (instant energy)
- Column 3: <e_tot> (average energy)
- Column 8: T (instant temperature)
- Column 9: <T> (average temperature)

**Practical tip:** Always check your data format before plotting! Use `head data.txt` to see the column headers.

**Hint:** For a clean plot, use column 9 (`<T>`) instead of column 8 (`T`). The averaged version is less noisy!

**Reference:** See Part 7️⃣ (Gnuplot Deep Dive) for ready-to-use templates.

**Reflection:** If you can answer all three of these questions, you've mastered the core concepts of Session 2A!

********************

---

## Part 9️⃣: Troubleshooting

### Temperature Skyrockets (> 1000 K)

**Causes:**

- Bad starting geometry
- Time step too large (atoms colliding)
- Thermostat coupling too weak

**Solutions:**

- Start from a clean optimized structure
- Reduce `-dt` (time step)
- Check force field compatibility

---

### Energy Blows Up (> 10000 kcal/mol)

**Often normal initially** but should stabilize. If it doesn't:

- Atoms too close? (bad initial structure)
- Wrong method? (UFF vs GFN-FF)
- Use `-md.dt 0.5` (smaller time step) instead of default

---

### Trajectory File is Huge

`.trj.xyz` files get large quickly. You can:

- Reduce `-dump` frequency (e.g., `-dump 10` instead of `-dump 1`)
- Only keep initial and final frames
- Compress with `gzip trj.xyz` or use https://brehm-research.de/bqb.php

---

## 🎓 You Completed Session 2A!

### Skills Acquired

- ✅ Understand MD theory (ensembles, thermostats)  
- ✅ Run Curcuma MD simulations  
- ✅ Extract data from MD output  
- ✅ Create publication-quality plots with Gnuplot  
- ✅ Analyze and interpret MD results  
- ✅ Compare different thermostat methods  

### What's Next

- **Session 2B:** Introduction to **Gromacs** (the industry standard)

  - More features, more control, better for large systems
  - File formats, workflows, input files

- **Session 3:** **AlphaFold and MD**

  - Get initial structures for proteins
  - What is the purpose of each simulation method?

---

## 📚 Additional Resources

- **Curcuma Documentation:** https://github.com/conradhuebler/curcuma
- **Gnuplot Manual:** http://www.gnuplot.info/
- **CSVR Paper:** Bussi, Donadio, Parrinello, *J. Chem. Phys.* **2007**, 126, 014101
- **MD Theory:** Frenkel & Smit, *Understanding Molecular Simulation*, 2nd ed.

---

## Questions for Live Discussion

1. Why is CSVR better than Berendsen for production MD?
2. What would happen if you ran MD at 5000 K instead of 298 K?
3. Can you predict the RMSD before running a simulation?
4. How long should MD run for a protein? (100 ps vs 1 ns vs 1 µs?)
5. What does temperature mean at the atomic scale?

---

## Funding

The development of this course material is funded by:

- TUBAFdigital
- European Union (Erasmus+ National Agency for Higher Education)

Funded by the European Union. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or Erasmus+ National Agency for Higher Education (German Academic Exchange Service). Neither the European Union nor the granting authority can be held responsible for them.

![Co-funded by the European Union](https://github.com/conradhuebler/TeacherTwinMolecularModelling/raw/main/EU.jpg)

---

*Session 2A — Molecular Dynamics*  
*Last updated: Decemver 02, 2025*  
*Course: Microcredential: Modeling interactions of high molecular weight compounds*
