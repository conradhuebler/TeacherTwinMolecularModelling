<!--
author: Molecular Modelling Course Team
language: en
narrator: US English Female
version: 1.0

Session 2A: Molecular Dynamics with Curcuma
Part of: Molecular Modelling and Quantum Chemistry (Master)
-->

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

Which ensemble is most realistic for a protein in aqueous solution?
- [ ] NVE (isolated system)
- [x] NVT (constant temperature, like in a water bath)
- [ ] NPT (constant pressure)
- [ ] None of the above

In NVE ensemble, what is conserved?
- [ ] Temperature
- [x] Total energy
- [ ] Volume
- [ ] Pressure

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

Berendsen thermostat rescales velocities. What does this accomplish?
- [x] Adjusts kinetic energy to match target temperature
- [ ] Changes potential energy
- [ ] Moves atoms to new positions
- [ ] Adds random forces

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

You're running a 10 ns production MD. Best thermostat?
- [ ] Berendsen (faster equilibration)
- [x] CSVR (proper NVT ensemble for correct statistics)
- [ ] Both simultaneously
- [ ] NVE (no thermostat needed)

---

## Part 3️⃣: Practical MD with Curcuma

### Basic Syntax

```bash
curcuma -md input.xyz -method methodname -thermostat thermostat_name [options]
```

**Example:**

```bash
# MD with Berendsen thermostat for 100 ps
curcuma -md structure.xyz -method gfnff -thermostat berendsen -md.maxtime 100

# MD with CSVR thermostat for 100 ps
curcuma -md structure.xyz -method gfnff -thermostat csvr -md.maxtime 100
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

What does the `<e_tot>` column represent?
- [ ] Total energy at a single snapshot
- [x] Running average of total energy (time-averaged)
- [ ] The final energy after MD
- [ ] The initial energy before MD

---

## Part 4️⃣: Exercise 1 — Compare Thermostats on Glucose

### Objective

Run the **optimized glucose** from Session 1 through two MD simulations:
1. **Berendsen thermostat** (fast, artificial)
2. **CSVR thermostat** (slower, realistic)

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
    -md.maxtime 100 -print 10 -dump 1 > md_berendsen.dat

# This creates:
# - glucose_start.trj.xyz (trajectory)
# - md_berendsen.dat (energy/temp data)
```

**Important flags:**
- `-md.maxtime 100` = Run for 100 picoseconds
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
    -md.maxtime 100 -print 10 -dump 1 > md_csvr.dat

# Creates:
# - glucose_start.trj.xyz (overwrites previous!)
# - md_csvr.dat (energy/temp data)
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

---

### ✅ Quick Check 5: Interpreting MD Results

If total energy drifts downward throughout the MD, this suggests:
- [ ] The thermostat is working perfectly
- [x] There might be numerical instabilities or poor parameters
- [ ] It's normal and expected
- [ ] The force field is wrong

Temperature fluctuations in CSVR should be:
- [ ] Zero (perfectly constant)
- [x] Non-zero but realistic (Boltzmann distribution)
- [ ] Larger than Berendsen
- [ ] Unpredictable

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
# Output: RMSD = [RMSD_GLUCOSE_MD] Angstrom
```

**Interpretation:**
- RMSD < 0.5 Å: Structure stable (good!)
- RMSD 0.5 - 2 Å: Significant change (flexible)
- RMSD > 2 Å: Major refolding or unfolding

💾 **Record:** `[RMSD_GLUCOSE_START_TO_END_BERENDSEN]` and `[RMSD_GLUCOSE_START_TO_END_CSVR]`

---

### ✅ Quick Check 6: RMSD Interpretation

High RMSD at the end of MD means:
- [ ] The simulation failed
- [x] The structure moved/changed significantly
- [ ] Better thermostat control
- [ ] Lower temperature

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
# Use optimized peptide from Session 1
cp ../session1_peptides/peptide_aaaaa.opt.xyz pep_aaaaa_start.xyz

# MD with CSVR
curcuma -md pep_aaaaa_start.xyz -method gfnff -thermostat csvr \
    -md.maxtime 100 -print 10 > md_pep_aaaaa.dat

mv peptide_aaaaa.trj.xyz pep_aaaaa_csvr.trj.xyz
```

💾 **Record:**
- Final temperature: `[T_AAAAA_MD]` K
- Final energy: `[E_AAAAA_MD]` kcal/mol

---

### Run MD on WRKLQ

```bash
cp ../session1_peptides/peptide_wrklq.opt.xyz pep_wrklq_start.xyz

curcuma -md pep_wrklq_start.xyz -method gfnff -thermostat csvr \
    -md.maxtime 100 -print 10 > md_pep_wrklq.dat

mv peptide_wrklq.trj.xyz pep_wrklq_csvr.trj.xyz
```

💾 **Record:**
- Final temperature: `[T_WRKLQ_MD]` K
- Final energy: `[E_WRKLQ_MD]` kcal/mol

---

### Compare

| Property | AAAAA | WRKLQ |
|----------|-------|-------|
| **Final Temperature** | [T_AAAAA_MD] K | [T_WRKLQ_MD] K |
| **Final Energy** | [E_AAAAA_MD] kcal/mol | [E_WRKLQ_MD] kcal/mol |
| **Energy per Atom** | [E_PER_ATOM_AAAAA_MD] | [E_PER_ATOM_WRKLQ_MD] |
| **RMSD (start to end)** | [RMSD_AAAAA_MD] Å | [RMSD_WRKLQ_MD] Å |

---

### Analysis Questions

1. Which peptide is more stable structurally (lower RMSD)?
2. Is AAAAA more stable because it's simpler?
3. Does heterogeneity (WRKLQ) lead to more flexibility?
4. Any relationship between sequence and dynamics?

### ✅ Quick Check 7: Peptide Dynamics

Why might AAAAA be more stable than WRKLQ?
- [x] Homogeneous sequences have fewer degrees of freedom
- [ ] CSVR doesn't work on diverse sequences
- [ ] Small peptides are always stable
- [ ] GFN-FF favors uniform structures

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

What does `plot 'data.txt' u 1:2 w l` mean?
- [ ] Use only column 1 for plotting
- [x] Plot column 1 vs column 2 with lines
- [ ] Plot data with points
- [ ] Use logarithmic scale

---

## Part 8️⃣: Synthesis & Summary

### What We Did

✅ **Concept:** NVT ensembles, thermostats, MD workflow  
✅ **Practice:** Run Curcuma MD with Berendsen and CSVR  
✅ **Analysis:** Extract data, analyze with AWK, plot with Gnuplot  
✅ **Comparison:** Glucose (Berendsen vs CSVR) and peptides (AAAAA vs WRKLQ)  

### Key Insights

1. **Berendsen is faster but artificial** — good for equilibration
2. **CSVR is slower but realistic** — use for production
3. **Temperature should fluctuate around target** — not constant!
4. **RMSD shows structural stability** — small = stable, large = flexible
5. **Visualization is essential** — plots reveal patterns that numbers hide

### ✅ Quick Check 9: Session Summary

After completing Session 2A, you understand:
- [x] How thermostats control temperature in MD simulations
- [x] The difference between Berendsen and CSVR
- [x] How to extract and plot MD data
- [ ] All of the above

---

## Part 9️⃣: Troubleshooting

### Temperature Skyrockets (> 1000 K)

**Causes:**
- Bad starting geometry
- Time step too large (atoms colliding)
- Thermostat coupling too weak

**Solutions:**
- Start from a clean optimized structure
- Reduce `-md.dt` (time step)
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
- Compress with `gzip trj.xyz`

---

## 🎓 You Completed Session 2A!

### Skills Acquired

✅ Understand MD theory (ensembles, thermostats)  
✅ Run Curcuma MD simulations  
✅ Extract data from MD output  
✅ Create publication-quality plots with Gnuplot  
✅ Analyze and interpret MD results  
✅ Compare different thermostat methods  

### What's Next

- **Session 2B:** Introduction to **Gromacs** (the industry standard)
  - More features, more control, better for large systems
  - File formats, workflows, input files

- **Session 3:** **AlphaFold vs. MD**
  - How do predictions compare to simulations?
  - When does each method succeed or fail?

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

*Session 2A — Molecular Dynamics (Curcuma Thermostats)*  
*Last updated: October 27, 2025*  
*Course: Molecular Modelling and Quantum Chemistry (Master)*
