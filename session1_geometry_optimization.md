<!--
author: Molecular Modelling Course Team
language: en
narrator: US English Female
version: 1.0

Session 1: Geometry Optimization
Part of: Molecular Modelling and Quantum Chemistry (Master)
-->

# Session 1: Geometry Optimization

## Your First Steps into Molecular Modelling

> **Welcome to Session 1!**
>
> In this session, you will learn what **geometry optimization** means, 
> why it matters, and how to perform it on real molecules.
>
> You will optimize **sugars** and **peptides** using different methods 
> and compare the results. All materials are self-contained — 
> you can work through this at your own pace.

---

## 🎯 Learning Objectives

By the end of this session, you will:

- ✅ Understand **geometry optimization** from first principles
- ✅ Know the difference between **force fields** (UFF vs. GFN-FF)
- ✅ Use **Curcuma** to optimize molecular structures
- ✅ Analyze results: **RMSD**, **energy changes**, **structural changes**
- ✅ Interpret optimization trajectories in **Avogadro**
- ✅ Compare optimization across **different methods**

---

## Part 1️⃣: Concept — Geometry Optimization from the Beginning

### What is Geometry Optimization?

**Geometry optimization** is the process of finding the **minimum energy configuration** of a molecule by adjusting atomic positions.

Think of it like this:

```
        Energy
           ▲
           │     ┌─────────────┐
           │    ╱               ╲
           │   ╱                 ╲
           │  ╱                   ╲
           │ ╱       ← LOCAL MIN    ╲
           │╱___OPTIMIZATION__________╲___► Atomic Positions
           
The algorithm searches for the valley (minimum)
```

---

### Why Do We Optimize Geometries?

1. **Experimental relevance:** Molecules naturally adopt low-energy conformations
2. **Accurate starting point:** For MD simulations, energy calculations, reactivity
3. **Structure determination:** We need realistic 3D coordinates
4. **Energy comparison:** Only meaningful when both structures are optimized

**Real-world example:** A protein in solution adopts a structure that minimizes its free energy. Before running a 100 ns MD simulation, we want to start from a realistic geometry.

---

### The Born-Oppenheimer (BO) Approximation

Geometry optimization relies on the **Born-Oppenheimer Approximation**:

> *Nuclei move slowly compared to electrons. We can calculate electronic energy at fixed nuclear positions, then use that energy to guide nuclear motion.*

**In practice:**
- Fix atomic positions → Calculate forces on each atom
- Move atoms in direction of forces (downhill in energy)
- Repeat until forces are near zero (converged)

**Key insight:** The energy surface is the **same for all methods** (UFF, GFN-FF, etc.), but they calculate it differently.

### ✅ Quick Check 1: BO Approximation

[[?]]
| What does the Born-Oppenheimer approximation assume?
| - [[X]] Nuclei move slowly compared to electrons, so we can decouple them
| - [[ ]] Electrons and nuclei move at the same speed
| - [[ ]] Nuclear motion doesn't affect electronic structure
| - [[ ]] The energy surface has no local minima

---

### Convergence Criteria

Optimization stops when these criteria are met:

| Criterion | Meaning | Typical Threshold |
|-----------|---------|-------------------|
| **Energy change** (dE) | Final energy doesn't change much | < 0.1 kcal/mol |
| **RMSD** (Atomic movement) | Atoms barely move between steps | < 0.01 Å |
| **Gradient norm** | Forces on atoms are nearly zero | < 0.001 kcal/(mol·Å) |

**All three must be satisfied** for a robust optimization.

### ✅ Quick Check 2: Convergence

[[?]]
| An optimization is converged when:
| - [[X]] Energy, atomic movement, and forces are all below thresholds
| - [[ ]] Energy stops changing completely
| - [[ ]] Only one criterion is satisfied
| - [[ ]] Maximum iterations is reached

---

## Part 2️⃣: Methods — UFF vs. GFN-FF

### Force Fields: What Are They?

A **force field** is a mathematical model that calculates molecular energy based on atomic positions. It's **fast** but **approximate**.

General form:
```
E_total = E_bonds + E_angles + E_dihedrals + E_nonbonded + ...
```

Each term describes how atoms interact.

---

### Method 1: UFF (Universal Force Field)

**UFF** = Universal Force Field

- **Type:** Classical force field
- **Based on:** Experimental bond lengths, angles, empirical parameters
- **Gradients:** Calculated **numerically** (slower but robust)
- **Speed:** Fast
- **Accuracy:** Reasonable for most organic molecules

**When to use UFF:**
- Large molecules (> 1000 atoms)
- Initial geometry guess (rough optimization)
- Conformational searching

**Limitations:**
- No quantum effects
- Fixed bond connectivity

---

### Method 2: GFN-FF (Geometry Function Nonbonding - Force Field)

**GFN-FF** = Semiempirical force field

- **Type:** Hybrid quantum/classical
- **Based on:** Quantum-derived parameters (xtb)
- **Gradients:** Analytical (faster than numerical)
- **Speed:** Faster than some QM, much faster than full QM
- **Accuracy:** Better than UFF for many systems (especially with heteroatoms, metals)

**When to use GFN-FF:**
- More accurate results needed
- Heteroatoms important (N, O, S, halides)
- Medium-sized molecules

**Advantages over UFF:**
- Includes dispersion (van der Waals)
- Better for polarized systems

---

### Quick Comparison

| Property | UFF | GFN-FF |
|----------|-----|--------|
| Speed | ⚡⚡⚡ | ⚡⚡ |
| Accuracy | ⭐⭐ | ⭐⭐⭐ |
| Quantum effects | ✗ | ✓ (partial) |
| Dispersion | ✗ | ✓ |
| Good for heteroatoms | ✗ | ✓ |
| Molecular size | Large OK | Medium best |

### ✅ Quick Check 3: When to Use Which Method

[[?]]
| You're optimizing a protein with 2000 atoms. Best choice?
| - [[X]] UFF (speed is critical for large systems)
| - [[ ]] GFN-FF (more accurate)
| - [[ ]] Both simultaneously
| - [[ ]] Neither works for proteins

[[?]]
| A molecule contains many sulfur atoms. Which method is preferable?
| - [[ ]] UFF (classical methods handle S well)
| - [[X]] GFN-FF (better for heteroatoms)
| - [[ ]] Both equally good
| - [[ ]] Neither can handle S

---

## Part 3️⃣: Practical Setup — Using Curcuma

### What is Curcuma?

**Curcuma** is an open-source molecular modelling tool developed by Conrad Huebler.

**It can do:**
- ✅ Geometry optimization (GeoOpt)
- ✅ Molecular dynamics (MD)
- ✅ Conformation searching (ConfScan)
- ✅ RMSD analysis and trajectory reordering
- ✅ Multiple force fields (UFF, GFN-FF, GFN1, GFN2)

**Command structure:**
```bash
curcuma -opt input.xyz -method methodname [options]
```

---

### Getting Input Structures

**Step 1: Find SMILES**

Go to **ChemSpider** (https://www.chemspider.com/):

1. Search for your molecule (e.g., "Glucose")
2. Find the correct entry
3. Copy the **SMILES string** (looks like: `C(C1C(C(C(C(O1)O)O)O)O)O`)

**Step 2: Convert SMILES to 3D Structure**

Use **Open Babel** (obabel):

```bash
obabel -ismi input.smi -O output.xyz -h --gen3d
```

- `-ismi` = input is SMILES
- `-O output.xyz` = output as XYZ file
- `-h` = add hydrogens
- `--gen3d` = generate 3D coordinates

**Example:** Glucose SMILES → 3D structure
```bash
echo "C(C1C(C(C(C(O1)O)O)O)O)O" > glucose.smi
obabel -ismi glucose.smi -O glucose.xyz -h --gen3d
```

Now you have `glucose.xyz` ready for optimization!

---

### Run Optimization with Curcuma

**Basic syntax:**

```bash
curcuma -opt input.xyz -method methodname
```

**Example 1: Optimize glucose with UFF**

```bash
curcuma -opt glucose.xyz -method uff
```

**Output files:**
- `glucose.opt.xyz` — Final optimized structure
- `glucose.trj.xyz` — Optimization trajectory (all steps)

**Example 2: Optimize with GFN-FF**

```bash
curcuma -opt glucose.xyz -method gfnff
```

**Example 3: Optimize peptide with GFN-FF (faster)**

```bash
curcuma -opt peptide.xyz -method gfnff -threads 4
```

---

### ✅ Quick Check 4: Curcuma Workflow

[[?]]
| What does the `.trj.xyz` file contain?
| - [[ ]] Only the final optimized structure
| - [[X]] The entire optimization trajectory (all intermediate steps)
| - [[ ]] The initial structure only
| - [[ ]] Energy values at each step

[[?]]
| Which command correctly optimizes a molecule?
| - [[ ]] `curcuma glucose.xyz`
| - [[X]] `curcuma -opt glucose.xyz -method gfnff`
| - [[ ]] `optimize glucose.xyz`
| - [[ ]] `curcuma glucose.xyz --optimize`

---

## Part 4️⃣: Exercise 1 — Optimize a Monosaccharide

### Task: Optimize Glucose with Two Methods

**Molecule:** Glucose (C₆H₁₂O₆)

**What you'll do:**
1. Get glucose structure (SMILES → 3D)
2. Optimize with **UFF**
3. Optimize with **GFN-FF**
4. Compare results
5. Analyze trajectory

---

### Step 1: Get the Structure

**Go to ChemSpider:**
- Search: "Glucose"
- Find: D-Glucose (or just Glucose)
- Copy SMILES: `C(C1C(C(C(C(O1)O)O)O)O)O`

**Convert to 3D:**

```bash
mkdir -p session1_glucose
cd session1_glucose

echo "C(C1C(C(C(C(O1)O)O)O)O)O" > glucose.smi
obabel -ismi glucose.smi -O glucose.xyz -h --gen3d
```

**Check your structure:**
```bash
cat glucose.xyz
```

You should see something like:
```
21
        0.0000    0.0000    0.0000
C        0.1234    0.5678    1.2345
O       -0.9876    ...
...
```

---

### Step 2: Optimize with UFF

```bash
curcuma -opt glucose.xyz -method uff
```

**This will create:**
- `glucose.opt.xyz` — Final structure
- `glucose.trj.xyz` — Trajectory of all steps

**Expected output (on screen):**
```
Step 1: E = -234.5 kcal/mol, Force norm = 0.450
Step 2: E = -238.2 kcal/mol, Force norm = 0.125
...
Converged: E = [EXPECTED_ENERGY_GLUCOSE_UFF]
```

💾 **Note your final energy:** `[EXPECTED_ENERGY_GLUCOSE_UFF]` kcal/mol

---

### Step 3: Optimize with GFN-FF

```bash
# Rename the previous opt file
mv glucose.opt.xyz glucose_uff.opt.xyz
mv glucose.trj.xyz glucose_uff.trj.xyz

# Optimize with GFN-FF
curcuma -opt glucose.xyz -method gfnff
```

**Final structure:** `glucose.opt.xyz`

💾 **Note your final energy:** `[EXPECTED_ENERGY_GLUCOSE_GFNFF]` kcal/mol

---

### Step 4: Calculate RMSD (Structural Difference)

Now compare the two optimized structures:

```bash
curcuma -rmsd glucose_uff.opt.xyz glucose.opt.xyz
```

**Output:** RMSD = `[RMSD_GLUCOSE_UFF_vs_GFNFF]` Å

**Interpretation:**
- RMSD < 0.1 Å: Very similar geometries
- RMSD 0.1 - 0.5 Å: Noticeable differences
- RMSD > 0.5 Å: Significantly different

### ✅ Quick Check 5: RMSD Meaning

[[?]]
| RMSD measures:
| - [[X]] The average distance atoms moved between two structures
| - [[ ]] The difference in energy between two structures
| - [[ ]] The number of bonds that changed
| - [[ ]] How long optimization took

---

### Step 5: Visualize the Trajectory

**Open Avogadro:**

```bash
avogadro glucose_uff.trj.xyz &
```

(The `&` runs it in the background)

**In Avogadro:**
1. File → Open → `glucose_uff.trj.xyz`
2. **Animate** (usually bottom toolbar: play button ▶)
3. Watch the atoms move during optimization!
4. You can **slow down/speed up** the animation

**What to observe:**
- Do atoms move a lot in early steps?
- Do movements get smaller as you progress?
- Any sudden jumps? (would indicate convergence issues)

**Repeat with GFN-FF trajectory:**
```bash
avogadro glucose_uff.trj.xyz &
avogadro glucose.trj.xyz &
# (Open both side-by-side in separate windows)
```

---

### Summary: Glucose Results

| Property | UFF | GFN-FF |
|----------|-----|--------|
| **Final Energy** | [EXPECTED_ENERGY_GLUCOSE_UFF] | [EXPECTED_ENERGY_GLUCOSE_GFNFF] |
| **Energy Difference** | [EXPECTED_DELTA_E_GLUCOSE] | |
| **Optimization Steps** | [EXPECTED_STEPS_GLUCOSE_UFF] | [EXPECTED_STEPS_GLUCOSE_GFNFF] |
| **RMSD to GFN-FF result** | [RMSD_GLUCOSE_UFF_vs_GFNFF] Å | — |

### ✅ Quick Check 6: Interpreting Results

[[?]]
| If GFN-FF gives much lower energy than UFF, this suggests:
| - [[X]] GFN-FF found a different (more stable) conformation
| - [[ ]] GFN-FF has a bug
| - [[ ]] UFF is always wrong
| - [[ ]] The structures are identical

---

## Part 5️⃣: Exercise 2 — Disaccharide (Sucrose)

### Task: Quick Comparison

**Molecule:** Sucrose (table sugar, C₁₂H₂₂O₁₁)

**What to do:**
1. Get SMILES from ChemSpider
2. Convert to 3D
3. Optimize with **GFN-FF only** (it's the better method)
4. Note the final energy
5. Compare trajectory length to glucose

---

### Step 1-2: Get Structure

```bash
mkdir -p session1_sucrose
cd session1_sucrose

# From ChemSpider: Sucrose SMILES
echo "[SUCROSE_SMILES_PLACEHOLDER]" > sucrose.smi
obabel -ismi sucrose.smi -O sucrose.xyz -h --gen3d
```

**If you can't find the exact SMILES, use this:**
```
C(C1C(C(C(C(O1)OC2(C(C(C(O2)CO)O)O)CO)O)O)O)O
```

---

### Step 3: Optimize

```bash
curcuma -opt sucrose.xyz -method gfnff
```

💾 **Final energy:** `[EXPECTED_ENERGY_SUCROSE_GFNFF]` kcal/mol

---

### Step 4: Analysis

**Questions to consider:**

- Is sucrose heavier than glucose? (Yes, it has more atoms)
- Is the final energy lower or higher per atom? (compute: Energy / number_of_atoms)
- How many optimization steps did it take?

### ✅ Quick Check 7: Molecular Size Effects

[[?]]
| Sucrose has more atoms than glucose. Its absolute energy should be:
| - [[ ]] Higher (same) as glucose
| - [[X]] Lower (more negative) than glucose
| - [[ ]] The same as glucose
| - [[ ]] Unpredictable

---

## Part 6️⃣: Exercise 3 — Peptides (Two Extremes)

### Why Compare Two Extreme Peptides?

Two very different amino acid sequences will show how **chemical diversity affects optimization**:

- **AAAAA** (Polyalanine) — Homogeneous, small side chains, simple
- **WRKLQ** (Heterogeneous) — Large/diverse: Trp (aromatic), Arg (charged), Lys (charged), Leu (hydrophobic), Gln (polar)

---

### Task: Build and Optimize Both

#### Peptide 1: AAAAA (Poly-Alanine)

**Get SMILES from ChemSpider or use:**
```
[AAAAA_PEPTIDE_SMILES_PLACEHOLDER]
```

**Build structure:**
```bash
mkdir -p session1_peptides
cd session1_peptides

echo "[AAAAA_PEPTIDE_SMILES_PLACEHOLDER]" > peptide_aaaaa.smi
obabel -ismi peptide_aaaaa.smi -O peptide_aaaaa.xyz -h --gen3d
```

**Optimize:**
```bash
curcuma -opt peptide_aaaaa.xyz -method gfnff
```

💾 **Final energy:** `[EXPECTED_ENERGY_AAAAA_GFNFF]` kcal/mol
💾 **Atoms:** Count from structure = `[N_ATOMS_AAAAA]`
💾 **Energy per atom:** `[E_PER_ATOM_AAAAA]` kcal/mol

---

#### Peptide 2: WRKLQ (Heterogeneous Sequence)

**Build structure:**
```bash
echo "[WRKLQ_PEPTIDE_SMILES_PLACEHOLDER]" > peptide_wrklq.smi
obabel -ismi peptide_wrklq.smi -O peptide_wrklq.xyz -h --gen3d
```

**Optimize:**
```bash
curcuma -opt peptide_wrklq.xyz -method gfnff
```

💾 **Final energy:** `[EXPECTED_ENERGY_WRKLQ_GFNFF]` kcal/mol
💾 **Atoms:** `[N_ATOMS_WRKLQ]`
💾 **Energy per atom:** `[E_PER_ATOM_WRKLQ]` kcal/mol

---

### Comparison Table

| Property | AAAAA | WRKLQ |
|----------|-------|-------|
| **Sequence** | Homogeneous | Heterogeneous |
| **Final Energy** | [EXPECTED_ENERGY_AAAAA_GFNFF] | [EXPECTED_ENERGY_WRKLQ_GFNFF] |
| **Number of Atoms** | [N_ATOMS_AAAAA] | [N_ATOMS_WRKLQ] |
| **Energy per Atom** | [E_PER_ATOM_AAAAA] | [E_PER_ATOM_WRKLQ] |
| **Optimization Steps** | [STEPS_AAAAA_GFNFF] | [STEPS_WRKLQ_GFNFF] |

---

### Analysis Questions

1. **Which peptide has lower energy per atom?**
   - Suggests which is more "stable" at this size scale

2. **Which took more optimization steps?**
   - Complex molecules often need more steps (more degrees of freedom)

3. **Visualize both trajectories (Avogadro):**
   ```bash
   avogadro peptide_aaaaa.trj.xyz &
   avogadro peptide_wrklq.trj.xyz &
   ```
   - Do the residues move differently?
   - Does WRKLQ seem "floppy" or "rigid"?

### ✅ Quick Check 8: Peptide Complexity

[[?]]
| WRKLQ should take more optimization steps than AAAAA because:
| - [[X]] More residue diversity → more degrees of freedom
| - [[ ]] It's always slower regardless of sequence
| - [[ ]] Larger residues always need more steps
| - [[ ]] GFN-FF is inefficient with long molecules

---

## Part 7️⃣: Optional — Fructose (Alternative Monosaccharide)

If you want to explore further:

**Fructose** is a structural isomer of glucose (same formula C₆H₁₂O₆, different structure).

**Prediction:** If two isomers have identical molecular formulas but different connectivity, their optimized energies should differ. By how much?

**Steps:**

```bash
# Get Fructose
echo "[FRUCTOSE_SMILES_PLACEHOLDER]" > fructose.smi
obabel -ismi fructose.smi -O fructose.xyz -h --gen3d

# Optimize
curcuma -opt fructose.xyz -method gfnff

# Compare to glucose
curcuma -rmsd glucose.opt.xyz fructose.opt.xyz
```

💾 **Fructose energy:** `[EXPECTED_ENERGY_FRUCTOSE_GFNFF]` kcal/mol
💾 **Energy difference (Fructose - Glucose):** `[DELTA_E_ISOMERS]` kcal/mol
💾 **Structural RMSD:** `[RMSD_GLUCOSE_vs_FRUCTOSE]` Å

---

## Part 8️⃣: Synthesis & Key Takeaways

### What We Learned

✅ **Geometry optimization** finds the lowest energy structure  
✅ **Force fields** (UFF, GFN-FF) are different approximations  
✅ **Convergence** requires energy, force, and RMSD criteria  
✅ **RMSD analysis** shows structural differences  
✅ **Trajectory visualization** reveals optimization dynamics  

---

### Key Concepts to Remember

1. **BO Approximation:** Nuclei move much slower than electrons
2. **Force fields have limitations:** No method is "best" for all cases
3. **Optimization is not guaranteed to find the global minimum:** Only local minima
4. **Energy differences matter:** Between methods and between structures
5. **Visualization is crucial:** Trajectories tell you if something went wrong

### ✅ Quick Check 9: BO Approximation Revisited

[[?]]
| Why is the Born-Oppenheimer approximation useful?
| - [[X]] It decouples nuclear and electronic motion, making calculations faster
| - [[ ]] It guarantees finding the global minimum
| - [[ ]] It eliminates the need for any numerical calculations
| - [[ ]] It makes all force fields equivalent

### ✅ Quick Check 10: When to Optimize

[[?]]
| Before running an MD simulation, why optimize the geometry first?
| - [[X]] To start from a realistic, low-energy conformation
| - [[ ]] MD doesn't need optimization
| - [[ ]] Optimization guarantees the MD will succeed
| - [[ ]] To speed up the MD calculation

---

## Part 9️⃣: Troubleshooting

### "Optimization didn't converge"

**Causes:**
- Starting geometry too distorted
- Method not suitable for your molecule
- Maximum iterations reached

**Solutions:**
- Visualize `.trj.xyz` — did atoms move too much?
- Try different method (UFF → GFN-FF or vice versa)
- Increase `MaxIter` in Curcuma settings

---

### "RMSD is huge (> 1 Å)"

**Possible reasons:**
- Two different conformations (different minima)
- Different protonation states
- Atom reordering messed up

**Check:**
- Visualize both structures in Avogadro
- Use `curcuma -rmsd struct1.xyz struct2.xyz -reorder`

---

### "Energy is very negative (< -1000 kcal/mol)"

**Not a bug!** Absolute energies depend on:
- Number of atoms
- Method (UFF vs GFN-FF)
- Unit conventions

**What matters:** Energy **differences** between structures

---

## 🎓 You Completed Session 1!

### Summary of Skills Acquired

✅ Structure preparation (SMILES → 3D with obabel)  
✅ Running Curcuma optimization  
✅ Analyzing results (energy, RMSD, trajectories)  
✅ Interpreting convergence  
✅ Comparing multiple methods  
✅ Using Avogadro for visualization  

### What's Next

- **Session 2:** Molecular Dynamics (MD)
  - How to simulate molecular motion over time
  - Thermostats (Berendsen vs. CSVR)
  - Temperature, energy, and trajectory analysis

- **Session 3:** AlphaFold vs. MD
  - Compare predicted structures to MD simulations
  - When does AlphaFold work? When does MD fail?

---

## 📚 Additional Resources

- **Curcuma GitHub:** https://github.com/conradhuebler/curcuma
- **Open Babel Documentation:** https://openbabel.org/
- **Avogadro Manual:** https://avogadro.cc/docs/
- **ChemSpider:** https://www.chemspider.com/
- **Force Field Theory:** Rappe et al. JACS 1992 (UFF paper)

---

## Questions for Discussion (in live Session)

1. Why did GFN-FF and UFF give different energies for glucose?
2. Did AAAAA optimize faster than WRKLQ? Why?
3. What does the trajectory tell us about the optimization process?
4. When would you use UFF instead of GFN-FF?
5. Can you think of molecules where geometry optimization might fail?

---

*Session 1 — Geometry Optimization*  
*Last updated: October 27, 2025*  
*Course: Molecular Modelling and Quantum Chemistry (Master)*
