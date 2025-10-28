<!--
author: Molecular Modelling Course Team
language: en
narrator: US English Female
version: 1.0

Session 3: AlphaFold and Molecular Dynamics
Part of: Molecular Modelling and Quantum Chemistry (Master)
-->

# Session 3: AlphaFold and Molecular Dynamics — Complementary Approaches

## Understanding Protein Structure AND Function

> **Sessions 1-2** taught you the tools: geometry optimization, MD simulations, Gromacs workflows.
>
> **Session 3** asks the bigger question: **How do we use these methods together?**
>
> AlphaFold and MD are not competitors — they answer different questions about proteins.
> This session teaches you how to combine them for deeper insights.

---

## 🎯 Learning Objectives

By the end of Session 3, you will:

- ✅ Understand what AlphaFold predicts and what it cannot
- ✅ Use MD to **validate** AlphaFold predictions
- ✅ Analyze MD trajectories to extract dynamics information
- ✅ Interpret per-residue flexibility (RMSF)
- ✅ Identify stable vs. flexible regions in proteins
- ✅ Combine structural predictions with functional dynamics
- ✅ Communicate scientific findings in a professional presentation

---

## Part 1️⃣: The Complementary Nature of AlphaFold and MD

### What Each Method Answers

**AlphaFold asks:** "What is the folded structure?"

- Input: Amino acid sequence
- Output: 3D structure coordinates + confidence scores
- Timescale: Instantaneous (static)
- Limitation: Single snapshot, no dynamics

**Molecular Dynamics asks:** "How does the structure behave?"

- Input: Starting structure + force field
- Output: Trajectory of all atomic positions over time
- Timescale: Nanoseconds to microseconds
- Limitation: Depends on starting structure and force field quality

**Together they ask:** "Is this structure stable? How does it move? Does it match what we know biologically?"

---

### Why They Complement Each Other

```
AlphaFold Result (Static)
      ↓
      Is it stable?
      ├─→ MD Simulation
      │   ├─ Energy landscape
      │   ├─ Atomic motion
      │   ├─ Flexibile regions
      │   └─ Ensemble
      └─→ MD Result (Dynamic)
          ├─ Confirms structure
          ├─ Reveals movement
          ├─ Shows biology
          └─ Validates prediction
```

### ✅ Quick Check 1: Complementarity

[[?]]
| What does AlphaFold reveal that MD cannot?
| - [[X]] The most likely native 3D fold from sequence alone
| - [[ ]] How atoms move over time
| - [[ ]] Protein dynamics at room temperature
| - [[ ]] Binding affinities

[[?]]
| What does MD reveal that AlphaFold cannot?
| - [[ ]] The native fold of unknown sequences
| - [[X]] Atomic motion, flexibility, and ensemble properties
| - [[ ]] Confidence in the prediction
| - [[ ]] Sequence homology information

---

## Part 2️⃣: Validation — Is the AlphaFold Structure Stable?

### Why Validate?

AlphaFold is trained on known structures, but:
- Proteins in solution are **dynamic**, not static
- AlphaFold may predict a **kinetically trapped state** (not thermodynamically stable)
- Some regions may be **artificially stable** in the training data
- Local forces may destabilize certain folds

**Validation with MD:** Run the AlphaFold structure through MD and check:
- Does it stay folded (low RMSD)?
- Does energy converge?
- Are there large rearrangements?

---

### Stability Criteria

| Metric | Interpretation | Good Value | Bad Value |
|--------|----------------|------------|-----------|
| **RMSD (vs AlphaFold)** | How much structure deviates | < 2 Å | > 5 Å |
| **Radius of Gyration** | Overall compactness | Stable, low fluctuations | Drifting up |
| **Energy (E_total)** | System stability | Converged after equilibration | Drifting/exploding |
| **Temperature (T)** | Thermal equilibration | Constant ~300 K | Fluctuating wildly |

### ✅ Quick Check 2: Stability

[[?]]
| If RMSD stays below 1 Å during a 500 ps MD, this suggests:
| - [[X]] The AlphaFold structure is stable under the force field
| - [[ ]] The force field is incorrect
| - [[ ]] The structure is too rigid
| - [[ ]] Nothing — RMSD is irrelevant

[[?]]
| A monotonically increasing RMSD over time indicates:
| - [[ ]] The simulation is working correctly
| - [[X]] The structure is unfolding or large rearrangements occur
| - [[ ]] Temperature is too low
| - [[ ]] Normal protein breathing

---

## Part 3️⃣: Analysis Methods — What to Extract from MD

### 1. Total Energy Over Time

**What it shows:** System equilibration and stability

**How to extract:**
```bash
gmx energy -f prod.edr -o energy.xvg
# Select "Total Energy" when prompted
```

**Interpretation:**
- First 50 ps: E might decrease (equilibration)
- After 50 ps: E should plateau around a mean value
- Large fluctuations: Normal (kT thermal energy)
- Monotonic drift: Problem (numerical instability, bad FF)

---

### 2. RMSD vs. Starting Structure

**What it shows:** How much the structure changed from AlphaFold input

**How to extract:**
```bash
gmx rms -f prod.xtc -s prod.tpr -o rmsd.xvg
# Select backbone C-alpha for protein
```

**Interpretation:**
- RMSD < 1 Å: Structure very stable (trustworthy)
- RMSD 1-3 Å: Some motion but overall stable (normal)
- RMSD > 3 Å: Large-scale rearrangement (investigate)
- Plateau: Good (converged)
- Linear increase: Unfolding (watch out!)

---

### 3. RMSF — Per-Residue Flexibility

**What it shows:** Which parts of the protein move most

**How to extract:**
```bash
gmx rmsf -f prod.xtc -s prod.tpr -o rmsf.xvg -res
# Gives flexibility for each residue
```

**Interpretation:**
- Low RMSF (< 1 Å): Rigid, probably important for structure
- High RMSF (2-5 Å): Flexible, likely loops or tails
- Very high (> 5 Å): Might indicate unfolded region or simulation artifact

**Biological meaning:**
- Secondary structure (α-helix, β-sheet): Usually low RMSF
- Loops connecting secondary structures: Usually high RMSF
- Binding sites: Often have intermediate RMSF (flexible for function)

### ✅ Quick Check 3: Analysis Methods

[[?]]
| What does RMSF (Root-Mean-Square Fluctuation) measure?
| - [[ ]] Average change between two structures
| - [[X]] How much each individual residue moves during MD
| - [[ ]] The total energy of the system
| - [[ ]] Temperature stability

[[?]]
| If a region has consistently high RMSF, it probably is:
| - [[ ]] Incorrectly folded
| - [[X]] Flexible, possibly functionally relevant (hinge, loop)
| - [[ ]] Unfolding
| - [[ ]] A simulation error

---

## Part 4️⃣: AlphaFold Input & Output

### Getting an AlphaFold Structure

**Option 1: AlphaFold Server (web-based)**

```
1. Go to: [ALPHAFOLD_LINK_PLACEHOLDER]
2. Paste protein sequence (FASTA format)
3. Wait for computation (~minutes to hours)
4. Download: model_final.pdb
```

**Option 2: Local Installation**

```bash
# (Setup provided by instructor, if available)
colabfold_search <input.fasta> <output_dir>
# Generates: ranked_0.pdb (best model)
```

**Option 3: Use PDB Database**

```
1. Go to: https://www.rcsb.org/
2. Search for protein (name or PDB ID)
3. Download: protein.pdb
4. (Use experimental structure, not AlphaFold prediction)
```

### AlphaFold Output Files

**Standard outputs:**
- `model_final.pdb` — Predicted structure (or ranked_0.pdb)
- Confidence scores (pLDDT per residue)
- PAE (Predicted Aligned Error) matrix

### Understanding Confidence Scores

**pLDDT (per-residue confidence):**
- **90+**: Very confident (light blue)
- **70-90**: Confident (blue)
- **50-70**: Moderate (yellow)
- **<50**: Low confidence (red)

**Interpretation:**
- High confidence regions are usually core structure
- Low confidence regions might be flexible or poorly trained
- Compare MD flexibility: do high-confidence = rigid?

### ✅ Quick Check 4: AlphaFold

[[?]]
| A region with pLDDT < 50 in AlphaFold prediction suggests:
| - [[X]] Uncertain prediction, possibly flexible or unstructured
| - [[ ]] The region is definitely incorrect
| - [[ ]] The entire protein is wrong
| - [[ ]] The region is very important

---

## Part 5️⃣: Exercise — Protein Analysis Workflow

### Your Task

You have:
- ✅ **AlphaFold prediction** (your protein, Session 2B input)
- ✅ **MD simulation results** (instructor ran production.tpr → prod.xtc, prod.edr)

Now analyze and compare them.

---

### Step 1: Get MD Results from Instructor

Files you'll receive:
```
prod.xtc        — Trajectory (all atomic positions)
prod.edr        — Energy/temperature data
prod.tpr        — Binary run file (needed for gmx commands)
prod.log        — Simulation log
```

**Copy to your working directory:**
```bash
mkdir -p session3_analysis
cd session3_analysis
# (Instructor will provide these files)
```

---

### Step 2: Extract Energy Data

```bash
gmx energy -f prod.edr -o energy.xvg
# When prompted, select:
# - "Total Energy"
# - Press Enter to confirm
# - Type "q" to quit
```

**Output:** `energy.xvg` (text file with energy over time)

**Examine:**
```bash
head -20 energy.xvg
tail -20 energy.xvg
```

**Extract just numbers (skip headers):**
```bash
grep -v "^@\|^#" energy.xvg > energy_clean.txt
awk '{print $1, $2}' energy_clean.txt > time_vs_energy.txt
```

---

### Step 3: Calculate RMSD vs. AlphaFold

First, align MD trajectory to AlphaFold structure:

```bash
# Convert AlphaFold PDB to Gromacs format (gro)
gmx editconf -f alphafold_structure.pdb -o af_structure.gro

# Calculate RMSD
gmx rms -f prod.xtc -s prod.tpr -o rmsd_vs_af.xvg
# When prompted, select "C-alpha" for protein
```

**Output:** `rmsd_vs_af.xvg`

**Extract:**
```bash
grep -v "^@\|^#" rmsd_vs_af.xvg > rmsd_clean.txt
awk '{print $1, $2}' rmsd_clean.txt > time_vs_rmsd.txt

# Check final RMSD
tail -1 rmsd_clean.txt
```

💾 **Record:** `[FINAL_RMSD_VS_ALPHAFOLD]` Å

---

### Step 4: Calculate RMSF (Per-Residue Flexibility)

```bash
gmx rmsf -f prod.xtc -s prod.tpr -o rmsf.xvg -res
```

**Output:** `rmsf.xvg` (one line per residue)

**Visualize:**
```bash
# Extract for plotting
grep -v "^@\|^#" rmsf.xvg > rmsf_clean.txt
awk '{print $1, $2}' rmsf_clean.txt > residue_vs_rmsf.txt
```

---

### Step 5: Create Plots

**Gnuplot script for analysis** (`plot_session3.gnu`):

```gnuplot
set terminal png size 1200,1200
set output 'session3_analysis.png'
set multiplot layout 2,2

# Plot 1: Energy
set title 'Total Energy During MD'
set xlabel 'Time (ps)'
set ylabel 'Energy (kJ/mol)'
set grid
plot 'time_vs_energy.txt' u 1:2 w l lw 2 title 'Total Energy'

# Plot 2: RMSD vs AlphaFold
set title 'RMSD vs AlphaFold Structure'
set xlabel 'Time (ps)'
set ylabel 'RMSD (Angstrom)'
set grid
plot 'time_vs_rmsd.txt' u 1:2 w l lw 2 title 'RMSD'

# Plot 3: RMSF per residue
set title 'Per-Residue Flexibility (RMSF)'
set xlabel 'Residue Number'
set ylabel 'RMSF (Angstrom)'
set grid
plot 'residue_vs_rmsf.txt' u 1:2 w l lw 2 title 'RMSF'

# Plot 4: Stability summary (combined)
set title 'Simulation Summary: Is the structure stable?'
set xlabel 'Time (ps)'
set ylabel 'Normalized Value'
set grid
# (Requires normalization, optional)

unset multiplot
```

**Run Gnuplot:**
```bash
gnuplot plot_session3.gnu
display session3_analysis.png
```

---

### Step 6: Interpret Your Results

**Ask yourself these questions:**

1. **Energy stability:**
   - Did energy converge after equilibration?
   - Are there unexplained spikes?
   - Is the mean energy reasonable?

2. **RMSD interpretation:**
   - Did RMSD plateau or keep increasing?
   - If plateau: what's the value? (<1 Å = very stable, 2-3 Å = normal motion)
   - If increasing: is the protein unfolding?

3. **Flexibility patterns:**
   - Which regions are rigid (low RMSF)?
   - Which are flexible (high RMSF)?
   - Do flexible regions match loop regions in structure?
   - Do they align with AlphaFold confidence (pLDDT)?

4. **Validation conclusion:**
   - Is the AlphaFold structure stable?
   - Does MD confirm or challenge the prediction?
   - Any concerning observations?

### ✅ Quick Check 5: Result Interpretation

[[?]]
| You observe: RMSD increases linearly from 0 to 5 Å over 500 ps. This means:
| - [[ ]] The simulation failed
| - [[X]] The structure is unfolding or rearranging significantly
| - [[ ]] The thermostat is too weak
| - [[ ]] Everything is normal

[[?]]
| Loops in the protein usually have:
| - [[X]] High RMSF (flexible)
| - [[ ]] Low RMSF (rigid)
| - [[ ]] Zero RMSF (frozen)
| - [[ ]] Unpredictable RMSF

---

## Part 6️⃣: Biological Interpretation — The Bigger Picture

### From Numbers to Biology

Your analysis gives you **numbers**. Now extract **biology**:

---

### Case 1: Stable Structure, Low Flexibility

**What you observe:**
- RMSD < 1.5 Å
- Mean RMSF < 1 Å overall
- Low pLDDT regions are also flexible (high RMSF)

**Interpretation:**
- AlphaFold prediction is **validated by MD**
- Structure is thermodynamically stable
- Confidence scores align with actual dynamics
- ✅ **Trust this prediction for binding studies, design, etc.**

---

### Case 2: Stable Core, Flexible Termini

**What you observe:**
- RMSD 1-2 Å overall
- N-terminus/C-terminus: High RMSF (5-10 Å)
- Core structure (secondary structure): Low RMSF (< 1 Å)

**Interpretation:**
- **Common in proteins** — termini are often intrinsically disordered
- AlphaFold may or may not capture this flexibility
- MD reveals terminal dynamics invisible in crystal structures
- ✅ **Expected and biologically meaningful**

---

### Case 3: Large Rearrangement, Increasing RMSD

**What you observe:**
- RMSD increases from 1 Å to 4 Å over 500 ps
- No plateau
- Secondary structure elements shift position

**Possibilities:**
1. **Force field artifact** — GFN-FF might not like this protein
2. **Incorrect starting structure** — AlphaFold mispredicted or chose wrong conformation
3. **Intrinsic flexibility** — Protein is designed to be mobile (hinge protein, signaling)
4. **Sampling problem** — 500 ps too short to reach equilibrium

**Next steps:**
- Compare multiple force fields (test with UFF in Session 1 context)
- Check literature for known conformations
- Run longer simulation if possible
- ⚠️ **Be cautious with this prediction**

---

### Case 4: High Confidence But High Flexibility

**What you observe:**
- AlphaFold pLDDT > 80 (high confidence)
- MD shows high RMSF (>2 Å) in same region

**Possibilities:**
1. **Functionally important flexibility** — Site for binding, enzyme active site
2. **Discrepancy between training data and real dynamics** — AlphaFold confident but MD shows movement
3. **Allosteric transition** — Region moves during function

**Interpretation:**
- High confidence doesn't always mean rigid
- MD reveals functional dynamics hidden in static predictions
- ✅ **These regions warrant functional investigation**

---

### ✅ Quick Check 6: Biological Interpretation

[[?]]
| A region has high pLDDT (AlphaFold confident) AND high RMSF (flexible in MD):
| - [[ ]] One of the analyses must be wrong
| - [[X]] The region might be functionally important and dynamically active
| - [[ ]] The protein is misfolded
| - [[ ]] This never happens

---

## Part 7️⃣: Presentation Guidelines

Your task is to present your findings in a **10-15 minute talk**.

### Presentation Structure

#### Slide 1: Title & Motivation
- Protein name, source, biological function
- Why study this protein?
- What are you investigating?

**Example:** "Ubiquitin (76 AA) — a key regulatory protein. How stable is the predicted structure? How flexible are functionally important regions?"

---

#### Slide 2: Methods
- AlphaFold prediction (link, confidence overview)
- MD setup (force field, duration, temperature)
- Analysis methods (RMSD, RMSF, energy)

**Example:**
```
- AlphaFold2 prediction: pLDDT 85 avg
- MD: GFN-FF, 500 ps, 300 K, NVT ensemble
- Analysis: RMSD (C-alpha), RMSF (all atoms), Energy
```

---

#### Slide 3: Results — Structure
- AlphaFold structure (screenshot from PyMOL)
- Color by pLDDT (confidence regions)
- Highlight interesting regions

---

#### Slide 4: Results — Energy
- Plot: Total Energy vs Time
- Is it converged? What's the mean?
- Any anomalies?

---

#### Slide 5: Results — RMSD
- Plot: RMSD vs AlphaFold over time
- Final value? Plateau or linear increase?
- Interpretation: Stable? Unfolding? Normal?

---

#### Slide 6: Results — Flexibility
- Plot: RMSF per residue
- Highlight rigid regions (secondary structure)
- Highlight flexible regions (loops, termini)
- Compare to AlphaFold pLDDT — do they match?

---

#### Slide 7: Comparison
- Overlay or RMSD matrix: AlphaFold vs. MD final frame
- Structural comparison: How different are they?
- Quantitative: alignment RMSD, etc.

---

#### Slide 8: Interpretation
- Is AlphaFold prediction validated by MD?
- What does MD add beyond static prediction?
- Are there surprising findings?

**Example:**
```
✓ Structure stable (RMSD < 1.5 Å)
✓ Confidence correlates with rigidity
✓ MD reveals loop dynamics invisible in AF
→ Prediction trustworthy for design studies
```

---

#### Slide 9: Biological Significance
- What does this mean for the protein's function?
- Which regions are important for stability?
- Which might be involved in binding/catalysis?
- How could this inform experiments or design?

**Example:**
```
- Flexible hinge region between domains
  → Likely important for conformational change
- Rigid core with high confidence
  → Suitable for structure-based drug design
```

---

#### Slide 10: Conclusions & Limitations
- Summary: What did you learn?
- How do AF and MD complement each other?
- Limitations: What didn't this study address?
- Future directions: What would you do next?

---

### Presentation Tips

✅ **Do:**
- Show your plots clearly (large fonts, good contrast)
- Explain what numbers mean (don't just show raw data)
- Relate findings back to biology
- Be honest about limitations
- Prepare for questions about methodology

❌ **Don't:**
- Say "AlphaFold is better" or "MD is better"
- Show every command you typed (show results, not process)
- Go too deep into technical details
- Assume audience knows all methods

### ✅ Quick Check 7: Presentation

[[?]]
| Your presentation should emphasize:
| - [[X]] How AlphaFold and MD answer different questions and complement each other
| - [[ ]] Which method is superior
| - [[ ]] The technical details of each simulation
| - [[ ]] Only the final conclusions, no data

---

## Part 8️⃣: Combining Methods — Advanced Perspective

### When to Use Each Method

**Use AlphaFold when:**
- You have a sequence but no structure
- You need a quick structural hypothesis
- You want to compare multiple homologs
- Timescale: minutes to hours on server

**Use MD when:**
- You need dynamics and flexibility information
- You want to validate a static prediction
- You need thermal stability assessment
- You're investigating functional motion
- Timescale: hours to days on cluster

**Use Both when:**
- Structure + dynamics matter for your question
- You want to validate predictions
- You're doing protein design (structure) + binding studies (dynamics)
- You want to understand allosteric mechanisms

---

### The Workflow: A Best Practice

```
1. Sequence → Predict with AlphaFold
           ↓
2. Check confidence (pLDDT, PAE)
           ↓
3. Run MD simulation from AF structure
           ↓
4. Validate: Is structure stable?
           ↓
        ├─ YES: Use prediction confidently
        │        ↓
        │   Structure-based applications
        │   (design, docking, etc.)
        │
        └─ NO or FLEXIBLE:
                ↓
           Investigate further
           • Try different force fields
           • Run longer MD
           • Check literature
           • Consider biology
```

---

### Machine Learning + Classical: The Future

**Emerging approaches:**
- **AF2 + Gromacs:** Systematic validation
- **Refinement:** Use MD to refine AF predictions
- **Ensemble:** Run multiple AFs, validate with MD
- **Dynamics in AF3:** AlphaFold now includes limited dynamics (future sessions!)

**Key insight:** ML predictions + classical simulations = **more trustworthy science**

---

## Part 9️⃣: Advanced Analysis (Optional, for Interested Students)

### Cluster Analysis

Group similar conformations in your trajectory:

```bash
gmx cluster -f prod.xtc -s prod.tpr -o rmsd-clust.xvg
# Groups structures by RMSD similarity
```

**Interpretation:** Protein samples multiple conformations? Single state?

---

### Principal Component Analysis (PCA)

Find major motions:

```bash
gmx covar -f prod.xtc -s prod.tpr -o covar.xvg -v eigenvec.trr
gmx anaeig -v eigenvec.trr -f prod.xtc -o proj.xvg
```

**Interpretation:** What are the dominant protein motions?

---

### Hydrogen Bonding Analysis

```bash
gmx hbond -f prod.xtc -s prod.tpr -o hbond.xvg
```

**Interpretation:** Intramolecular H-bonds stable?

---

### Secondary Structure Evolution

```bash
gmx dssp -f prod.xtc -s prod.tpr -o ss.xvg
```

**Interpretation:** Does secondary structure persist?

---

### ✅ Quick Check 8: Advanced Analysis

[[?]]
| Cluster analysis on a trajectory might reveal:
| - [[X]] Whether the protein samples one conformation or multiple states
| - [[ ]] If the simulation has bugs
| - [[ ]] Which atoms are most important
| - [[ ]] The force field accuracy

---

## Part 🔟: Summary & Key Takeaways

### What You've Learned

✅ **Conceptually:**
- AlphaFold and MD are complementary, not competitive
- AF predicts structure; MD simulates dynamics
- Each method has strengths and limitations
- Combined use gives the most insight

✅ **Practically:**
- How to extract energy, RMSD, RMSF from MD
- How to interpret these quantities biologically
- How to validate AF predictions with MD
- How to communicate findings professionally

✅ **Scientifically:**
- What questions each method answers
- How to use results for biological insight
- How ML predictions and classical simulations combine
- The importance of validation and uncertainty

---

### Key Principles

1. **Validate predictions** — Don't trust a single method
2. **Question outliers** — If MD disagrees with AF, investigate
3. **Interpret biologically** — Numbers must mean something for the organism
4. **Acknowledge limitations** — Every method has assumptions
5. **Combine approaches** — Complementarity is powerful

---

### The Bigger Picture: Molecular Modelling in Biology

```
Sequence
   ↓ [AlphaFold]
Structure (prediction)
   ↓ [MD]
Dynamics (simulation)
   ↓ [Analysis]
Mechanism (hypothesis)
   ↓ [Experiment]
Validation (truth)
```

This cycle — prediction → simulation → experiment — is **modern structural biology**.

---

### ✅ Final Quiz: Integration

[[?]]
| The primary value of combining AlphaFold and MD is:
| - [[X]] Validation: Structure + Dynamics → Robust understanding
| - [[ ]] Finding which method is "best"
| - [[ ]] Faster computation
| - [[ ]] Better visualization

[[?]]
| If MD shows large RMSD from AF structure, it could mean:
| - [[X]] Structure unstable, or AF mispredicted, or force field artifact
| - [[ ]] MD definitely failed
| - [[ ]] AF is always wrong for this protein
| - [[ ]] Nothing concerning

[[?]]
| Per-residue flexibility (RMSF) is most directly useful for:
| - [[X]] Identifying flexible regions that might be functionally important
| - [[ ]] Calculating protein size
| - [[ ]] Determining sequence
| - [[ ]] Nothing; it's just a number

[[?]]
| A region with low AlphaFold confidence AND high MD RMSF:
| - [[X]] Warrants investigation as potentially functionally relevant
| - [[ ]] Is certainly misfolded
| - [[ ]] Should be ignored
| - [[ ]] Means the protein is unstable

---

## 🎓 Your Final Deliverable: The Presentation

### What to Submit

**Online Presentation (10-15 minutes):**
1. Slides (PowerPoint, PDF, or Google Slides)
2. All figures (energy, RMSD, RMSF plots)
3. Structural comparison (PyMOL/Avogadro screenshot)
4. Brief written summary (1-2 pages)

**Grading Criteria:**
- ✓ Correct methodology (you understand the tools)
- ✓ Accurate interpretation (numbers → biology)
- ✓ Clear communication (understandable to peers)
- ✓ Honest assessment (acknowledged limitations)
- ✓ Biological insight (why does this matter?)

---

## 📚 Resources for This Session

### Analysis Tools

| Tool | Purpose | Command |
|------|---------|---------|
| `gmx energy` | Extract energy data | `gmx energy -f prod.edr -o E.xvg` |
| `gmx rms` | Calculate RMSD | `gmx rms -f prod.xtc -s prod.tpr` |
| `gmx rmsf` | Per-residue flexibility | `gmx rmsf -f prod.xtc -s prod.tpr` |
| `gmx cluster` | Cluster analysis | `gmx cluster -f prod.xtc -s prod.tpr` |
| Gnuplot | Data plotting | `gnuplot script.gnu` |
| PyMOL | Structure visualization | `pymol structure.pdb` |

### References

- **Lemkul (2024):** Analysis sections of the Gromacs tutorial
- **Jumper et al. (2021):** "Highly accurate protein structure prediction with AlphaFold2" *Nature*
- **Allen & Tildesley (1987):** *Computer Simulation of Liquids* (MD theory)

---

## Questions & Discussion

This is a **live session** — bring questions!

- Did your MD agree or disagree with AF? Why?
- What surprised you about your results?
- How would you improve the analysis?
- What other proteins would you like to study?
- How might this inform experimental design?

---

*Session 3 — AlphaFold and Molecular Dynamics: Complementary Tools*  
*Last updated: October 27, 2025*  
*Course: Molecular Modelling and Quantum Chemistry (Master)*
