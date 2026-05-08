<!--
author: Molecular Modelling Course Team
language: en
narrator: US English Female
version: 1.0

[Open this course in LiaScript](https://liascript.github.io/course/?https://raw.githubusercontent.com/conradhuebler/TeacherTwinMolecularModelling/main/BashRPG.md)
-->

> ![Co-funded by the European Union](https://github.com/conradhuebler/TeacherTwinMolecularModelling/raw/main/EU.jpg)
> *Funded by TUBAFdigital and the European Union (Erasmus+).*

# Pre-Course: Linux Console Basics
## Your Gateway to Molecular Modelling

> **Welcome, Future Molecular Modellers!** 
>
> Before we dive into proteins and simulations, you need to befriend the **terminal**. 
> This pre-course takes ~2 hours and teaches you the essentials to navigate your VirtualBox without fear.
>
> Think of your console as a **quest interface** — we'll unlock it step by step. 🗺️

---

## 🎯 Learning Objectives

By the end of this pre-course, you will:

- ✅ Navigate the filesystem confidently (`cd`, `pwd`, `ls`)
- ✅ Create, copy, and manage files and folders
- ✅ Edit text files with `nano`
- ✅ Use pipes and wildcards to combine commands
- ✅ Search files with `grep` and find
- ✅ Understand file paths (absolute vs. relative)
- ✅ Be comfortable running commands on your VirtualBox

**Time investment:** ~2 hours (can be split across days)

---

## Part 1️⃣: Welcome to the Terminal

### What is the Terminal/Console?

The **terminal** (also called *console* or *shell*) is a text-based interface to your computer. Instead of clicking buttons, you type commands. 

**Why use it?**

- ✨ Powerful automation (run 100 commands at once)
- ✨ Precise control (no ambiguity, always exact)
- ✨ Scientific computing standard (Gromacs, Avogadro, your tool all use it)
- ✨ Reproducibility (commands can be written down and re-run)

### First Steps: Open Your Terminal

On your **OpenSUSE + KDE** (VirtualBox):

1. Click the **Application Menu** (bottom left)
2. Search for **"Konsole"** or **"Terminal"**
3. Click to open

You should see something like:

```
user@hostname:~$ 
```

**Congratulations!** You're in. 🎉

---

### The Terminal Anatomy

```
user@hostname:~/some/path$ command -option argument
│     │         │          │       │       │
└─────┴─────────┴──────────┤       │       │
   User, computer, where   │       │       │
   you are right now       │       │       │
                  ┌────────┘       │       │
                  │  The command   │       │
                  │  (what to do)  │       │
                  └────────────────┘       │
                    Flags/options  ────┐   │
                    (modify behavior)  │   │
                    Arguments ─────────┴───┘
                    (input data)
```

### ✅ Quick Check 1: Terminal Anatomy

In the command `ls -la /home`, what is "-la"?

- [[X]] Flags/options that modify the command
- [[ ]] The file to list
- [[ ]] The output
- [[ ]] A path

---

## Part 2️⃣: Navigation — Finding Your Way

### Command: `pwd` — Where Am I?

**`pwd`** = "**Print Working Directory**"

It tells you exactly where you are in the filesystem.

**Try it:**

```bash
pwd
```

Expected output:

```
/home/user
```

The `/` at the start means **root** of the filesystem. Everything starts there.

---

### The Filesystem Tree

Your computer's files are organized like this:

```
/                          ← ROOT (everything starts here)
├── home/                  ← User home folders
│   └── user/              ← YOUR home folder (~)
│       ├── Documents/
│       ├── Downloads/
│       └── Desktop/
├── usr/                   ← System programs & files
├── etc/                   ← Configuration files
├── tmp/                   ← Temporary files
└── opt/                   ← Optional software (Gromacs, etc.)
```

When you open the terminal, you **start in** `/home/user` (aka `~` = "home shortcut")

### ✅ Quick Check 2: Filesystem Structure

What does the `~` symbol represent?

- [[ ]] The root of the filesystem
- [[X]] Your home directory
- [[ ]] The parent folder
- [[ ]] A temporary folder

An absolute path always starts with:

- [[ ]] `~`
- [[X]] `/`
- [[ ]] `..`
- [[ ]] A letter (like `C:` on Windows)

---

### Command: `ls` — List Files

**`ls`** = "**LiSt**"

Shows what's in your current folder.

**Try it:**

```bash
ls
```

Expected output (depends on your system):

```
Desktop    Documents    Downloads    Pictures    Videos
```

**Make it fancier:**

```bash
ls -la
```

- `-l` = **Long format** (details: permissions, size, date)
- `-a` = **All** (including hidden files, those start with `.`)

Example output:

```
drwxr-xr-x 10 user user 4096 Oct 27 10:30 .
drwxr-xr-x  3 root root 4096 Oct 20 15:00 ..
drwxr-xr-x  2 user user 4096 Sep 15 12:00 Desktop
drwxr-xr-x  2 user user 4096 Sep 15 12:00 Documents
-rw-r--r--  1 user user  245 Oct 26 09:00 .bashrc
```

**Read it:**

- `d` = directory (folder), `-` = regular file
- `rwx` = read, write, execute permissions
- `user user` = owner, group
- `4096` = file size in bytes
- `Oct 27 10:30` = last modified

---

### Quiz: Understanding `ls -la`

What does the `d` at the start of `drwxr-xr-x` mean?

- [[X]] It's a directory (folder)
- [[ ]] It's a device file
- [[ ]] It's a deleted file

---

### Command: `cd` — Change Directory

**`cd`** = "**Change Directory**"

Move to a different folder.

**Try it:**

```bash
cd Documents
pwd
```

Now you should see:

```
/home/user/Documents
```

---

### Path Types: Absolute vs. Relative

**Absolute path:** Starts with `/`, always works from anywhere

```bash
cd /home/user/Documents    # Works from everywhere
```

**Relative path:** Doesn't start with `/`, relative to where you are NOW

```bash
cd Documents               # Only works if Documents is in current folder
```

**Special shortcuts:**

| Command | Meaning |
|---------|---------|
| `cd ~` | Go to home folder |
| `cd` | Go to home folder (same as above) |
| `cd ..` | Go up one level (parent folder) |
| `cd -` | Go back to previous folder |
| `cd .` | Current folder (useful later) |

**Try the journey:**

```bash
pwd                   # See where you are
cd ..                 # Go up one level
pwd                   # Check location
cd -                  # Go back
pwd                   # Should be in Documents again
```

---

### Challenge 1️⃣: Navigate to `/usr/bin` and Back

**Your task:**

1. Navigate to `/usr/bin` (absolute path)
2. List the files
3. Go back home
4. Print your working directory to confirm

**Solution (reveal when ready):**

```bash
cd /usr/bin
ls
cd ~
pwd
# Should print: /home/user
```

---

## Part 3️⃣: Working with Files & Folders

### Command: `mkdir` — Make Directory

Create new folders.

```bash
mkdir my_project
ls
```

Output:

```
my_project
```

**Create nested folders:**

```bash
mkdir -p my_project/data/raw
# -p = "parents" (creates all needed folders)
```

---

### Command: `touch` — Create Empty File

Create an empty file (or update timestamp of existing file).

```bash
cd my_project
touch README.txt
ls
```

Output:

```
README.txt
```

---

### Command: `cp` — Copy Files

Copy files or folders.

**Copy a file:**

```bash
cp README.txt README_backup.txt
ls
```

Output:

```
README.txt    README_backup.txt
```

**Copy a folder (recursively):**

```bash
cp -r data data_copy
# -r = recursive (includes all contents)
```

---

### Command: `mv` — Move or Rename

Move/rename files and folders.

**Rename:**

```bash
mv README_backup.txt README.bak
ls
```

**Move to subfolder:**

```bash
mv README.bak data/
ls data/
```

Output:

```
README.bak
```

---

### Command: `rm` — Remove Files

**Delete a file:**

```bash
rm README.bak
```

⚠️ **CAREFUL:** No trash/undo in terminal! Use with caution.

**Delete a folder (recursively):**

```bash
rm -r data_copy
# -r = recursive
```

---

### Quiz: File Commands

What command removes a file permanently?

- [[ ]] `delete file.txt`
- [[X]] `rm file.txt`
- [[ ]] `erase file.txt`

---

### Challenge 2️⃣: Build a Project Structure

**Your task:**

1. Create a folder `protein_sim` in your home
2. Inside, create subfolders: `inputs`, `outputs`, `analysis`
3. Create a file `log.txt` in the main folder
4. Copy `log.txt` to `outputs/log_backup.txt`
5. List everything to confirm

**Solution:**

```bash
cd ~
mkdir -p protein_sim/{inputs,outputs,analysis}
cd protein_sim
touch log.txt
cp log.txt outputs/log_backup.txt
ls -la
ls -la outputs/
```

---

## Part 4️⃣: Editing Text Files

### Command: `nano` — Simple Text Editor

`nano` is a beginner-friendly text editor. Open any file:

```bash
nano log.txt
```

You now see:

```
  GNU nano 7.2                                          New Buffer
  
  
  
  
^G Get Help     ^O Write Out    ^W Where Is     ^K Cut Text     ^J Justify
^X Exit         ^R Read File    ^\ Replace      ^U Paste Text   ^T To Spell
```

---

### Nano Basics

| Key | Action |
|-----|--------|
| `Ctrl+X` | Exit (prompted to save) |
| `Ctrl+O` | Save (Write Out) |
| `Ctrl+W` | Find text |
| `Ctrl+K` | Cut line |
| `Ctrl+U` | Paste |
| `Ctrl+G` | Help |

**Try it:**

1. Type some text:

```
This is my first simulation project!
```

2. Press `Ctrl+O` → Press `Enter` to save
3. Press `Ctrl+X` to exit
4. Open it again: `nano log.txt`
5. You should see your text!

---

### View File Content (without editing)

**`cat`** — Print entire file to screen

```bash
cat log.txt
```

Output:

```
This is my first simulation project!
```

**`less`** — View large files (page by page)

```bash
less log.txt
# Press 'q' to quit
```

---

## Part 5️⃣: Powerful Combinations - Pipes & Wildcards

### Pipes: `|` — Chain Commands

A **pipe** `|` sends output of one command as input to the next.

**Syntax:**

```bash
command1 | command2
```

**Example — Find long output easier:**

```bash
ls -la | grep "txt"
# Shows only lines containing "txt"
```

**Example — Count files:**

```bash
ls | wc -l
# wc = "word count" (-l = lines, i.e., count them)
```

**Real use case from molecular modelling:**

```bash
gmx energy -f ener.edr | tail -20
# Shows last 20 lines of energy output
```

---

### Wildcards: `*` and `?` — Pattern Matching

**`*`** = any characters (0 or more)

```bash
ls *.txt          # All files ending in .txt
ls protein_*      # All files starting with protein_
ls *.pdb *.gro    # All .pdb or .gro files
```

**`?`** = exactly one character

```bash
ls file?.txt      # file1.txt, fileA.txt, but not file12.txt
```

**Real use case:**

```bash
cp *.xyz /backup/     # Copy all .xyz files to backup
rm energy_*.edr       # Remove all energy files
```

---

### Challenge 3️⃣: Create and Organize Files

**Your task:**

1. Go to `protein_sim/inputs`
2. Create 3 files: `mol1.xyz`, `mol2.xyz`, `opt.inp`
3. List only `.xyz` files using wildcard
4. Count how many files are in this folder

**Solution:**

```bash
cd ~/protein_sim/inputs
touch mol1.xyz mol2.xyz opt.inp
ls *.xyz
ls | wc -l      # Should be 3
```

---

## Part 6️⃣: Searching & Finding

### Command: `grep` — Find Text in Files

**`grep`** searches for text patterns in files.

**Syntax:**

```bash
grep "pattern" filename
```

**Example:**

Create a file with content:

```bash
cat > simulation.log << EOF
Step 1000, Energy = -245.5 kcal/mol
Step 2000, Energy = -250.3 kcal/mol
Step 3000, Temperature = 298.15 K
EOF
```

Now search:

```bash
grep "Energy" simulation.log
```

Output:

```
Step 1000, Energy = -245.5 kcal/mol
Step 2000, Energy = -250.3 kcal/mol
```

**Case-insensitive search:**

```bash
grep -i "energy" simulation.log
# Same result, even if you wrote "ENERGY"
```

**Count matches:**

```bash
grep -c "Energy" simulation.log
# Output: 2
```

---

### Real Use Case: Extract Energy Values

After a Gromacs simulation, you have many lines. Extract only energy:

```bash
grep "Potential" energy.log | tail -1
# Gets last energy value
```

---

### Command: `find` — Locate Files

Find files by name, size, date, etc.

**Basic syntax:**

```bash
find /path -name "pattern"
```

**Examples:**

```bash
find . -name "*.xyz"              # All .xyz in current folder and subfolders
find ~ -name "protein*"           # All files starting with "protein"
find /tmp -name "*.tmp"           # All .tmp files in /tmp
```

---

### Quiz: Pipes and Wildcards

What does this command do?

```bash
ls *.pdb | wc -l
```

- [[X]] Counts how many `.pdb` files are in the current folder
- [[ ]] Lists all `.pdb` files and their full paths
- [[ ]] Removes all `.pdb` files

---

## Part 7️⃣: The Big Quest — Let's Practice!

Now let's combine everything in a realistic scenario:

**Scenario:** You're setting up for a protein simulation campaign. You need to:

1. Create a project structure
2. Add some mock input files
3. Find and organize them
4. Create a summary file

---

### Quest Steps:

**Step 1: Create structure**

```bash
cd ~
mkdir -p MD_campaign/{proteins,simulations,analysis}
cd MD_campaign
```

**Step 2: Add some protein files**

```bash
cd proteins
touch protein_1A.pdb protein_2B.pdb protein_3C.pdb
touch ligand_A.mol2 ligand_B.mol2
touch README.txt
cd ..
ls proteins/
```

**Step 3: Create simulation setup**

```bash
cd simulations
touch run1.mdp run2.mdp run3.mdp
touch equilibration.gro production.gro
cd ..
```

**Step 4: Find all `.pdb` files**

```bash
find . -name "*.pdb"
```

Expected output:

```
./proteins/protein_1A.pdb
./proteins/protein_2B.pdb
./proteins/protein_3C.pdb
```

**Step 5: Count them**

```bash
find . -name "*.pdb" | wc -l
# Output: 3
```

**Step 6: Create a summary**

```bash
nano SIMULATION_LOG.txt
```

Type:

```
My First MD Campaign
=====================
Date: Oct 27, 2025
Status: Setup Phase

Files prepared:
- 3 protein structures
- 2 ligands
- 3 simulation input files
```

Save (`Ctrl+O` → Enter → `Ctrl+X`)

**Step 7: Verify**

```bash
ls -la
cat SIMULATION_LOG.txt
```

---

### Final Challenge 🏆: Your Own Project

**Do this:**

1. Create a folder called `my_molecule` in your home
2. Create subfolders: `structures`, `calculations`, `results`
3. Create 3 files: `molecule.xyz`, `settings.txt`, `notes.md`
4. Edit `settings.txt` with your favorite method parameters
5. List everything with `ls -la`
6. Count all files: `find my_molecule | wc -l`

**Show your work:**

```bash
cd ~/my_molecule
ls -laR    # R = recursive, shows everything
```

✅ **You just completed a realistic workflow!**

---

## Part 8️⃣: Common Traps & Tips

### Don't Do This ❌

```bash
rm -rf /              # This deletes EVERYTHING (don't do!)
rm -r *               # In wrong folder = disaster
```

### Do This Instead ✅

```bash
rm -i filename        # -i = ask for confirmation
ls -la first          # Check before deleting
rm one_specific_file  # Be specific
```

---

### Useful Shortcuts

| Shortcut | Action |
|----------|--------|
| `Tab` | Auto-complete (press twice for options) |
| `↑` / `↓` | Navigate command history |
| `Ctrl+C` | Stop current command |
| `Ctrl+L` | Clear screen |
| `Ctrl+A` | Go to start of line |
| `Ctrl+E` | Go to end of line |

**Tab completion example:**

```bash
cd Doc[TAB]           # Auto-completes to "Documents"
ls *.py[TAB][TAB]     # Shows all matching files
```

### ✅ Quick Check 11: Shortcuts and Safety

How do you stop a running command?

- [[ ]] Press Enter
- [[X]] Press Ctrl+C
- [[ ]] Type "stop"
- [[ ]] Wait forever

What does the `Tab` key do in the terminal?

- [[ ]] Creates indentation
- [[X]] Auto-completes file/folder names
- [[ ]] Switches between open windows
- [[ ]] Does nothing

---

### The Holy Trinity of Help

**1. Manual pages (detailed):**

```bash
man ls              # Read full documentation (press 'q' to exit)
man grep
```

**2. Help flag (quick):**

```bash
ls --help           # Short summary
grep --help
```

**3. Search online:**

```bash
# Google: "how to grep in linux"
# Stack Overflow is your friend
```

---

## Part 9️⃣: Glossary & Quick Reference

| Term | Meaning |
|------|---------|
| **~** | Home folder shortcut |
| **/** | Root of filesystem / path separator |
| **.** | Current folder |
| **..** | Parent folder (one level up) |
| **-** | Previous folder |
| **Absolute path** | Starts with `/`, works from anywhere |
| **Relative path** | No `/`, relative to current location |
| **Flag/Option** | Modifies command behavior (starts with `-`) |
| **Pipe** | `\|` chains commands together |
| **Wildcard** | `*` (any chars) or `?` (one char) |

---

## Commands At A Glance

```bash
# Navigation
pwd                 # Print working directory
cd path             # Change directory
ls                  # List files
ls -la              # List with details + hidden

# Folders
mkdir foldername    # Create folder
mkdir -p a/b/c      # Create nested folders
cd ..               # Go up one level

# Files
touch filename      # Create empty file
cp file1 file2      # Copy file
mv file1 file2      # Move or rename
rm file             # Delete file
rm -r folder        # Delete folder

# Viewing
cat file            # Show file contents
less file           # View large file (press q to exit)
nano file           # Edit file

# Searching
grep "text" file    # Search in file
find /path -name "*pattern*"  # Find files

# Pipes & Wildcards
command1 | command2 # Pipe (chain commands)
ls *.txt            # Wildcard: all .txt files

# Counting
wc -l file          # Count lines
ls | wc -l          # Count files in folder
```

---

## 🎓 You Did It!

You've completed the **Linux Basics Pre-Course**. 

### What You Learned:

- ✅ Navigate the filesystem confidently  
- ✅ Create, copy, move, delete files  
- ✅ Edit text with `nano`  
- ✅ Use pipes and wildcards  
- ✅ Search and find files  
- ✅ Understand paths and structure  

### What's Next:

- **Session 1:** Geometry optimization with [TOOL_NAME]
- **Session 2:** Molecular Dynamics with Gromacs
- **Session 3:** AlphaFold vs. MD comparison

---

## 🏆 Final Mastery Quiz

Time to test everything you've learned!

What does `pwd` stand for?

- [[ ]] Print Working Data
- [[X]] Print Working Directory
- [[ ]] Package Working Directory
- [[ ]] Print Web Directory

---

How would you create a folder `analysis` inside an existing folder `project`?

- [[ ]] `mkdir project analysis`
- [[X]] `mkdir project/analysis` or `cd project` then `mkdir analysis`
- [[ ]] `create project/analysis`
- [[ ]] `touch project/analysis`

---

What does `ls -la` show that `ls` doesn't?

- [[X]] Hidden files, permissions, owner, size, and modification date
- [[ ]] Only folder names
- [[ ]] Deleted files
- [[ ]] File contents

---

If you're in `/home/user` and type `cd ./Desktop`, where do you go?

- [[ ]] To the root directory
- [[X]] To `/home/user/Desktop`
- [[ ]] To `./Desktop` at the root level
- [[ ]] The command is invalid

---

How many files match the pattern `data_?.csv`?

- [[ ]] Any number
- [[X]] Only files like data\_1.csv, data\_A.csv, data\_X.csv (exactly one character between underscore and dot)
- [[ ]] Files like data\_1.csv, data\_12.csv, data\_ABC.csv
- [[ ]] Only files named exactly `data_?.csv`

---

What does `grep -i "energy" log.txt` do?

- [[ ]] Counts lines in the file
- [[X]] Searches for "energy" in log.txt, ignoring uppercase/lowercase
- [[ ]] Deletes lines containing "energy"
- [[ ]] Copies lines with "energy" to a new file

---

The command `cp -r folder1 folder2` does what?

- [[ ]] Copies only files from folder1 to folder2
- [[X]] Copies the entire folder1 and all its contents as folder2
- [[ ]] Moves folder1 to folder2
- [[ ]] Creates a link instead of copying

---

What's the safest way to delete a file you're unsure about?

- [[ ]] `rm filename`
- [[X]] `rm -i filename` (asks for confirmation)
- [[ ]] `mv filename /trash` (no such directory)
- [[ ]] There's no safe way in the terminal

---

How would you list only `.pdb` files and count them?

- [[ ]] `ls *.pdb count`
- [[X]] `ls *.pdb | wc -l`
- [[ ]] `count *.pdb`
- [[ ]] `ls -count *.pdb`

---

What does `find . -name "*.xyz"` do?

- [[ ]] Lists all .xyz files only in the current folder
- [[X]] Finds all .xyz files in the current folder and all subfolders
- [[ ]] Creates new .xyz files
- [[ ]] Searches the entire filesystem for .xyz files

---

## 📚 Additional Resources

- **GNU Bash Manual:** https://www.gnu.org/software/bash/manual/
- **Linux Command Line Tutorial:** https://linuxcommand.org/
- **Nano Tutorial:** https://www.nano-editor.org/

---

## Questions?

If something is unclear, **ask your instructors in Session 1!**

The terminal might seem intimidating, but:
- 🎯 Everyone starts as a beginner
- 🎯 Practice makes perfect
- 🎯 You'll be comfortable in a few days

**Good luck with your molecular modelling journey!** 🧬🚀

---

## Funding

The development of this course material is funded by:

- TUBAFdigital
- European Union (Erasmus+ National Agency for Higher Education)

Funded by the European Union. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or Erasmus+ National Agency for Higher Education (German Academic Exchange Service). Neither the European Union nor the granting authority can be held responsible for them.

![Co-funded by the European Union](https://github.com/conradhuebler/TeacherTwinMolecularModelling/raw/main/EU.jpg)

---

*Last updated: October 27, 2025*
*Course: Molecular Modelling and Quantum Chemistry (Master)*
