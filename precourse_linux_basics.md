<!--
author: Molecular Modelling Course Team
language: en
narrator: US English Female
version: 1.1

This is a LiaScript course. View it at:
https://liascript.github.io/course/?<URL>
-->

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
   User, computer, where    │       │       │
   you are right now         │       │       │
                    ┌────────┘       │       │
                    │  The command   │       │
                    │  (what to do)   │       │
                    └────────────────┘       │
                       Flags/options ────┐   │
                       (modify behavior)  │   │
                       Arguments ─────────┴───┘
                       (input data)
```

### ✅ Quick Check 1: Terminal Anatomy

[[?]]
| In the command `ls -la /home`, what is "-la"?
| - [[X]] Flags/options that modify the command
| - [[ ]] The file to list
| - [[ ]] The output
| - [[ ]] A path

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

[[?]]
| What does the `~` symbol represent?
| - [[ ]] The root of the filesystem
| - [[X]] Your home directory
| - [[ ]] The parent folder
| - [[ ]] A temporary folder

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

### ✅ Quick Check 3: Understanding `ls -la`

[[?]]
| What does the `d` at the start of `drwxr-xr-x` mean?
| - [[X]] It's a directory (folder)
| - [[ ]] It's a device file
| - [[ ]] It's a deleted file
| - [[ ]] It indicates the owner

[[?]]
| Which flag shows hidden files (starting with a dot)?
| - [[ ]] `-l`
| - [[X]] `-a`
| - [[ ]] `-f`
| - [[ ]] `-h`

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

### ✅ Quick Check 4: Paths and Navigation

[[?]]
| What does `cd ~` do?
| - [[X]] Takes you to your home folder
| - [[ ]] Takes you to the root folder
| - [[ ]] Goes up one level
| - [[ ]] Shows your current directory

[[?]]
| An absolute path always starts with:
| - [[ ]] `~`
| - [[X]] `/`
| - [[ ]] `..`
| - [[ ]] A letter (like `C:` on Windows)

[[?]]
| If you're in `/home/user/Documents` and type `cd ..`, where do you end up?
| - [[ ]] `/home/user/Documents` (no change)
| - [[X]] `/home/user`
| - [[ ]] `/home`
| - [[ ]] `/`

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

### ✅ Quick Check 5: File Operations

[[?]]
| What flag should you use to copy a folder with all its contents?
| - [[X]] `-r` (recursive)
| - [[ ]] `-f` (force)
| - [[ ]] `-a` (all)
| - [[ ]] `-c` (copy)

[[?]]
| How can you rename a file with the `mv` command?
| - [[ ]] `mv rename file.txt newname.txt`
| - [[X]] `mv file.txt newname.txt`
| - [[ ]] `ren file.txt newname.txt`
| - [[ ]] You need a special command for renaming

[[?]]
| What happens if you run `rm file.txt`?
| - [[ ]] It asks for confirmation first
| - [[X]] It deletes the file permanently (no undo!)
| - [[ ]] It moves it to trash
| - [[ ]] It makes a backup first

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

### ✅ Quick Check 6: Nano Editor

[[?]]
| How do you save a file and exit nano?
| - [[ ]] `Ctrl+S` then `Ctrl+Q`
| - [[X]] `Ctrl+X` (then save when prompted)
| - [[ ]] Just close the terminal
| - [[ ]] `Ctrl+Z`

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

### ✅ Quick Check 7: Viewing Files

[[?]]
| What's the difference between `cat` and `less`?
| - [[X]] `cat` shows everything at once; `less` shows it page-by-page
| - [[ ]] `less` is faster than `cat`
| - [[ ]] `cat` creates files; `less` views them
| - [[ ]] They're identical

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

### ✅ Quick Check 8: Pipes and Wildcards

[[?]]
| What does `ls *.pdb | wc -l` do?
| - [[X]] Counts how many `.pdb` files are in the current folder
| - [[ ]] Lists all `.pdb` files and their full paths
| - [[ ]] Removes all `.pdb` files
| - [[ ]] Shows the first `.pdb` file

[[?]]
| What does the `*` wildcard match?
| - [[ ]] Exactly one character
| - [[X]] Any number of characters (0 or more)
| - [[ ]] Only file extensions
| - [[ ]] Hidden files

[[?]]
| In the command `rm mol_*.xyz`, the wildcard matches:
| - [[ ]] Only one file
| - [[X]] All files like mol_1.xyz, mol_A.xyz, mol_test.xyz, etc.
| - [[ ]] All files in the directory
| - [[ ]] Nothing (invalid syntax)

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

### ✅ Quick Check 9: Grep Command

[[?]]
| What does `grep -c "pattern" file.txt` do?
| - [[ ]] Deletes lines matching the pattern
| - [[X]] Counts how many lines contain the pattern
| - [[ ]] Shows the first line with the pattern
| - [[ ]] Copies lines with the pattern to a new file

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

### ✅ Quick Check 10: Complex Commands

[[?]]
| What command would you use to find all `.gro` files in your home directory and subdirectories?
| - [[X]] `find ~ -name "*.gro"`
| - [[ ]] `ls ~ -name "*.gro"`
| - [[ ]] `grep -r "*.gro" ~`
| - [[ ]] `find *.gro`

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

[[?]]
| How do you stop a running command?
| - [[ ]] Press Enter
| - [[X]] Press Ctrl+C
| - [[ ]] Type "stop"
| - [[ ]] Wait forever

[[?]]
| What does the `Tab` key do in the terminal?
| - [[ ]] Creates indentation
| - [[X]] Auto-completes file/folder names
| - [[ ]] Switches between open windows
| - [[ ]] Does nothing

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

## 🏆 Final Mastery Quiz

Time to test everything you've learned!

[[?]]
| What does `pwd` stand for?
| - [[ ]] Print Working Data
| - [[X]] Print Working Directory
| - [[ ]] Package Working Directory
| - [[ ]] Print Web Directory

[[?]]
| How would you create a folder `analysis` inside an existing folder `project`?
| - [[ ]] `mkdir project analysis`
| - [[X]] `mkdir project/analysis` or `cd project` then `mkdir analysis`
| - [[ ]] `create project/analysis`
| - [[ ]] `touch project/analysis`

[[?]]
| What does `ls -la` show that `ls` doesn't?
| - [[X]] Hidden files, permissions, owner, size, and modification date
| - [[ ]] Only folder names
| - [[ ]] Deleted files
| - [[ ]] File contents

[[?]]
| If you're in `/home/user` and type `cd ./Desktop`, where do you go?
| - [[ ]] To the root directory
| - [[X]] To `/home/user/Desktop`
| - [[ ]] To `./Desktop` at the root level
| - [[ ]] The command is invalid

[[?]]
| How many files match the pattern `data_?.csv`?
| - [[ ]] Any number
| - [[X]] Only files like data_1.csv, data_A.csv, data_X.csv (exactly one character between underscore and dot)
| - [[ ]] Files like data_1.csv, data_12.csv, data_ABC.csv
| - [[ ]] Only files named exactly `data_?.csv`

[[?]]
| What does `grep -i "energy" log.txt` do?
| - [[ ]] Counts lines in the file
| - [[X]] Searches for "energy" in log.txt, ignoring uppercase/lowercase
| - [[ ]] Deletes lines containing "energy"
| - [[ ]] Copies lines with "energy" to a new file

[[?]]
| The command `cp -r folder1 folder2` does what?
| - [[ ]] Copies only files from folder1 to folder2
| - [[X]] Copies the entire folder1 and all its contents as folder2
| - [[ ]] Moves folder1 to folder2
| - [[ ]] Creates a link instead of copying

[[?]]
| What's the safest way to delete a file you're unsure about?
| - [[ ]] `rm filename`
| - [[X]] `rm -i filename` (asks for confirmation)
| - [[ ]] `mv filename /trash` (no such directory)
| - [[ ]] There's no safe way in the terminal

[[?]]
| How would you list only `.pdb` files and count them?
| - [[ ]] `ls *.pdb count`
| - [[X]] `ls *.pdb | wc -l`
| - [[ ]] `count *.pdb`
| - [[ ]] `ls -count *.pdb`

[[?]]
| What does `find . -name "*.xyz"` do?
| - [[ ]] Lists all .xyz files only in the current folder
| - [[X]] Finds all .xyz files in the current folder and all subfolders
| - [[ ]] Creates new .xyz files
| - [[ ]] Searches the entire filesystem for .xyz files

---

## 🎓 You Did It!

You've completed the **Linux Basics Pre-Course**. 

### What You Learned:

✅ Navigate the filesystem confidently  
✅ Create, copy, move, delete files  
✅ Edit text with `nano`  
✅ Use pipes and wildcards  
✅ Search and find files  
✅ Understand paths and structure  
✅ **11 Quick Checks** throughout the course  
✅ **10 Mastery Quiz Questions** at the end  

### What's Next:

- **Session 1:** Geometry optimization with [TOOL_NAME]
- **Session 2:** Molecular Dynamics with Gromacs
- **Session 3:** AlphaFold vs. MD comparison

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

*Last updated: October 27, 2025*
*Course: Molecular Modelling and Quantum Chemistry (Master)*
