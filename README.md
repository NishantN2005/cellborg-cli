
# Cellborg CLI

This tool provides a small, easy-to-use dashboard for managing single-cell
projects on your computer. It is intended for researchers who want a simple
way to organize data folders and run basic QC steps without writing code.

If this tool helps your work, please star the repository — it really helps!


## Get the code — step-by-step
1. Click the green "Code" button and copy the HTTPS URL (it looks like
   `https://github.com/NishantN2005/cellborg-cli.git`).

   [Screenshot placeholder: GitHub Code button]

2. Open the Terminal application on your computer and run these commands:

```bash
git clone https://github.com/NishantN2005/cellborg-cli.git
cd cellborg-cli
```

If you prefer a graphical option, you can select "Open with GitHub Desktop"
from the same menu and follow the prompts in GitHub Desktop.

   [Screenshot placeholder: Terminal showing git clone]

## Easy setup (three steps)

Follow these three steps to prepare and run the dashboard locally. After each
subsection in Step 3 there is a placeholder where you can add a screenshot.

Step 1 — Prepare Python

1. Install Python 3.10 or newer if you do not already have it. On macOS you
   can download the installer from python.org or use Homebrew.

   [Screenshot placeholder: Python download page]

2. Create and activate a Python virtual environment, then install packages:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

   [Screenshot placeholder: Terminal showing virtualenv commands]

Step 2 — (Optional) macOS system package

On macOS, the folder-picker dialog uses Tcl/Tk (tkinter). If the Add Project
dialog does not open, you may need to install Tcl/Tk via Homebrew:

```bash
brew install tcl-tk
```

   [Screenshot placeholder: Homebrew install output]

Step 3 — Start and use the dashboard (user steps)

3a. Start the dashboard server

```bash
python main.py
```

The app runs a small local web server. Open your browser to
`http://localhost:8080` to see the dashboard.

   [Screenshot placeholder: Browser open to localhost:8080]

3b. Add a project (Import)

- Click the "Add New Project" button in the dashboard.
- A folder picker will appear — choose the folder on your computer that holds
  the project (for example, your 10x output folder). The tool will copy that
  folder into `./projects/` so the dashboard can manage it.

   [Screenshot placeholder: Add New Project dialog]

3c. View and manage

- After the copy completes, the project appears in the project list. Click it
  to view metadata and use the QC helpers.

   [Screenshot placeholder: Project selected in dashboard]

## Notes & troubleshooting (non-technical)

- The app copies your selected folder — it does not delete your original data.
  If a folder with the same name already exists in `projects/`, a number will
  be added (for example `myproj-1`).
- The QC and plotting features use scientific Python packages (Scanpy, AnnData
  and their dependencies). These are listed in `requirements.txt` but may need
  extra system libraries. You can still use the dashboard to manage and view
  projects without running the full QC.

## Want help or improvements?

- If something doesn't work, tell me which step failed (which command and any
  error message) and I can help troubleshoot. The GitHub repository's
  "Issues" tab is a good place to report problems.

If you find this useful, please star the repo — it helps the project grow.

