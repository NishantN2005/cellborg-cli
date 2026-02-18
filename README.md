<img width="1728" height="1117" alt="Screenshot 2026-02-17 at 7 00 06 PM" src="https://github.com/user-attachments/assets/e2b0e85c-c61a-444f-b654-37d527486b54" />
# Cellborg CLI

This tool provides a small, easy-to-use dashboard for managing single-cell
projects on your computer. It is intended for researchers who want a simple
way to organize data folders and run basic QC steps without writing code.

If this tool helps your work, please star the repository — it really helps!


## Get the code — step-by-step
1. Click the green "Code" button and copy the HTTPS URL (it looks like
   `https://github.com/NishantN2005/cellborg-cli.git`).

   <img width="1703" height="786" alt="Screenshot 2026-02-17 at 6 48 47 PM" src="https://github.com/user-attachments/assets/7e8678f0-f026-4789-8ef7-cbfbe1ef1579" />


2. Open the Terminal application on your computer and go to your desired location, then run these commands:

```bash
git clone https://github.com/NishantN2005/cellborg-cli.git
cd cellborg-cli
```

## Easy setup (three steps)

Follow these three steps to prepare and run the dashboard locally. After each
subsection in Step 3 there is a placeholder where you can add a screenshot.

Step 1 — Prepare Python

1. Install Python 3.10 or newer if you do not already have it. On macOS you
   can download the installer from python.org or use Homebrew.

   <img width="1709" height="882" alt="Screenshot 2026-02-17 at 6 50 18 PM" src="https://github.com/user-attachments/assets/9963df2b-14c3-4fdf-85d3-9568a0f09ecc" />


2. Create and activate a Python virtual environment, then install packages:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```


Step 2 — (Optional) macOS system package

On macOS, the folder-picker dialog uses Tcl/Tk (tkinter). If the Add Project
dialog does not open, you may need to install Tcl/Tk via Homebrew:

```bash
brew install tcl-tk
```


Step 3 — Start and use the dashboard (user steps)

3a. Start the dashboard server

```bash
python main.py
```

The app runs a small local web server which should open up automatically. If not try going to google and typing in localhost:8080 into the search bar.

   <img width="1728" height="1117" alt="Screenshot 2026-02-17 at 6 52 14 PM" src="https://github.com/user-attachments/assets/1f1b133b-9987-4302-a1d4-3a05cf55e7a4" />

3b. Add a project (Import)

- Click the "Add New Project" button in the dashboard.
- A folder picker will appear — choose any file within the folder that you want to upload. NOTE: The entire folder will get copied NOT just the file selected. This has to do more with the compatibility issues regarding the packages used than anything else significant.
<img width="1728" height="1117" alt="Screenshot 2026-02-17 at 6 55 14 PM" src="https://github.com/user-attachments/assets/81fe0525-e691-4f5c-a349-52a277f9cfe9" />

   
- Give the project a title and description.

   <img width="1728" height="1117" alt="Screenshot 2026-02-17 at 6 55 42 PM" src="https://github.com/user-attachments/assets/a239c087-2ca4-460c-a1aa-72dfe778574a" />


3c. View and manage

- After the copy completes, the project appears in the project list. Click it
  to view metadata and use the QC helpers.

   <img width="1728" height="1117" alt="Screenshot 2026-02-17 at 6 56 06 PM" src="https://github.com/user-attachments/assets/4b299225-5534-40c1-8042-bcad355edc40" />

QC
- Click 'Run QC', the screen will freeze for a little while some background processes run. If it is succcessful, you will be directed to the qc page.
  
<img width="1728" height="1117" alt="Screenshot 2026-02-17 at 6 57 33 PM" src="https://github.com/user-attachments/assets/5cfa06a2-a334-4af3-a79b-56961eb1b659" />

- Once you filter your data and click 'Save', click 'Go Back to Projects'
  <img width="1728" height="1117" alt="Screenshot 2026-02-17 at 6 58 37 PM" src="https://github.com/user-attachments/assets/1937659f-437c-451b-9c83-5b38bfa348b0" />

- Notice QC changes to PA (Processing & Annotations), click on that (same deal with the pausing so only click once).

<img width="1728" height="1117" alt="Screenshot 2026-02-17 at 6 59 23 PM" src="https://github.com/user-attachments/assets/e82822fa-a1ba-464d-a9e6-3cda3020c125" />

- Select your desired Leiden resolution and click 'Run Clustering'
<img width="1728" height="1117" alt="Screenshot 2026-02-17 at 6 59 53 PM" src="https://github.com/user-attachments/assets/4594a0ef-a7e2-48e0-a95d-60f234252d05" />

- Click 'Next Page' and Annotate your clusters, switch between genes using the dropdown on the top-left of the graph
<img width="1728" height="1117" alt="Screenshot 2026-02-17 at 7 00 57 PM" src="https://github.com/user-attachments/assets/784009a2-f7af-4250-9a8b-d8a2c2fd1e90" />

- Once everything looks good and you saved your annotations, click 'Next Page'
<img width="1728" height="1117" alt="Screenshot 2026-02-17 at 7 01 46 PM" src="https://github.com/user-attachments/assets/08327963-10d0-403a-ba94-7ab155edc9b4" />

- Click 'Run Analysis' and you will get directed to the analysis dashboard
   
- <img width="1728" height="1117" alt="Screenshot 2026-02-17 at 7 02 01 PM" src="https://github.com/user-attachments/assets/ce00c61c-dcb6-4cd8-81ad-40085cdc08ae" />

- All of the analysis features are implemented and their use is very intuitive. In the future, I will put all of these instructions in wiki's so it is formatted better, but for now have fun!

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

