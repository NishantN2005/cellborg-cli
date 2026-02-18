# Cellborg CLI — Simple project dashboard for biologists

This tool provides a small, easy-to-use dashboard for managing single-cell
projects on your computer. It is designed for researchers who want a simple
way to organize data folders and run basic QC steps without writing code.

If this tool helps your work, please star the repository — it really helps!

If you find this useful, please star the repo — it helps the project grow.

## What it does (in plain words)

- Shows a list of projects found in the `projects/` folder.
- Lets you import (copy) a project folder from your computer using a Finder
  / file-picker dialog. The folder will be copied into `./projects/` so the app
  can keep a local copy.
- Provides helpers for running QC and generating plots (advanced features).

## Easy setup (3 steps)

1. Install Python 3.10 or newer. On macOS you can use the official installer
   from python.org or Homebrew.

2. Create and activate a Python virtual environment, then install packages:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

3. Start the dashboard and open it in your web browser:

```bash
python main.py
```

The dashboard will pop-up with the projects list and right side of the screen being empty.

Click add projects to copy over the desired folder, once you do this enter a title and description. The project should now be visible in the list.


When you click the project options should appear on the right to do QC.

When you click 'Run QC' the page will stay still for a little, this is the processing happening in the background (I haven't implemented loaders yet). You will know it is finished if you get directed to the next page, else you will see an error log in the console.


Once you finish QC, you will return to the dashboard and the QC button will change into 'Run PA'. PA stands for 'Processing & Annotations'. 

Click PA


Once PA is done, you will return to the dashboard. Click 'Run Analysis'.

## Notes & troubleshooting (non-technical)
- The app copies your selected folder — it does not delete your original data.
  If a folder with the same name already exists in `projects/`, a number will
  be added (for example `myproj-1`).
- The QC/plot features use scientific Python packages (Scanpy, AnnData). These
  are included in `requirements.txt` but may need extra system libraries. You
  can still use the dashboard to manage and view projects without running QC.

## Want help or improvements?

- If something doesn't work, tell me which step failed (which command and the error message) and I can help troubleshoot by clicking 'Issues' at the top of the github repo and creating a new issue.

