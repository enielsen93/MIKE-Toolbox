# MIKE Toolbox

The **MIKE Toolbox** is a collection of ArcGIS Python toolboxes designed to support hydraulic modeling workflows using MIKE Urban, MIKE+ and related tools. This is the **first public release**, containing a comprehensive suite of tools aimed at improving productivity, QA/QC, and integration with GIS environments.

### 🧰 Included Toolboxes

1. **Display MIKE Urban** – Visualize MIKE Urban database contents in ArcGIS.
2. **Display MIKE+** – Visualize MIKE+ contents in ArcGIS
3. **Dandas Toolbox** – Visualize and import Dandas contents
4. **Compare MIKE Databases** – Compare elements across MIKE Urban databases.
5. **Validate MIKE Urban** – Run a suite of consistency checks on MIKE Urban models.
6. **Catchment Toolbox** – Tools for analyzing and editing catchments
7. **Set GroundLevel From DHM** – Assign terrain elevations based on DHM raster data using API.
8. **Pipe Design Toolbox** – Automatic Pipe Sizing and misc. tools for Detailed Pipe Design
9. **MIKE+ Field Calculator** – Field Calculator for MIKE+ Databases (SQLITE)
10. **Display MIKE Results** – Import and visualize simulation results.
11. **MIKE Flood Toolbox** – Support tools for MIKE FLOOD workflows.
12. **Filter Map to Selection** – Filter map content to the current MIKE selection.
13. **Export MIKE Urban To CAD or DDS** – Export functionality for CAD or Dandas DDS formats.
14. **Compress Folder** – Utility to zip model folders for archiving or sharing.
15. **Run MOUSE in Batch** – Batch-execution of MIKE Urban simulations using MOUSE engine.

## Getting Started
These instructions will guide you through installing and using MIKE Urban Tools.

---

## Installation

### Step 1: Download the Tools
1. [Click here to download the ZIP file](https://github.com/enielsen93/MIKE-Urban-Tools/archive/refs/heads/main.zip).
2. Extract the contents of the ZIP file to a local folder on your computer.

### Step 2: Load the Tools into ArcGIS
1. Open **ArcMap** or **ArcGIS Pro**.
2. In the **Catalog** pane, browse to the folder where you extracted the tools.
3. The tools are now available for use!

---

## Installing Python Requirements in ArcGIS Pro
If you are using ArcGIS Pro, follow these steps to safely install dependencies `mikeio` and `mikeio1d`.

### 1. Clone your Python environment
Open **ArcGIS Pro → Project → Python → Manage Environments**, select your default environment, and click **Clone**.  
> This ensures the default Python environment remains intact.

### 2. Activate the cloned environment
Select your cloned environment and click **Activate**.

### 3. Open a Command Prompt
Press `Win + R`, type `cmd`, and hit **Enter**.

### 4. Install packages
Copy and paste these commands (adjust the path if your ArcGIS Pro is installed elsewhere):

```cmd
"C:\Progra~1\ArcGIS\Pro\bin\Python\scripts\propy.bat" -m pip install --user mikeio mikeio1d mikecore
```
---

## Need Help?
If you encounter any issues or have questions, feel free to [open an issue](https://github.com/enielsen93/MIKE-Toolbox/issues) on this repository. We’re here to help!
