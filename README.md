# BellaVista

BellaVista is a visualization tool utilizing the [napari](https://napari.org/) viewer for interactive exploration of imaging-based spatial transcriptomic data. This BellaVista release is designed for data already processed through the KoLab-MERFISH pipeline. BellaVista can be used to visualize the 4 key dataset components of each KoLab-MERFISH dataset: (1) Cell boundary & nuclear images (WGA & DAPI), (2) Cell segmentation boundaries, (3) Transcript locations, (4) Cell network connectivity graphs. Note: BellaVista is purely a visualization tool - it does not perform data processing or analysis :)

## Quick Start (with sample data)

This short demo will load a sample FOV from the TAC mouse heart. BellaVista is installed and run via the command line. Run the steps below in the `Terminal` (macOS/Linux) or `PowerShell` (Windows) application. 

> **Note:** BellaVista requires a monitor display to render images in the napari viewer.

### 1. Install [uv](https://docs.astral.sh/uv/getting-started/installation/)
```
# macOS/Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# Windows
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"
```
> [!IMPORTANT]
> Restart the shell (close and reopen the terminal) to complete the installation. 

### 2. Launch BellaVista 
```
uvx -p 3.12 bellavista
```

<br>

> **Linux note**: If you see the error: `Could not load the Qt platform plugin "xcb" in "" even though it was found.` Install the required system libraries: `sudo apt install libxcb-cursor0 libxkbcommon-x11-0 libxcb-xinerama0`. Then re-run step 2!

<br>

> [!NOTE]
> It will take a few minutes to download and create the required data files. The terminal will print updates & display progress bars for time consuming steps. The sample data is a single FOV from the TAC mouse heart and is ~17MB. 

After successfully loading BellaVista, you should see the message `Data Loaded!` in the terminal. A napari window should appear displaying the sample data similar to the image below:

<p align="middle">
<img src="https://raw.githubusercontent.com/pkosurilab/BellaVista-KoLab-MERFISH/main/images/BellaVista_demo_launch_screen.png" alt="BellaVista demo sample TAC dataset initial screen" width="800" />
</p>

Now, you can interactively move around the napari canvas to explore the data.
Try zooming in & out, plotting cell-type-specific transcripts, cell boundaries, and cell networks!

<p align="middle">
<img src="https://raw.githubusercontent.com/pkosurilab/BellaVista-KoLab-MERFISH/main/images/BellaVista_sample_demo.png" alt="BellaVista demo sample TAC dataset" width="800" />
</p>

## BellaVista Widget Menu

BellaVista uses the napari interface and features a widget located on the right side of the napari window to plot each dataset feature. The widget has four components:


<!-- 
<p align="left">
  
  <h3 align="left">Image Widget</h3>
  <img src="https://raw.githubusercontent.com/pkosurilab/BellaVista-KoLab-MERFISH/main/images/image_widget.png"
  alt="BellaVista image widget" width="200" />

  The image widget features a dropdown menu with the images that are available to plot. The demo dataset contains WGA and DAPI images. 
</p> -->


<p align="left">
  
  <h3 align="left">Gene Widget</h3>
  <img src="https://raw.githubusercontent.com/pkosurilab/BellaVista-KoLab-MERFISH/main/images/gene_widget.png"
  alt="BellaVista gene widget" width="200" />
  
  The gene widget features a dropdown menu with the genes in the gene panel, and a dropdown menu to select cell-type-specific transcripts. 

* "all transcripts" = all transcript molecules (un-partitioned) -- this includes all transcripts, both transcripts assigned to cells and unassigned transcripts
* "CM transcripts" = cardiomyocyte transcripts
* "EC transcripts" = endothelial cell transcripts
* "IC transcripts" = immune cell transcripts
* "FB transcripts" = fibroblast transcripts
* "CM 1-2-3 transcripts" = cardiomyocyte transcripts from CM1, CM2 & CM3 subclusters (spatially-scrubbed)

To plot the transcripts from a gene, select the gene from the dropdown gene list, select your cell-type of interest, then press `Plot`! Each individual transcript molecule will be plotted as a single point. The layer color is random, and can be changed using the "Layer Color" text input & button. The input color value can be provided as a hexcode value e.g. `#FF00FF` or by name e.g. `Magenta`.

</p>

<p align="left">
  
  <h3 align="left">Segmentation Widget</h3>
  <img src="https://raw.githubusercontent.com/pkosurilab/BellaVista-KoLab-MERFISH/main/images/segmentation_widget.png"
  alt="BellaVista segmentation widget" width="200" />

  The segmentation widget features a dropdown menu with cell-type-specific cell segmentation boundaries. 

* "all boundaries" = all segmentation mask boundaries
* "CM boundaries" = cardiomyocyte boundaries
* "EC boundaries" = endothelial cell boundaries
* "IC boundaries" = immune cell boundaries
* "FB boundaries" = fibroblast boundaries
* "CM 1-2-3 boundaries" = cardiomyocyte boundaries from CM1, CM2 & CM3 subclusters (spatially-scrubbed)

The boundaries for each cell type will be colored as follows: CM: pink, EC: green, IC: blue, FB: yellow. To change the color of the cell boundaries, use the "colormap" option on the layer control menu. For the best visualization, you should select a colormap starting with "single-hue"!

</p>

<p align="left">
  
  <h3 align="left">Network Widget</h3>
  <img src="https://raw.githubusercontent.com/pkosurilab/BellaVista-KoLab-MERFISH/main/images/network_widget.png"
  alt="BellaVista network widget" width="200" />

  The network widget can be used to plot the cell connectivity networks. The centroid of each cell is plotted as a node (napari point layer), and is colored by its corresponding cell type: CM: pink, EC: green, IC: blue, FB: yellow. Cardiomyocyte->cell-type-specific edges are plotted, and are colored by the cell-type identity of the corresponding cardiomyocyte's neighbor. Note: this network is cardiomyocyte centric, meaning if two cells are connected, but neither is a cardiomyocyte, then the connection will not be shown.
  </p>

<p align="left">
  
  <h3 align="left">Location Widget</h3>
  <img src="https://raw.githubusercontent.com/pkosurilab/BellaVista-KoLab-MERFISH/main/images/location_widget.png"
  alt="BellaVista location widget" width="200" />

  The location widget can be used to save and move to marked camera coordinates. Additionally, the "Load CSV" & "Export CSV" buttons can be used to export your saved locations or load previously saved locations from a previous session. 
  </p>

> [!NOTE]  
> If the input files for a feature was not provided or an error occurred during data processing, the corresponding widget may not be available. A detailed log file `error_log.log` can be found in the `BellaVista_outputs` subfolder inside the data folder. 

> **That's it for the quick start!** 

# Visualizing full MERFISH datasets

We will share the commands to visualize the Sham and TAC datasets from the mouse heart, including the private url-links. To visualize a dataset hosted on the web, use the following single-line command:

```
uvx -p 3.12 bellavista --dataset-url "url-link-to-dataset"
```

> [!NOTE]
> We recommend having at least 16GB of RAM and 10GB of disk space available to download and visualize the full Sham and TAC datasets. The files will be downloaded and extracted in the folder you're currently in (working directory of your terminal). 
>
> Each dataset contains WGA & DAPI images, tens-of-thousands of cells, and hundreds-of-millions of transcripts. For visualization, these data will be converted to visualization files that will also require approximately the same amount of space as the raw datasets. So please keep this in mind when downloading the data!
>
> It will take a few minutes to download and create the required data files. The terminal will print updates & display progress bars for time consuming steps.


After successfully launching BellaVista, you should see the message `Data Loaded!` in the terminal. A napari window should appear displaying the data similar to the image below (TAC dataset shown here):

<p align="middle">
<img src="https://raw.githubusercontent.com/pkosurilab/BellaVista-KoLab-MERFISH/main/images/TAC_WGA.png" alt="BellaVista TAC dataset zoom out" width="800" />
</p>

Now, you can interactively move around the napari canvas to explore the data.
Try zooming in & out, plotting cell-type-specific transcripts, cell boundaries, and cell networks!

<p align="middle">
<img src="https://raw.githubusercontent.com/pkosurilab/BellaVista-KoLab-MERFISH/main/images/BellaVista_demo.png" alt="BellaVista TAC dataset" width="800" />
</p>

> **That's it for the example Sham and TAC datasets!** 
