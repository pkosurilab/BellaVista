<p align="right">
  <img alt="Bella Vista logo" src="https://github.com/pkosurilab/BellaVista/blob/pypi-documentation/images/bellavista_logo_graphic.png?raw=true" width="150">
</p>

# Bella Vista: Open-Source Visualization for Imaging-Based Spatial Transcriptomics

Bella Vista is an open-source Python package developed for 10x Genomics Xenium, Vizgen MERSCOPE, 
and custom (home-built) MERFISH datasets utilizing [napari](https://napari.org/) for interactive data exploration. 
We developed Bella Vista to help the spatial transcriptomics community explore their data and create reproducible paper-ready figures. See our [work-in-progress documentation website](https://baby-stringbean.readthedocs.io/) for further documentation and usage instructions. 
> **Note**
> 
> This project is currently under active development. 

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://github.com/pkosurilab/BellaVista/blob/main/images/bellavista_figure_darkmode.png?raw=true" width="900">
    <img alt="Bella Vista workflow" src="https://github.com/pkosurilab/BellaVista/blob/main/images/bellavista_figure.png?raw=true" width="900">
  </picture>
</p>
<p align="center">


## Installation
The following instructions require that you have [Anaconda](https://www.anaconda.com/) installed.
- It is recommended to create an Anaconda virtual environment to prevent conflicting package dependencies. 
<!-- - The package can be installed from PyPI via [pip](https://pypi.org/project/pip/) (recommended) or from the [GitHub repository](https://github.com/pkosurilab/BellaVista). -->
- Bella Vista requires python 3.9 or above.

Create and activate a new virtual environment:

```
conda create -n bellavista_env python=3.9
conda activate bellavista_env
```

<!-- ### Installation via pip:
```
pip install bellavista
```
--- -->
### Installation from GitHub repository:

```
conda install git
git clone https://github.com/pkosurilab/BellaVista
pip install -e BellaVista
```

---
## Getting Started (with sample data)
Download sample data: [Xenium mouse brain dataset (Replicate 3)](https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard)

<img src="https://github.com/pkosurilab/BellaVista/blob/main/images/xenium_testdata_location.png?raw=true" alt="Xenium sample data website location" width="600" />
https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard
<br/><br/>

- Copy and save contents below into a new JSON file called `xenium_example.json`
- Replace the paths in the `data_folder` and `bella_vista_output_folder` properties
- **Note**: JSON files cannot interpret the blackslash character (\\) instead you should use a forward slash (/)

```
{ 
    "system": "xenium", 
    "data_folder": "/path/to/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs",
    "bella_vista_output_folder": "/path/to/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs/bellavista_outs",
    "create_bellavista_inputs": true,

    "visualization_parameters": {
        "plot_image": true,
        "plot_transcripts": true,
        "plot_allgenes": true,
        "plot_cell_seg": false,
        "plot_nuclear_seg": false,
        "transcript_point_size": 0.75,
        "contrast_limits": [0, 5000],
        "rotate_angle": 180
    },

    "input_files": {
        "images": "morphology_mip.ome.tif",
        "z_plane": 5,
        "transcript_filename": "transcripts.parquet",
        "cell_segmentation": "cell_boundaries.parquet",
        "nuclear_segmentation": "nucleus_boundaries.parquet"
    }
}
```
Run Bella Vista with Xenium sample data:
```
bellavista xenium_example.json
```
**Note**: It will take a few minutes to create the required data files.\
The terminal will print updates & have progress bars for time consuming steps.

Once successfully loaded, you should see the message `Data Loaded!` in the terminal.\
A napari window should appear displaying the data similar to the image below:

<img src="https://github.com/pkosurilab/BellaVista/blob/main/images/xenium_brain_position0_allgenes.png?raw=true" alt="Initial napari load page" width="600" />
<br/><br/>

Now, you can interactively move around the napari canvas to explore the data!\
Try zooming in & out, toggling layers on & off to see different spatial patterns:

<img src="https://github.com/pkosurilab/BellaVista/blob/pypi-documentation/images/xenium_brain_position0_selectgenes.png?raw=true" alt="Zoom out of napari" width="600" /> 
<br/><br/>
<img src="https://github.com/pkosurilab/BellaVista/blob/main/images/xenium_brain_position1.png?raw=true" alt="Zoom out of napari with selected genes visible" width="600" />

The example JSON file can also be found on the Bella Vista GitHub repository: https://github.com/pkosurilab/BellaVista/tree/main/sample_json/xenium_example.json

<!-- ## To run Bella Vista with your own data, see the [Bella Vista GitHub](https://github.com/pkosurilab/BellaVista) for documentation and usage instructions -->
