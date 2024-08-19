<p align="right">
  <picture><img alt="Bella Vista logo." src="./bellavista_logo.png" width="150">
  </picture>
</p>
<p align="right">

# Bella Vista: Open-Source Visualization for Imaging-Based Spatial Transcriptomics
Bella Vista is an open-source Python package developed for 10x Genomics Xenium, Vizgen MERSCOPE, 
and custom (home-built) MERFISH datasets utilizing [napari](https://napari.org/) for interactive data exploration. 
We developed Bella Vista to help the spatial transcriptomics community explore their data and create reproducible paper-ready figures. See the [Bella Vista GitHub](https://github.com/pkosurilab/BellaVista) for documentation and usage instructions. 

## Installation
The following instructions require that you have [Anaconda](https://www.anaconda.com/) installed.
- It is recommended to create an Anaconda virtual environment to prevent conflicting package dependencies. 
- The package can be installed from PyPI via [pip](https://pypi.org/project/pip/) (recommended) or from the [GitHub repository](https://github.com/pkosurilab/BellaVista).

Create and activate a new virtual environment:

```
conda create -n bellavista_env python=3.9
conda activate bellavista_env
```

### Installation via pip:
```
pip install bellavista
```
---
### Installation from GitHub repository:

```
git clone https://github.com/pkosurilab/BellaVista
```
```
cd BellaVista
```
```
pip install .
```
---
## Getting Started (with sample data)
Download sample data: [Xenium mouse brain dataset (replicate 3)](https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard)
- Copy and save contents below into a new JSON file called `xenium_example.json`
- Replace the paths in the `data_folder` and `bella_vista_output_folder` properties
```
{ 
    "system": "xenium", 
    "data_folder": "/path/to/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs",
    "bella_vista_output_folder": "/path/to/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs/bellavista_outs",
    "create_bellavista_inputs": true,

    "parameters": {
        "plot_image": true,
        "plot_transcripts": true,
        "plot_allgenes": true,
        "plot_cell_seg": false,
        "plot_nuclear_seg": false,
        "transcript_point_size": 0.75,
        "contrast_limits": [800,6000]
    },

    "input_files": {
        "images": "morphology_mip.ome.tif",
        "z_plane": 0,
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

<img src="./xenium_example/xenium_position0.png" width="500" />

Now, you can interactively move around the napari canvas to explore the data!\
Try zooming in & out, toggling layers on & off to see different spatial patterns:

<img src="./xenium_example/xenium_position2_ALL.png" width="500" /> <img src="./xenium_example/xenium_position2_select.png" width="500" />

The example JSON file can also be found on the Bella Vista GitHub repository: https://github.com/pkosurilab/BellaVista/tree/main/sample_json/xenium_brain_rep3.json

## To run Bella Vista with your own data, see the [Bella Vista GitHub](https://github.com/pkosurilab/BellaVista) for documentation and usage instructions