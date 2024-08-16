# Bella Vista: Open-Source Visualization for Imaging-Based Spatial Transcriptomics
Bella Vista is an open-source Python package developed for 10x Genomics Xenium, Vizgen MERSCOPE, 
and custom (home-built) MERFISH datasets utilizing [napari](https://napari.org/) for interactive data exploration. 
We developed Bella Vista to help the spatial transcriptomics community explore their data and create reproducible paper-ready figures.

## Installation
The following instructions require that you have [Anaconda](https://www.anaconda.com/) installed.\
It is recommended to create an Anaconda virtual environment to prevent conflicting package dependencies. \
Bella Vista requires python 3.9

Create and activate a new virtual environment:

```
conda create -n bellavista_env python=3.9
conda activate bellavista_env
```
Clone the Bella Vista package:

```
git clone https://github.com/pkosurilab/BellaVista
cd BellaVista
```
Install required packages:
```
pip install -r requirements.txt
```
Download sample data: [Xenium mouse brain dataset (replicate 3)](https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard)\
Paste local path to sample data in the `data_folder` property in `BellaVista/sample_json/xenium_brain_rep3.json`\
\
Load Bella Vista with sample data:
```
cd bellavista
python bellavista.py ../sample_json/xenium_brain_rep3.json
```
**Note**: It will take a few minutes to create the required data files.\
The terminal will print updates & have progress bars for time consuming steps.

Once successfully loaded, you should see the message `Data Loaded!` in the terminal.\
A napari window should appear displaying the data similar to the image below.

<p align="center">
  <img src="./test/xenium_position0.png" width="800" />
</p>

Now, you can interactively move around the napari canvas to explore the data!\
Try zooming in & out, toggling layers on & off to see different spatial patterns!


<img src="./test/xenium_position2_ALL.png" width="500" /> <img src="./test/xenium_position2_select.png" width="500" />

To visualize your own data, you will need to create a input JSON file as explained below.\
With your JSON file, you can load your data with the command:

```
python bellavista.py /path/to/your_dataset.json
```

--- 
## Getting Started

The *only* input file required to visualize data in Bella Vista is a JSON file containing dataset-specific parameters. 
The input parameters for Bella Vista are the same for all technologies, however, the required input files are different depending on which technology you are using. 
Input file documentation and sample JSON files for [Xenium](#xenium-input-files), [MERSCOPE](#merscope-input-files), and [MERFISH-MERlin](#merfish-merlin-input-files) can be found below. 

Input JSON file organization:
```
  {
    "system": "spatial_technology",
    "data_folder": "/path/to/dataset",
    "bella_vista_output_folder": "path/to/dataset/bellavista_outs",
    "create_bellavista_inputs": true,
    
    "parameters": {
        "plot_image": true,
        "plot_transcripts": true,
        "plot_cell_seg": true
    },
  
    "input_files": {
        "images": "image_filename.tif",
        "transcript_filename": "transcript_filename",
        "cell_segmentation": "cell_seg_filename"
    }
  } 
```
## Input properties for JSON file

**system**: *string*
- "Xenium", "MERSCOPE" or "MERlin"

**plot_image**: *string*
- Path to folder containing dataset outputs
  
**bella_vista_output_folder**: *string*
- Path to save & load Bella Vista input files from
  
**create_bellavista_inputs**: *boolean*
- Create required input files for Bella Vista from dataset outputs
- Must be True when first loading data
- Can be False in subsequent runs (since files have already been created)
  
---
## Input parameters for JSON file

**plot_image**: *boolean, default=False*
- Display image(s)

**plot_transcripts**: *boolean, default=False*
- Plot gene transcript spatial coordinates

**plot_allgenes**: *boolean, default=True*
- Plot transcripts for all gene IDs.
- If False, only gene IDs in `selected_genes` will be plotted

**selected_genes**: *1D array of strings, default=None*
- Plot transcripts only for specified gene IDs

**plot_cell_seg**: *boolean, default=False*
- Plot cell segmentation

**plot_nuclear_seg**: *boolean, default=False*
- Plot nuclear segmentation

**transcript_point_size**: *float, default=1.0*
- Point size for plotting transcript coordinates

**contrast_limits**: *tuple array of integers, default=None*
- Values in the range [0, 65535].
- Contrast limits for a displayed image
--- 
## Xenium input files

**transcript_filename**: *string*
- relative path to a Parquet or CSV file containing transcript spatial locations.
- If None, no transcripts will be prepared

**images**: *string or 1D array of strings*
- relative path to image file(s).
- Must be an OME-TIFF or TIFF file.
- If None, no images will be displayed

**z_plane**: *integer, default=0*
- z-plane of image to be used

**cell_segmentation**: *string*
- relative path to Parquet or Zarr file containing cell segmentations.
- If None, no cell segmentations will be prepared

**nuclear_segmentation**: *string*
- relative path to Parquet or Zarr file containing nuclear segmentations.
- If None, no nuclear segmentations will be prepared
---
#### Xenium example
Example JSON file for publicly-available [Xenium mouse brain dataset (replicate 3)](https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard)
```
    { 
    "system": "xenium", 
    "data_folder": "/Users/kosuri_lab/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs",
    "bella_vista_output_folder": "/Users/kosuri_lab/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs/bellavista_outs",
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
This can also be found here: [BellaVista/sample_json/xenium_brain_rep3.json](./sample_json/xenium_brain_rep3.json)


## MERSCOPE input files

work in progress...

## MERFISH-MERlin input files

work in progress...

