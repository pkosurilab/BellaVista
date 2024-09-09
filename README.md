# Bella Vista

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://github.com/pkosurilab/BellaVista/blob/updates/misc-changes/images/bellavista_figure_darkmode.png?raw=true" width="900">
    <img alt="Bella Vista workflow" src="https://github.com/pkosurilab/BellaVista/blob/updates/misc-changes/images/bellavista_figure.png?raw=true" width="900">
  </picture>
</p>
<p align="center">

Bella Vista enables visualization of imaging-based spatial transcriptomics data. It is an open-source Python package currently supporting 10x Genomics Xenium, Vizgen MERSCOPE, and custom (home-built) MERFISH datasets, utilizing [napari](https://napari.org/) for interactive data exploration. We developed Bella Vista to help the spatial transcriptomics community explore their data and create reproducible paper-ready figures. For more information, see our [documentation website](https://baby-stringbean.readthedocs.io/).

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

Download sample data: Xenium mouse brain dataset (Replicate 3). To download the dataset, 10x Genomics may ask you to fill out a questionnaire.

[https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard](https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard)

<img src="https://github.com/pkosurilab/BellaVista/blob/main/images/xenium_testdata_location.png?raw=true" alt="Xenium sample data website location" width="600" />


1. Copy and save contents below into a new JSON file called `xenium_sample.json`
      - This sample JSON can also be found in the GitHub repository: [BellaVista/sample_json/xenium_sample.json](https://github.com/pkosurilab/BellaVista/blob/updates/misc-changes/sample_json/xenium_sample.json)
        
2. Replace the paths in `data_folder` and `bella_vista_output_folder` parameters
      - JSON files cannot interpret the blackslash character (\\), instead you should use a forward slash (/)

```
{ 
      "system": "xenium", 
      "data_folder": "/path/to/xenium_brain_rep3",
      "bella_vista_output_folder": "/path/to/xenium_brain_rep3/bellavista_outs",
      "create_bellavista_inputs": true,

      "visualization_parameters": {
          "plot_image": true,
          "plot_transcripts": true,
          "plot_allgenes": true,
          "genes_visible_on_startup": false,
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
3. In the terminal, run Bella Vista with the Xenium sample JSON:
    - The JSON file argument should contain the file path to the JSON file.
```
bellavista xenium_example.json
```
> [!NOTE]  
> It will take a few minutes to create the required data files.\
> The terminal will print updates & have progress bars for time consuming steps.

> [!IMPORTANT]  
> This is a large dataset, so if the program crashes or encounters a memory-related error, you may need to visualize a smaller subset of the data.
> For more information, see [What should I do if the program crashes?](https://baby-stringbean.readthedocs.io/en/latest/faq.html#frequently-asked-questions) in our FAQ on the documentation website.

Once successfully loaded, you should see the message `Data Loaded!` in the terminal.\
A napari window should appear displaying the data similar to the image below:

<img src="https://github.com/pkosurilab/BellaVista/blob/updates/misc-changes/images/xenium_initial.png?raw=true" alt="Initial napari load page"/>
<br/><br/>

Now, you can interactively move around the napari canvas to explore the data!\
Try zooming in & out, toggling layers on & off to see different spatial patterns:

<div style="position: relative; width: 100%; display: flex; justify-content: space-between; align-items: flex-end;">
  <img src="https://github.com/pkosurilab/BellaVista/blob/updates/misc-changes/images/xenium_brain_position_0_select.png?raw=true" alt="zoom out screenshot" style="width: 49%">
  <img src="https://github.com/pkosurilab/BellaVista/blob/updates/misc-changes/images/xenium_brain_position_1.png?raw=true" alt="zoom in screenshot" style="width: 49%">
</div>

<br/>

> [!TIP] 
> To visualize a single layer, and hide all other layers, `Option/Alt-click` on the visibility button (the eye, to the left of the layer name). 
>
> Check out [our FAQ](https://baby-stringbean.readthedocs.io/en/latest/faq.html#helpful-napari-tips) on the documentation website for more tips!


> [!NOTE]  
> Gene colors are assigned randomly every time Bella Vista is launched. So, the gene colors displayed in your window will be different from the image above. Refer to [our FAQ](https://baby-stringbean.readthedocs.io/en/latest/faq.html#helpful-napari-tips) on the documentation website for information on how to configure gene colors and other customizable visualization options.
>
> To reproduce the same colors every time you launch Bella Vista, refer to the [figure guide](https://baby-stringbean.readthedocs.io/en/latest/figure_guide.html) on the documentation website.


For an exact reproduction of the two screenshots above, please refer to the figure guide: [Reproducing sample figures (Xenium)](https://baby-stringbean.readthedocs.io/en/latest/figure_guide.html#reproducing-sample-figures-xenium) on the documentation website.


## To run Bella Vista with your own data, refer to the tutorials on the [documentation website](https://baby-stringbean.readthedocs.io). 