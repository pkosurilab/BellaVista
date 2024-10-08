# Bella Vista

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://github.com/pkosurilab/BellaVista/blob/main/images/bellavista_figure_darkmode.png?raw=true" width="900">
    <img alt="Bella Vista workflow" src="https://github.com/pkosurilab/BellaVista/blob/main/images/bellavista_figure.png?raw=true" width="900">
  </picture>
</p>
<p align="center">

Bella Vista enables visualization of imaging-based spatial transcriptomics data. It is an open-source Python package currently supporting 10x Genomics Xenium, Vizgen MERSCOPE, and custom (home-built) MERFISH datasets, utilizing [napari](https://napari.org/) for interactive data exploration. We developed Bella Vista to help the spatial transcriptomics community explore their data and create reproducible paper-ready figures. For more information, see our [documentation website](https://bellavista.readthedocs.io/en/latest/).

## Installation
The following instructions require that you have [Anaconda](https://www.anaconda.com/) installed.
- In MacOS, run the following commands from the Terminal.
- In Windows, run the following commands from the Anaconda Prompt.
- Bella Vista requires Python 3.9 or above and is dependent on GPU for rendering. 

Create and activate a new virtual environment:

```
conda create -n bellavista_env python
conda activate bellavista_env
```

Installation via pip:
```
pip install bellavista
```

Alternatively, you can install from GitHub:

```
conda install git
git clone https://github.com/pkosurilab/BellaVista
pip install -e BellaVista
```

---
## Quickstart (with sample data)

Below is a short tutorial for loading Bella Vista with sample Xenium data.

1. Download sample data: [Xenium mouse brain dataset (Replicate 3)](https://www.10xgenomics.com/datasets/fresh-frozen-mouse-brain-replicates-1-standard)

      - To download the dataset, 10x Genomics may ask you to fill out a questionnaire.
      - Unzip the downloaded zip file. This will create a folder named "Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs".
      - Take note of your local path to this folder, as you will need this path when running Bella Vista.

<img src="https://github.com/pkosurilab/BellaVista/blob/main/images/xenium_testdata_location.png?raw=true" alt="Xenium sample data website location" width="600" />

2. Run Bella Vista from the command line with the Xenium sample data:

      - Note: Before running this command, replace "/path/to/" with the actual path to the Xenium sample data folder.

```
bellavista --xenium-sample /path/to/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs
```
<br/>

> [!NOTE]  
> It will take a few minutes to create the required data files.

<br/>
Once successfully loaded, you should see the message `Data Loaded!` in the terminal. 

A napari window should appear displaying the data similar to the image below:

<img src="https://github.com/pkosurilab/BellaVista/blob/updates/misc-changes/images/xenium_initial.png?raw=true" alt="Initial napari load page"/>

<br/>

> [!TIP]
> This is a large dataset, so if the program encounters a memory-related error, try visualizing a smaller subset of the data:
> ```
> bellavista --xenium-sample-lite /path/to/Xenium_V1_FF_Mouse_Brain_MultiSection_3_outs
>```
> 
Now, you can interactively move around the napari canvas to explore the data!\
Try zooming in & out, toggling layers on & off to see different spatial patterns:


<p align="left">
  <img src="https://github.com/pkosurilab/BellaVista/blob/main/images/xenium_brain_position_0_select.png?raw=true" alt="zoom out screenshot" />
  <img src="https://github.com/pkosurilab/BellaVista/blob/main/images/xenium_brain_position_1.png?raw=true" alt="zoom in screenshot" />
  <img src="https://github.com/pkosurilab/BellaVista/blob/updates/misc-changes/images/xenium_brain_position_2_cellbounds.png?raw=true" alt="zoom in cellbounds screenshot" />
</p>

<br/>

> [!TIP] 
> To visualize a single layer, and hide all other layers, `Option/Alt-click` on the visibility button (the eye, to the left of the layer name). 
>
> Check out our [FAQ](https://bellavista.readthedocs.io/en/latest/faq.html#helpful-napari-tips) on the documentation website for more tips!


> [!NOTE]  
> Gene colors are assigned randomly every time Bella Vista is launched. So, the gene colors displayed in your window will be different from the image above. Refer to our [FAQ](https://bellavista.readthedocs.io/en/latest/faq.html#helpful-napari-tips) on the documentation website for information on how to configure gene colors and other customizable visualization options.
>
> To reproduce the same colors every time you launch Bella Vista, refer to the [figure guide](https://bellavista.readthedocs.io/en/latest/figure_guide.html) on the documentation website.
<br/>

For an exact reproduction of the screenshots above, please refer to the figure guide: [Reproducing sample figures (Xenium)](https://bellavista.readthedocs.io/en/latest/figure_guide.html#reproducing-sample-figures-xenium) on the documentation website.


### To run Bella Vista with your own data, refer to the tutorials on the [documentation website](https://bellavista.readthedocs.io/en/latest/tutorials.html).
