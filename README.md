[![Manuscript](https://img.shields.io/badge/DOI-10.1016/j.bpj.2024.11.3199-orange.svg?logo=doi)](https://doi.org/10.1016/j.bpj.2024.11.3199)
[![pypi](https://img.shields.io/badge/pypi-bellavista-blue.svg?logo=pypi)](https://pypi.org/project/bellavista)
[![bioconda](https://img.shields.io/badge/bioconda-bellavista-blue.svg?logo=anaconda)](https://anaconda.org/bioconda/bellavista)
[![Docker Repository on Quay](https://img.shields.io/badge/container-bellavista-blue?logo=docker)](https://quay.io/repository/bgruening/bellavista)
[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/root?tool_id=interactive_tool_bellavista)

# BellaVista

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://github.com/pkosurilab/BellaVista/blob/main/images/bellavista_figure_darkmode.png?raw=true" width="900">
    <img alt="BellaVista workflow" src="https://github.com/pkosurilab/BellaVista/blob/main/images/bellavista_figure.png?raw=true" width="900">
  </picture>
</p>
<p align="center">

BellaVista enables visualization of imaging-based spatial transcriptomics data. It is an open-source Python package currently supporting 10x Genomics Xenium, Vizgen MERSCOPE, and custom (home-built) MERFISH datasets, utilizing [napari](https://napari.org/) for interactive data exploration. We developed BellaVista to help the spatial transcriptomics community explore their data and create reproducible paper-ready figures. For more information, see our [documentation website](https://bellavista.readthedocs.io/en/latest/).

## Installation
The following instructions require that you have [Anaconda](https://www.anaconda.com/) installed.
- In MacOS, run the following commands from the Terminal.
- In Windows, run the following commands from the Anaconda Prompt.
- BellaVista requires Python 3.9 or above and is dependent on GPU for rendering. 

Create and activate a new virtual environment:

```
conda create -n bellavista_env python=3.12
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

BellaVista can also be installed directly from bioconda:

```
conda create -n bellavista_env -c bioconda bellavista
conda activate bellavista_env
```

A Docker image of the BellaVista Tool is also provided that can be deployed locally, in compute clusters, or in the cloud. The container includes all dependencies required by BellaVista to run in an isolated container.

The container can be pulled via:

```
docker pull quay.io/bgruening/bellavista:latest
```

For more information about the container, please refer to the [docker-bellavista repository](https://github.com/usegalaxy-eu/docker-bellavista).

---
## Quickstart (with sample data)

Below is a short tutorial for loading BellaVista with sample Xenium data.

1. Download sample data from Zenodo: [Xenium mouse brain dataset (Replicate 3)](https://zenodo.org/records/14279832)

      - Unzip the downloaded zip file. This will create a folder named "xenium_mouse_brain_rep3".
      - Take note of your local path to this folder, as you will need this path when running BellaVista.

<img src="https://github.com/pkosurilab/BellaVista/blob/updates/zenodo-tutorial/images/zenodo_download.png?raw=true" alt="Xenium sample data zenodo" width="600" />

2. Run BellaVista from the command line with the Xenium sample data:

      - Note: Before running this command, replace "/path/to/" with the actual path to the Xenium sample data folder.

```
bellavista --xenium-sample /path/to/xenium_mouse_brain_rep3
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
> bellavista --xenium-sample-lite /path/to/xenium_mouse_brain_rep3
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
> Gene colors are assigned randomly every time BellaVista is launched. So, the gene colors displayed in your window will be different from the image above. Refer to our [FAQ](https://bellavista.readthedocs.io/en/latest/faq.html#helpful-napari-tips) on the documentation website for information on how to configure gene colors and other customizable visualization options.
>
> To reproduce the same colors every time you launch BellaVista, refer to the [figure guide](https://bellavista.readthedocs.io/en/latest/figure_guide.html) on the documentation website.
<br/>

For an exact reproduction of the screenshots above, please refer to the figure guide: [Reproducing sample figures (Xenium)](https://bellavista.readthedocs.io/en/latest/figure_guide.html#reproducing-sample-figures-xenium) on the documentation website.


### To run BellaVista with your own data, refer to the tutorials on the [documentation website](https://bellavista.readthedocs.io/en/latest/tutorials.html).

## Using BellaVista on Galaxy Platform

Galaxy is a free, open-source system for analyzing data, authoring workflows, training and education, publishing tools, managing infrastructure, and more. 
You can easily run Bellavista as an Interactive tool on Galaxy.

[Try it out here!](https://usegalaxy.eu/?tool_id=interactive_tool_bellavista&version=latest)

You can check "[Galaxy Basics for everyone](https://training.galaxyproject.org/training-material/topics/introduction/tutorials/galaxy-intro-101-everyone/tutorial.html)" training to learn more about data analysis using Galaxy!