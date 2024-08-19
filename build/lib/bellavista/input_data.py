#Author: Annabelle Coles
from . import input_data_xenium
from . import input_data_merscope
from . import input_data_merlin

from ngff_writer.array_utils import to_tczyx
from ngff_writer.writer import open_ngff_zarr
from dask_image.imread import imread
import dask.array as da
import dask
import pickle
import os  
import json
import math
from datetime import datetime
from typing import Dict, List


def create_inputs(json_file: Dict):

    """
        Create required input files for bellavista.py

        The list of transformation is used to write registered images and 
        the transformation list is archived.

        Args:
            json_file: A dictionary storing input file parameters parsed from user's input JSON file.

    """

    data_folder = json_file.get("data_folder")
    if not os.path.exists(data_folder):
        raise FileNotFoundError(f"Directory {data_folder} does not exist.")

    json_file_input_files = json_file.get("input_files")

    bellavista_output_folder = json_file.get("bella_vista_output_folder")
    if not os.path.exists(bellavista_output_folder):
        print(f'Directory {bellavista_output_folder} does not exist -- creating the directory!')
        os.makedirs(bellavista_output_folder)

    # create a dictionary to store any exceptions caused by invalid/missing input files
    exceptions = {}

    system = json_file.get("system")

    if (system.lower() == "xenium"):
        print("Creating Bella Vista input files for 10x Genomics Xenium")
        exceptions = create_ome_zarr(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        exceptions = input_data_xenium.create_micron_pixel(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        exceptions = input_data_xenium.create_transcripts(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        exceptions = input_data_xenium.create_segmentations(data_folder, bellavista_output_folder, json_file_input_files, exceptions)

        with open(os.path.join(bellavista_output_folder, 'exceptions.json'), 'w') as f:
            json.dump(exceptions, f)


    elif (system.lower() == "merscope"):
        print("Creating Bella Vista input files for Vizgen MERSCOPE")
        exceptions = create_ome_zarr(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        exceptions = input_data_merscope.create_micron_pixel(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        exceptions = input_data_merscope.create_transcripts(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        exceptions = input_data_merscope.create_segmentations(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        
        with open(os.path.join(bellavista_output_folder, 'exceptions.json'), 'w') as f:
            json.dump(exceptions, f)

    else:
        print("Creating Bella Vista input files for Custom MERFISH (MERlin)")
        exceptions = create_ome_zarr(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        exceptions = input_data_merlin.create_micron_pixel(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        exceptions = input_data_merlin.create_transcripts(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        exceptions = input_data_merlin.create_segmentations(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        
        with open(os.path.join(bellavista_output_folder, 'exceptions.json'), 'w') as f:
            json.dump(exceptions, f)

    print("Bella Vista input files created!")


def create_ome_zarr(data_folder: str, bellavista_output_folder: str, json_file_input_files: Dict, exceptions: Dict):
    
    """
        Converts TIFF to multiscale pyramidal OME-Zarr images.

        Args:
            data_folder: Path to folder containing dataset files.
            bellavista_output_folder: Path to output data folder to store input files for bellavista.py
            json_file_input_files: A dictionary storing input file parameters parsed from user's input JSON file, containing image paths.
            exceptions: An empty dictionary.

        Returns:
            exceptions: An updated version of exceptions dictionary parsed if an error occured.
    """

    # if images have been processed previously, exit
    if os.path.exists(os.path.join(bellavista_output_folder, 'OMEzarrImages')):
        return exceptions

    # open GL max texture size
    napari_size_limit = 16384

    print("Creating OME-zarr image", datetime.now())
    images = json_file_input_files.get("images")
    z_plane = json_file_input_files.get("z_plane", 0)
    
    dask.config.set({"array.slicing.split_large_chunks": False})

    try:
        # images passed is a single image
        if isinstance(images, str):
            image_file_names = os.path.basename(images.replace('\\', '/')).split('.')[0]
            try:
                img_data = imread(os.path.join(data_folder, images))[z_plane]
                # calculate number of pyramidal scales required to prevent napari downscaling image
                max_dimension = max(img_data.shape)
                min_n_scales = math.ceil(math.log(max_dimension / napari_size_limit, 2)) + 1
                min_n_scales = max(min_n_scales, 1)

                with open_ngff_zarr(store=bellavista_output_folder, dimension_separator="/", overwrite=True) as f:
                    collection = f.add_collection(name="OMEzarrImages")
                    collection.add_image(image_name="Images",
                                            array=to_tczyx(img_data, axes_names=( "y", "x")), n_scales=min_n_scales)
            except:
                img_data = imread(os.path.join(data_folder, images))[0]
                max_dimension = max(img_data.shape)
                min_n_scales = math.ceil(math.log(max_dimension / napari_size_limit, 2)) + 1
                min_n_scales = max(min_n_scales, 1)

                with open_ngff_zarr(store = bellavista_output_folder, dimension_separator="/", overwrite=True) as f:
                    collection = f.add_collection(name = "OMEzarrImages")
                    collection.add_image(image_name = "Images", array=to_tczyx(img_data, axes_names=("y", "x")), n_scales = min_n_scales)

        # images passed is a list of images
        elif isinstance(images, List):

            img_shapes = []
            image_list = []
            # check if any of the images passed are 3D
            # assumption that image dimensions are (z,y,x)
            for image in images: 

                img_data = imread(os.path.join(data_folder, image))
                img_shape = img_data.shape
                img_shapes.append(img_shape)

                if(img_shape[0] > 1):
                    # extract image from specified z plane
                    try: 
                        img_data_2d = img_data[z_plane]
                        z_plane = 0
                    except: 
                        img_data_2d = img_data[z_plane]

                image_list.append(da.expand_dims(img_data_2d, axis=0))

            
            img_x_shapes = [shape[2] for shape in img_shapes]
            img_y_shapes = [shape[1] for shape in img_shapes]
            max_dimension = max(max(img_x_shapes), max(img_y_shapes))
            min_n_scales = math.ceil(math.log(max_dimension / napari_size_limit, 2)) + 1
            min_n_scales = max(min_n_scales, 1)

            image_file_names = [os.path.basename(image.replace('\\', '/')).split('.')[0] for image in images]
            with open_ngff_zarr(store = bellavista_output_folder, dimension_separator="/", overwrite=True) as f:
                collection = f.add_collection(name = "OMEzarrImages")
                collection.add_image(image_name = "Images", array=to_tczyx(da.concatenate(image_list), 
                                        axes_names=("c", "z", "y", "x")), n_scales = min_n_scales)
        
        # save image names, these will be used as layer names when loading napari
        f = open(os.path.join(bellavista_output_folder, "image_file_names.pkl"),"wb")
        pickle.dump(image_file_names, f)
        f.close()
        print("OME-zarr image saved!",datetime.now())

    except:
        # update exceptions dictionary to document an error
        if not os.path.exists(os.path.join(bellavista_output_folder, 'OMEzarrImages')):
            exceptions['valid_image'] = False

    return exceptions