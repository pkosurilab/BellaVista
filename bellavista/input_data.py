#Author: Annabelle Coles
from .utils.ngff_writer.array_utils import to_tczyx
from .utils.ngff_writer.writer import open_ngff_zarr
from dask_image.imread import imread
import dask.array as da
import dask
import pickle
import os  
import json
import math
from datetime import datetime
from typing import Dict, List
import logging

def setup_logger(bellavista_output_folder: str):
    # Remove any handlers that were set up earlier
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
        
    # Set a log file in the output folder
    logging.basicConfig(filename=os.path.join(bellavista_output_folder, 'error_log.log'),  # Absolute path
                        level=logging.WARNING,
                        format='%(asctime)s - %(levelname)s - %(message)s')

def create_inputs(json_file: Dict):

    '''
        Create required input files for bellavista.py

        The list of transformation is used to write registered images and 
        the transformation list is archived.

        Args:
            json_file: A dictionary storing input file parameters parsed from user's input JSON file.

    '''

    data_folder = json_file.get('data_folder')
    data_folder = data_folder.replace('\\', '/')
    if not os.path.exists(data_folder):
        raise FileNotFoundError(f'Directory {data_folder} does not exist.')

    json_file_input_files = json_file.get('input_files')

    bellavista_output_folder = os.path.join(data_folder, "BellaVista_output")
    if not os.path.exists(bellavista_output_folder):
        print(f'Directory {bellavista_output_folder} does not exist -- creating the directory!')
        os.makedirs(bellavista_output_folder)
    else: 
        print(f'Directory {bellavista_output_folder} exists, processing new input files')
    
    setup_logger(bellavista_output_folder)

    # add local imports here so errors will be logged to error_log.log
    from . import input_data_xenium
    from . import input_data_merscope
    from . import input_data_merlin
    
    # create a dictionary to store any exceptions caused by invalid/missing input files
    exceptions = {}

    system = json_file.get('system')

    if (system.lower() == 'xenium'):
        print('Creating Bella Vista input files for 10x Genomics Xenium')
        exceptions = create_ome_zarr(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        exceptions = input_data_xenium.create_micron_pixel(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        exceptions = input_data_xenium.create_transcripts(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        exceptions = input_data_xenium.create_segmentations(data_folder, bellavista_output_folder, json_file_input_files, exceptions)

        with open(os.path.join(bellavista_output_folder, 'exceptions.json'), 'w') as f:
            json.dump(exceptions, f)


    elif (system.lower() == 'merscope'):
        print('Creating Bella Vista input files for Vizgen MERSCOPE')
        exceptions = create_ome_zarr(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        exceptions = input_data_merscope.create_micron_pixel(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        exceptions = input_data_merscope.create_transcripts(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        exceptions = input_data_merscope.create_segmentations(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        
        with open(os.path.join(bellavista_output_folder, 'exceptions.json'), 'w') as f:
            json.dump(exceptions, f)

    else:
        print('Creating Bella Vista input files for Custom MERFISH (MERlin)')
        exceptions = create_ome_zarr(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        exceptions = input_data_merlin.create_micron_pixel(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        exceptions = input_data_merlin.create_transcripts(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        exceptions = input_data_merlin.create_segmentations(data_folder, bellavista_output_folder, json_file_input_files, exceptions)
        
        with open(os.path.join(bellavista_output_folder, 'exceptions.json'), 'w') as f:
            json.dump(exceptions, f)

    print('Bella Vista input files created!', end='\n\n')


def create_ome_zarr(data_folder: str, bellavista_output_folder: str, json_file_input_files: Dict, exceptions: Dict):
    
    '''
        Converts TIFF to multiscale pyramidal OME-Zarr images.

        Args:
            data_folder: Path to folder containing dataset files.
            bellavista_output_folder: Path to output data folder to store input files for bellavista.py
            json_file_input_files: A dictionary storing input file parameters parsed from user's input JSON file, containing image paths.
            exceptions: An empty dictionary.

        Returns:
            exceptions: An updated version of exceptions dictionary parsed if an error occured.
    '''
    try:
        # if images have been processed previously, exit early
        ome_zarr_path = os.path.join(bellavista_output_folder, 'OMEzarrImages')
        image_names_path = os.path.join(bellavista_output_folder, 'image_file_names.pkl')

        if os.path.exists(ome_zarr_path) and os.path.exists(image_names_path):
            print('OME-Zarr images have been processed previously. Skipping reprocessing.')
            return exceptions

        images = json_file_input_files.get('images')
        z_plane = json_file_input_files.get('z_plane', 0)

        if images is None:
            print('No input image files provided, skipping processing. `plot_image` will remain false until a valid file is given.')
            logging.warning('MISSING INPUT FILE: No input image files provided --> cannot process images. `plot_image` will remain false until a valid file is given.')
            exceptions['valid_image'] = False
            return exceptions
        
        print(f'Creating OME-Zarr image at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')

        # OpenGL max texture size - napari will downscale image if OME-Zarr tile exceeds this shape
        napari_size_limit = 16384
        dask.config.set({'array.slicing.split_large_chunks': False})

        def process_image(image_path: str, z_plane: int):
            '''Helper function to process single images.'''
            try:
                img_data = imread(os.path.join(data_folder, image_path))[z_plane]
            except IndexError:  # if z_plane is out of bounds
                img_data = imread(os.path.join(data_folder, image_path))
                if len(img_data) > 1:
                    print(f'z-plane {z_plane} is out of bounds, using plane 0 instead.')
                    img_data = imread(os.path.join(data_folder, image_path))[0]
                else: 
                    img_data = imread(os.path.join(data_folder, image_path))[0]
            return img_data

        def calculate_scales(img_shape):
            '''Calculate the number of scales for pyramidal images.'''
            max_dimension = max(img_shape)
            min_n_scales = math.ceil(math.log(max_dimension / napari_size_limit, 2)) + 1
            return max(min_n_scales, 1)

        if isinstance(images, str):  # single image case
            print(f'Processing single image: {images}')
            image_file_names = os.path.basename(images.replace('\\', '/')).split('.')[0]
            img_data = process_image(images, z_plane)
            n_scales = calculate_scales(img_data.shape)

            with open_ngff_zarr(store=bellavista_output_folder, dimension_separator='/', overwrite=True) as f:
                collection = f.add_collection(name='OMEzarrImages')
                collection.add_image(image_name='Images', array=to_tczyx(img_data, axes_names=('y', 'x')), n_scales=n_scales)

        elif isinstance(images, list):  # multiple images case
            print(f'Processing image list: {images}')
            image_list = []
            img_shapes = []

            for image in images:
                img_data = imread(os.path.join(data_folder, image))
                img_shape = img_data.shape
                img_shapes.append(img_shape)

                if img_shape[0] > 1:  # 3D image case
                    img_data_2d = process_image(image, z_plane)
                    image_list.append(da.expand_dims(img_data_2d, axis=0))
                else:
                    image_list.append(img_data)

            image_list = [da.expand_dims(img, axis=0) for img in image_list]
            max_dimension = max(max(shape[1] for shape in img_shapes), max(shape[2] for shape in img_shapes))
            n_scales = calculate_scales((max_dimension, max_dimension))

            image_file_names = [os.path.basename(image.replace('\\', '/')).split('.')[0] for image in images]
            with open_ngff_zarr(store=bellavista_output_folder, dimension_separator='/', overwrite=True) as f:
                collection = f.add_collection(name='OMEzarrImages')
                collection.add_image(image_name='Images', array=to_tczyx(da.concatenate(image_list), axes_names=('c', 'z', 'y', 'x')), n_scales=n_scales)

        # save image names, these will be used as layer names when loading napari
        with open(image_names_path, 'wb') as f:
            pickle.dump(image_file_names, f)

        print(f'OME-Zarr image saved successfully at {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')

    except Exception as e:
        print(f'An error occurred in create_ome_zarr. `plot_image` will remain false until the error is resolved.')
        print(f'Please check the log file for details: {os.path.join(bellavista_output_folder, "error_log.log")}', end='\n\n')
        logging.error(f'Error in create_ome_zarr: {e}', exc_info=True)
        exceptions['valid_image'] = False

    return exceptions