#Author: Annabelle Coles
import pickle
import pandas as pd
import os
import numpy as np
from tqdm import tqdm
import tifffile
import zarr
import re
import xml.etree.ElementTree as ET
from typing import Dict
import logging

def get_namespace(element):
    match = re.match(r'\{.*\}', element.tag)
    return match.group(0) if match else ''

def create_micron_pixel(data_folder: str, bellavista_output_folder: str, json_file_input_files: Dict, exceptions: Dict):

    '''
        Calculates scale factors and translations to align images to transcripts. Transformations are stored as a pickled dictionary. 
        This will be used when loading the image in napari. 

        Args:
            data_folder: Path to folder containing dataset files.
            bellavista_output_folder: Path to output data folder to store input files for bellavista.py
            json_file_input_files: A dictionary storing input file parameters parsed from user's input JSON file, containing image paths.
            exceptions: A dictionary documenting errors. 

        Returns:
            exceptions: An updated version of exceptions dictionary parsed if an error occured.
    '''

    # if transformations have been processed previously, exit
    if os.path.exists(os.path.join(bellavista_output_folder, 'um_to_px_transforms.pkl')):
        print("Micron to pixel transforms have been processed previously. Skipping reprocessing.")
        return exceptions
    
    images = json_file_input_files.get('images')
    
    if images is None: 
        print('No images provided')
        exceptions['valid_image'] = False
        return exceptions

    try:
        print('Calculating micron to pixel transform')
        if isinstance(images, list):
            image = images[0]
        else:
            image = images

        # Access image metadata
        with tifffile.TiffFile(os.path.join(data_folder, image)) as tif:
            ome_metadata = tif.ome_metadata

        # Extract transformations needed for image-transcript alignment
        root = ET.fromstring(ome_metadata)
        namespace = get_namespace(root)
        ns = {'ome': namespace[1:-1]}
        plate = root.find('ome:Plate', ns)

        well_origin_x = float(plate.attrib['WellOriginX'])
        well_origin_y = float(plate.attrib['WellOriginY'])

        image = root.find('ome:Image', ns)

        pixels = image.find('ome:Pixels', ns)
        pixels_size_x = float(pixels.attrib['PhysicalSizeX'])
        pixels_size_y = float(pixels.attrib['PhysicalSizeY'])

        um_to_px_transform_dict = {'um_per_pixel_x': pixels_size_x, 'um_per_pixel_y': pixels_size_y,
                                'x_shift': well_origin_x, 'y_shift': well_origin_y}

        with open(os.path.join(bellavista_output_folder, 'um_to_px_transforms.pkl'), 'wb') as f:
            pickle.dump(um_to_px_transform_dict, f)

        print('Micron to pixel transform saved!')

    except Exception as e: 
        # Log the exception with traceback
        print(f'An error occurred in create_micron_pixel. Please check the log file for details: {os.path.join(bellavista_output_folder, 'error_log.log')}', end='\n\n')
        logging.error(f'Error in create_micron_pixel: {e}', exc_info=True)
        # update exceptions dictionary to document an error
        exceptions['valid_image'] = False
        
    return exceptions


def create_transcripts(data_folder: str, bellavista_output_folder: str, json_file_input_files: Dict, exceptions: Dict):

    '''
        Saves pickled dictionary mapping transcript coordinates for each gene. This will be used when loading transcripts in napari. 

        Args:
            data_folder: Path to folder containing dataset files.
            bellavista_output_folder: Path to output data folder to store input files for bellavista.py
            json_file_input_files: A dictionary storing input file parameters parsed from user's input JSON file, containing path to transcript file.
            exceptions: A dictionary documenting errors. 

        Returns:
            exceptions: An updated version of exceptions dictionary parsed if an error occured.
    '''
    # if transcripts have been processed previously, exit
    if os.path.exists(os.path.join(bellavista_output_folder, 'gene_dict.pkl')):
        print("Transcripts have been processed previously. Skipping reprocessing.")
        return exceptions
    
    transcript_filename = json_file_input_files.get('transcript_filename')

    if transcript_filename is None: 
        print('No transcript file provided')
        exceptions['valid_txs'] = False
        return exceptions

    try:
        print('Creating transcripts')
        if transcript_filename.endswith('.csv.gz'):
            txs_locations_df = pd.read_csv(os.path.join(data_folder, transcript_filename), delimiter=',',
                                        compression='gzip')
        elif (transcript_filename.endswith('.csv')):
            txs_locations_df = pd.read_csv(os.path.join(data_folder, transcript_filename), delimiter=',')
        elif (transcript_filename.endswith('.parquet')):
            txs_locations_df = pd.read_parquet(os.path.join(data_folder, transcript_filename))
        else:
            error_message = f'Error in create_transcripts: "{transcript_filename}" is not a valid file type. Must be a CSV or Parquet file'
            print(error_message, end='\n\n')
            logging.error(error_message)
            exceptions['valid_txs'] = False
            return exceptions

        # convert gene names from byte to str datatype
        txs_locations_df['feature_name'] = txs_locations_df['feature_name'].apply(lambda x: x.decode('utf-8') if isinstance(x, bytes) else x)
        
        # Map all transcript locations to genes in a dictionary
        gene_dict = {}
        for gene, group in tqdm(txs_locations_df.groupby('feature_name'),
                                desc=f'Processing genes',
                                total=txs_locations_df['feature_name'].nunique()):
            coords = group[['y_location', 'x_location']].values
            gene_dict[gene] = coords

        with open(os.path.join(bellavista_output_folder, 'gene_dict.pkl'), 'wb') as f:
            pickle.dump(gene_dict, f)

        print('Transcripts saved!')
        
    except Exception as e: 
        # Log the exception with traceback
        print(f'An error occurred in create_transcripts. Please check the log file for details: {os.path.join(bellavista_output_folder, 'error_log.log')}', end='\n\n')
        logging.error(f'Error in create_transcripts: {e}', exc_info=True)
        exceptions['valid_txs'] = False
    return exceptions

def convert_to_tuples(coords: np.ndarray, counter: int) -> np.ndarray:
    coords = [(coords[i], coords[i+1]) for i in range(0, len(coords), 2)]
    flipped_coords = np.asarray(coords)[:, ::-1]
    coords = np.asarray([tuple(row) for row in flipped_coords])
    zeros_array = np.tile([counter, 0], (coords.shape[0], 1))
    result_array = np.hstack((zeros_array, coords))
    
    return result_array


def process_segmentations(data_folder: str, bellavista_output_folder: str, segmentation_file: str, seg_type: str):

    '''
        Extracts cell or nuclear boundaries.

        Args:
            data_folder: Path to folder containing dataset files.
            bellavista_output_folder: Path to output data folder to store input files for bellavista.py
            segmentation_file: Path to segmentation data, must be Parquet or CSV file.
            seg_type: 'cell' or 'nuclear' depending on type of segmentation. 

        Returns:
            False if an error occured.
    '''
    if segmentation_file is None:
        print(f'No {seg_type} segmentation file provided')
        return False

    try:
        print(f'Creating {seg_type} boundaries')
        file_type = None
        if (segmentation_file.endswith('.parquet')):
            cell_df = pd.read_parquet(os.path.join(data_folder, segmentation_file))
            file_type = 'parquet'
        elif (segmentation_file.endswith('.csv.gz')):
            cell_df = pd.read_csv(os.path.join(data_folder, segmentation_file), delimiter=',', compression='gzip')
            file_type = 'csv'
        elif (segmentation_file.endswith('.csv')):
            cell_df = pd.read_csv(os.path.join(data_folder, segmentation_file), delimiter=',')
            file_type = 'csv'

        if (file_type == 'parquet' or file_type == 'csv'):
            allcellbounds_array = []

            for i, (cell_id, group) in tqdm(enumerate(cell_df.groupby('cell_id')),
                                            desc=f'Processing {seg_type} boundaries',
                                            total=cell_df['cell_id'].nunique()):
                coords = group[['vertex_y', 'vertex_x']].values
                # create 4D numpy array storing segmentation coordinates
                zeros_array = np.column_stack((np.full(coords.shape[0], i), np.zeros(coords.shape[0])))
                allcellbounds_array.extend(np.hstack((zeros_array, coords)))

        else:
            file_type = 'zarr'
            # indexing based on Xenium documentation as of Aug 2024: https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/advanced/xoa-output-zarr
            if seg_type == 'cell':
                seg_index = 1
            elif seg_type == 'nuclear':
                seg_index = 0

            cells = zarr.open(os.path.join(data_folder, segmentation_file), mode='r')[seg_index]
            # pairs x,y coordinates from x,y-1D arrays
            coords = [convert_to_tuples(coords, idx) for idx, coords in enumerate(cells, start=0)]
            allcellbounds_array = [item for sublist in coords for item in sublist]
        
        with open(os.path.join(bellavista_output_folder, f'{seg_type}_boundary_coords.pkl'), 'wb') as f:
            pickle.dump(allcellbounds_array, f)

        print(f'{seg_type} boundaries saved!')
        return True

    except Exception as e: 
        # Log the exception with traceback
        print(f'An error occurred in process_segmentations. Please check the log file for details: {os.path.join(bellavista_output_folder, 'error_log.log')}', end='\n\n')
        logging.error(f'Error in process_segmentations ({seg_type}) {e}', exc_info=True)
        return False 

def create_segmentations(data_folder: str, bellavista_output_folder: str, json_file_input_files: Dict, exceptions: Dict):
    
    '''
        Saves extracted cell and nuclear boundaries as a pickle. This will be used when loading segmentations in napari. 

        Args:
            data_folder: Path to folder containing dataset files.
            bellavista_output_folder: Path to output data folder to store input files for bellavista.py
            json_file_input_files: A dictionary storing input file parameters parsed from user's input JSON file, containng segmentation paths.
            exceptions: An empty dictionary.

        Returns:
            exceptions: An updated version of exceptions dictionary parsed if an error occured.
    '''

    cell_segmentation_file = json_file_input_files.get('cell_segmentation')
    nuclear_segmentation_file = json_file_input_files.get('nuclear_segmentation')
    
    # process segmentations
    if os.path.exists(os.path.join(bellavista_output_folder, f'cell_boundary_coords.pkl')):
        print(f"Cell segmentations have been processed previously. Skipping reprocessing.")
    else:
        if not process_segmentations(data_folder, bellavista_output_folder, segmentation_file=cell_segmentation_file,
                                seg_type='cell'):
            # update exceptions dictionary to document an error
            exceptions[f'valid_cell_seg'] = False

    if os.path.exists(os.path.join(bellavista_output_folder, f'nuclear_boundary_coords.pkl')):
        print(f"Nuclear segmentations have been processed previously. Skipping reprocessing.")
    else:
        if not process_segmentations(data_folder, bellavista_output_folder, segmentation_file=nuclear_segmentation_file,
                                    seg_type='nuclear'):
            # update exceptions dictionary to document an error
            exceptions[f'valid_nuclear_seg'] = False
    return exceptions