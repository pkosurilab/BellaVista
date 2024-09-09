#Author: Annabelle Coles
import pickle
import pandas as pd
import os
import numpy as np
import h5py
from tqdm import tqdm
import shapely
from typing import Dict
import logging

def create_micron_pixel(data_folder: str, bellavista_output_folder: str, json_file_input_files: Dict, exceptions: Dict):

    """
        Extracts scale factors and translations to align images to transcripts. Transformations are stored as a pickled dictionary. 
        This will be used when loading the image in napari. 

        Args:
            data_folder: Path to folder containing dataset files.
            bellavista_output_folder: Path to output data folder to store input files for bellavista.py
            json_file_input_files: A dictionary storing input file parameters parsed from user's input JSON file, containing path to transformation file.
            exceptions: A dictionary documenting errors. 

        Returns:
            exceptions: An updated version of exceptions dictionary parsed if an error occured.
    """

    # if transformations have been processed previously, exit
    if os.path.exists(os.path.join(bellavista_output_folder, "um_to_px_transforms.pkl")):
        print("Micron to pixel transforms have been processed previously. Skipping reprocessing.")
        return exceptions
    
    um_to_px_transform = json_file_input_files.get("um_to_px_transform")
    
    if um_to_px_transform is None: 
        print("No um_to_px_transform file provided")
        exceptions['valid_image'] = False
        return exceptions
    
    try:
        print('Calculating micron to pixel transform')
        um_to_px_transform = pd.read_csv(os.path.join(data_folder, um_to_px_transform), header=None, delimiter=r"\s+")

        scaling_x = float(1/um_to_px_transform.iloc[0][0])
        scaling_y = float(1/um_to_px_transform.iloc[1][1])
        translation_x = - (float(um_to_px_transform.iloc[0][2]) * scaling_x)
        translation_y = - (float(um_to_px_transform.iloc[1][2]) * scaling_y)

        um_to_px_transform_dict = {'um_per_pixel_x': scaling_x, 'um_per_pixel_y': scaling_y, 'x_shift': translation_x, 'y_shift': translation_y}

        with open(os.path.join(bellavista_output_folder, 'um_to_px_transforms.pkl'), 'wb') as f:
            pickle.dump(um_to_px_transform_dict, f)
            
        print('Micron to pixel transform saved!')
    
    except Exception as e: 
        # Log the exception with traceback
        print(f'An error occurred in create_micron_pixel. Please check the log file for details: {os.path.join(bellavista_output_folder, 'error_log.log')}', end='\n\n')
        logging.error(f"Error in create_micron_pixel: {e}", exc_info=True)
        exceptions['valid_image'] = False

    return exceptions

def create_transcripts(data_folder: str, bellavista_output_folder: str, json_file_input_files: Dict, exceptions: Dict):

    """
        Saves pickled dictionary mapping transcript coordinates for each gene. This will be used when loading transcripts in napari. 

        Args:
            data_folder: Path to folder containing dataset files.
            bellavista_output_folder: Path to output data folder to store input files for bellavista.py
            json_file_input_files: A dictionary storing input file parameters parsed from user's input JSON file, containg path to transcript file.
            exceptions: A dictionary documenting errors. 

        Returns:
            exceptions: An updated version of exceptions dictionary parsed if an error occured.
    """

    # if transcripts have been processed previously, exit
    if os.path.exists(os.path.join(bellavista_output_folder, "gene_dict.pkl")):
        print("Transcripts have been processed previously. Skipping reprocessing.")
        return exceptions
    
    transcript_filename = json_file_input_files.get("transcript_filename")
    
    if transcript_filename is None: 
        print("No transcript file provided")
        exceptions['valid_txs'] = False
        return exceptions

    try:
        print('Creating transcripts')
        txs_locations_df = pd.read_csv(os.path.join(data_folder, transcript_filename), delimiter=',')

        # Map all transcript locations to genes in a dictionary
        gene_dict = {}
        for gene, group in tqdm(txs_locations_df.groupby('gene'), 
                                                    desc=f'Processing genes', 
                                                    total=txs_locations_df['gene'].nunique()):
            coords = group[['global_y', 'global_x']].values
            gene_dict[gene] = coords

        with open(os.path.join(bellavista_output_folder, 'gene_dict.pkl'), 'wb') as f:
            pickle.dump(gene_dict, f)

        print("Transcripts saved!")

    except Exception as e: 
        # Log the exception with traceback
        print(f'An error occurred in create_transcripts. Please check the log file for details: {os.path.join(bellavista_output_folder, 'error_log.log')}', end='\n\n')
        logging.error(f"Error in create_transcripts: {e}", exc_info=True)
        exceptions['valid_txs'] = False
    return exceptions

def process_segmentations(data_folder: str, bellavista_output_folder: str, segmentation_file: str, seg_type: str, z_plane: int):
     
    """
        Extracts cell or nuclear boundaries.

        Args:
            data_folder: Path to folder containing dataset files.
            bellavista_output_folder: Path to output data folder to store input files for bellavista.py
            segmentation_file: Path to segmentation data, must be Parquet or folder containing HDF5 files.
            seg_type: "cell" or "nuclear" depending on type of segmentation. 
            z_plane: Integer index to extract segmentations from.

        Returns:
            False if an error occured.
    """
    
    if segmentation_file is None:
        print(f"No {seg_type} segmentation file provided")
        return False
    
    try: 
        print(f'Creating {seg_type} boundaries')
        if(segmentation_file.endswith(".parquet")): 
            print(f'Processing boundaries from Parquet file')
            cell_df = pd.read_parquet(os.path.join(data_folder, segmentation_file))
            cell_df_z_plane = cell_df[cell_df['ZIndex'] == int(z_plane)]
            counter = 0
            allcellbounds_array = []
            for ind, row in tqdm(cell_df_z_plane.iterrows(), desc=f'Processing {seg_type} boundaries', total=len(cell_df_z_plane)):
                # Convert binary to Shapely geometry object, extract boundary coordinates of geometry
                multi_poly = shapely.from_wkb(row['Geometry'])
                for poly in multi_poly.geoms:
                    coords = np.asarray(poly.exterior.coords)
                    # create 4D numpy array storing segmentation coordinates
                    zeros_array = np.tile([counter, 0], (coords.shape[0], 1))
                    result_array = np.hstack((zeros_array, coords))
                    allcellbounds_array.extend(result_array)
                    counter += 1

        # segmentations are stored in HDF5 files
        else: 
            print(f'Processing boundaries from HDF5 files')
            allcellbounds_array = []
            files = os.listdir(os.path.join(data_folder, segmentation_file))
            counter = 0
            # extract boundaries from each file
            for file in tqdm(files, desc=f'processing {seg_type} boundaries', total=len(files)):
                if file.endswith(".hdf5"):
                    with h5py.File(os.path.join(data_folder, segmentation_file, file), "r") as f:
                    
                        a_group_key = list(f.keys())[0]
                        cells = list(f[a_group_key])
                        for i, key in enumerate(cells):
                            try:
                                zGroup = f[a_group_key][key]['zIndex_' + str(z_plane)]
                                for geom in zGroup:
                                    geom_coords = list(zGroup[geom]['coordinates'])[0]
                                    # Flip x and y values to match napari coordinate system
                                    flipped_coords = geom_coords[:, ::-1]
                                    # Convert back to a list of tuples
                                    geom_coords = np.asarray([tuple(row) for row in flipped_coords])
                                    # create 4D numpy array storing segmentation coordinates
                                    zeros_array = np.tile([counter, 0], (geom_coords.shape[0], 1))
                                    result_array = np.hstack((zeros_array, geom_coords))
                                    allcellbounds_array.extend(result_array)
                                    counter += 1
                            except: pass 

        with open(os.path.join(bellavista_output_folder, f'{seg_type}_boundary_coords.pkl'), 'wb') as f:
            pickle.dump(allcellbounds_array, f)

        print(f"{seg_type} boundaries saved!")
        return True
        
    except Exception as e: 
        # Log the exception with traceback
        print(f'An error occurred in process_segmentations. Please check the log file for details: {os.path.join(bellavista_output_folder, 'error_log.log')}', end='\n\n')
        logging.error(f"Error in process_segmentations ({seg_type}): {e}", exc_info=True)
        print()
        return False

def create_segmentations(data_folder: str, bellavista_output_folder: str, json_file_input_files: Dict, exceptions: Dict):

    """
        Saves extracted cell and nuclear boundaries as a pickle. This will be used when loading segmentations in napari. 

        Args:
            data_folder: Path to folder containing dataset files.
            bellavista_output_folder: Path to output data folder to store input files for bellavista.py
            json_file_input_files: A dictionary storing input file parameters parsed from user's input JSON file, containing segmentation paths.
            exceptions: An empty dictionary.

        Returns:
            exceptions: An updated version of exceptions dictionary parsed if an error occured.
    """

    cell_segmentation_file = json_file_input_files.get("cell_segmentation")
    nuclear_segmentation_file = json_file_input_files.get("nuclear_segmentation")
    z_plane = json_file_input_files.get("z_plane", 0)

    # process segmentations
    if os.path.exists(os.path.join(bellavista_output_folder, f'cell_boundary_coords.pkl')):
        print(f"Cell segmentations have been processed previously. Skipping reprocessing.")
    else:
        if not process_segmentations(data_folder, bellavista_output_folder, segmentation_file=cell_segmentation_file,
                            seg_type="cell", z_plane=z_plane):
            # update exceptions dictionary to document an error
            exceptions[f'valid_cell_seg'] = False

    if os.path.exists(os.path.join(bellavista_output_folder, f'nuclear_boundary_coords.pkl')):
        print(f"Nuclear segmentations have been processed previously. Skipping reprocessing.")
    else:
        if not process_segmentations(data_folder, bellavista_output_folder, segmentation_file=nuclear_segmentation_file,
                              seg_type="nuclear", z_plane=z_plane):
            # update exceptions dictionary to document an error
            exceptions[f'valid_nuclear_seg'] = False

    return exceptions