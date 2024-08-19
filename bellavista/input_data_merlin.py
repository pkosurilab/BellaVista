#Author: Annabelle Coles
import pickle
import pandas as pd
import os
import numpy as np
import h5py
from tqdm import tqdm
from typing import Dict
from json import load


def create_micron_pixel(data_folder: str, bellavista_output_folder: str, json_file_input_files: Dict, exceptions: Dict):

    """
        Extracts scale factors and calculates translations to align images to transcripts. Transformations are stored as a pickled dictionary. 
        This will be used when loading the image in napari. 

        Args:
            data_folder: Path to folder containing dataset files.
            bellavista_output_folder: Path to output data folder to store input files for bellavista.py
            json_file_input_files: A dictionary storing input file parameters parsed from user's input JSON file, contianing microscope parameter and position list paths.
            exceptions: A dictionary documenting errors. 

        Returns:
            exceptions: An updated version of exceptions dictionary parsed if an error occured.
    """

    # if transformations have been processed previously, exit
    if os.path.exists(os.path.join(bellavista_output_folder, "um_to_px_transforms.pkl")):
        return exceptions
    try:
        print('Calculating micron to pixel transform')

        microscope_parameters = json_file_input_files.get("microscope_parameters")
        positions_list = json_file_input_files.get("positions_list")
        with open(os.path.join(data_folder, microscope_parameters), 'r') as f:
            microscope_parameters = load(f)

        um_per_pixel = float(microscope_parameters.get("microns_per_pixel"))
        # find coordinates of left-most, top coordinate of the imageing area 
        positions = pd.read_csv(os.path.join(data_folder, positions_list), header=None)
        minx, miny = float(min(positions[0])), float(min(positions[1]))
        um_to_px_transform_dict = {'um_per_pixel_x': um_per_pixel, 'um_per_pixel_y': um_per_pixel, 'x_shift': minx, 'y_shift': miny}
        f = open(os.path.join(bellavista_output_folder, "um_to_px_transforms.pkl"),"wb")
        pickle.dump(um_to_px_transform_dict, f)
        f.close()
        print('Micron to pixel transform saved!')

    except:
        # update exceptions dictionary to document an error
        if not os.path.exists(os.path.join(bellavista_output_folder, "um_to_px_transforms.pkl")):
            exceptions['valid_image'] = False
    return exceptions


def create_transcripts(data_folder: str, bellavista_output_folder: str, json_file_input_files: Dict, exceptions: Dict):

    """
        Saves pickled dictionary mapping transcript coordinates for each gene. This will be used when loading transcripts in napari. 

        Args:
            data_folder: Path to folder containing dataset files.
            bellavista_output_folder: Path to output data folder to store input files for bellavista.py
            json_file_input_files: A dictionary storing input file parameters parsed from user's input JSON file containg path to transcript file.
            exceptions: A dictionary documenting errors. 

        Returns:
            exceptions: An updated version of exceptions dictionary parsed if an error occured.
    """
    
    # if transcripts have been processed previously, exit
    if os.path.exists(os.path.join(bellavista_output_folder, "gene_dict.pkl")):
        return exceptions
    try:
        print('Creating transcripts')
        transcript_filename = json_file_input_files.get("transcript_filename")
        codebook = json_file_input_files.get("codebook")
        txs_locations_df = pd.read_csv(os.path.join(data_folder, transcript_filename), delimiter=',')
        #map codebook to barcode ids 
        codebook_df = pd.read_csv(os.path.join(data_folder, codebook), delimiter=',')
        name_map = codebook_df['name'].to_dict() #map barcode_id to gene IDs 
        txs_locations_df['gene'] = txs_locations_df['barcode_id'].map(name_map)

        gene_dict = {}
        for gene, group in tqdm(txs_locations_df.groupby('gene'), 
                                                    desc=f'Processing genes', 
                                                    total=txs_locations_df['gene'].nunique()):
            coords = group[['global_y', 'global_x']].values
            gene_dict[gene] = coords

        f = open(os.path.join(bellavista_output_folder, "gene_dict.pkl"),"wb")
        pickle.dump(gene_dict, f)
        f.close()
        print("Transcripts saved!")

    except:
        # update exceptions dictionary to document an error
        if not os.path.exists(os.path.join(bellavista_output_folder, "gene_dict.pkl")):
            exceptions['valid_txs'] = False
    return exceptions

def process_segmentations(data_folder: str, bellavista_output_folder: str, segmentation_folder: str, seg_type: str, z_plane: int):
    """
        Extracts cell or nuclear boundaries.

        Args:
            data_folder: Path to folder containing dataset files.
            bellavista_output_folder: Path to output data folder to store input files for bellavista.py
            segmentation_folder: Path to folder containing HDF5 files.
            seg_type: "cell" or "nuclear" depending on type of segmentation. 
            z_plane: Integer index to extract segmentations from.

        Returns:
            False if an error occured.
    """

    try:
        print(f'Creating {seg_type} boundaries')
        allcellbounds_array = []
        files = os.listdir(os.path.join(data_folder, segmentation_folder))
        counter = 0
        for file in tqdm(files, desc=f'processing {seg_type} boundaries', total=len(files)):
            if file.endswith(".hdf5"):
                with h5py.File(os.path.join(data_folder, segmentation_folder, file), "r") as f:
                
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
                                zeros_array = np.tile([counter, 0], (geom_coords.shape[0], 1))
                                result_array = np.hstack((zeros_array, geom_coords))
                                allcellbounds_array.extend(result_array)
                                counter += 1
                        except: pass 

        f = open(os.path.join(bellavista_output_folder,f"{seg_type}_boundary_coords.pkl"),"wb")
        pickle.dump(allcellbounds_array, f)
        f.close()

        print(f"{seg_type} boundaries saved!")
    except:
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
    # if segmentations have been processed previously, exit
    if os.path.exists(os.path.join(bellavista_output_folder, f"cell_boundary_coords.pkl")):
        return exceptions
    if os.path.exists(os.path.join(bellavista_output_folder, f"nuclear_boundary_coords.pkl")):
        return exceptions
    
    cell_segmentation_folder = json_file_input_files.get("cell_segmentation")
    nuclear_segmentation_folder = json_file_input_files.get("nuclear_segmentation")
    z_plane = json_file_input_files.get("z_plane", 0)

    if not process_segmentations(data_folder, bellavista_output_folder, segmentation_folder=cell_segmentation_folder,
                            seg_type="cell", z_plane=z_plane):
        # update exceptions dictionary to document an error
        if not os.path.exists(os.path.join(bellavista_output_folder, f"cell_boundary_coords.pkl")):
            exceptions[f'valid_cell_seg'] = False

    if not process_segmentations(data_folder, bellavista_output_folder, segmentation_folder=nuclear_segmentation_folder,
                              seg_type="nuclear", z_plane=z_plane):
        # update exceptions dictionary to document an error
        if not os.path.exists(os.path.join(bellavista_output_folder, f"nuclear_boundary_coords.pkl")):
            exceptions[f'valid_nuclear_seg'] = False
    return exceptions
