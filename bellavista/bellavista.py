from . import input_data
import argparse
from json import load
import importlib.resources as pkg_resources
import sys
import os
import pickle 
import numpy as np
import napari
import random
import matplotlib.pyplot as plt
from napari.utils.colormaps import Colormap
from napari.utils.colormaps import AVAILABLE_COLORMAPS
from matplotlib import colors as mpl_colors
from tqdm import tqdm

def rotation_matrix(rotate_angle):
    # convert angle from degrees to radians
    angle_radians = np.radians(rotate_angle)
    # compute cosine and sine of the angle
    cos_theta = np.cos(angle_radians)
    sin_theta = np.sin(angle_radians)
    # construct the rotation matrix
    rotation_mat = np.array([
        [cos_theta, -sin_theta],
        [sin_theta, cos_theta]
    ])
    return rotation_mat

def bellavista(
        bellavista_output_folder,
        plot_image=False,
        image_colormap=None,
        plot_transcripts=False,
        plot_allgenes=True,
        genes_visible_on_startup=False,
        selected_genes=None,
        plot_cell_seg=False,
        plot_nuclear_seg=False,
        transcript_point_size=1,
        contrast_limits=None,
        rotate_angle=0
):

    print('Loading Bella Vista:')
    txs_colors = plt.cm.gist_rainbow(np.linspace(0, 1, 1025))

    # check for any exceptions caused while creating input data
    with open(os.path.join(bellavista_output_folder, 'exceptions.json'), 'r') as f:
        exceptions_json_file = load(f)
    if exceptions_json_file.get('valid_image') is not None:
        plot_image = False
    if exceptions_json_file.get('valid_txs') is not None:
        plot_transcripts = False
    if exceptions_json_file.get('valid_cell_seg') is not None:
        plot_cell_seg = False
    if exceptions_json_file.get('valid_nuclear_seg') is not None:
        plot_nuclear_seg = False

    if not (plot_image or plot_transcripts or plot_cell_seg or plot_nuclear_seg):
        print('No data to plot. Exiting.')
        exit()

    viewer = napari.Viewer()

    if(plot_image): 

        image_file_names = pickle.load(open(os.path.join(bellavista_output_folder,'image_file_names.pkl'),'rb'))
        # load transformations to align images with transcripts & segmentations
        um_per_pixel_dict = pickle.load(open(os.path.join(bellavista_output_folder,'um_to_px_transforms.pkl'),'rb'))
        um_per_pixel_x = um_per_pixel_dict['um_per_pixel_x']
        um_per_pixel_y = um_per_pixel_dict['um_per_pixel_y']
        x_shift = um_per_pixel_dict['x_shift']
        y_shift = um_per_pixel_dict['y_shift']

        rotation_mat = rotation_matrix(rotate_angle)
        rotated_yshift = (rotation_mat[0][0] * y_shift) + (rotation_mat[0][1] * x_shift)
        rotated_xshift = (rotation_mat[1][0] * y_shift) + (rotation_mat[1][1] * x_shift)

        print('Loading image...')
        viewer.open(os.path.join(bellavista_output_folder, 'OMEzarrImages/Images'), plugin='napari-ome-zarr', scale = (1, um_per_pixel_y, um_per_pixel_x), \
                        translate = (0, rotated_yshift, rotated_xshift), blending = 'additive', colormap=image_colormap, contrast_limits = contrast_limits, \
                            name = image_file_names, channel_axis = 1, rotate=rotate_angle) #create napari image layer 
    
    if(plot_transcripts):
        gene_dict = pickle.load(open(os.path.join(bellavista_output_folder,'gene_dict.pkl'),'rb'))
        sorted_gene_names = sorted(gene_dict.keys(), reverse=True) #reverse the order of the genes so that the genes are in alphabetical order
    
        # load subset of user-specified genes
        if not plot_allgenes and selected_genes is not None and len(selected_genes) > 0:
            for gene in tqdm(selected_genes, desc = 'Loading gene transcripts...', total = len(selected_genes)):
                color = random.choice(txs_colors)
                try:
                    viewer.add_points(gene_dict[gene], size=transcript_point_size, border_width=0, name=gene, face_color=color, border_color=color, visible=genes_visible_on_startup, rotate=rotate_angle) #create napari point layer for each selected gene 
                except KeyError:
                    try:
                        lowercase_geneID = gene.lower()
                        geneID = next((key for key in gene_dict if key.lower() == lowercase_geneID), None)
                        if geneID is None:
                            raise KeyError(f'{gene} not found in dataset. Please check the spelling. Skipping {gene}')
                        txs_points = gene_dict.get(geneID)
                        viewer.add_points(txs_points, size=transcript_point_size, border_width=0, name=geneID, face_color=color, border_color=color, visible=genes_visible_on_startup, rotate=rotate_angle) #create napari point layer for each selected gene 
                    except KeyError as e:
                        print(e)
        # load all genes
        else:
            for gene in tqdm(sorted_gene_names, desc = 'Loading gene transcripts...', total = len(sorted_gene_names)):
                color = random.choice(txs_colors)
                viewer.add_points(gene_dict[gene], size=transcript_point_size, border_width=0, name=gene, face_color=color, border_color=color, visible=genes_visible_on_startup, rotate=rotate_angle) #create napari point layer for each gene 

    if(plot_cell_seg):
        print('Loading cell segmentation...')
        # layer color can be changed to valid matplotlib colors
        color = list(mpl_colors.to_rgba('white')) 
        colors = np.array([color, color])
        # create a napari Colormap for layer
        custom_cmap = Colormap(colors, name='cell_cmap')
        AVAILABLE_COLORMAPS['cell_cmap'] = custom_cmap

        coords = pickle.load(open(os.path.join(bellavista_output_folder,'cell_boundary_coords.pkl'),'rb'))
        viewer.add_tracks(coords, name='cell boundaries', colormap='cell_cmap', visible=False, blending = 'opaque', rotate=rotate_angle)

    if(plot_nuclear_seg):
        print('Loading nuclear segmentation...')
        color = list(mpl_colors.to_rgba('gray'))
        colors = np.array([color, color])
        custom_cmap = Colormap(colors, name='nuclear_cmap')
        AVAILABLE_COLORMAPS['nuclear_cmap'] = custom_cmap

        coords = pickle.load(open(os.path.join(bellavista_output_folder,'nuclear_boundary_coords.pkl'),'rb'))
        viewer.add_tracks(coords, name='nuclear boundaries', colormap='nuclear_cmap', visible=False, blending = 'opaque', rotate=rotate_angle)
    
    viewer.reset_view()
    print('Data loaded!')
    viewer.scale_bar.visible = True
    viewer.scale_bar.unit = 'um'
    napari.run()

def main():

    
    parser = argparse.ArgumentParser(description='Process input file for Bellavista.')
    # user provided JSON file
    parser.add_argument('positional_input_file', type=str, nargs='?', help='Path to the input JSON file')
    parser.add_argument('-i', '--input-file', type=str)
    # sample JSON files
    parser.add_argument('--xenium-sample', type=str)
    parser.add_argument('--xenium-sample-lite', type=str)
    args = parser.parse_args()

    # "Quick Start" tutorial sample JSONs
    if args.xenium_sample:
        sample_data_folder = args.xenium_sample
        if not sample_data_folder:
            print('Error: No sample data folder path provided. Usage: "bellavista --xenium-sample /path/to/xenium_sample_data"')
            sys.exit(1)
        with pkg_resources.open_text('bellavista.utils.quickstart', 'quickstart_xenium_sample.json') as f:
            json_file = load(f)
            json_file['data_folder'] = sample_data_folder
    
    elif args.xenium_sample_lite:
        sample_data_folder = args.xenium_sample_lite
        if not sample_data_folder:
            print('Error: No sample data folder path provided. Usage: "bellavista --xenium-sample-lite /path/to/xenium_sample_data"')
            sys.exit(1)
        with pkg_resources.open_text('bellavista.utils.quickstart', 'quickstart_xenium_sample_lite.json') as f:
            json_file = load(f)
            json_file['data_folder'] = sample_data_folder
    
    # user defined JSON
    else:
        input_file = args.input_file if args.input_file else args.positional_input_file
        if not input_file:
            print('Error: No input JSON file provided. You must provide an input file either as the first argument or with the -i/--input_file option. If following the "Quick Start" tutorial, use the "--xenium-sample" option')
            parser.print_help()
            sys.exit(1)

        # load dataset-specific JSON (first argument)
        with open(input_file, 'r') as f:
            json_file = load(f)
    
    data_folder = json_file.get('data_folder')
    data_folder = data_folder.replace('\\', '/')
    json_file_param = json_file.get('visualization_parameters')
    create_bellavista_inputs = json_file.get('create_bellavista_inputs', True)
        
    if(create_bellavista_inputs == True):
        input_data.create_inputs(json_file)

    bellavista_output_folder = os.path.join(data_folder, "BellaVista_output")

    bellavista(
        bellavista_output_folder=bellavista_output_folder,
        plot_image=json_file_param.get('plot_image', False),
        image_colormap=json_file_param.get('image_colormap'),
        plot_transcripts=json_file_param.get('plot_transcripts', False),
        plot_allgenes=json_file_param.get('plot_allgenes', True),
        genes_visible_on_startup=json_file_param.get('genes_visible_on_startup', False),
        selected_genes=json_file_param.get('selected_genes'),
        plot_cell_seg=json_file_param.get('plot_cell_seg', False),
        plot_nuclear_seg=json_file_param.get('plot_nuclear_seg', False),
        transcript_point_size=json_file_param.get('transcript_point_size',1),
        contrast_limits=json_file_param.get('contrast_limits'),
        rotate_angle=json_file_param.get('rotate_angle', 0)
    )

if __name__ == '__main__':
    main()