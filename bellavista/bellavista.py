from json import load
import sys
import os
import pickle 
import numpy as np
import napari
import random
import matplotlib.pyplot as plt
import input_data
from napari.utils.colormaps import Colormap
from napari.utils.colormaps import AVAILABLE_COLORMAPS
from matplotlib import colors as mpl_colors
from tqdm import tqdm

#test git in vs code

def bellavista(
        bella_vista_output_folder,
        plot_image=False,
        plot_transcripts=False,
        plot_allgenes=True,
        selected_genes=None,
        plot_cell_seg=False,
        plot_nuclear_seg=False,
        transcript_point_size=None,
        contrast_limits=None
):

    print("Loading Bella Vista")
    txs_colors = plt.cm.gist_rainbow(np.linspace(0, 1, 1025))

    # check for any exceptions caused while creating input data
    with open(os.path.join(bella_vista_output_folder, "exceptions.json"), 'r') as f:
        exceptions_json_file = load(f)
    if exceptions_json_file.get("valid_image") is not None:
        plot_image = False
    if exceptions_json_file.get("valid_txs") is not None:
        plot_transcripts = False
    if exceptions_json_file.get("valid_cell_seg") is not None:
        plot_cell_seg = False
    if exceptions_json_file.get("valid_nuclear_seg") is not None:
        plot_nuclear_seg = False


    with napari.gui_qt():
        viewer = napari.Viewer()

        if(plot_image): 
            image_file_names = pickle.load(open(os.path.join(bella_vista_output_folder,"image_file_names.pkl"),"rb"))
            # load transformations to align images with transcripts & segmentations
            um_per_pixel_dict = pickle.load(open(os.path.join(bella_vista_output_folder,"um_to_px_transforms.pkl"),"rb"))
            um_per_pixel_x = um_per_pixel_dict['um_per_pixel_x']
            um_per_pixel_y = um_per_pixel_dict['um_per_pixel_y']
            x_shift = um_per_pixel_dict['x_shift']
            y_shift = um_per_pixel_dict['y_shift']

            print("Loading image...")
            if contrast_limits is None:
                viewer.open(os.path.join(bella_vista_output_folder, "OMEzarrImages/Images"), plugin='napari-ome-zarr', scale = (1, um_per_pixel_y, um_per_pixel_x), \
                            translate = (0, y_shift, x_shift), blending = 'additive', name = image_file_names, channel_axis= 1) #create napari image layer 
            else:
                viewer.open(os.path.join(bella_vista_output_folder, "OMEzarrImages/Images"), plugin='napari-ome-zarr', scale = (1, um_per_pixel_y, um_per_pixel_x), \
                            translate = (0, y_shift, x_shift), blending = 'additive', contrast_limits = contrast_limits, name = image_file_names, channel_axis= 1) #create napari image layer 

        if(plot_transcripts):
            if transcript_point_size is None:
                transcript_point_size = 1.0
            else: 
                transcript_point_size = float(transcript_point_size)
            gene_dict = pickle.load(open(os.path.join(bella_vista_output_folder,"gene_dict.pkl"),"rb"))
            sorted_gene_names = sorted(gene_dict.keys(), reverse=True) #reverse the order of the genes so that the genes are in alphabetical order

            # load subset of user-specified genes
            if not plot_allgenes and selected_genes is not None and len(selected_genes) > 0:
                for gene in tqdm(selected_genes, desc = 'Loading gene transcripts...', total = len(selected_genes)):
                    color = random.choice(txs_colors)
                    viewer.add_points(gene_dict[gene], size=transcript_point_size, border_width=0, name=gene, face_color=color, border_color=color, visible=True) #create napari point layer for each selected gene 
            # load all genes
            else:
                for gene in tqdm(sorted_gene_names, desc = 'Loading gene transcripts...', total = len(sorted_gene_names)):
                    color = random.choice(txs_colors)
                    viewer.add_points(gene_dict[gene], size=transcript_point_size, border_width=0, name=gene, face_color=color, border_color=color, visible=True) #create napari point layer for each gene 

        if(plot_cell_seg):
            print("Loading cell segmentation...")
            # layer color can be changed to valid matplotlib colors
            color = list(mpl_colors.to_rgba('white')) 
            colors = np.array([color, color])
            # create a napari Colormap for layer
            custom_cmap = Colormap(colors, name='cell_cmap')
            AVAILABLE_COLORMAPS['cell_cmap'] = custom_cmap

            coords = pickle.load(open(os.path.join(bella_vista_output_folder,"cell_boundary_coords.pkl"),"rb"))
            viewer.add_tracks(coords, name='cell boundaries', colormap='cell_cmap', visible=True, blending = 'opaque')

        if(plot_nuclear_seg):
            print("Loading nuclear segmentation...")
            color = list(mpl_colors.to_rgba('gray'))
            colors = np.array([color, color])
            custom_cmap = Colormap(colors, name='nuclear_cmap')
            AVAILABLE_COLORMAPS['nuclear_cmap'] = custom_cmap

            coords = pickle.load(open(os.path.join(bella_vista_output_folder,"nuclear_boundary_coords.pkl"),"rb"))
            viewer.add_tracks(coords, name='nuclear boundaries', colormap='nuclear_cmap', visible=True, blending = 'opaque')
        
        print("Data loaded!")
        viewer.scale_bar.visible = True
        viewer.scale_bar.unit = 'um'
        viewer.camera.zoom = 0.6
        napari.run()

def main():

    # load dataset-specific JSON (first argument)
    with open(sys.argv[1], 'r') as f:
        json_file = load(f)
        json_file_param = json_file.get("parameters")
        create_bellavista_inputs = json_file.get("create_bellavista_inputs", False)
    
    if(create_bellavista_inputs == True):
        input_data.create_inputs(json_file)

    bellavista(
        bella_vista_output_folder=json_file.get("bella_vista_output_folder"),
        plot_image=json_file_param.get("plot_image"),
        plot_transcripts=json_file_param.get("plot_transcripts"),
        plot_allgenes=json_file_param.get("plot_allgenes"),
        selected_genes=json_file_param.get("selected_genes"),
        plot_cell_seg=json_file_param.get("plot_cell_seg"),
        plot_nuclear_seg=json_file_param.get("plot_nuclear_seg"),
        transcript_point_size=json_file_param.get("transcript_point_size"),
        contrast_limits=json_file_param.get("contrast_limits")
    )

if __name__ == '__main__':
    main()