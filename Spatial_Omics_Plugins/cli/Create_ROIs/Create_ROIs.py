"""

Creating spatial --omics ROIs based on selection

"""

import os
import sys
sys.path.append('..')
import json

import numpy as np
import pandas as pd

import girder_client
from ctk_cli import CLIArgumentParser
from shapely.geometry import Polygon, Point, shape
import geopandas as gpd
from wsi_annotations_kit import wsi_annotations_kit as wak
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

from tqdm import tqdm

import anndata as ad

from io import BytesIO
from typing_extensions import Union
import requests


class SpatAnnotation:
    """
    gc = girder_client.GirderClient authenticated object
    omics_item_id = itemId in Girder for file containing omics information
    image_item_id = itemId for image to post SpatAnnotations to
    roi_name = name of ROIs to use (e.g. Spots, Cells, etc.)
    mpp = microns-per-pixel for slide
    definitions_file = csv file with columns Main_Types, Sub_Types, and Cell_States (not mandatory) for grouping together cell subtypes
    omics_args = additional arguments for genes to include in annotations if genes are the main omics 
    roi_shape_args = either exterior coordinates or radius (microns) to use to create rois
    
    """
    def __init__(
            self,
            gc,
            omics_item_id:str,
            image_item_id:str,
            roi_name: str,
            mpp:Union[float,None],
            definitions_file = None,
            omics_args = None,
            roi_shape_args = None
    ):
        self.gc = gc
        self.omics_item_id = omics_item_id
        self.image_item_id = image_item_id
        self.omics_args = omics_args
        self.definitions_file = definitions_file

        self.roi_name = roi_name
        self.mpp = mpp

        self.roi_shape_args = roi_shape_args

        self.user_token = self.gc.get('token/session')['token']

        if not self.definitions_file is None:
            # Creating main_cell_types, cell_subtypes, states, etc. table
            self.definitions = pd.read_csv(
                BytesIO(
                    requests.get(
                        f'{self.gc.urlBase}/item/{self.definitions_file}/download?token={self.user_token}'
                    ).content
                )
            )
        else:
            self.definitions = None

        # Getting format of omics file
        omics_item_info = self.gc.get(f'/item/{self.omics_item_id}')
        omics_item_name = omics_item_info['name']
        file_extension = omics_item_name.split('.')[-1]

        if file_extension == 'h5ad': 
            self.process_h5ad(omics_item_name)
        elif file_extension == 'rds': 
            self.process_rds()
        elif file_extension=='csv':
            self.omics = pd.read_csv(omics_item_name)
        

        if not self.definitions is None:
            # Processing omics data according to main-sub-state breakdown provided
            processed_omics = self.process_omics()
        else:
            if self.omics_args is None:
                processed_omics = {'Gene Counts': self.omics}
            elif 'columns' in self.omics_args:
                include_columns = [i.strip() for i in self.omics_args['columns'].split(',')]
                processed_omics = {
                    i: self.omics[i]
                    for i in include_columns
                }
            else:
                processed_omics = {'Cell_Subtypes': self.omics}

        if self.mpp is None and not 'geojson' in self.roi_shape_args:
            self.find_mpp()

        spat_annotations = self.create_shapes(processed_omics)

        self.save(spat_annotations)

    def find_mpp(self):
        # Finding minimum distance spots from first spot and using that to determine MPP
        # spot centroids are 100um apart and spot diameters are 55um
        spot_x_coords = self.coordinates['imagecol'].tolist()
        spot_y_coords = self.coordinates['imagerow'].tolist()

        # Test spot is the first one
        test_spot = [spot_x_coords[0],spot_y_coords[0]]

        # Spot distances
        spot_dists = np.array([self.distance(test_spot, [x, y]) for x, y in zip(spot_x_coords, spot_y_coords)])
        spot_dists = np.unique(spot_dists[spot_dists > 0])
        min_spot_dist = np.min(spot_dists)

        # Minimum distance between the test spot and another spot = 100um (same as doing 100/min_spot_dist)
        self.mpp = 1/(min_spot_dist/100)

    def distance(self, point1, point2):
        # Distance between 2 points
        return (((point1[0]-point2[0])**2)+((point1[1]-point2[1])**2))**0.5

    def process_h5ad(self,filename:str):

        # Downloading item locally
        self.gc.downloadItem(
            itemId = self.omics_item_id,
            dest = os.getcwd()+'/'
        )
        ann_data_object = ad.read_h5ad(
            os.getcwd()+'/'+filename
        )

        # Coordinates are stored in obsm['spatial']
        #TODO: This should be adaptable to other types of "spatial" measurements
        self.coordinates = pd.DataFrame(
            data = ann_data_object.obsm['spatial'],
            index = ann_data_object.obs_names,
            columns = ['imagecol','imagerow']
        )

        # For extracting gene expression data
        if not self.omics_args is None:
            if self.omics_args['method'] == 'highest_mean':
                # Getting the genes with the highest mean value
                mean_vals = ann_data_object.var['mean'].sort_values(ascending=False)
                highest_n_mean_genes = list(mean_vals.iloc[0:self.omics_args['n']].index)

                subset_ann_data = ann_data_object[:,highest_n_mean_genes]

            elif self.omics_args['method'] == 'specific_list':
                # Passing a list of specific genes to include
                subset_ann_data = ann_data_object[:,self.omics_args['list']]

            elif self.omics_args['method'] == 'highly_variable':

                # Checking the "flavor"
                # See: https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html
                flav = ann_data_object.uns['hvg']['flavor']
                if flav=='seurat':
                    # Default value, uses dispersion (normalized)
                    # See: https://satijalab.org/seurat/reference/findvariablefeatures (mean.var.plot)
                    dispersions = ann_data_object.var['dispersions_norm'].sort_values(ascending=False)
                    highest_dispersed_genes = list(dispersions.iloc[0:self.omics_args['n']].index)

                    subset_ann_data = ann_data_object[:, highest_dispersed_genes]

                elif flav=='seurat_v3' or flav=='seurat_v3_paper':
                    # Same link as above but using "vst"
                    hv_rank = ann_data_object.var['highly_variable_rank'].sort_values(ascending=False)
                    highest_hv_rank = list(hv_rank.iloc[0:self.omics_args['n']].index)

                    subset_ann_data = ann_data_object[:, highest_hv_rank]

                else:
                    # FLAVA FLAV
                    print(f'Flavor: {flav} not supported')
                    raise NotImplementedError

            self.omics = pd.DataFrame(
                data = subset_ann_data.X,
                index = subset_ann_data.obs_names,
                columns = subset_ann_data.var_names
            )

    def process_rds(self):

        # Processing RDS file containing spot coordinates and omics data
        robjects.r('library(Seurat)')
        robjects.r('library(stringr)')

        robjects.r('''
                    # Function to extract normalized cell type fractions from RDS file
                    get_cell_norms <- function(rds_file){
                        print(rds_file)

                        read_rds_file <- readRDS(url(rds_file))

                        if ("predsubclassl2" %in% names(read_rds_file@assays)){
                            cell_type_fract <- GetAssayData(read_rds_file@assays[["predsubclassl2"]])
                        } else if ("pred_subclass_l2" %in% names(read_rds_file@assays)){
                            cell_type_fract <- GetAssayData(read_rds_file@assays[["pred_subclass_l2"]])
                        }

                        # Normalizing so that the columns sum to 1
                        print('Generating normalized cell type fractions')
                        cell_type_fract <- cell_type_fract[1:nrow(cell_type_fract)-1,]
                        cell_type_norm <- sweep(cell_type_fract,2,colSums(cell_type_fract),'/')
                        cell_type_norm[is.na(cell_type_norm)] = 0

                        return(as.data.frame(cell_type_norm))
                    }
                    
                    ''')

        robjects.r('''
                   # Function to get spot coordinates
                   get_spot_coords <- function(rds_file){
                        print('Getting spot coordinates')

                        read_rds_file <- readRDS(url(rds_file))
                        spot_coords <- read_rds_file@images[["slice1"]]@coordinates
                   
                        return(as.data.frame(spot_coords))
                   }
                   ''')
        
        robjects.r('''
                   # Function to get UMI Counts
                   get_umi_counts <- function(rds_file){
                        print('Getting UMI counts')
                   
                        read_rds_file <- readRDS(url(rds_file))
                        umi_count <- data.frame(read_rds_file@meta.data[["nCount_Spatial"]],row.names=rownames(read_rds_file@images[["slice1"]]@coordinates))
                        colnames(umi_count)<- 'UMI Counts'
                   
                        return(umi_count)
                   }
                   
                   ''')

        get_cell_norms = robjects.globalenv['get_cell_norms']
        get_spot_coords = robjects.globalenv['get_spot_coords']
        get_umi_counts = robjects.globalenv['get_umi_counts']

        rds_file_address = f'{self.gc.urlBase.replace("//","/").replace("http:","http:/")}item/{self.rds_file}/download?token={self.user_token}'

        print(f'rds_file_address: {rds_file_address}')
        cell_norm_output = get_cell_norms(rds_file_address)
        spot_coords_output = get_spot_coords(rds_file_address)
        #umi_counts_output = get_umi_counts(rds_file_address)

        # Converting R dataframes to pandas dataframes
        with (robjects.default_converter + pandas2ri.converter).context():
            cell_type_norms = robjects.conversion.get_conversion().rpy2py(cell_norm_output)
            spot_coordinates = robjects.conversion.get_conversion().rpy2py(spot_coords_output)
            #umi_counts = robjects.conversion.get_conversion().rpy2py(umi_counts_output)

        self.omics = cell_type_norms
        self.coordinates = spot_coordinates
        #self.umi_counts = umi_counts

    def process_omics(self):

        sub_types_list = self.definitions['Sub_Types'].tolist()
        main_types_list = self.definitions['Main_Types'].tolist()

        if 'Cell_States' in self.definitions:
            cell_states_list = self.definitions['Cell_States'].tolist()
        else:
            cell_states_list = []

        main_cell_types = pd.DataFrame()
        cell_states = {}

        for m_idx,m in enumerate(main_types_list):
            m_subtypes = sub_types_list[m_idx].split('.')

            # Getting columns of self.omics to aggregate for this main cell type
            agg_subs = [i for i in self.omics.columns.tolist() if i in m_subtypes]

            if len(agg_subs)>0:
                # This should be a dataframe with rows having the "barcode" or "cell_id" and columns having the main cell type
                agg_main_sub = self.omics.iloc[:,[agg_subs]].sum(axis=1).to_frame()
                agg_main_sub.columns = [m]
                agg_main_sub.index = list(self.omics.index)

                if main_cell_types.empty:
                    main_cell_types = agg_main_sub
                else:
                    main_cell_types = pd.concat([main_cell_types,agg_main_sub],axis=1,ignore_index=True)

                if not len(cell_states)==0:
                    m_states = cell_states_list[m_idx].split('.')
                    unique_m_states = np.unique(m_states).tolist()
                    non_normed_dict = pd.DataFrame()
                    for cs in unique_m_states:
                        st_states = [m_subtypes[i] for i in range(len(m_subtypes)) if m_subtypes[i] in agg_subs and m_states[i]==cs]

                        if len(st_states)>0:
                            if non_normed_dict.empty:
                                non_normed_dict = self.omics.iloc[:,[st_states]].sum(axis=1).to_frame()
                            else:
                                non_normed_dict = pd.concat([non_normed_dict,self.omics.iloc[:,[st_states]].sum(axis=1).to_frame()],axis=1,ignore_index=True)
                    
                    non_normed_dict.columns = unique_m_states

                    cell_states[m] = (non_normed_dict/non_normed_dict.sum(axis=1)).fillna(0)

        # Normalizing main cell types
        main_cell_types = (main_cell_types/main_cell_types.sum(axis=1)).fillna(0)

        processed_omics = {
            'Cell_Subtypes': self.omics,
            'Main_Cell_Types': main_cell_types,
            'Cell_States': cell_states
        }

        return processed_omics

    def create_shapes(self, processed_omics:dict):

        roi_annotations = wak.Annotations(mpp = self.mpp)
        roi_annotations.add_names([self.roi_name])

        # Iterating through barcodes or cell ids to create
        if 'diameter' in self.roi_shape_args or 'radius' in self.roi_shape_args:
            # Creating circle ROIs
            if 'diameter' in self.roi_shape_args:
                pixel_diameter = int((1/self.mpp)*self.roi_shape_args['diameter'])
                pixel_radius = int(pixel_diameter/2)
            elif 'radius' in self.roi_shape_args:
                pixel_radius = int((1/self.mpp)*self.roi_shape_args['radius'])
        else:
            pixel_radius = None

        if 'geojson' in self.roi_shape_args:
            self.gc.downloadItem(
                itemId = self.roi_shape_args['geojson'],
                dest = os.getcwd()+'/'
            )
            geojson_name = self.gc.get(f'/item/{self.roi_shape_args["geojson"]}')['name']
            with open(os.getcwd()+'/'+geojson_name,'r') as f:
                roi_shapes = json.load(f)

                f.close()
            
            roi_df = gpd.GeoDataFrame.from_features(roi_shapes['features'])
        else:
            roi_df = None


        for b_idx,b in tqdm(enumerate(list(self.omics.index)),total=self.omics.shape[0]):
            
            if not pixel_radius is None:
                b_coords = self.coordinates.loc[b]
                b_x = b_coords['imagecol']
                b_y = b_coords['imagerow']

                spat_poly = Point(b_x,b_y).buffer(pixel_radius)
            elif not roi_df is None:
                spat_poly = shape(roi_df.iloc[b_idx,:]['geometry'])
            else:
                print('Shape argument not implemented')
                print(self.roi_shape_args)
                raise NotImplementedError
            
            spat_properties = {}

            if 'Main_Cell_Types' in processed_omics:
                # This is the only weird one since cell states are stored per main cell type
                spat_properties['Main_Cell_Types'] = processed_omics['Main_Cell_Types'].loc[b].to_dict()
                spat_properties['Cell_Subtypes'] = processed_omics['Cell_Subtypes'].loc[b].to_dict()
                spat_properties['Cell_States'] = {
                    i: processed_omics['Cell_States'][i].loc[b].to_dict()
                    for i in list(processed_omics['Cell_States'].keys())
                }
            else:
                spat_properties = {
                    i: processed_omics[i].loc[b]
                    for i in list(processed_omics.keys())
                }

            roi_annotations.add_shape(
                poly = spat_poly,
                box_crs = [0,0],
                structure = self.roi_name,
                name = b,
                properties = spat_properties
            )
        
        return roi_annotations
    
    def save(self, annotations):

        annotations = wak.Histomics(annotations).json
        _ = self.gc.post(
            f'/annotation/item/{self.image_item_id}',
            data = json.dumps(annotations),
            headers = {
                'X-HTTP-Method': "POST",
                'Content-Type': 'application/json'
            }
        )



def main(args):

    for a in vars(args):
        print(f'{a}: {getattr(args,a)}')

    gc = girder_client.GirderClient(apiUrl=args.girderApiUrl)
    gc.setToken(args.girderToken)

    image_file_id = args.image_file_id
    file_info = gc.get(f'/file/{image_file_id}')
    image_item_id = file_info['itemId']

    image_metadata = gc.get(f'/item/{image_item_id}/tiles')

    mpp = image_metadata['mm_x']

    omics_file_id = args.omics_file_id
    file_info = gc.get(f'/file/{omics_file_id}')
    omics_item_id = file_info['itemId']

    if args.definitions_file == '':
        definitions_file = None
    else:
        definitions_file = args.definitions_file

    if len(list(args.omics_args.keys()))>0:
        omics_args = args.omics_args
    else:
        omics_args = None


    SpatAnnotation(
        gc = gc,
        image_item_id = image_item_id,
        omics_item_id = omics_item_id,
        roi_name = args.roi_name,
        mpp = mpp,
        definitions_file = definitions_file,
        omics_args = omics_args,
        roi_shape_args = args.roi_shape_args
    )



if __name__=='__main__':
    main(CLIArgumentParser().parse_args())
