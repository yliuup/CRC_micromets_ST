
######
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import anndata
import smurf as su
import pickle
import copy
import gzip
import os
from PIL import Image
from tqdm import tqdm
import sys

param = sys.argv[1]

save_path = '。/2.visium_HD_smurf/'+param+'/'
# Check if the directory exists
if not os.path.exists(save_path):
    # If it doesn't exist, create it
    os.makedirs(save_path)

path  =  '。/Visium_HD/alignment/'+param+'/outs/binned_outputs/square_002um/'


so = su.prepare_dataframe_image(path+'spatial/tissue_positions.parquet',
                                '。/Visium_HD/image/high/'+param+'.tif',
                               'HE')   

with gzip.GzipFile(save_path + 'image_to_segmentation.npy.gz', 'w') as f:
    np.save(f, so.image_temp())


####---stardist-----#######

import os
import sys
import gzip
import numpy as np
from tqdm import tqdm
from stardist.models import StarDist2D
from csbdeep.utils import normalize
from skimage.measure import label

param = sys.argv[1]
save_path = '/rsrch5/home/genomic_med/yliu47/rsrch9/Analysis/CRC_revision/2.visium_HD_smurf/'+param+'/'

os.makedirs(save_path, exist_ok=True)

with gzip.GzipFile(save_path + 'image_to_segmentation.npy.gz', 'r') as f:
    image = np.load(f)

def stardist_patch(img, model, nms_thresh=0.3, prob_thresh=None):
   
    axis_norm = (0,1,2) if img.ndim==3 else (0,1)
    img_n = normalize(img, 2, 98, axis=axis_norm)
    labels, _ = model.predict_instances(
        img_n,
        prob_thresh=prob_thresh,
        nms_thresh=nms_thresh,
        overlap_label=None
    )
    return labels

def segment_large_image_stardist(
    image,
    model,
    loop=4700,
    gap=80,
    nms_thresh=0.3,
    prob_thresh=None,
    verbose=True
):
    i_max, j_max = image.shape[:2]
    final_mask = np.zeros((i_max, j_max), dtype=np.int32)
    current_label = 1
    for i in tqdm(range(0, i_max, loop - gap), disable=not verbose):
        for j in range(0, j_max, loop - gap):
            patch = image[i:min(i+loop, i_max), j:min(j+loop, j_max)]
            masks = stardist_patch(patch, model, nms_thresh, prob_thresh)
            if masks is None or masks.max()==0:
                continue
            masks = (masks > 0) * (masks + current_label)
            mask_slice = final_mask[i:i+masks.shape[0], j:j+masks.shape[1]]
            update = (mask_slice == 0) & (masks > 0)
            mask_slice[update] = masks[update]
            final_mask[i:i+masks.shape[0], j:j+masks.shape[1]] = mask_slice
            current_label = final_mask.max() + 1
    final_mask = label(final_mask > 0, connectivity=1)
    return final_mask


model = StarDist2D.from_pretrained("2D_versatile_he")  


mask = segment_large_image_stardist(
    image,
    model,
    loop=4700,
    gap=80,
    nms_thresh=0.01,      
    prob_thresh=0.1,   
    verbose=True
)


np.save(save_path + 'segmentation_stardist.npy', mask)
print(f"StarDist segmentation save {save_path + 'segmentation_stardist.npy'}")


#######smurf

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import anndata
import smurf as su
import pickle
import copy
import gzip
import os
import sys
param = sys.argv[1]


save_path = './2.visium_HD_smurf/'+param+'/'


# Check if the directory exists
if not os.path.exists(save_path):
    # If it doesn't exist, create it
    os.makedirs(save_path)


load_path = './2.visium_HD_smurf/'+param+'/'


path  =  './Visium_HD/alignment/'+param+'/outs/binned_outputs/square_002um/'


so = su.prepare_dataframe_image(path+'spatial/tissue_positions.parquet',
                                './Visium_HD/image/high/'+param+'.tif',
                               'HE')   


so.segmentation_final = np.load(load_path+'segmentation_stardist.npy')
so.generate_cell_spots_information()


adata = sc.read_10x_mtx(path +'filtered_feature_bc_matrix')
adata = copy.deepcopy(adata[so.df[so.df.in_tissue == 1]['barcode']])
adata.obs

sc.pp.filter_genes(adata, min_counts=500)

su.nuclei_rna(adata,so)
adata_sc = copy.deepcopy(so.final_nuclei)
adata_sc

sc.pp.filter_cells(adata_sc, min_counts=5)
adata_raw = copy.deepcopy(adata_sc)
adata_raw


os.chdir(save_path)

import pickle

adata_sc.write("adata_sc.h5ad")
adata_raw.write("adata_raw.h5ad")
adata.write("adata.h5ad")

with open("extra_data.pkl", "wb") as f:
    pickle.dump({"so": so}, f)


import pickle
import scanpy as sc


adata_sc = su.singlecellanalysis(adata_sc,resolution=1)
su.itering_arragement(adata_sc, adata_raw, adata, so, resolution=1, save_folder = save_path, show = True, keep_previous = False)


###
adatas_final = sc.read_h5ad(save_path +'adatas.h5ad')

with open(save_path +'cells_final.pkl', 'rb') as file:
    cells_final = pickle.load(file)

with open(save_path +'weights_record.pkl', 'rb') as file:
    weights_record = pickle.load(file)

pct_toml_dic, spots_X_dic, celltypes_dic, cells_X_plus_dic, nonzero_indices_dic, nonzero_indices_toml, cells_before_ml, cells_before_ml_x, groups_combined, spots_id_dic,spots_id_dic_prop = su.make_preparation(cells_final,
                                                                                                                so, adatas_final, adata, weights_record, maximum_cells = 6000)

import torch
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
device

spot_cell_dic = su.start_optimization(spots_X_dic, celltypes_dic, cells_X_plus_dic, nonzero_indices_toml, device,
                                      num_epochs=1000, learning_rate=0.1, print_each=100, epsilon=0.00001,)
with open(save_path +'spot_cell_dic.pkl', 'wb') as f:
    pickle.dump(spot_cell_dic, f)



import scanpy as sc


path  =  './Visium_HD/alignment/'+param+'/outs/binned_outputs/square_002um/'

adata1 = sc.read_10x_mtx( './Visium_HD/alignment/'+param+'/outs/binned_outputs/square_002um/filtered_feature_bc_matrix')

in_tissue_barcodes = so.df[so.df.in_tissue == 1]['barcode']
adata1 = adata1[in_tissue_barcodes]

from scipy.sparse import csr_matrix, lil_matrix
from numpy.linalg import norm
import numpy as np
import tqdm

def calculate_weight_to_celltype2(adatas_final, adata, cells_final, so):
    cell_ids = list(adatas_final.obs.index.astype(float))
    data_temp = csr_matrix(adata.X)
    final_X = lil_matrix(np.zeros([len(adatas_final.obs), data_temp.shape[1]]))
    for i in tqdm.tqdm(range(len(adatas_final.obs))):
        set_cell_index = [
            so.set_toindex_data[key]
            for key in cells_final[cell_ids[int(i)]]
            if key in so.set_toindex_data
        ]
        final_X[int(i)] = lil_matrix(
            data_temp[set_cell_index, :].sum(axis=0).reshape(1, -1)
        )
    unique_clusters = adatas_final.obs.leiden.astype(str).unique()
    weight_to_celltype = np.zeros((len(unique_clusters), adata.shape[1]))
    for i, cluster in enumerate(unique_clusters):
        mask = (adatas_final.obs.leiden.astype(str) == cluster).values
        weight_to_celltype[i] = final_X[mask].mean(axis=0).A1  # Convert to 1D array
    weight_to_celltype = weight_to_celltype / norm(weight_to_celltype, axis=1).reshape(-1, 1)
    return weight_to_celltype


weight_to_celltype = calculate_weight_to_celltype2(
    adatas_final, adata1, cells_final, so
)

adata_sc_final2 = su.get_finaldata(
    adata1,
    adatas_final,
    spot_cell_dic,
    weight_to_celltype,
    cells_before_ml,
    groups_combined,
    pct_toml_dic,
    nonzero_indices_dic,
    spots_X_dic=None,
    nonzero_indices_toml=nonzero_indices_toml,
    cells_before_ml_x=None
)

adata_sc_final2 = su.get_finaldata(adata, adatas_final, spot_cell_dic, weights_record, cells_before_ml, groups_combined,
                               pct_toml_dic,  nonzero_indices_dic, spots_X_dic, cells_before_ml_x = cells_before_ml_x, so = so)
adata_sc_final2

adata_sc_final2.write(save_path+"adata_sc_final.h5ad")

su.make_pixels_cells(so, adata, cells_before_ml, spot_cell_dic, spots_id_dic, spots_id_dic_prop, nonzero_indices_dic)

with open("final_so_data.pkl", "wb") as f:
    pickle.dump({"so": so}, f)

su.plot_results(so.image_temp(), so.pixels_cells, dpi = 500, save = save_path + 'final_results.png')
