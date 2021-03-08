from sklearn.preprocessing import MinMaxScaler
import networkx as nx
from sklearn.metrics import pairwise_distances

import numpy as np 
import pandas as pd
from pandas import DataFrame
import metrics as me

def run_silhouette_metrics(adata, batch_key='orig.ident', cell_label='paper.cell.type', embed='X_pca'):
    # global silhouette coefficient
    sil_global = me.silhouette(adata, group_key='paper.cell.type', embed=embed, metric='euclidean')
    # silhouette coefficient per batch
    _, sil_clus = me.silhouette_batch(adata, batch_key=batch_key, group_key=cell_label,
                                       embed=embed, metric='euclidean', verbose=False)
    sil_clus = sil_clus['silhouette_score'].mean()
    il_score_sil = me.isolated_labels(adata, label_key=cell_label, batch_key=batch_key,
                                       cluster=False, n=1, verbose=False)
    il_score_clus = me.isolated_labels(adata, label_key=cell_label, batch_key=batch_key,
                                cluster=True, n=1, verbose=False)
    
    return(sil_global, sil_clus, il_score_clus, il_score_sil)

def run_silhouette_metrics_methods(adata_list, df,
                                  batch_key='orig.ident', cell_label='paper.cell.type', embed='X_pca', extra_label=''):
    silhouette_scores = [[], [], [], []]
    for method in adata_list:
        score1, score2, score3, score4 = run_silhouette_metrics(method,
                                                                batch_key=batch_key,
                                                                cell_label=cell_label,
                                                                embed=embed)
        silhouette_scores[0].append(score1)
        silhouette_scores[1].append(score2)
        silhouette_scores[2].append(score3)
        silhouette_scores[3].append(score4)
        
    df['ASW_label'+extra_label] = silhouette_scores[0]
    df['ASW_label/batch'+extra_label] = silhouette_scores[1]
    df['isolated_label_F1'+extra_label] = silhouette_scores[2]
    df['isolated_label_silhouette'+extra_label] = silhouette_scores[3]




def accuracy_paired_omics(adata, omic_layer, variable, cell_name=None, percent=False):
    """
    will match cell barcode from paired measurement for 2 layers. 
    I will return the ratio of cells for which the RNA and ATAC barcode end up in the same cluster.
    
    adata : coembed multiomic object
    variable : cell clustering obs variable
    cell_name : obs variable containing the matching barcodes for the 2 omic layers
    omic_layer : obs variable containing the batch/omic layer of origin
    percent=True  return percentage. if false, return ratio
    
    return
    ------
    
    accuracy: float ratio of cells for which the barcodes end up in the same barcodes
    
    
    """
    
    # extract important informations from the adata.obs
    df = adata.obs[[omic_layer, variable]]
    if cell_name != None:
        df.index = adata.obs[cell_name]
    else:
        df.index = adata.index
        
    
    # split RNA and ATAC cells in 2 dataframes
    omic_layer_variable = list(set(df[omic_layer]))
    df_atac = df[df[omic_layer]==omic_layer_variable[0]]
    df_rna = df[df[omic_layer]==omic_layer_variable[1]]

    
    # only keep cells that are present in the 2 tables
    atac_cells = []
    for name in df_atac.index.tolist():
        if name in df_rna.index:
            atac_cells.append(True)
        else:
            atac_cells.append(False)
        
    rna_cells = []
    for name in df_rna.index.tolist():
        if name in df_atac.index:
            rna_cells.append(True)
        else:
            rna_cells.append(False)
        
    df_rna = df_rna[rna_cells].sort_index()
    df_atac = df_atac[atac_cells].sort_index()
    
    
    # get the accuracy
    x =(df_rna[variable] == df_atac[variable]).tolist()
    Correct_values = x.count(True)
    False_values = x.count(False)
    accuracy = Correct_values/(Correct_values+False_values)
    if percent==True:
        accuracy=accuracy*100
    return(accuracy)

def accuracy_paired_omics_per_cell_type(adata, omic_layer, variable, cell_type, cell_name=None, percent=False):
    """
    will match cell barcode from paired measurement for 2 layers. 
    I will return the dict of ratio of cells for which the RNA and ATAC barcode end up in the same cluster.
    But the ratios are per cell types
    
    adata : coembed multiomic object
    variable : cell clustering obs variable
    cell_name : obs variable containing the matching barcodes for the 2 omic layers
    omic_layer : obs variable containing the batch/omic layer of origin
    cell_type : obs variable containing the ground truth cell type
    percent=True  return percentage. if false, return ratio
    
    return
    ------
    
    accuracy: dict of float ratio of cells for which the barcodes end up in the same barcodes per cell type
    
    
    """
    
    # extract important informations from the adata.obs
    df = adata.obs[[omic_layer, variable, cell_type]]
    if cell_name != None:
        df.index = adata.obs[cell_name]
    else:
        df.index = adata.index

    cell_type_dict = {}
    for current_cell_type in sorted(set(df[cell_type])):
        df_cell_type = df[df[cell_type]==current_cell_type]
        
        
        # split RNA and ATAC cells in 2 dataframes
        omic_layer_variable = list(set(df[omic_layer]))
        df_atac = df_cell_type[df_cell_type[omic_layer]==omic_layer_variable[0]]
        df_rna = df_cell_type[df_cell_type[omic_layer]==omic_layer_variable[1]]
        
        # only keep cells that are present in the 2 tables
        atac_cells = []
        for name in df_atac.index.tolist():
            if name in df_rna.index:
                atac_cells.append(True)
            else:
                atac_cells.append(False)
        
        rna_cells = []
        for name in df_rna.index.tolist():
            if name in df_atac.index:
                rna_cells.append(True)
            else:
                rna_cells.append(False)
        
        df_rna = df_rna[rna_cells].sort_index()
        df_atac = df_atac[atac_cells].sort_index()
    
        # get the accuracy
        x =(df_rna[variable] == df_atac[variable]).tolist()
        Correct_values = x.count(True)
        False_values = x.count(False)
        accuracy = Correct_values/(Correct_values+False_values)
        if percent==True:
            accuracy=accuracy*100
        cell_type_dict[current_cell_type] = accuracy
        
    return(cell_type_dict)



def get_pairs(adata, cell_names=None):
    # extract graph tuples for paired cells
    df = adata.obs.copy()
    df

    g = nx.Graph(adata.obsp['connectivities'])
    node_list = []
    for node in g.nodes:
        node_list.append(node)
    df['nodes'] = node_list

    #cell_names = [x.split('_')[0] for x in df.index.tolist()]
    if cell_names == None:
        df['org_cell_name'] = cell_names
    
    df['org_cell_name'] = df[cell_names]
    dico_tuple = {}
    tuples = []
    index = 0
    for name in df['org_cell_name']:
        if name not in dico_tuple.keys():
            dico_tuple[name] = [df['nodes'][index]]
        else:
            dico_tuple[name].append(df['nodes'][index])
            tuples.append(dico_tuple[name])
        index += 1

    return(tuples)


def distance_between_matching_barcodes(adata, cell_names, absolute=True):
    distances = pairwise_distances(adata.X)
    distances_mean = distances.mean(axis=0)
    pairs = get_pairs(adata, cell_names=cell_names)
    list_distance_barcodes = [np.nan]*len(adata.obs_names)

    if absolute == True:
        for position in pairs:
            list_distance_barcodes[position[0]] = distances[position[0], position[1]]
            list_distance_barcodes[position[1]] = distances[position[0], position[1]]
    else:
        for position in pairs:
            list_distance_barcodes[position[0]] = distances[position[0], position[1]]/distances_mean[position[0]]
            list_distance_barcodes[position[1]] = distances[position[0], position[1]]/distances_mean[position[0]]
            #return(list_distance_barcodes)
    adata.obs['euclidean_pairwise_distance_between_matching_barcodes'] = list_distance_barcodes
    
def average_distance_between_matching_barcodes(adata, cell_names=None, cell_type=None, absolute=True):
    distance_between_matching_barcodes(adata, cell_names, absolute)
    
    if cell_type == None:
        average_metric = np.mean(adata.obs['euclidean_pairwise_distance_between_matching_barcodes'])
        return(average_metric)
    else:
        average_dist_per_cluster = {}
        for elem in list(set(adata.obs[cell_type])):
            average_dist_per_cluster[elem] = \
            np.mean(adata[adata.obs[cell_type] == elem,:].obs['euclidean_pairwise_distance_between_matching_barcodes'])
        return(average_dist_per_cluster)
    
# normalising ASW by 0.5x + 0.5
def norm_metrics(df, metric=['ASW_label/batch', 'ASW_label', 'cLisi'], copy=True):
    if copy ==True:
        df_tmp = df.copy()
    else:
        df_tmp = df
    if type(metric) == str:
        metric = [metric]
    for n in metric:
        df_tmp[n] = [0.5*x + 0.5 for x in df_tmp[n]]
    if copy ==True:
        return(df_tmp)
    
def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def linear_scaling(df, metric='pairwise_distance', copy=False):
    if copy ==True:
        df_tmp = df.copy()
    else:
        df_tmp = df
        
    if type(metric) == str:
        metric = [metric]
        
    for n in metric:
        df_tmp[n] = NormalizeData(df_tmp[n])
        
    if copy ==True:
        return(df_tmp)
