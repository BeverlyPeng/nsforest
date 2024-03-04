
import argparse
import numpy as np
import pandas as pd
import time
from tqdm import tqdm
import glob

def _parse_args(args, **kwargs): 
    """\
    Simple argument parsing
    """
    parser = argparse.ArgumentParser(**kwargs, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-a", metavar = "arguments", type = str, 
                        help="all input for NSForest algorithm in csv format")
    parser.add_argument("-c", metavar = "cluster", type = str, 
                        help="parallelizing NSForest per cluster")
    return parser.parse_args()

def get_kwargs(filename): 
    """\
    Returning dictionary version of csv, type casting values

    Parameters
    ----------
    filename
        csv representing all arguments
    
    Returns
    -------
    kwargs: dictionary
    """
    df = pd.read_csv(filename).astype(str)
    df = df[df["val"] != "nan"]
    kwargs = dict(zip(df["argument"], df["val"]))
    for key in kwargs.keys(): 
        if "true" in kwargs[key].lower(): 
            kwargs[key] = True
        elif "false" in kwargs[key].lower(): 
            kwargs[key] = False
        elif "n_" == key[:2]: 
            kwargs[key] = int(kwargs[key])
        elif key in ["beta"]: 
            kwargs[key] = float(kwargs[key])
        elif "_list" in key: 
            kwargs[key] = str_to_list(kwargs[key])
        else: 
            kwargs[key] = str(kwargs[key])
    # Creating default output_folder and outputfilename names
    if "outputfilename" not in kwargs: 
        kwargs["outputfilename"] = kwargs["filename"].split("/")[-1].replace(".h5ad", "")
    if "output_folder" not in kwargs: 
        kwargs["output_folder"] = "outputs/"
    elif kwargs["output_folder"][-1] != "/": 
        kwargs["output_folder"] = kwargs["output_folder"] + "/"
    return kwargs

def str_to_list(values): 
    """\
    Converting str representation of list to list

    Parameters
    ----------
    values: str representation of list
    
    Returns
    -------
    values: list
    """
    values = [val.replace("[", "").replace("]", "").replace(", ", ",").replace("'", "").replace('"', "").split(",") for val in values]
    return values

def get_medians(adata, cluster_header):
    """\
    Calculating the median expression per gene for each cluster

    Parameters
    ----------
    adata
        AnnData. Annotated data matrix.
    cluster_header
        Column in `adata`'s `.obs` storing cell annotation.
    
    Returns
    -------
    cluster_medians: gene-by-cluster dataframe
    """
    cluster_medians = pd.DataFrame()
    for cl in tqdm(sorted(set(adata.obs[cluster_header].astype(str))), desc="Calculating medians per cluster"):
        adata_cl = adata[adata.obs[cluster_header]==cl,]
        medians_cl = adata_cl.to_df().median()
        cluster_medians = pd.concat([cluster_medians, pd.DataFrame({cl: medians_cl})], axis=1) #gene-by-cluster
    return cluster_medians

def preprocessing_medians(adata, cluster_header, positive_genes_only = True):
    """\
    Calculating the median expression and filtering genes

    Parameters
    ----------
    adata
        AnnData. Annotated data matrix.
    cluster_header
        Column in `adata`'s `.obs` storing cell annotation.
    
    Returns
    -------
    adata: anndata with cluster_medians in adata.varm and subset by genes_selected
    """
    print("Calculating medians...")
    start_time = time.time()
    ## get medians
    cluster_medians = get_medians(adata, cluster_header) #gene-by-cluster
    ## attach calculated medians to adata
    print("Saving calculated medians as... adata.varm.medians_" + cluster_header)
    adata.varm['medians_' + cluster_header] = cluster_medians #gene-by-cluster
    print("--- %s seconds ---" % (time.time() - start_time))

    if positive_genes_only:
        ## select only genes with median > 0
        genes_selected = cluster_medians.index[cluster_medians.sum(axis=1)>0].to_list()
        print(f"Only positive genes selected. {len(genes_selected)} positive genes out of {adata.n_vars} total genes")
        ## subset data with only positive genes
        adata = adata[:,genes_selected].copy()
    return adata

def preprocessing_binary(adata, cluster_header, medians_header):
    """\
    Calculating the binary scores

    Parameters
    ----------
    adata
        AnnData. Annotated data matrix.
    cluster_header
        Column in `adata`'s `.obs` storing cell annotation.
    medians_header
        Key in `adata`'s `.varm` storing median expression matrix. 
    
    Returns
    -------
    adata: anndata with binary_scores in adata.varm
    """
    print("Calculating binary scores...")
    start_time = time.time()
    ## get medians
    cluster_medians = adata.varm[medians_header].transpose() #cluster-by-gene
    n_total_clusters = cluster_medians.shape[0]
    
    ## calculate binary scores based on cluster_medians for all genes per cluster
    binary_scores = []
    for cl in tqdm(cluster_medians.index, desc="Calculating binary scores per cluster"):
        ## get binary scores for all genes in each row (cluster) in df
        binary_scores_cl = [sum(np.maximum(0,1-cluster_medians[i]/cluster_medians.loc[cl,i]))/(n_total_clusters-1) for i in cluster_medians.columns]
        binary_scores.append(binary_scores_cl)

    ## binary scores matrix and handle nan
    binary_scores = pd.DataFrame(binary_scores, index=cluster_medians.index, columns=cluster_medians.columns).fillna(0) #cluster-by-gene
    ## attach pre-calculated binary scores to adata
    print("Saving calculated binary scores as... adata.varm.binary_scores_" + cluster_header)
    adata.varm['binary_scores_' + cluster_header] = binary_scores.transpose() #gene-by-cluster
    print("--- %s seconds ---" % (time.time() - start_time))
    print("median:", binary_scores.stack().median())
    print("mean:", binary_scores.stack().mean())
    print("std:", binary_scores.stack().std())
    
    return adata

def prepare_markers(filename, col_cluster, col_marker, output_folder = "", outputfilename = ""): 
    """\
    Converting marker_genes_csv to dictionary (clusterName: marker_genes)

    Parameters
    ----------
    filename
        csv with markers per cluster. 
    col_cluster
        Column name in file storing cell annotation.
    col_marker
        Column name in file storing markers.
    output_folder
        Output folder for missing genes. 
    outputfilename
        Prefix for missing genes file. 
    
    Returns
    -------
    marker_genes_dict: dictionary (cluster: list of markers)

    Creates files
    -------------
    {output_folder}{outputfilename}_not_found.csv: list of markers not found in anndata
    """
    marker_genes_csv = pd.read_csv(filename) 
    not_found = []
    marker_genes_dict = {}
    # if markers are represented as list
    if "[" in list(marker_genes_csv[col_marker])[0]: # "['NTNG1', 'EYA4']"
        marker_genes_dict = dict(zip(marker_genes_csv[col_cluster], str_to_list(marker_genes_csv[col_marker])))
    else: # 'NTNG1'\n'EYA4'
        for cluster, marker in zip(marker_genes_csv[col_cluster], marker_genes_csv[col_marker]): 
            if cluster not in marker_genes_dict: 
                marker_genes_dict[cluster] = [marker]
            else: 
                marker_genes_dict[cluster].append(marker)
    # print out too
    if len(not_found) > 1: 
        df = pd.DataFrame()
        df[col_marker] = not_found
        df.to_csv(output_folder + outputfilename + "_not_found.csv", index = False, header = False)

    return marker_genes_dict 

def subset_kwargs(function, kwargs): 
    """\
    Subsetting the arguments per function call

    Parameters
    ----------
    function
        Relevant function call. 
    kwargs
        dictionary of arguments. 
    
    Returns
    -------
    subset_kwargs: list of arguments
    """
    subset_kwargs = {}
    dictionary = {"nsf.NSF": ["cluster_header", "medians_header", "binary_scores_header", 
                              "cluster_list", "gene_selection", 
                              "n_trees", "n_jobs", "beta", "n_top_genes", "n_binary_genes", "n_genes_eval", 
                              "output_folder", "outputfilename"], 
                  "dtwml.DT": ["cluster_header", "medians_header", "binary_scores_header", "marker_genes_dict", 
                              "beta", "exact_genes_eval", "output_folder", "outputfilename"]
                 }
    for key in kwargs: 
        if key in dictionary[function]: 
            subset_kwargs[key] = kwargs[key]
    return subset_kwargs

def combine_results(input_folder, output_folder = "", outputfilename = ""): 
    """\
    Creating single files from each cluster's results

    Parameters
    ----------
    input_folder
        Folder containing all output files for all clusters. 
    output_folder
        Output folder. 
    outputfilename
        Prefix for all output files. 

    Creates files
    -------------
    {output_folder}{outputfilename}_results.csv
    {output_folder}{outputfilename}_diagonal.csv
    {output_folder}{outputfilename}_markers.csv
    {output_folder}{outputfilename}_supplementary.csv
    """
    for suffix in ["results", "diagonal", "markers", "supplementary"]: 
        files = glob.glob(input_folder + "*_" + suffix + ".csv")
        df = pd.DataFrame()
        for file in files: 
            data = pd.read_csv(file)
            df = pd.concat([df, data])
        filename = output_folder + outputfilename + "_" + suffix + ".csv"
        print(filename)
        df.to_csv(filename, index = False)
    return
