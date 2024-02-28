
_epilog = """
Examples

python3 nsforest -a arguments.csv -p
python3 nsforest -a arguments_layer1.csv

"""

import sys
import numpy as np
import pandas as pd
import scanpy as sc
import utils
import nsforest
import decisiontreewithmarkerlist
import plotting
import calculate_diagonals
import os

def main():

    # Preparing kwargs
    args = utils._parse_args(sys.argv[1:], epilog=_epilog, description=__doc__)
    kwargs = utils.get_kwargs(args.a)

    # If parallelizing NSForest
    if args.c: 
        kwargs["cluster_list"] = [args.c.replace("\r", "").replace('"', "").replace("'", "")]
        if not os.path.exists(kwargs["output_folder"]): 
            os.makedirs(kwargs["output_folder"])
        kwargs["output_folder"] = kwargs["output_folder"] + "nsforest_percluster/"
        kwargs["outputfilename"] = kwargs["cluster_list"][0].replace(" ", "_").replace("/", "_")
    
    # Printing inputs
    print("Input values...")
    print(kwargs)

    # Setting up folder directories
    if not os.path.exists(kwargs["output_folder"]):
        os.makedirs(kwargs["output_folder"])

    # Loading h5ad file
    adata = sc.read_h5ad(kwargs["h5ad"])
    print("adata...")
    print(adata)
    
    # Calculating median gene expression matrix if no header given
    to_save = False
    if "medians_header" not in kwargs: 
        adata = utils.preprocessing_medians(adata, kwargs["cluster_header"])
        kwargs["medians_header"] = 'medians_' + kwargs["cluster_header"]
        to_save = True
    # Calculating binary score gene expression matrix if no header given
    if "binary_scores_header" not in kwargs: 
        adata = utils.preprocessing_binary(adata, kwargs["cluster_header"], kwargs["medians_header"])
        kwargs["binary_scores_header"] = 'binary_scores_' + kwargs["cluster_header"]
        to_save = True
    # Saving new h5ad
    if to_save: 
        filename = kwargs["h5ad"].replace(".h5ad", "_" + kwargs["cluster_header"] + "_precalculated.h5ad")
        print("Saving new anndata object as...\n", filename)
        adata.write_h5ad(filename)
    else: 
        print(f'Using pre-calculated cluster_medians and binary_scores in {kwargs["h5ad"]} from {kwargs["medians_header"]} {kwargs["binary_scores_header"]}')
    print(f'\tanndata.varm[{kwargs["medians_header"]}] and anndata.varm[{kwargs["binary_scores_header"]}]')

    # Saving clusters in file
    if "save_clusters" in kwargs and kwargs["save_clusters"]: 
        df = pd.DataFrame()
        df["cluster"] = np.unique(adata.obs[kwargs["cluster_header"]])
        filename = f'{kwargs["input_folder"]}/clusters_{kwargs["cluster_header"]}.csv'
        print("Saving cluster list as...\n", filename)
        df.to_csv(filename, index = False, header = False)

    if "compute_decision_trees" in kwargs and kwargs["compute_decision_trees"]: 
        # Using nsforest to generate marker genes
        if "marker_genes_csv" not in kwargs: 
            nsf_kwargs = utils.subset_kwargs("nsf.NSF", kwargs)
            nsforest.NSForest(adata, **nsf_kwargs)
        # Inputting marker genes file to create decision trees
        else: 
            kwargs["marker_genes_dict"] = utils.prepare_markers(kwargs["marker_genes_csv"], 
                                                                kwargs["marker_genes_csv_clustercol"], 
                                                                kwargs["marker_genes_csv_markercol"])
            dt_kwargs = utils.subset_kwargs("dtwml.DT", kwargs)
            decisiontreewithmarkerlist.DecisionTree(adata, exact = True, **dt_kwargs)
    
    if "combine_results" in kwargs and kwargs["combine_results"]: 
        # input_folder, output_folder, outputfilename
        utils.combine_results(kwargs["output_folder"], 
                              kwargs["output_folder"].replace("nsforest_percluster/", ""), 
                              kwargs["outputfilename"])

    # Plotting results
    if "plot_results" in kwargs and kwargs["plot_results"]: 
        nsf_results_df = pd.read_csv(kwargs["output_folder"] + kwargs["outputfilename"] + "_results.csv")
        plotting.plot(nsf_results_df, kwargs["output_folder"], kwargs["outputfilename"])
        plotting.plot_scanpy(adata, kwargs["cluster_header"], nsf_results_df)

    # Getting output marker list
    if "calculate_diagonals" in kwargs and kwargs["calculate_diagonals"]: 
        nsf_markers_df = pd.read_csv(kwargs["output_folder"] + kwargs["outputfilename"] + "_markers.csv")
        if "subclade_list" in kwargs: 
            subclade_list = utils.str_to_list(kwargs["subclade_list"])
            nsf_markers_df = nsf_markers_df[nsf_markers_df['clusterName'].isin(subclade_list)].reset_index(drop=True)
        cluster_medians = adata.varm[kwargs["medians_header"]].transpose()
        calculate_diagonals.on_target_ratio(nsf_markers_df, cluster_medians, kwargs["output_folder"], kwargs["outputfilename"])
        # TODO: add diagonal heatmap

main()