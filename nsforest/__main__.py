
_epilog = """
Examples

python3 nsforest -a arguments.csv
python3 nsforest -a arguments_layer1.csv
python3 nsforest -a arguments_hlca_core.csv

"""

import sys
import numpy as np
import pandas as pd
import scanpy as sc
import utils
import nsforest
import decisiontreewithmarkerlist
import calculate_diagonals
# import os
# import psutil
# ASCTplusB tables print("MEMORY CHECK", psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2, "MiB")
def main():

    # Preparing kwargs
    args = utils._parse_args(sys.argv[1:], epilog=_epilog, description=__doc__)
    kwargs = utils.get_kwargs(args.a)
    if args.c: 
        kwargs["cluster_list"] = [args.c.replace("\r", "")]
        # Adding cluster folder to output_folder
        if "output_folder" in kwargs: 
            kwargs["output_folder"] = kwargs["output_folder"] + args.c.replace("\r", "").replace(" ", "_") + "/"
        # Re-assigning outputfilename to be cluster specific
        kwargs["outputfilename"] = args.c.replace("\r", "").replace(" ", "_")
    
    print("Input values...")
    print(kwargs)

    # Loading h5ad file
    adata = sc.read_h5ad(kwargs["filename"])
    print("adata...")
    print(adata)
    
    # Calculating median gene expression matrix if not pre-calculated
    to_save = False
    if "medians_header" not in kwargs: 
        adata = utils.preprocessing_medians(adata, kwargs["cluster_header"])
        kwargs["medians_header"] = 'medians_' + kwargs["cluster_header"]
        to_save = True
    # Calculating binary score gene expression matrix if not pre-calculated
    if "binary_scores_header" not in kwargs: 
        adata = utils.preprocessing_binary(adata, kwargs["cluster_header"], 'medians_' + kwargs["cluster_header"])
        kwargs["binary_scores_header"] = 'binary_scores_' + kwargs["cluster_header"]
        to_save = True
    # Saving new h5ad
    if to_save: 
        print("Saving new anndata object as " + kwargs["filename"].replace(".h5ad", "_precalculated.h5ad"))
        adata.write_h5ad(kwargs["filename"].replace(".h5ad", "_precalculated.h5ad"))
    else: 
        print(f'Using pre-calculated cluster_medians and binary_scores in {kwargs["filename"]}')
    print(f'\tanndata.varm[{kwargs["medians_header"]}] and anndata.varm[{kwargs["binary_scores_header"]}]')

    # Return if p parameter specified
    if args.p: 
        # Saving clusters as csv
        df = pd.DataFrame()
        df["cluster"] = np.unique(adata.obs[kwargs["cluster_header"]].astype(str))
        name = kwargs["output_folder"].replace("_outputs", "_inputs") + "clusters_" + kwargs["outputfilename"] + ".csv"
        df.to_csv(name, index = False, header = False)
        print("Saving clusters as...")
        print(name)
        return

    if "marker_genes_csv" not in kwargs: 
        nsf_kwargs = utils.subset_kwargs("nsf.NSF", kwargs)
        nsforest.NSForest(adata, **nsf_kwargs)
    else: 
        # Getting custom marker genes list from csv
        kwargs["marker_genes_dict"] = utils.prepare_markers(kwargs["marker_genes_csv"])
        dt_kwargs = utils.subset_kwargs("dtwml.DT", kwargs)
        decisiontreewithmarkerlist.DecisionTree(adata, **dt_kwargs)

    # Getting output marker list
    nsf_markers_df = pd.read_csv(kwargs["output_folder"] + kwargs["outputfilename"] + "_markers.csv")
    if "subclade_list" in kwargs: 
        subclade_list = utils.str_to_list(kwargs["subclade_list"])
        nsf_markers_df = nsf_markers_df[nsf_markers_df['clusterName'].isin(subclade_list)].reset_index(drop=True)
    cluster_medians = adata.varm[kwargs["medians_header"]].transpose()
    calculate_diagonals.calculate_diagonals(nsf_markers_df, cluster_medians, kwargs["output_folder"], kwargs["outputfilename"])

main()