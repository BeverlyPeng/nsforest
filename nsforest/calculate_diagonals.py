
import pandas as pd
import utils
import scanpy as sc

# 
def violin_plot(adata, nsf_results_df, cluster_header): 

    adata = sc.read_h5ad("demo_data/hlca_core_ann_finest_level_precalculated.h5ad")
    nsf_results_df = pd.read_csv("outputs_ann_finest_level/ann_finest_level_results.csv")

    # Getting results
    nsf_results_df = nsf_results_df.dropna()
    # nsf_results_df['NSForest_markers'] = utils.str_to_list(nsf_results_df['NSForest_markers'])
    nsf_results_df['NSForest_markers'] = str_to_list(nsf_results_df['NSForest_markers'])

    # reorder the clusters based on the dendrogram order
    dend_header = "dendrogram_" + cluster_header
    if dend_header not in adata.uns: 
        sc.tl.dendrogram(adata, groupby=cluster_header)
    ## obtain dendrogram cluster order
    dend_order = adata.uns[dend_header].get('categories_ordered')
    df_dend_order = pd.DataFrame({'clusterName': dend_order})
    nsf_results_df.index = nsf_results_df['clusterName']
    nsf_results_df = nsf_results_df.reindex(dend_order)

    # store cluster-marker information in a dictionary
    markers_dict = dict(zip(nsf_results_df["clusterName"], nsf_results_df["NSForest_markers"]))
    markers = []
    for key in sorted(markers_dict.keys()): 
        markers.extend(markers_dict[key])
    print(len(markers))
    print(markers)

    # use scanpy plot function to plot the marker genes
    sc.pl.stacked_violin(adata, markers, cluster_header, standard_scale = "var", dendrogram = True, save = "_2.png") # change dendrogram=False to dendrogram=True to include dendrogram if you made it!

    # sc.pl.dotplot(adata, markers, cluster_header, standard_scale = "var", save = "_2.png")

    return

# Calculate ratio of diagonal/total expression
def on_target_ratio(nsf_markers_df, cluster_medians, output_folder, outputfilename): 
    

    marker_list = nsf_markers_df['markerGene']
    clusters = pd.unique(nsf_markers_df['clusterName'])

    df = pd.DataFrame()
    total_target_exp = 0
    total_subclade_exp = 0
    for m in range(len(marker_list)):
        # get median expression of marker in its target cluster
        marker = marker_list[m]
        target_cluster = nsf_markers_df.loc[m, 'clusterName']
        target_exp = cluster_medians.loc[target_cluster, marker]
        total_target_exp += target_exp

        # get median expression values of this marker m in all other clusters
        for cluster in clusters:
            if cluster == target_cluster: # don't include expression in target cluster! (so this ratio could be  > 1!)
                continue
            else:
                total_subclade_exp += cluster_medians.loc[cluster, marker]

        df_marker = pd.DataFrame({'marker': [marker],
                                'total_target_exp': [total_target_exp],
                                'total_subclade_exp': [total_subclade_exp],
                                'difference': [total_subclade_exp - total_target_exp],
                                'ratio': [total_target_exp/total_subclade_exp]})
        df = pd.concat([df, df_marker]).reset_index(drop=True)

    df.to_csv(output_folder + outputfilename + "_diagonal.csv", index = False)
    print("Saving diagonals in " + output_folder + outputfilename + "_diagonal.csv")

    return df
