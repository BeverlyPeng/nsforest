
import pandas as pd

# Calculate ratio of diagonal/total expression
def calculate_diagonals(nsf_markers_df, cluster_medians, output_folder, outputfilename): 
    
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
