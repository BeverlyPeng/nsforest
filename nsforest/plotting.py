
import pandas as pd
import plotly.express as px
import scanpy as sc
import utils

def plot(NSForest_results = pd.DataFrame(), output_folder = None, outputfilename = None): 

    fig = px.box(NSForest_results, y='f_score', points='all', range_y=[-.05,1.05],
                title=f"F-beta score median = {round(NSForest_results['f_score'].median(),3)}",
                width=400, height=500, hover_name='clusterName')
    fig.write_html(output_folder + outputfilename + "_boxplot_fscore.html")

    fig = px.scatter(NSForest_results, x='clusterSize', y='f_score', range_y=[-.05,1.05],
                 width=700, height=500, hover_name='clusterName')
    fig.write_html(output_folder + outputfilename + "_scatter_fscore.html")

    fig = px.box(NSForest_results, y='PPV', points='all', range_y=[-.05,1.05],
             title=f"Positive predictive value median = {round(NSForest_results['PPV'].median(),3)}",
             width=400, height=500, hover_name='clusterName')
    fig.write_html(output_folder + outputfilename + "_boxplot_ppv.html")

    fig = px.scatter(NSForest_results, x='clusterSize', y='PPV', range_y=[-.05,1.05],
                 width=700, height=500, hover_name='clusterName')
    fig.write_html(output_folder + outputfilename + "_scatter_ppv.html")

    return

def plot_scanpy(adata, cluster_header, NSForest_results = pd.DataFrame()): 
    
    markers_dict = dict(zip(NSForest_results["clusterName"], NSForest_results["NSForest_markers"]))
    markers = []
    for key in sorted(markers_dict.keys()): 
        markers.extend(markers_dict[key])
    print(len(markers))
    sc.pl.dotplot(adata, markers, cluster_header, standard_scale = "var", save = "dotplot_broad.png")

    return
