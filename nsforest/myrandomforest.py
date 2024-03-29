
import pandas as pd
import scanpy as sc
from sklearn.ensemble import RandomForestClassifier

## run Random Forest on the binary dummy variables ==> outputs all genes ranked by Gini impurity
def myRandomForest(adata, df_dummies, cl, n_trees, n_jobs, n_top_genes, binary_dummies):
    """\
    Returning top genes sorted by gini index from sklearn.ensemble's RandomForestClassifier. 

    Parameters
    ----------
    adata: AnnData
        Annotated data matrix.
    df_dummies: pd.DataFrame
        Dummy dataframe for one vs all model. 
    cl: str
        Specified cluster. 
    n_trees: int
        Number of `n_estimators` in sklearn.ensemble's RandomForestClassifier.
    n_jobs: int
        Number of `n_jobs` in sklearn.ensemble's RandomForestClassifier. 
    n_top_genes: int
        Taking the top `n_top_genes` ranked by sklearn.ensemble's RandomForestClassifier for sklearn.tree's DecisionTreeClassifier. 
    binary_dummies: pd.DataFrame
        Dataframe of binary scores filtered by `gene_selection`. 

    Returns
    -------
    top_rf_genes: top `n_top_genes` genes ranked by gini index
    """
    x_train = adata.to_df()
    y_train = df_dummies[cl]
    
    ## pre-select genes based on gene_selection criterium
    ind_genes_selected = binary_dummies.loc[cl] == 1
    x_train = x_train.loc[:,ind_genes_selected]   # subset x_train
    print("\t", "Pre-selected", x_train.shape[1], "genes to feed into Random Forest.")

    rf_clf = RandomForestClassifier(n_estimators=n_trees, n_jobs=n_jobs, random_state=123456) #<===== criterion="gini", by default
    rf_clf.fit(x_train, y_train)
    ## get feature importance and rank/subset top genes
    top_rf_genes = pd.Series(rf_clf.feature_importances_, index=x_train.columns).sort_values(ascending=False)[:n_top_genes]
    
    return top_rf_genes  
