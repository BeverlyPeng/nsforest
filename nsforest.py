
### Libraries ###
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
import time
from tqdm import tqdm
import os
import myrandomforest
import mydecisiontreeevaluation
# import utils
# import psutil

# v4.0 includes new parameter, "gene_selection," which determines whether BinaryFirst is used or not and its cutoff value
def NSForest(lmao):

    """
    Performs NSForest algorithm to find an optimal list of marker genes. 

    :param lmao: a string to test with 
    :type lmao: str
    :return: a list
    :rtype: list[str]
    """
    
    return ["test1", "test2", "test3"]
