import sys, os.path
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from collections import defaultdict
import glob
from scipy.stats import iqr
from scipy.spatial.distance import pdist
import numpy as np
#from skbio import DistanceMatrix
#from skbio.stats.distance import mantel
#from scipy.cluster.hierarchy import ward, complete, average, dendrogram, fcluster, linkage
import re
import glob
import time
import math

#annotation dictionaries /Users/tdangelo/Desktop/projects/kelp/gorg_classifier/diffab_results/annotation_info
ko = "/Users/tdangelo/Desktop/projects/kelp/gorg_classifier/diffab_results/annotation_info/kegg_annotation_dictionary.txt"
kegg = pd.read_csv(ko, sep = "\t", index_col = 0)
pf = "/Users/tdangelo/Desktop/projects/kelp/gorg_classifier/diffab_results/annotation_info/pfam_gorg_combined.txt"
pfam = pd.read_csv(pf, sep = "\t", index_col = 0)
og = "/Users/tdangelo/Desktop/projects/kelp/gorg_classifier/diffab_results/annotation_info/nog_annotation_information_plus.txt"
nog = pd.read_csv(og, sep = "\t", index_col = 0)

# "Check the order of your groups in the output file. Positive effect size means greater abundance in the second group and negative in the first one."
# negative diff.btw value means it is more enriched in the FIRST rab.win category
# positive diff.btw value means it is more enriched in the SECOND rab.win category

outdir = "./aldex2/filtered/"
aldex_files = glob.glob(os.path.join("aldex2/", "ko*pruned10*"))

for f in aldex_files:
    print(f)
    (dir, file) = os.path.split(f)
    sname = file.replace(".txt", "_filtered_annotated.txt")
    oname = os.path.join(outdir, sname) 
    df = pd.read_csv(f, sep ="\t", index_col = 0)
    
    positive = df.columns.values.tolist()[2].rsplit('.',3)[2] ##double check these results
    negative = df.columns.values.tolist()[1].rsplit('.',3)[2]
    
    df["effect_abs"] = abs(df["effect"])
    df = df[df["we.eBH"] < 0.05]
    df["Enriched_Class"] = df["diff.btw"].apply(lambda x: positive if x > 0 else negative)
    df = df.reset_index()
    df = df.rename(columns={"index" : "KO"})
    df_defs = kegg.merge(df, how='right', on='KO')
    df_defs.to_csv(oname, sep="\t", header = True)

outdir = "./aldex2/filtered/"
aldex_files = glob.glob(os.path.join("aldex2/", "OG*pruned10*"))

for f in aldex_files:
    print(f)
    (dir, file) = os.path.split(f)
    sname = file.replace(".txt", "_filtered_annotated.txt")
    oname = os.path.join(outdir, sname) 
    df = pd.read_csv(f, sep ="\t", index_col = 0)
    
    positive = df.columns.values.tolist()[2].rsplit('.',3)[2] ##double check these results
    negative = df.columns.values.tolist()[1].rsplit('.',3)[2]
    
    df["effect_abs"] = abs(df["effect"])
    df = df[df["we.eBH"] < 0.05]
    df["Enriched_Class"] = df["diff.btw"].apply(lambda x: positive if x > 0 else negative)
    df = df.reset_index()
    df = df.rename(columns={"index" : "eggNOG"})
    df_defs = nog.merge(df, how='right', on='eggNOG')
    df_defs.to_csv(oname, sep="\t", header = True)

outdir = "./aldex2/filtered/"
aldex_files = glob.glob(os.path.join("aldex2/", "pfam*pruned10*"))

for f in aldex_files:
    print(f)
    (dir, file) = os.path.split(f)
    sname = file.replace(".txt", "_filtered_annotated.txt")
    oname = os.path.join(outdir, sname) 
    df = pd.read_csv(f, sep ="\t", index_col = 0)
    
    positive = df.columns.values.tolist()[2].rsplit('.',3)[2] ##double check these results
    negative = df.columns.values.tolist()[1].rsplit('.',3)[2]
    
    df["effect_abs"] = abs(df["effect"])
    df = df[df["we.eBH"] < 0.05]
    df["Enriched_Class"] = df["diff.btw"].apply(lambda x: positive if x > 0 else negative)
    df = df.reset_index()
    df = df.rename(columns={"index" : "Pfam"})
    df_defs = pfam.merge(df, how='right', on='Pfam')
    df_defs.to_csv(oname, sep="\t", header = True)
    
outdir = "./ancombc2/filtered/"
ancom_files = glob.glob(os.path.join("ancombc2/", "ko*pruned10*"))

for f in ancom_files:
    
    print(f)
    (dir, file) = os.path.split(f)
    sname = file.replace(".txt", "_filtered_annotated.txt")
    oname = os.path.join(outdir, sname) 
    df = pd.read_csv(f, sep ="\t", index_col = 0)
    
    df = df[df.iloc[:,-1]==True]
    df = df[df.iloc[:,-3]==True]
    df = df[df.iloc[:,-5] < 0.05]
   
    df = df.rename(columns={"taxon" : "KO"})
    df_defs = kegg.merge(df, how='right', on='KO')
    df_defs.to_csv(oname, sep="\t", header = True)

outdir = "./ancombc2/filtered/"
ancom_files = glob.glob(os.path.join("ancombc2/", "pfam*pruned10*"))

for f in ancom_files:

    print(f)
    (dir, file) = os.path.split(f)
    sname = file.replace(".txt", "_filtered_annotated.txt")
    oname = os.path.join(outdir, sname) 
    df = pd.read_csv(f, sep ="\t", index_col = 0)
    
    df = df[df.iloc[:,-1]==True]
    df = df[df.iloc[:,-3]==True]
    df = df[df.iloc[:,-5] < 0.05]
   
    df = df.rename(columns={"taxon" : "Pfam"})
    df_defs = pfam.merge(df, how='right', on='Pfam')
    df_defs.to_csv(oname, sep="\t", header = True)

outdir = "./ancombc2/filtered/"
ancom_files = glob.glob(os.path.join("ancombc2/", "OG*pruned10*"))

for f in ancom_files:
    
    print(f)
    (dir, file) = os.path.split(f)
    sname = file.replace(".txt", "_filtered_annotated.txt")
    oname = os.path.join(outdir, sname) 
    df = pd.read_csv(f, sep ="\t", index_col = 0)
    
    df = df[df.iloc[:,-1]==True]
    df = df[df.iloc[:,-3]==True]
    df = df[df.iloc[:,-5] < 0.05]
   
    df = df.rename(columns={"taxon" : "eggNOG"})
    df_defs = nog.merge(df, how='right', on='eggNOG')
    df_defs.to_csv(oname, sep="\t", header = True)
    
    
    
## write a loop to filter aldex2 results by ancombc2 results
## read in aldex2.lower(), read in ancombc2 of same comaprison, set up to do the one liner below
## aldex2 as the subsetted, ancombc2 after the if statement

#wo_dir = os.path.join(abs_path, o_dir) 
abs_path = os.path.abspath(os.getcwd()) 
wo_dir = os.path.join(abs_path, "composition_intersection/")
if not os.path.exists(wo_dir):
    os.makedirs(wo_dir)

comp_df = pd.DataFrame(columns=['test', 'aldex2', 'ancombc2', 'intersection'])

aldex_files = glob.glob(os.path.join("aldex2/filtered/", "*"))

for f in aldex_files:
    
    (dir, ald_file) = os.path.split(f)
    if ald_file.startswith("ko"):
        print(ald_file)
    
        df1 = pd.read_csv(f, sep ="\t", index_col = 0)
    
        (dir, ald_file) = os.path.split(f)
        sname = ald_file.replace("_filtered_annotated.txt", "_intersection.txt")
        outname = os.path.join(wo_dir, sname)
    
        stat_test1 = ald_file.rsplit('_aldex2_',1)[0].lower()
        print(stat_test1)
    
        ancom_files = glob.glob(os.path.join("ancombc2/filtered/", "*"))
    
        for f in ancom_files:
        
            df2 = pd.read_csv(f, sep ="\t", index_col = 0)
        
            (dir, anc_file) = os.path.split(f)
            stat_test2 = anc_file.rsplit('_ancom_',1)[0].lower()
        
            if stat_test1 == stat_test2:
                print(stat_test2)
                int_df = df1[df1["KO"].isin([value for value in df1.KO.to_list() if value in df2.KO.to_list()])]
                int_df.to_csv(outname, sep="\t", header = True)
                comp_df.loc[len(comp_df.index)] = [stat_test1, len(df1.KO.to_list()), len(df2.KO.to_list()), len(int_df.KO.to_list())]    
            else:
                continue
    
    elif ald_file.startswith("pfam"):
    
        df1 = pd.read_csv(f, sep ="\t", index_col = 0)
    
        (dir, ald_file) = os.path.split(f)
        sname = ald_file.replace("_filtered_annotated.txt", "_intersection.txt")
        outname = os.path.join(wo_dir, sname)
    
        stat_test1 = ald_file.rsplit('_aldex2_',1)[0].lower()
    
        ancom_files = glob.glob(os.path.join("ancombc2/filtered/", "*"))
    
        for f in ancom_files:
        
            df2 = pd.read_csv(f, sep ="\t", index_col = 0)
        
            (dir, anc_file) = os.path.split(f)
            stat_test2 = anc_file.rsplit('_ancom_',1)[0].lower()
        
            if stat_test1 == stat_test2:
                int_df = df1[df1["Pfam"].isin([value for value in df1.Pfam.to_list() if value in df2.Pfam.to_list()])]
                int_df.to_csv(outname, sep="\t", header = True)
                comp_df.loc[len(comp_df.index)] = [stat_test1, len(df1.Pfam.to_list()), len(df2.Pfam.to_list()), len(int_df.Pfam.to_list())]
            else:
                continue
    
    
    else:
        
        df1 = pd.read_csv(f, sep ="\t", index_col = 0)
    
        (dir, ald_file) = os.path.split(f)
        sname = ald_file.replace("_filtered_annotated.txt", "_intersection.txt")
        outname = os.path.join(wo_dir, sname)
    
        stat_test1 = ald_file.rsplit('_aldex2_',1)[0].lower()
    
        ancom_files = glob.glob(os.path.join("ancombc2/filtered/", "*"))
    
        for f in ancom_files:
        
            df2 = pd.read_csv(f, sep ="\t", index_col = 0)
        
            (dir, anc_file) = os.path.split(f)
            stat_test2 = anc_file.rsplit('_ancom_',1)[0].lower()
        
            if stat_test1 == stat_test2:
                int_df = df1[df1["eggNOG"].isin([value for value in df1.eggNOG.to_list() if value in df2.eggNOG.to_list()])]
                int_df.to_csv(outname, sep="\t", header = True)
                comp_df.loc[len(comp_df.index)] = [stat_test1, len(df1.eggNOG.to_list()), len(df2.eggNOG.to_list()), len(int_df.eggNOG.to_list())]
            else:
                continue
