#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter
from adjustText import adjust_text
# import mplcursors
#import prince
#from pca import pca
from scipy import interpolate
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


# get_ipython().run_line_magic('matplotlib', 'inline')
# **First part consist of parsing the two excel tables containg OTUs counts and their taxonomy and generate one file containing all the data**
def pars_tax(file):
    """
    This function takes a .taxonomy file which should be saved under .txt file
    and split the taxonomic levels in different columns, also eliminates the bootstraps values
    """
    global df4
    tax = pd.read_table(file)
    df1 = pd.DataFrame(tax)
    df2 = df1["Taxonomy"].str.split(";", expand=True)
    df2 = df2.drop([6], axis=1)
    column_names = ["kingdom", "phylum", "class", "order", "family", "genus"]
    df2.columns = column_names
    df3 = df2.replace(to_replace=r"\(.*\)", value="", regex=True)
    otu_names = df1["OTU"]
    df4 = pd.concat([otu_names, df3], axis=1)
    return df4

def pars_otu(file):
    """
    This function takes a .subsampled.shared file which should be saved under .txt file
    and transposing it in order to fit to the taxonomy file
    """
    global df1
    otu = pd.read_table(file)
    df1 = pd.DataFrame(otu)
    df1 = df1.drop(["numOtus", "label"], axis=1)
    df1 = df1.set_index("Group").transpose()
    return df1

# pars_tax("TS.pick.opti_mcc.0.03.cons.txt")
# pars_otu("TS.pick.opti_mcc.0.03.subsample.txt")

def concat_tables(otu, tax):
    taxonomy = tax.set_index("OTU")
    otu.rename(index={'Group': 'OTU'}, inplace=True)
    concat_table = pd.concat([taxonomy, otu], axis=1, sort=True)
    return concat_table.dropna()


def keep_more_than(tax_level, perc, df):
    """this function keep the chosen taxonomic level above the threshold"""
    column_names = ["kingdom", "phylum", "class", "order", "family", "genus"]
    sample_names = set(df.columns.to_list()) - set(column_names)
    for sample in sample_names:
        df["perc_{}".format(sample)] = df[sample] / df[sample].sum()
    df = df.groupby(tax_level).sum() 
    fractions = set(df.columns.to_list()) - set((sample_names)) 
    fractions = list(fractions)   
    try: 
        mask = (df[fractions[0]] > perc) | (df[fractions[1]] > perc) | (df[fractions[2]] > perc) | (df[fractions[3]] > perc) | (df[fractions[4]] > perc) | (df[fractions[5]] > perc) | (df[fractions[6]] > perc) | (df[fractions[7]] > perc) | (df[fractions[8]] > perc) | (df[fractions[9]] > perc) | (df[fractions[10]] > perc) | (df[fractions[11]] > perc)               
    except: IndexError       
    try: 
        mask = (df[fractions[0]] > perc) | (df[fractions[1]] > perc) | (df[fractions[2]] > perc) |(df[fractions[3]] > perc) | (df[fractions[4]] > perc) | (df[fractions[5]] > perc) |(df[fractions[6]] > perc) | (df[fractions[7]] > perc) | (df[fractions[8]] > perc) | (df[fractions[9]] > perc) | (df[fractions[10]] > perc)                
    except: IndexError      
    try: 
        mask = (df[fractions[0]] > perc) | (df[fractions[1]] > perc) | (df[fractions[2]] > perc) |(df[fractions[3]] > perc) | (df[fractions[4]] > perc) | (df[fractions[5]] > perc) |(df[fractions[6]] > perc) | (df[fractions[7]] > perc) | (df[fractions[8]] > perc) | (df[fractions[9]] > perc)               
    except: IndexError     
    try: 
        mask = (df[fractions[0]] > perc) | (df[fractions[1]] > perc) | (df[fractions[2]] > perc) |(df[fractions[3]] > perc) | (df[fractions[4]] > perc) | (df[fractions[5]] > perc) |(df[fractions[6]] > perc) | (df[fractions[7]] > perc) | (df[fractions[8]] > perc)
    except: IndexError     
    try: 
        mask = (df[fractions[0]] > perc) | (df[fractions[1]] > perc) | (df[fractions[2]] > perc) |(df[fractions[3]] > perc) | (df[fractions[4]] > perc) | (df[fractions[5]] > perc) |(df[fractions[6]] > perc) | (df[fractions[7]] > perc)               
    except: IndexError       
    try: 
        mask = (df[fractions[0]] > perc) | (df[fractions[1]] > perc) | (df[fractions[2]] > perc) |(df[fractions[3]] > perc) | (df[fractions[4]] > perc) | (df[fractions[5]] > perc) |(df[fractions[6]] > perc)                
    except: IndexError    
    try: 
        mask = (df[fractions[0]] > perc) | (df[fractions[1]] > perc) | (df[fractions[2]] > perc) |(df[fractions[3]] > perc) | (df[fractions[4]] > perc) | (df[fractions[5]] > perc)
    except: IndexError    
    try: 
       mask = (df[fractions[0]] > perc) | (df[fractions[1]] > perc) | (df[fractions[2]] > perc) |(df[fractions[3]] > perc) | (df[fractions[4]] > perc)
    except: IndexError     
    try: 
        mask = (df[fractions[0]] > perc) | (df[fractions[1]] > perc) | (df[fractions[2]] > perc) |(df[fractions[3]] > perc)                
    except: IndexError   
    try: 
        mask = (df[fractions[0]] > perc) | (df[fractions[1]] > perc) | (df[fractions[2]] > perc)                
    except: IndexError     
    try: 
        mask = (df[fractions[0]] > perc) | (df[fractions[1]] > perc)                
    except: IndexError    
    try: 
        mask = (df[fractions[0]] > perc)                 
    except: IndexError

    df_new = df[mask]
    df_to_keep = df[~mask] 

    my_dict = df_to_keep.sum().to_dict()
    my_dict[tax_level] = 'Less than {}%'.format(perc)
    df_new.reset_index(inplace=True)

    row  = pd.DataFrame([my_dict], columns=my_dict.keys())
    df_new = pd.concat([df_new, row], axis =0, sort=True)
    return df_new


def barplot(df):
    samples_to_name = list(df.columns)[1:len(df.columns)//2+1]
    samples_to_plot = list(df.columns)[len(df.columns)//2+1:len(df.columns)]
    df = df.set_index(df.columns[0])
    plt.style.use('ggplot')
    b = df[samples_to_plot].transpose().plot(kind='bar', stacked=True, alpha=0.8)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.63)) 
    plt.xlabel("samples", fontsize="large")
    plt.ylabel("Relative abundance", fontsize="large")
    plt.show()
    fig = b.get_figure()
    return fig.savefig("Output_1.jpg")


def heatmap_plot(df):
    df = df.set_index(df.columns[0])
    samples_to_plot = list(df.columns)[len(df.columns)//2:len(df.columns)+1]
    h = sns.heatmap(df[samples_to_plot], cmap="BuPu", annot=True, linewidths=2, linecolor="#B7C3D0", fmt=".2f")
    plt.xlabel("Samples")
    plt.show()
    fig = h.get_figure()
    return fig.savefig("Output.jpg")
#plt.savefig('heatmap_figure.png')


def info(df, combo):
    #distinct taxa above the threshold
    dist_otus = len((df[combo].to_list())) - 1
    # samples names
    sampl_names = list(df.columns[1:round(len(df.columns)/2)])
    return ("this field will contain explanations about the results \nAbout this file:\nThe samples names are {}\n Number of distinct taxa above the treshhold is {}".format(sampl_names, dist_otus))


def PCA_f(taxa, num, data, df):
    scaler = StandardScaler()
    sampl_names = list(df.columns[1:round(len(df.columns)/2)])
    d_pca = data.groupby(taxa).sum()[sampl_names]
    std_d_pca = scaler.fit_transform(d_pca)
    # d_pca.to_csv("Y.csv")
    pca = PCA()
    pc = pca.fit_transform(std_d_pca)
    d_pca['PC 1'] = pc[:,0]
    d_pca['PC 2'] = pc[:,1]
    sns.set_style("white")
    x = pc[:,0]
    y = pc[:,1]
    n = d_pca.index[:num]
    plt.figure(figsize=(10, 5))
    p1 = sns.scatterplot(x="PC 1", y="PC 2", data=d_pca, s=100)
    texts = []
    for x, y, s in zip(x, y, n):
        texts.append(plt.text(x, y, s))

    f = interpolate.interp1d(pc[:,0], pc[:,1])
    x = np.arange(min(pc[:,0]), max(pc[:,0]), 0.0005)
    y = f(x)

    plt.xlabel("PC1 ({}%)".format((pca.explained_variance_ratio_[0] * 100).round()), fontsize="xx-large")
    plt.ylabel("PC2 ({}%)".format((pca.explained_variance_ratio_[1] * 100).round()), fontsize="xx-large")
    plt.title(("Principal component analysis plot of 16S rRNA gene amplicon sequencing"), fontsize="xx-large")
    adjust_text(texts, x=x, y=y, autoalign='y',
                only_move={'points':'y', 'text':'y'}, force_points=0.35,
                arrowprops=dict(arrowstyle="->", color='b', lw=0.88))

    d_pca['PC 1'] = pc[:,0]
    d_pca['PC 2'] = pc[:,1]

    x = pc[:,0]
    y = pc[:,1]

    n = d_pca.index[:num]

    coeff = np.transpose(pca.components_[0:2, :])
    for i, j in zip(range(coeff.shape[0]), d_pca.columns):
        plt.plot(coeff[i,0], coeff[i,1],color = 'g',alpha = 0.5, marker='o', markersize=20)
        plt.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, j, color = 'black', ha = 'left', va = 'top', fontsize=14)
    fig = p1.get_figure()
    return fig.savefig("Output_2.jpg")

