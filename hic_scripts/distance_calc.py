import sys
import pandas as pd
import string
import math
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as colors
import seaborn as sns
from itertools import chain
import os
from matplotlib.backends.backend_pdf import PdfPages
    
def parse_args(args):
    parser = argparse.ArgumentParser(description='Check the help flag')
    parser.add_argument('-g',
                        '--gene_list',
                        dest='genes',
                        help='Required: List of genes in gff format.',
                        required=True)
    parser.add_argument('-p1',
                        '--pdb1',
                        dest='pdb1',
                        help='REQUIRED: PDB file for conditional sample.',
                        required=True)
    parser.add_argument('-p2',
                        '--pdb2',
                        dest='pdb2',
                        help='REQUIRED: PDB file for control sample.',
                        required=True)
    parser.add_argument('-n',
                        '--na_bins',
                        dest='na_bins',
                        help='REQUIRED: List of missing (NA) bins required.',
                        required=True)
    parser.add_argument('-o',
                        '--output',
                        dest='output',
                        help='Output directory.',
                        required=False)
    return parser.parse_args()


def input_params(args):
    if args.pdb1 is None or args.pdb2 is None:
        print("PDB files are required.")
        exit(1)
    else:
        pdb1 = pd.read_csv(args.pdb1, sep=r"\s+", usecols=[4,5,6,7,8], header=None)
        pdb2 = pd.read_csv(args.pdb2, sep=r"\s+", usecols=[4,5,6,7,8], header=None)
    if args.na_bins is None:
        print("List of bins not found in dataset is required")
        exit(1)
    else:
        na_bins = pd.read_csv(args.na_bins, sep=r"\s+", header=None)
    if args.output is None:
        print("Output will be saved to same directory as first PDB file.")
        out = os.path.dirname(args.pdb1)
    else:
        out = args.output
    return pdb1, pdb2, na_bins, out


def gene_list(gl, bins, cl, res):
    gl = pd.read_csv(gl, sep='\t', header=None)
    gl = gl[gl[2] == "protein_coding_gene"]
    chr_breaks = gl[0].value_counts().sort_index().tolist()
    tot_breaks = 0
    for i in range(len(chr_breaks)):
        chr_breaks[i] += tot_breaks
        tot_breaks = chr_breaks[i]
    gl[9] = ((gl[4]-gl[3])/2+gl[3]).floordiv(res).add(1).astype(int)
    gl[8] = gl[8].str.split(";", expand = True)[0].str.split("=", expand = True)[1]
    for i in gl[0].unique():
        chain = cl[int(bins.loc[bins[0] == i].index.values)]
        gl.replace(i, chain, inplace=True)
    gl[1] = gl[0].map(str) + gl[9].map(str)
    gl = gl[[1,8]]
    gl.columns = ["bin","gene"]
    return gl, chr_breaks

def read_na_bins(na):
    gl = pd.read_csv(gl, sep='\t', header=None)

def gene_pdb(pdb, gl):
    pdb = pd.read_csv(pdb, sep=r"\s+", usecols=[4,5,6,7,8], header=None)
    pdb[1] = pdb[4].map(str) + pdb[5].map(str)
    pdb = pdb[[1,4,6,7,8]]
    pdb.columns = ["bin","chr","x","y","z"]
    print(gl)
    #pdb = pdb[[1,6,7,8]]
    #pdb.columns = ["bin","x","y","z"]
    gpdb = pd.merge(gl, pdb, how="left", on="bin")
    #gpdb.drop("bin", axis=1, inplace=True)
    return gpdb

def telomere_pdb(pdb):
    pdb = pd.read_csv(pdb, sep=r"\s+", usecols=[1,4,5,6,7,8], header=None)
    pdb.columns = ['bin', 'chr', 'chr_bin', 'x', 'y', 'z']
    idx_max = pdb.groupby('chr')['chr_bin'].idxmax()
    idx_min = pdb.groupby('chr')['chr_bin'].idxmin()
    pdb_max = pdb.loc[idx_max]
    pdb_min = pdb.loc[idx_min]
    tel_pdb = pd.concat([pdb_max, pdb_min]).sort_values(by=['chr', 'chr_bin'])
    tel_pdb = tel_pdb[['bin', 'x', 'y', 'z']]
    return tel_pdb

def dist_array(gpdb, tpdb):
    d_arr = []
    t_arr = []
    for idx,row in gpdb.iterrows():
        d_list = []
        c1 = [row["x"],row["y"],row["z"]]
        for idx2,row2 in gpdb.iterrows():
            c2 = [row2["x"],row2["y"],row2["z"]]
            if row["chr"] != row2["chr"]:
                d = math.dist(c1,c2)
                d_list.append(d)
        d_arr.append(d_list)
    return np.array(d_arr)
"""    for idx,row in tpdb.iterrows():
        t_list = []
        c1 = [row["x"],row["y"],row["z"]]
        for idx2,row2 in tpdb.iterrows():
            c2 = [row2["x"],row2["y"],row2["z"]]
            t = math.dist(c1,c2)
            t_list.append(t)
        t_arr.append(t_list)
    return np.array(d_arr), np.array(t_arr)
"""
def plot_dist(d1, d2, t1, t2, gl, chr_breaks, chroms):#, pdf):
    if d2 is not None:
        d = np.subtract(d1, d2)
        d[np.isnan(d)] = 0
        np.savetxt("/mnt/f/Deitsch_collab/HiC/output_files/D2_distances.txt",d,delimiter="\t")
        d = np.rot90(d)
        max_val = np.max(np.abs(d))
        sns.set(font_scale=d.shape[1]/20)
        plt.figure(figsize=(d.shape[1],d.shape[0]))
        hm = sns.heatmap(d, vmax=max_val, vmin=-max_val,
                            linewidths=0.5, linecolor='black',
                            cmap='RdBu_r', center=0, square=True,
                            xticklabels=gl.gene.tolist(), 
                            yticklabels=list(reversed(gl.gene.tolist())),
                            cbar_kws={"shrink": .8})
        hm.axhline(y=0, color='black', linewidth=d.shape[1]/10)
        hm.axhline(y=d.shape[1], color='black', linewidth=d.shape[1]/10)
        hm.axvline(x=0, color='black', linewidth=d.shape[0]/10)
        hm.axvline(x=d.shape[0], color='black', linewidth=d.shape[0]/10)
        for i in chr_breaks:
            j = d.shape[1] - i
            hm.axhline(y=j, color='black', linewidth=d.shape[1]/10)
        for i in chr_breaks:
            hm.axvline(x=i, color='black', linewidth=d.shape[0]/10)
        #pdf.savefig()
        plt.close()

    if t2 is not None:
        t = np.subtract(t1, t2)
        t[np.isnan(t)] = 0
        t = np.rot90(t)
        max_val = np.max(np.abs(t))
        sns.set(font_scale=t.shape[1]/10)
        plt.figure(figsize=(t.shape[1],t.shape[0]))
        hm = sns.heatmap(t, vmax=max_val, vmin=-max_val,
                                cmap='RdBu_r', center=0, square=True,
                                cbar_kws={"shrink": .8})
        hm.set_yticks(range(1,len(chroms)*2,2))
        hm.set_yticklabels(list(reversed(chroms)), rotation=0)
        hm.set_xticks(range(1,len(chroms)*2,2))
        hm.set_xticklabels(chroms, rotation=90)
        hm.axhline(y=0, color='black', linewidth=t.shape[1]/10)
        hm.axhline(y=t.shape[1], color='black', linewidth=t.shape[1]/10)
        hm.axvline(x=0, color='black', linewidth=t.shape[0]/10)
        hm.axvline(x=t.shape[0], color='black', linewidth=t.shape[0]/10)
        for i in range(2,len(chroms)*2,2):
            hm.axhline(y=i, color='black', linewidth=0.5)
            hm.axvline(x=i, color='black', linewidth=0.5)
        #pdf.savefig()
        plt.close()

def main():
    args = parse_args(sys.argv[1:])
    pdb1, pdb2,  = input_params(args)

    chain_lst = list(string.ascii_uppercase)
    bins = sizes2bins(sys.argv[2], res)

    gl, chr_breaks = gene_list(sys.argv[1], chain_lst, res)

    gene_pdb1 = gene_pdb(sys.argv[3], gl)
    gene_pdb2 = gene_pdb(sys.argv[4], gl)

    tel_pdb1 = telomere_pdb(sys.argv[3])
    tel_pdb2 = telomere_pdb(sys.argv[4])

    tel_pdb1 = tel_pdb1.merge(tel_pdb2, how="inner", on="bin")[["bin","x_x","x_y","z_x"]]
    tel_pdb1.columns = ["bin","x","y","z"]
    tel_pdb2 = tel_pdb2.merge(tel_pdb1, how="inner", on="bin")[["bin","x_x","x_y","z_x"]]
    tel_pdb2.columns = ["bin","x","y","z"]

    tel_list = bins[[2,3]].values.tolist()
    tel_list = list(map(int, chain.from_iterable(tel_list)))
    missing_tel = list(set(tel_list) - set(tel_pdb1["bin"]))
    missing_tel_df = pd.DataFrame({"bin":missing_tel, "x":0, "y":0, "z":0})
    
    tel_pdb1 = pd.concat([tel_pdb1, missing_tel_df], axis=0).sort_values("bin")
    tel_pdb2 = pd.concat([tel_pdb2, missing_tel_df], axis=0).sort_values("bin")
    print(gene_pdb1)
    print(gene_pdb2)
    exit(1)
    d1 = dist_array(gene_pdb1, tel_pdb1)
    d2 = dist_array(gene_pdb2, tel_pdb2)
    print(np.nanmean(d1))
    print(np.nanmean(d2))
    print(np.nanmean(t1))
    print(np.nanmean(t2))
    #out = os.path.dirname(os.path.abspath(sys.argv[3]))
    #pdf = PdfPages("%s/%s" % (out, os.path.basename(str(sys.argv[3]).replace(".pdb", "_diff_distance.pdf"))))
    #plot_dist(d1, d2, t1, t2, gl, chr_breaks, bins[0].tolist())#, pdf)
    #pdf.close()

if __name__ == "__main__":
    main()