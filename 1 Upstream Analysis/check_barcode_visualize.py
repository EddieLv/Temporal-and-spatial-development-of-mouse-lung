# Description: Give a stat report of the barcode file depending on a ref-barcode and visualize
# Author: genger
# Data: 2022-08-18-Thursday

import pandas as pd
import numpy as np
from tqdm import *
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

import os
import sys
import getopt

## change the path 
path = "./"
os.chdir(path)

## usage
def usage():
    print("")
    print("Stat and visualize your barcode file.")
    print("uage: python %s -option <argument>" %sys.argv[0])
    print(" -h/--help ")
    print(" --bc=<STRING> barcode input file.")
    print(" --refA=<STRING> bc_A_url barcode file.")
    print(" --refB=<STRING> bc_B_url barcode file.")
    print(" --fig=<STRING> figure output file.")
    print(" --map_stat=<STRING> mapped pixel stat file.")
    print(" --unmap_stat=<STRING> unmapped pixel stat file.")
    print(" --barcodeNum=<INTEGER> num of your barcode.")

## deal with options
try: 
    opts, args = getopt.getopt(sys.argv[1:], "h", ["help","bc=","refA=","refB=","fig=","map_stat=","unmap_stat=","barcodeNum="])
except getopt.GetoptError:
    print("ERROR: Get option error.\nYou can contact the author through wechat 13958598285.")
    usage()
    sys.exit(2)

for opt, val in opts:
    if opt in ("-h", "--help"):
        usage()
        sys.exit(1)
    else:
        if opt in ("--bc"):
            bc_url = val 
        if opt in ("--refA"):
            bc_A_url = val
        if opt in ("--refB"):
            bc_B_url = val
        if opt in ("--fig"):
            fig_url = val
        if opt in ("--map_stat"):
            map_stat_url = val
        if opt in ("--unmap_stat"):
            unmap_stat_url = val
        if opt in ("--barcodeNum"):
            num = int(val)

try:
    print("Your bc is:",bc_url)
    print("Your refA is:",bc_A_url)
    print("Your refB is:",bc_B_url)
    print("Your fig is:",fig_url)
    print("Your map_stat is:",map_stat_url)
    print("Your unmap_stat is:",unmap_stat_url)
    print("Your barcodeNum is:",num)
except:
    print("ERROR: Missing option.")
    usage()
    sys.exit(2)

bc_A=pd.read_table(bc_A_url,header=None)[0].values
bc_B=pd.read_table(bc_B_url,header=None)[0].values
bc=pd.read_table(bc_url,header=None)
bc_stat=bc.value_counts()

bc_out=pd.DataFrame({"iB":np.repeat(range(1,num+1),num,axis=0),"bc_B":np.repeat(bc_B,num,axis=0),"iA":list(range(1,num+1))*num,"bc_A":list(bc_A)*num})
bc_out["pixel"]=bc_out[["bc_B","bc_A"]].apply(lambda x: "_".join(x),axis=1)
bc_out["raw_count"]=bc_out["pixel"].apply(lambda x: bc_stat[x] if x in bc_stat.index else 0)
bc_out.to_csv(map_stat_url,index=None)

bc_stat=pd.DataFrame(bc_stat)
bc_stat.columns=["count"]
ambiguous_pixel=[i for i in bc_stat.index if i not in bc_out["pixel"].values]
mis_stat=bc_stat.loc[ambiguous_pixel]
mis_stat.to_csv(unmap_stat_url)

plt.figure(figsize=(8,8))
normlizer=mpl.colors.Normalize(vmin=np.min(bc_out["raw_count"]), vmax=np.quantile(bc_out["raw_count"],0.99))
sns.scatterplot(data=bc_out,x="iB",y="iA",marker="s",hue="raw_count",hue_norm=normlizer,palette="rainbow",alpha=0.9)
plt.xlim(0,num+1)
plt.ylim(0,num+1)
if num == 50:
	plt.xticks(ticks=[1,10,20,30,40,50], fontsize=8, rotation=60)
	plt.yticks(ticks=[1,10,20,30,40,50], fontsize=8, rotation=60)
	plt.gca().invert_xaxis()
	plt.gca().invert_yaxis()
if num == 96:
	plt.xticks(ticks=[1,20,40,60,80,96], fontsize=8, rotation=60)
	plt.yticks(ticks=[1,20,40,60,80,96], fontsize=8, rotation=60)
#plt.xticks(ticks=np.append(np.arange(1, num+1, 5),[num]), labels=bc_B[np.append(np.arange(0, num, 5),[num-1])], fontsize=8, rotation=60)
#plt.yticks(ticks=np.append(np.arange(1, num+1, 5),[num]), labels=bc_A[np.append(np.arange(0, num, 5),[num-1])], fontsize=8, rotation=60)
plt.legend(loc="upper right", bbox_to_anchor=(1.2,1))
plt.title("raw_Reads heatmap")
#plt.axis("off")
plt.savefig(fig_url,dpi=500,bbox_inches="tight")



