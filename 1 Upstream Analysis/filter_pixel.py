# Description: Crop pixels based on two svgs(one is original, the other one is cropped)
# Author: genger
# Data: 2022-05-12-Thursdays

import seaborn as sns
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import *

import os
import sys
import getopt

## change the path 
path = "./"
os.chdir(path)

## usage
def usage():
    print("")
    print("Crop pixels based on two svgs(one is original, the other one is cropped).")
    print("uage: python %s -option <argument>" %sys.argv[0])
    print(" -h/--help ")
    print(" --coorR=<STRING> raw svg coordinate file.")
    print(" --coorC=<STRING> cropped svg coordinate file.")
    print(" --output=<STRING> output pixel file.")

## deal with options
try: 
    opts, args = getopt.getopt(sys.argv[1:], "h", ["help","coorR=","coorC=","output="])
except getopt.GetoptError:
    print("ERROR: Get option error.\nYou can contact the author through wechat 13958598285.")
    usage()
    sys.exit(2)

for opt, val in opts:
    if opt in ("-h", "--help"):
        usage()
        sys.exit(1)
    else:
        if opt in ("--coorR"):
            coord_ori = val
        if opt in ("--coorC"):
            coord_crop = val
        if opt in ("--output"):
            output_url = val

try:
    print("Your coorR is:",coord_ori)
    print("Your coorC is:",coord_crop)
    print("Your output is:",output_url)
except:
    print("ERROR: Missing option.")
    usage()
    sys.exit(2)

def coord2df(url):
    df=pd.read_table(url,header=None,sep=" ")
    #df["row"]=df[0].apply(lambda x: (math.ceil(x*10) / 10)).astype("str")
    #df["col"]=df[1].apply(lambda x: (math.ceil(x*10) / 10)).astype("str")
    #df["row"]=df[0].apply(lambda x: round(x,1)).astype("str")
    #df["col"]=df[1].apply(lambda x: round(x,1)).astype("str")
    #df["coord"]=df["col"].str.cat(df["row"],sep="_")
    return df
    
def judge_coord(x,y,df):
    n=0
    df["x_match"]=df[0].apply(lambda xx: abs(xx-x) < 0.1)
    df["y_match"]=df[1].apply(lambda yy: abs(yy-y) < 0.1)
    df["match"]=df["x_match"] & df["y_match"]
    n=df.loc[df["match"]==True].shape[0]
    if n == 0:
        return False
    if n == 1:
        return True
    if n > 1:
        return "Multi"
        
ori_df=coord2df(coord_ori)
crop_df=coord2df(coord_crop)

if np.sqrt(ori_df.shape[0]) == 96:
    num=96
    ori_df["iB"]=np.repeat(range(1,num+1),num,axis=0)
    ori_df["iB"]=ori_df["iB"].astype("str")
    ori_df["iA"]=list(range(num,0,-1))*num
    ori_df["iA"]=ori_df["iA"].astype("str")
    ori_df["BA"]=ori_df["iB"].str.cat(ori_df["iA"],sep="x")

if np.sqrt(ori_df.shape[0]) == 50:
    num=50
    ori_df["iB"]=np.repeat(range(num,0,-1),num,axis=0)
    ori_df["iB"]=ori_df["iB"].astype("str")
    ori_df["iA"]=list(range(1,num+1))*num
    ori_df["iA"]=ori_df["iA"].astype("str")
    ori_df["BA"]=ori_df["iB"].str.cat(ori_df["iA"],sep="x")
    
match=[]
for i in tqdm(range(ori_df.shape[0])):
    match.append(judge_coord(ori_df.iloc[i,0],ori_df.iloc[i,1],crop_df))
ori_df["match"]=match

print(ori_df["match"].value_counts())
print("if [Multi] exists, please contact with the author!")
ori_df=ori_df.loc[ori_df["match"]==False]
ori_df["BA"].to_csv(output_url,sep="\t",index=None,header=None)



