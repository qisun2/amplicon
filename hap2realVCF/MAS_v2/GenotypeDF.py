import os
import sys
import pandas as pd
import numpy as np

import numpy as np
from scipy.spatial.distance import pdist
import allel

class GenotypeDF:
    def __init__(self, hap_df):
        self.dataframe = hap_df
        self.samples,self.alleles,self.gt, self.depth=self.__parse()
     
    def __parse(self):
        df=self.dataframe
        alleles=df["Locus"].tolist()
        samples=df.columns.tolist()[2:]
        gt=df[samples].applymap(GenotypeDF.get_gt)
        gt.index=alleles
        depth=df[samples].applymap(GenotypeDF.get_dp)
        #depth.index.name = alleles
        depth.index = alleles
        return samples, alleles, gt, depth

    def get_gt(x):
        #print(x)
        a=x.split(":")[0]
        if a=="./.":
            return np.nan
        else:
            return a.split("/")

    def get_dp(x):
        a=x.split(":")[1]
        if a=="0":
            return [0,0]
        else:
            return [int(i) for i in a.split(",")]

    def dfun(u,v):
        d=len(set(u).symmetric_difference(set(v)))/(len(u)+len(v)) 
        return d

    #def mdist():
        

    def indv_missing(self):
        print("missing",self.gt.shape)
        return self.gt.isna().sum()/self.gt.shape[0]

    def _format_desire(cell_value, column_name, key_dic):
        if column_name not in key_dic:
            print("marker is not in key_dic, error")
            
        outL=[]
        for i in cell_value:
            if i in key_dic[column_name]:
                outL.append(1)
             
        return outL
        
    def binary_desire(self,key_dic):
        self.desire2_gt=self.gt.apply(lambda x: x.apply(_format_desire,args=(x.name,key_dic)))
    
    def get_major(df):
        keydic={}
        for (columnName, columnData) in df.iteritems():
            print('Colunm Name : ', columnName)
            print('Column Contents : ', columnData.values)
            d={}
            for tlist in columnData.values:
                for x in tlist:
                    if x not in tdic:
                        d[x]=0
                    d[x]+=1
            s = [(k, d[k]) for k in sorted(d, key=d.get, reverse=True)]
            keydic[columnName]=s[0][0]
        return keydic

    def binary_major(self):
        key_dic=get_major(self.gt)       

        self.major2_gt=self.gt.apply(lambda x: x.apply(_format_desire,args=(x.name,key_dic)))

    def keep_pos(self,keep_list):
        out_df=self.dataframe.loc[keep_list,:]
        return GenotypeDF(out_df)

    def rm_pos(self,rm_list):
        out_df=self.dataframe.drop([rm_list],axis=0)
        return GenotypeDF(sum_df)
        
    def keep_indv(self,keep_list):
        keep_list.insert(0,"Haplotypes")
        keep_list.insert(0,"Locus")
        out_df=self.dataframe.reindex(columns = keep_list)
        return GenotypeDF(out_df)

    def rm_indv(self,rm_list):
        out_df=self.dataframe.drop([rm_list],axis=1)
        return GenotypeDF(sum_df)

    def __add__(self, other):
        
        #self.dataframe other.dataframe
        sum_df0=self.dataframe.merge(other.dataframe, left_on='Locus', right_on='Locus',how="outer",suffixes=('', '_right'))
        print(sum_df0.shape)
        sum_df=sum_df0.drop(['Haplotypes_right'], axis=1)
        sum_df1=sum_df.fillna("./.:0")
        #sum_df2 =sum_df1.replace(r'^\s*$', '1(1.00);', regex=True)
        return GenotypeDF(sum_df1)




