#!/usr/bin/env python
import os
import sys
mdir = os.path.dirname(os.path.realpath(sys.argv[0])) + "/MAS_v2"
sys.path.append(mdir)
import pandas as pd
import GenotypeDF

import screed
from screed import ScreedDB


### 

inp=sys.argv[1]
alleleFasta=sys.argv[2]
outDir=sys.argv[3]


isExist = os.path.exists(outDir)
if not isExist:
    os.makedirs(outDir)

df=pd.read_csv(inp,sep="\t")
p1 = GenotypeDF.GenotypeDF(df)

df1=p1.gt
df2=p1.depth
alleles=p1.alleles

screed.read_fasta_sequences(alleleFasta)
fadb = ScreedDB(alleleFasta)

for s in p1.samples:
    oup=open(os.path.join(outDir,"%s.fastq" % s), "w")
    L1=df1.loc[:,s].tolist()
    L2=df2.loc[:,s].tolist()
    L3=[]
    for x in L2:
		#print(x)
        if len(x)==1:
            L3.append([int(x[0]/2),int(x[0]/2)])

        else:
            L3.append(x)

    for x in range(len(L1)):
		#print(x)
		#print(L1)
        for y in [0,1]:
		    #if not pd.isna(L1[x]):
			#print(type( L1[x] ))
            if not isinstance(L1[x], float):
				#print(L1[x])
                allele = alleles[x] + "#" + str(L1[x][y])
                if allele in fadb:
                    seqStr = fadb[allele].sequence
                    seqId=f"{allele}_{y}"
                    qualStr = "I" * len(seqStr)
					#seq=">%s_%s\n%s\n" % (allele, str(y), fadb[allele].sequence)
                    oup.write(f"@{seqId}\n{seqStr}\n+\n{qualStr}\n")

    oup.close()

	
	
