#! /usr/bin/env python

import csv
import re
import argparse
import glob
import os
import collections
import itertools
from collections import OrderedDict


parser = argparse.ArgumentParser()
parser.add_argument('--variantfiles', required=False, default="./Called/*_postCentrifuge.gatk_position_variants_cf4_cr4_fr75_ph4_outmode000.tab", help='Folder containing the variant tab files (default is "Called")')
parser.add_argument('--goi', required=False, default="/home/avujkovic/ESX/ID_all_genes.txt", help='list of all locustag IDs of ESX genes of interest')
args = parser.parse_args()

files=glob.glob(args.variantfiles)
#a function to read files with csv.Dictreader
def readFile(file):
    return csv.DictReader(open(file,mode='r'), delimiter='\t')

#read the file ID_esx which contains the list of locustags of interest
DictIDesx = readFile(args.goi)
#create a list of all the locustags of interest
listlocustags=[]
#parse the ID_esx file to list all the locustags
for row in DictIDesx:
    listlocustags.append(row['ID'])



#make a dictionary of all the variants
#rowDict = {'#Pos':row['#Pos'],'Ref':row['Ref'],'Type':row['Type'],'Allel':row['Allel'],'Subst':row['Subst'],'Gene':row['Gene'],'GeneName':row['GeneName'],'Product':row['Product'],"Lineage":lineage}

for gene in listlocustags:
    GlobalDictVariant=[]
    allsamples = []
    uniquepositions = []
    print(gene)
    for file in files:
        filename =re.sub("./Called/","",file)
        filename =re.sub("_postCentrifuge.gatk_position_variants_cf4_cr4_fr75_ph4_outmode000.tab","",filename)
        allsamples.append(filename)
        #print(filename)
        DictVariantfile = readFile(file)

        for row in DictVariantfile:
            sampleDictvariant={}
            if row["Gene"] == gene and row["Type"] == "SNP":
                sampleDictvariant["sample"]=filename
                sampleDictvariant["gene"]= row["Gene"]
                sampleDictvariant["pos"]= row['#Pos']
                sampleDictvariant["ref"]= row['Ref']
                sampleDictvariant["allel"]= row['Allel']
                GlobalDictVariant.append(sampleDictvariant)
                if row['#Pos'] not in uniquepositions:
                    uniquepositions.append(row['#Pos'])

    results=[]
    for pos in uniquepositions:
        result={"pos":pos}
        for elem in allsamples:
            result[elem]=""
        WT=""
        for dict in GlobalDictVariant:
            if dict["pos"] == result["pos"]:
                result[dict["sample"]]= dict["allel"]
                WT=dict["ref"]
        #print(WT+ "dit is het WT")

        for key in result:
            if bool(result[key]) == False:
                #print(key)
                result[key]= WT
            else:
                continue
        results.append(result)
    fieldnames=['pos']+allsamples
    #sort all SNP's according position
    sorted_results=sorted(results,key = lambda i: i['pos'])
    #open file and wright results
    with open("./SNPtableallgenes/"+gene+"_snp.fasta",'w', newline='') as f:
        w = csv.DictWriter(f,fieldnames, delimiter=',')
        w.writeheader()
        for result in sorted_results:
            w.writerow(result)
