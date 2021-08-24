from __future__ import division
import statistics
import random
import numpy.random as npr
from search import runSearch
import time
import numpy as np
import matplotlib.pyplot as plt
import pickle
import math


# Attempt to calculate P-value by generating random feature vector from feature distribution
# Note: change distance range in parseFASTA for dist optmization

t0 = time.time()
def parseFASTA(file, dct):

    over = 0
    seqID =[]
    seqs=[]
    lengths = []
    tempIDs =[]
    total = 0
    temp=""
    with open(file) as f:
        line = f.readline()
        seqID.append(line.strip())
        count = 1
        while line:
            line = f.readline()
            if line.strip():  # ignore blank lines
                if (count % 2 == 0): #ID
                    seqID.append(line.strip())
                    temp = line.strip()
                else: # seq
                    seq = line.strip()
                    # print(seq)
                    length = len(seq)
                    if length > 70 and length <100 and "X"not in seq and "Z" not in seq:
                        tempIDs.append(temp)
                        over+=1

                    lengths.append(length)
                    # print(list(seq))
                    for i in seq:
                        if i != "X" and i!="Z":
                            dct[i] += 1
                            total+=1
                    # print(length)
                    seqs.append(seq)
            count +=1
        print("over: ",over)
        return seqID, seqs, lengths, dct, total, tempIDs


aa = 'ARNDCQEGHILKMFPSTWYVU'
keys = list(aa)
aaFreq = {key: 0 for key in keys}

seqID, seqs , lengths, dct, total, ids= parseFASTA('mobidb_all.fasta', aaFreq)
print(len(ids))
print(ids[0])

arr_df= pickle.load(open('mobiALL_db.pkl', 'rb'))
p = pickle.load(open('mobiALL_hnsw.pkl','rb'))
IDRs_names = pickle.load(open('mobiALL_IDRs_names.pkl', 'rb'))

names = {}

temp_df = arr_df.loc[ids,:]
arr = [0] *1000
i = 0
while i <1000:
    # print(i)
    col_names = list(arr_df.columns)
    rand_seq = {feature: 0.0 for feature in col_names}
    while len(col_names) >0:
    
        random_feature = random.choice(col_names)
        col_names.remove(random_feature)
        random_row = temp_df.sample()
        row_name = list(random_row.index)[0]
        value = arr_df.loc[row_name, random_feature]
        rand_seq[random_feature] = value


    random_feature_vector = np.asarray(list(rand_seq.values()))
    labels, distances = p.knn_query(random_feature_vector, k =1)
    dists = [math.sqrt(distances[0][i]) for i in range(len(distances[0]))]
    top_matches = [IDRs_names[labels[0][i]] for i in range(len(labels[0]))]
    
    # if top_matches[0] in names:
    #     names[top_matches[0]] +=1
    # else:
    #     names[top_matches[0]] =0

    arr[i] = dists[0]
    i+=1


# 'MAELSEQVQNLSINDNNENGYVPPHLRGKPRSARNNSSNYNNNNGGYNGGRGGGSFFSNNRRGGYGNGGFFGGNNGGSRSNGRSGGRW'

# fus MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNYGQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRGGRGRGGSGGGGGGGGGGYNRSSGGYEPRGRGGGRGGRGGMGGSDRGGFNKFGGPRDQGSRHDSEQDNSD
# ddx3 MSHVAVENALGLDQQFAGLDLNSSDNQSGGSTASKGRYIPPHLRNREATKGFYDKDSSGWSSSKDKDAYSSFGSRSDSRGKSSFFSDRGSGSRGRFDDRGRSDYDGIGSRGDRSGFGKFERGGNSRWCDKSDEDDWSKP
# ddx4 MGDEDWEAEINPHMSSYVPIFEKDRYSGENGDNFNRTPASSSEMDDGPSRRDHFMKSGFASGRNFGNRDAGECNKRDNTSTMGGFGVGKSFGNRGFSNSRFEDGDSSGFWRESSNDCEDNPTRNRGFSKRGGYRDGNNSEASGPYRRGGRGSFRGCRGGFGLGSPNNDLDPDECMQRTGGLFGSRRPVLSGTGNGDTSQSRSGSGSERGGYKGLNEEVITGSGKNSWKSEAEGGESSDTQGP
query_dist , matches= runSearch("MAELSEQVQNLSINDNNENGYVPPHLRGKPRSARNNSSNYNNNNGGYNGGRGGGSFFSNNRRGGYGNGGFFGGNNGGSRSNGRSGGRW" ,10)
print(query_dist)
print(matches)

# hold = []
# for i in matches:
#     if i in names:
#         # if names[i] >0:
#         hold.append((i,names[i]))


counter_lst=[]
for hit in query_dist:
    counter = 0
    for i in arr:
        # print(i)
        if i <= hit:
            counter+=1
    # print(counter)
    counter_lst.append(counter)
# print(counter_lst)
p_values = [c/len(arr) for c in counter_lst]
# p_value = counter/len(arr) 
print(p_values)


plt.boxplot(arr,showfliers=False)
plt.show()

