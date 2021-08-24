from __future__ import division
import statistics
import random
import numpy.random as npr
from search import runSearch

import time
import numpy as np
import matplotlib.pyplot as plt


# Attempt to calculate P-values from generating random sequences from distance distribution

t0 = time.time()
def parseFASTA(file, dct):
    over = 0
    seqID =[]
    seqs=[]
    lengths = []
    total = 0
    with open(file) as f:
        line = f.readline()
        seqID.append(line.strip())
        count = 1
        while line:
            line = f.readline()
            if line.strip():  # ignore blank lines
                if (count % 2 == 0): #ID
                    seqID.append(line.strip())
                else: # seq
                    seq = line.strip()
                    length = len(seq)
                    # if length > 270 and length <300:
                    #     over+=1
                    lengths.append(length)
                    for i in seq:
                        if i != "X" and i!="Z":
                            dct[i] += 1
                            total+=1
                    seqs.append(seq)
            count +=1
        print("over: ",over)
        return seqID, seqs, lengths, dct, total

def generate_seq(new_count, keys):

    randStr = ''
    for count,letter in zip(new_count, keys):
        for i in range(count):
            randStr += letter
    randLst = list(randStr)
    random.shuffle(randLst)
    # print(randLst)
    res = ''.join(randLst)
    # print(res)
    return res

def random_sequence(length):  
    sample_string = 'ARNDCQEGHILKMFPSTWYV' # define the specific string   
    # random.seed(10)
    result = ''.join((random.choice(sample_string)) for x in range(length))  
    return result




t = time.time()

aa = 'ARNDCQEGHILKMFPSTWYVU'
keys = list(aa)
aaFreq = {key: 0 for key in keys}

seqID, seqs , lengths, dct, total= parseFASTA('mobidb_all.fasta', aaFreq)
prob = [(c/total) for c in dct.values()]

t2 = time.time()
def calculate_distribution(length,arr_size):
    arr = [0] *arr_size
    i = 0
    while i < arr_size:
        new_count = npr.multinomial(length, prob)
        seq = generate_seq(new_count, keys)
        # seq = random_sequence(int(statistics.median(lengths)))
        d = runSearch(seq,10, top_one =True)
        arr[i] = d
        i+=1
    return arr


print("start")

data = []
# len_arr = [30,50,100,150,200,250,300]
len_arr = [100]
medians = []
npa = []
for i in len_arr:
    arr = calculate_distribution(i, 100)
    medians.append(statistics.median(arr))
    npa = np.asarray(arr, dtype=np.float32)
    data.append(npa)

# Attempt to calculate P-values using shuffling
# arr = [0]*500
# names = {}
# j = 0
# while j<500:
#     s ='MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNYGQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRGGRGRGGSGGGGGGGGGGYNRSSGGYEPRGRGGGRGGRGGMGGSDRGGFNKFGGPRDQGSRHDSEQDNSD'
#     l = list(s)
#     random.shuffle(l)
#     result = ''.join(l)
#     d, top_matches= runSearch(result, 1,top_one =True)
#     arr[j] = d

#     if top_matches[0] in names:
#         names[top_matches[0]] +=1
#     else:
#         names[top_matches[0]] =0
#     j+=1

t1 = time.time()
print(t1-t0)

query_dist, matches= runSearch("MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNYGQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRGGRGRGGSGGGGGGGGGGYNRSSGGYEPRGRGGGRGGRGGMGGSDRGGFNKFGGPRDQGSRHDSEQDNSD", 20)
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

# plt.hist(npa, bins=50)

plt.show()


#plot dist graph

plt.boxplot(data,showfliers=False, positions=len_arr, widths = 10)
x = np.linspace(0,300)
line_fit = np.polyfit(len_arr,medians,1)
len_arr_np = np.asarray(len_arr)
medians_np = np.asarray(medians)
print(line_fit)
plt.plot(len_arr_np, line_fit[0]* len_arr_np+line_fit[1])
plt.gca().set(title='Distance distribtuon vs length (all species)', ylabel='Minimum distance distribution', xlabel='sequence length (n=1000)')
plt.show()
