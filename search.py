import os
import pandas as pd
from helpers import*
from run_features import*
import time
import pickle
import hnswlib
import numpy as np
import re
import math

def runSearch(seq, threshold= 10, top_one = False):
    
    #if feature vector database does not exist
    if not os.path.isfile("mobiALL_db.pkl"):
        db_mol_feats= run_feats('mobidb_all.fasta')
        all_IDRs = parseFeaturesDict(db_mol_feats)
        all_IDRs_names = all_IDRs.index
        pickle.dump(all_IDRs_names, open( "mobiALL_IDRs_names.pkl", "wb" ))
        df_all_IDRs_standardized, transformation = z_score(all_IDRs)
        pickle.dump(transformation,open( 'mobiALL_transformation.pkl', 'wb'))
        pickle.dump(df_all_IDRs_standardized,open("mobiALL_db.pkl", "wb"))


    if not os.path.isfile("mobiALL_hnsw.pkl"):
        dim = 94
        num_elements = 325415
        arr_df= pickle.load(open('mobiALL_db.pkl', 'rb'))
        arr_db = arr_df.to_numpy()

        #  data
        data = np.float32(arr_db)
        data_labels = np.arange(num_elements)

        # Declaring index
        p = hnswlib.Index(space = 'l2', dim = dim) # possible options are l2, cosine or ip

        # Initializing index - the maximum number of elements should be known beforehand
        p.init_index(max_elements = num_elements, ef_construction = 50, M = 64)

        # Element insertion (can be called several times):
        p.add_items(data, data_labels)

        # Controlling the recall by setting ef:
        p.set_ef(300) # ef should always be > k

        pickle.dump(p, open('mobiALL_hnsw.pkl', 'wb'))


    #download all of Uniprot annotations
    
    # if not os.path.isfile("uniprot_Annotations.pkl"):
    #     names = pickle.load(open('mobiALL_IDRs_names.pkl', 'rb'))
    #     temp_names = [re.sub('[0-9]$|[_>]', '', i) for i in list(names)]
    #     clean_names = set(temp_names)
    #     uniprot_output = ''
    #     for val in clean_names:
    #         uniprot_output = uniprot_output + get_protein_sequences([val])
    #     # uniprot_output = get_protein_sequences(clean_names)
    #     uniprot_Regex = re.findall(">(.*?)PE=", uniprot_output)
    #     annotations= {key:[] for key in clean_names}
    #     for i in uniprot_Regex:
    #         clean = re.findall('\|(.*)\|', i)
    #         word = i.split("|")
    #         annotations[clean[0]].append(word[-1].partition(" ")[-1])
    #     links_list = get_uniprotlink(clean_names)
    #     for i in range(len(links_list)):
    #         annotations[clean_names[i]].append(links_list[i])
    #     pickle.dump(annotations, open('uniprot_Annotations.pkl', 'wb'))
        
   
    #calculate feature vector for query seqeucnce
    query_mol_feats = run_feats(seq)
    query = parseFeaturesDict(query_mol_feats)
    transformation = pickle.load(open('mobiALL_transformation.pkl', 'rb'))

    #normalize feature vector
    query_standardized = transform(query, transformation)
    arr_query = query_standardized.to_numpy()
    
    #find closest k vectors
    p = pickle.load(open('mobiALL_hnsw.pkl','rb'))
    labels, distances = p.knn_query(arr_query, k =threshold)

    IDRs_names = pickle.load(open('mobiALL_IDRs_names.pkl', 'rb'))

    top_matches = [IDRs_names[labels[0][i]] for i in range(len(labels[0]))]
    # dists = [math.sqrt(distances[0][i]) for i in range(len(distances[0]))]

    clean_matches = [re.sub('[0-9]$|[_>]', '', match) for match in top_matches]
    values = {key:[] for key in clean_matches}
      
    #get Uniprot annotations
    uniprot_seq = get_protein_sequences(clean_matches)
    uniprot_seqRegex = re.findall(">(.*?)PE=", uniprot_seq)    

    for i in uniprot_seqRegex:
        clean = re.findall('\|(.*)\|', i)
        word = i.split("|")
        values[clean[0]].append(word[-1].partition(" ")[-1])
    
    #get Uniprot links
    links_list = get_uniprotlink(clean_matches)
    for i in range(len(links_list)):
             values[clean_matches[i]].append(links_list[i])
  
    return values


# output = runSearch('MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNYGQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRGGRGRGGSGGGGGGGGGGYNRSSGGYEPRGRGGGRGGRGGMGGSDRGGFNKFGGPRDQGSRHDSEQDNSD',10)
# for i in output.values():
#     print(i)