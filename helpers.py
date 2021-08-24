import math
import pandas as pd
from collections  import OrderedDict
import urllib.request
import ssl
import certifi


def parseFeaturesDict(featdict):
    '''
    Takes in a dictionary with key as protein ID (if exist) 
    nd values as a map of feature names to feature values

    Returns a dataframe of feature vectors

    '''
    IDs=[key for key in featdict.keys()]
    first_entry=IDs[0]
    featureNames = [entry for entry in featdict[first_entry]]
    all_values =[]
    for seqid,feats in featdict.items():
        values = []
        for featid,featval in feats.items():
            values.append(featval)
        all_values.append(values)
    df = pd.DataFrame(all_values, columns=featureNames,index = IDs)
    return df
  

def eucdis(v1, v2):
    '''
    Takes in 2 feature vectors
    Returns the euclidean distance

    '''
    dist = [(a - b)**2 for a, b in zip(v1, v2)]
    dist = math.sqrt(sum(dist))
    return dist

def rank_features(v1,v2,featNames):
    
    feats = [abs((a - b)) for a, b in zip(v1, v2)]
    feats_dict = dict(zip(featNames,feats))
    return feats_dict

def get_similar_features(data, n=5, order=True):
    """Get top n features by difference. 

    Returns a dictionary or an `OrderedDict` if `order` is true.
    """ 
    top = sorted(data.items(), key=lambda x: x[1])[:n]
    if order:
        return OrderedDict(top)
    return dict(top)


def z_score(df):
    '''
    Takes in a dataframe of feature vectors and returns:
    1: normalized dataframe
    2: mean and std dev used for transformation
    '''
    # copy the dataframe
    # df_std = df.copy()

    transformation = []
    # apply the z-score method
    for column in df.columns:
        transformation.append((df[column].mean(), df[column].std()))
        if df[column].std() != 0: #don't divide by 0
            df[column] = (df[column] - df[column].mean()) / df[column].std()
    
    return df, transformation


def transform(df, transformation, idx =[]):
    '''
    Takes in dataframe of one sequence and tranformation vector (mean,std dev)
    Returns normalized dataframe
    '''
    # copy the dataframe
    # df_std = df.copy()

    ## search for sepcific features
    if idx != []:
        # print("HERE")
        i= idx.pop(0)
        for column in df.columns:
            # print(i)
            if transformation[i][1] != 0: #don't divide by 0
                df[column] = (df[column] - transformation[i][0]) / transformation[i][1]
            if len(idx) != 0:
                i = idx.pop(0)
            else:
                break

    else:
        i= 0
        for column in df.columns:
            if transformation[i][1] != 0: #don't divide by 0
                df[column] = (df[column] - transformation[i][0]) / transformation[i][1]
            i+=1
    return df


def calculate_dist(arr_db, arr_query):
    '''
    Takes in 2 feature vectors
    Returns list with euclidean distances
    '''
    all_dist = []
    for i in range(len(arr_db)):
        dist = eucdis(arr_query[0],arr_db[i])
        all_dist.append(dist)
    return all_dist


def get_protein_sequences(uniprot_list):
    """Retrieves the sequences from the UniProt database based on the list of
    UniProt ids.
    In general, 
        1. Compose your query here with the advanced search tool:
    https://www.uniprot.org/uniprot/?query=id%3Ap40925+OR+id%3Ap40926+OR+id%3Ao43175&sort=score
        2. Replace `&sort=score` with `&format=fasta`
        3. Edit this function as necessary
    Returns:
        string of fasta annotation
    """
    # This makes it so we match only the ENTRY field
    uniprot_list = ['id%3A'+id for id in uniprot_list]
    line = '+OR+'.join(uniprot_list)
    url = f'https://www.uniprot.org/uniprot/?query={line}&format=fasta'
    with urllib.request.urlopen(url, context=ssl.create_default_context(cafile=certifi.where())) as f:
        fasta = f.read().decode('utf-8').strip()
    return fasta

def get_uniprotlink(uniprot_list):
    """
    """
    uniprot_list = [id for id in uniprot_list]
    links_list = [f'https://www.uniprot.org/uniprot/{id}' for id in uniprot_list]
    return links_list

def parseFASTA(file):
    '''
    Takes in a fasta file
    Returns  2 lists
    List 1: protein seq ID
    List 2: protein sequences

    '''
    seqID =[]
    seqs=[]
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
                    seqs.append(line.strip())
            count +=1
        return seqID, seqs

def parseFeatures(file) -> list:
    '''
    Takes in a feature vector text file and returns a list of 3 lists
    List 1: feature values (float)
    List 2: feature names
    List 3: protein IDs
    '''
    with open(file) as f:
        IDs = []
        featureNames = []
        values_vector = []
        all_vectors =[]
        line = f.readline() #first line (features names)
        line_clean = line.split("\t")
        temp = list(filter(None,  line_clean))
        featureNames.extend(temp)
        while line:
            line = f.readline()
            if line.strip():  # ignore blank lines
                str = line.split("\t")
                ID = str.pop(0) #remove ID 
                IDs.append(ID[1:]) #remove "<"
                nums = [float(numeric_str) for numeric_str in str]
                values_vector.append(nums)
        all_vectors.append(values_vector)
        all_vectors.append(featureNames)
        all_vectors.append(IDs)
        return all_vectors