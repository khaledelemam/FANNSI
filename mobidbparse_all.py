import re
import time

t0 = time.time()

# Parse MobiDB files

#restrictions: grab only th_50 and idr length > 20

def getIDRS(seq, nums):
    counter = 0
    idr = ""
    allIdrs = []
    for i, j in zip(range(len(seq)), range(len(nums))):
        if int(nums[j]) == 1:
            idr += seq[i]
            counter+=1
        elif int(nums[j]) == 0:
            if counter > 20:
                # print(j)
                # print(counter)
                allIdrs.append(idr)
                idr = ""
                counter = 0
            else:
                idr = ""
                counter = 0
    if counter > 20:
        allIdrs.append(idr)
    return allIdrs


count_c = 0
count_b = 1
all ={}
with open('all.fasta') as f:
    line = f.readline()
    # print(line)
    regexID= re.findall("^(.*?)\|", line)
    tempID = regexID[0]
    seq = f.readline().strip()
    # allSeqs[tempID] = {"seq": seq, "num": ""}
    all[tempID] = []
    # seq = f.readline()
    # print(regexID)
    flag = True
    while line:
        line = f.readline()
        regex50 = re.findall("th_50", line)
        regexID= re.findall("^(.*?)\|", line)
        # print(regexID)

        #found binary sequence
        if regex50 == ["th_50"] and regexID[0] == tempID:
            count_c+=1
            nums = f.readline().strip()
            # print(nums)
            # allSeqs[regexID[0]]["num"] = nums
            # print(line)
            # print(getIDRS(seq,nums))
            idrs = getIDRS(seq,nums)
            # if len of idr > 20
            # if len(idrs) > 0:
            all[tempID] = idrs
            flag = True
            # else:
            #     del all[tempID]
    
        #new id
        elif len(regexID) > 0:
            # print (all)
            if regexID[0] != tempID:
            
                count_b+=1

                # if th_50 does not exist or no idrs
                if len(all[tempID])  == 0:
                    del all[tempID]         
                tempID = regexID[0]
                seq = f.readline().strip()
                all[tempID] = []   

if len(all[tempID])  == 0:
     del all[tempID] 

# print(allSeqs)
# print(all)
# ["lol", "kay"]
# all[">P31946"] = ["lol", "kay"]
t11 = time.time()

tt = t11 - t0
print(tt)
count = 0
with open("mobidb_all.fasta",'w+') as f:
    for id, seq in all.items():
        count+=1
        if len(seq) >= 1:
            for i in range(len(seq)):
                idnew = id + "_" + str(i+1)
                f.write(idnew + '\n')
                f.write(seq[i] + '\n')
        # else:
        #     f.write(id + '\n')
        #     f.write(seq[0] + '\n')

       

t1 = time.time()
t = t1 - t0
# tt = t11 - t0
print(t)
# print(tt)
print("chosen IDRs: ",count)
print("all IDRs: ", count_b)
print("th_50 absent: ",count_c)








#fus
# 1- 286
# 368 - 432
# 437 - 526
#['MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNYGQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRGGRGRGGSGGGGGGGGGGYNRSSGGYEPRGRGGGRGGRGGMGGSDRGGFNKFGGPRDQGSRHDSEQDNSDNNT', 'ATRRADFNRGGGNGRGGRGRGGPMGRGGYGGGGSGGGGRGGFPSGGGGGGGQQRAGDWKCPNPT', 'FSWRNECNQCKAPKPDGPGGGPGGSHMGGNYGDDRRGGRGGYDRGGYRGRGGDRGGFRGGRGGGDRGGFGPGKMDSRGEHRQDRRERPY']

#uniprot
# 1 - 286
# 375 – 424
# 444 – 526

#current db
# >P35637_1
#'MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNYGQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRGGRGRGGSGGGGGGGGGGYNRSSGGYEPRGRGGGRGGRGGMGGSDRGGFNKFGGPRDQGSRHDSEQDNSDNNT'
# myresults: MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNYGQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRGGRGRGGSGGGGGGGGGGYNRSSGGYEPRGRGGGRGGRGGMGGSDRGGFNKFGGPRDQGSRHDSEQDNSD
# >P35637_2
# ADFNRGGGNGRGGRGRGGPMGRGGYGGGGSGGGGRGGFPSGGGGGGGQQRAGDWKCPNPTCENMNFSWRNECNQCKAPKPDGPGGGPGGSHMGGNYGDDRRGGRGGYDRGGYRGRGGDRGGFRGGRGGGDRGGFGPGKMDSRGEHRQDRRERPY
# myresults: ATRRADFNRGGGNGRGGRGRGGPMGRGGYGGGGSGGGGRGGFPSGGGGGGGQQRAGDWKCPNPT   'FSWRNECNQCKAPKPDGPGGGPGGSHMGGNYGDDRRGGRGGYDRGGYRGRGGDRGGFRGGRGGGDRGGFGPGKMDSRGEHRQDRRERPY

# test
#1-237
#267-294
#321-592 

#uniprot taf-15
#1 – 237	
#324 – 356	
#373 – 592	

# >Q92804_1
# MSDSGSYGQSGGEQQSYSTYGNPGSQGYGQASQSYSGYGQTTDSSYGQNYSGYSSYGQSQSGYSQSYGGYENQKQSSYSQQPYNNQGQQQNMESSGSQGGRAPSYDQPDYGQQDSYDQQSGYDQHQGSYDEQSNYDQQHDSYSQNQQSYHSQRENYSHHTQDDRRDVSRYGEDNRGYGGSQGGGRGRGGYDKDGRGPMTGSSGGDRGGFKNFGGHRDYGPRTDADSESDNSDNNTIFVQGLGEGVSTDQVGEFFKQIGIIKTNKKTGKPMINLYTDKDTGKPKGEATVSFDDPPSAKAAIDWFDGKEFHGNI
# >Q92804_2
# RPEFMRGGGSGGGRRGRGGYRGRGGFQGRGGDPKSGDWVCPNPSCGNMNFARRNSCNQCNEPRPEDSRPSGGDFRGRGYGGERGYRGRGGRGGDRGGYGGDRSGGGYGGDRSSGGGYSGDRSGGGYGGDRSGGGYGGDRGGGYGGDRGGGYGGDRGGGYGGDRGGYGGDRGGGYGGDRGGYGGDRGGYGGDRGGYGGDRGGYGGDRSRGGYGGDRGGGSGYGGDRSGGYGGDRSGGGYGGDRGGGYGGDRGGYGGKMGGRNDYRNDQRNRPY