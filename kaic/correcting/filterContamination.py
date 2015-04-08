'''
Created on Apr 7, 2015

@author: kkruse1
'''


def filterContamination(samSample, samContaminant, output):
    # collect all ID's of mapped contaminant
    ids = set();
    with open(samContaminant, 'r') as s:
        for x in s:
            x = x.rstrip()
            fields = x.split("\t")
            if len(fields) > 0 and not fields[0].startswith("@"):
                readId = fields[0]
                ids.add(readId)
                
    print "Number of unique IDs in contaminant: ", len(ids)
    
    nOriginal = 0
    nFiltered = 0
    with open(samSample, 'r') as s:
        with open(output, 'w') as o:
            for x in s:
                #x = x.rstrip();
                nOriginal+=1
                fields = x.split("\t")
                if len(fields) > 0 and fields[0] not in ids:
                    o.write(x)
                    nFiltered+=1
    
    print "Kept %d of %d reads " % (nFiltered, nOriginal);