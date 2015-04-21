'''
Created on Apr 7, 2015

@author: kkruse1
'''
from __future__ import division

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
    
    
def filterContaminationLowMem(samSample, samContaminant, output):
    
    nOriginal = 0
    nFiltered = 0
    # collect all ID's of mapped contaminant
    with open(samContaminant, 'r') as sc:
        with open(samSample, 'r') as s:
            with open(output, 'w') as o:
                # skip headers
                line1 = s.readline()
                while line1 != '' and line1.startswith("@"):
                    o.write(line1)
                    line1 = s.readline()
                line2 = sc.readline()
                while line2 != '' and line2.startswith("@"):
                    line2 = sc.readline()
                
                
                while line1 != '' and line2 != '':
                    id1 = line1.split("\t")[0]
                    id2 = line2.split("\t")[0]
                    if id1 < id2:
                        o.write(line1)
                        line1 = s.readline()
                        nOriginal += 1
                        nFiltered += 1
                    elif id1 > id2:
                        line2 = sc.readline()
                    else:
                        line1 = s.readline()
                        nOriginal += 1
                        line2 = sc.readline()
                
                # rest must be uncontaminated
                while line1 != '':
                    o.write(line1)
                    line1 = s.readline()
                    nOriginal += 1
                    nFiltered += 1
                    
    print "Kept %d of %d reads (%.2f%%)" % (nFiltered, nOriginal, nFiltered/nOriginal*100);