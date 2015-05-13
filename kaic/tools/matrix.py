'''
Created on May 11, 2015

@author: kkruse1
'''

def readMatrixFromFile(file_name, delim="\t"):
    with open(file_name, 'r') as f:
        M = [ map(float,line.split(delim)) for line in f ]
    return M



def writeMatrixToFile(M,file_name,delim="\t",row_names=None,col_names=None):
    # check if M is 2D matrix
    try:
        n_rows = len(M)
        n_cols = len(M[0])
    except IndexError:
        raise IndexError("Input must be 2D matrix")

    with open(file_name, 'w') as o:
        for name in col_names:
            o.write(name + delim)
        o.write("\n")
        
        for i in n_rows:
            if row_names:
                o.write(row_names[i] + delim)
            
            for j in n_cols:
                o.write("%.6E%s" + (M[i,j],delim))

    print "Done writing to file."