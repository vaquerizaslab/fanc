'''
Created on Apr 7, 2015

@author: kkruse1
'''

import os
import subprocess

def simpleMap(fastq, 
              indexPath,
              outputSam,
              filterDuplicates=True,
              qualityCutoff=30):
    
    rs = ''.join(random.SystemRandom().choice(string.uppercase + string.digits) for _ in xrange(6))  # @UndefinedVariable
    tmpFilename = outputSam + '.' + rs + '.tmp'
    
    bowtieExecutablePath = subprocess.Popen("which bowtie2", shell=True, stdout=subprocess.PIPE).stdout.read().rstrip();
    
    subprocess.call([bowtieExecutablePath, '--very-sensitive', '--no-unal', '-x', indexPath, '-f', fastq, '-S', tmpFilename], shell=True)
    
    process = subprocess.Popen("head -n 500 " + tmpFilename + " | grep \"^@\" | wc -l", stdout=subprocess.PIPE, shell=True);
    
    nHeaderLines = process.communicate();
    
    if not filterDuplicates:
        #subprocess.call('{ head -n %d %s & tail -n +%d %s | sort -k1,1 -k5,5rn; } > %s' % (nHeaderLines, tmpFilename, nHeaderLines, tmpFilename, outputSam), shell=True);
        subprocess.call("{ head -n %d %s & tail -n +%d %s | awk \'{ if($5 > %d) print $0; }\' | sort -k1,1 -k5,5rn; } > %s" % (nHeaderLines, tmpFilename, nHeaderLines+1, tmpFilename, qualityCutoff, outputSam), shell=True);
    else:
        subprocess.call("{ head -n %d %s & tail -n +%d %s | awk \'{ duplicate=0; for(i=12; i<=NF; ++i) if($i ~ /^YS/) duplicate=1; if(duplicate == 0 && $1 !~ /^@/ && $2 != 4 && $5 > %d) print $0; }\' | sort -k1,1 -k5,5rn; } > %s" % (nHeaderLines, tmpFilename, nHeaderLines+1, tmpFilename, qualityCutoff, outputSam), shell=True);
        
    os.unlink(tmpFilename);
    