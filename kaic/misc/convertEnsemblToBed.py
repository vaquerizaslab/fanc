'''
Created on Aug 28, 2015

@author: kkruse1
'''
import argparse
import os.path

if __name__ == '__main__':
    parser = argparse.ArgumentParser();
    
    parser.add_argument(
        'input',
        help='''Ensembl tab-separated data file'''
    );
    
    parser.add_argument(
        'output',
        help='''Folder to generate the individual BED files in'''
    );
    
    args = parser.parse_args()
    in_file = os.path.expanduser(args.input)
    out_folder = os.path.expanduser(args.output)
    
    try:
        os.makedirs(out_folder)
    except OSError:
        pass
    
    with open(in_file, 'r') as f:
        o = None
        mode = 0
        is_first_line = False
        previous_feature = ''
        previous_line = None
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('seqname'):
                mode = 1
                header = line.split('\t')
            elif line.startswith('SeqRegion'):
                mode = 2
                header = line.split('\t')
                file_name = "_".join(previous_line.split())
                o = open(out_folder + '/' + file_name + '.bed', 'w')
                if previous_line.lower().startswith('gene'):
                    mode = 3
                    out_header = 'chrom\tstart\tend\tgene\tscore\tstrand\ttype'
                else:
                    out_header = 'chrom\tstart\tend\tscore\tstrand\tname'
                ignore = set([0,1,2])
                for idx, field in enumerate(header):
                    if idx not in ignore and field.lower() != 'name':
                        out_header += "\t%s" % ''.join(field.split())
                o.write(out_header+ '\n')
            elif line == '':
                o.close()
                mode = 0
            else:
                fields = line.split('\t')
                if mode == 1:
                    if fields[2] != previous_feature:
                        o = open(out_folder + '/' + fields[2] + '.bed', 'w')
                        previous_feature = fields[2]
                        out_header = 'chrom\tstart\tend\tscore\tstrand\tname'
                        ignore = set([0,3,4,5,6,8])
                        for idx, field in enumerate(header):
                            if idx not in ignore:
                                out_header += "\t%s" % field
                        o.write(out_header+ '\n')
                    chromosome = fields[0]
                    start = fields[3]
                    end = fields[4]
                    score = fields[5]
                    strand = fields[6]
                    name = fields[8]
                    
                    out_line = "%s\t%s\t%s\t%s\t%s\t%s" % (chromosome, start, end, score, strand, name)
                    
                    ignore = set([0,3,4,5,6,8])
                    for idx, field in enumerate(fields):
                        if idx not in ignore:
                            out_line += "\t%s" % field
                    
                    o.write(out_line + "\n")
                elif mode == 2:
                    chromosome = fields[0]
                    start = fields[1]
                    end = fields[2]
                    score = 1
                    strand = '+'
                    
                    out_line = "%s\t%s\t%s\t%s\t%s" % (chromosome, start, end, score, strand)
                    
                    end_line = ''
                    ignore = set([0,1,2])
                    for idx, field in enumerate(fields):
                        if header[idx].lower() == 'name':
                            out_line += '\t%s' % field
                        elif idx not in ignore:
                            end_line += "\t%s" % field
                    
                    o.write(out_line + end_line + "\n")
                elif mode == 3:
                    chromosome = fields[0]
                    start = fields[1]
                    end = fields[2]
                    score = 1
                    strand = 1

                    out_line = "%s\t%s\t%s" % (chromosome, start, end)

                    end_line = '\t%s\t%s\texon' % (score, strand)
                    ignore = set([0,1,2])
                    for idx, field in enumerate(fields):
                        if header[idx].lower() == 'name':
                            out_line += '\t%s' % field
                        elif idx not in ignore:
                            end_line += "\t%s" % field

                    o.write(out_line + end_line + "\n")
            previous_line = line
        
        o.close()
