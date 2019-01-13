from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import random

hubs = ['NC_000913.3', 'NC_016810.1', 'NC_006155.1']

def overlap(start1, end1, start2, end2):
    """
    Does the range (start1, end1) overlap with (start2, end2)?
    """
    return end1 >= start2 and end2 >= start1

random_igrs = []

for hub in hubs:
    with open('/home/hdutcher/iwa/global_entero/tree_gbs/' + hub + '.gb', 'r') as infile:
        print('parsing genbank ', hub)
        igrs = []
        record = SeqIO.read(infile, 'genbank')
        genes = [f for f in record.features if f.type == 'gene']
        for i, g in enumerate(genes):
            if i != (len(genes) - 1): #exclude the last gene
                start1 = int(genes[i].location.start)
                end1 = int(genes[i].location.end)
                start2 = int(genes[i+1].location.start)
                end2 = int(genes[i+1].location.end)
                if overlap(start1, end1, start2, end2):
                    continue
                else:
                    if start2 - end1 >= 116:
                        flen = start2 - end1
                        offset = int(flen/2) - 58
                        igr_id = record.id + '/' + str(end1+offset) + '-' + str(start2-offset)
                        igr_seq = record.seq[end1+offset:start2-offset]
                        igr = SeqRecord(igr_seq, id=igr_id, description="")
                        igrs.append(igr)
                    else:
                        continue
        rand100 = random.sample(igrs, 100)
        random_igrs.extend(rand100)

ct = 0
for r in random_igrs:
    ct += 1
    with open(r.id.split('/')[0] + '_' + str(ct) + '.fa', 'w') as outfile:
        SeqIO.write(r, outfile, 'fasta')