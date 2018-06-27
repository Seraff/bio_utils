#!/usr/bin/env python3
# Author: Serafim Nenarokov, 2018

from optparse import OptionParser
import subprocess

import intervals as I
from Bio import SeqIO
import progressbar

BLAST_FORMAT = 'qseqid qlen sseqid slen length evalue pident bitscore mismatch gaps qstart qend sstart send'
BLAST_DELIMITER = ','

def parse_options():
    parser = OptionParser()
    parser.add_option("-f", "--fasta", help="Input fasta file (.fasta)")
    parser.add_option("-b", "--hits", help="Input BLAST hit file (.csv)")
    parser.add_option("-o", "--output", default="filtered.fasta", help="Filtered contigs file (.fasta)")
    parser.add_option("--threshold", default=15, help="Filter if hits coverage less than threshold (in %)")
    parser.add_option("-t", "--target", default="query", help="subject or query (query by default)")
    return parser.parse_args()[0]


def load_intervals_from_blast(path, target='query'):
    id_key = 'qseqid' if target == 'query' else 'sseqid'
    start_key = 'qstart' if target == 'query' else 'sstart'
    end_key = 'qend' if target == 'query' else 'send'

    blast_dict = {}
    with open(path) as f:
        for line in f:
            entry = parse_blast_entry(line)
            seqid = entry[id_key]
            interval = I.closed(int(entry[start_key]), int(entry[end_key]))

            if seqid not in blast_dict:
                blast_dict[seqid] = interval
            else:
                blast_dict[seqid] = blast_dict[seqid] | interval

    return blast_dict


def parse_blast_entry(entry):
    return dict(zip(BLAST_FORMAT.split(' '), entry.split(BLAST_DELIMITER)))


def calc_intervals_sum_len(intervals):
    len = 0
    for i in intervals:
        len += i.upper - i.lower + 1
    return len



if __name__ == '__main__':
    options = parse_options()

    cnt = int(subprocess.check_output('grep -c ">" %s' % options.fasta, shell=True))
    pb = progressbar.ProgressBar(max_value=cnt)

    intervals = load_intervals_from_blast(options.hits, target=options.target)
    fasta_sequences = SeqIO.parse(open(options.fasta), 'fasta')

    with open(options.output, 'w') as out_f:
        for i, record in enumerate(fasta_sequences):
            hits_len = 0

            if record.id in intervals:
                hits_len = calc_intervals_sum_len(intervals[record.id])

            pct = float(hits_len)/float(len(record)) * 100

            if pct <= float(options.threshold):
                SeqIO.write(record, out_f, "fasta")

            pb.update(i+1)
