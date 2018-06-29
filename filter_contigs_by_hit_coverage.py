#!/usr/bin/env python3
# Author: Serafim Nenarokov, 2018

import re, textwrap
import argparse
import subprocess

import intervals as I
from Bio import SeqIO
import progressbar

BLAST_DEFAULT_FORMAT = 'qseqid qlen sseqid slen length evalue pident bitscore mismatch gaps qstart qend sstart send'
BLAST_DEFAULT_DELIMITER = ','


def parse_options():
    usage = "./filter_contigs_by_hit_coverage.py"
    description = """
    The script filters contigs in fasta file by the BLAST 
    hits coverage and provided criteria: threshold and conditions.
    \n\n
    Example: ./filter_contigs_by_hit_coverage.py -f ./test.fasta -b --threshold 20 ./test.csv --condition le
    This command will calculate hit coverage of contigs from ./test.fasta 
    and extract them if the following condition is satisfied: hit_coverage <= 20%
    \n\n
    Example: ./filter_contigs_by_hit_coverage.py -f ./test.fasta -b --threshold 50 ./test.csv --condition b
    Condition for contigs extracted: hit_coverage > 20%
    """

    description = re.sub(r"(?<=\w)\n(?=\w)", "\n", description)

    parser = argparse.ArgumentParser(usage, description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-f", "--fasta", help="Input fasta file (.fasta)")
    parser.add_argument("-b", "--hits", help="Input BLAST hit file (.csv)")
    parser.add_argument("-o", "--output", default="filtered.fasta", help="Filtered contigs file (.fasta)")
    parser.add_argument("--threshold", help="Filter if hits coverage according to threshold (in %%)")
    parser.add_argument("--condition", help="""Condition for the threshold. Possible values:
                           b - bigger (>)
                           l - less (<)
                           e - equal (==)
                           be - bigger or equal (>=)
                           le - less or equal (<=)""")
    parser.add_argument("-t", "--target", default="query", help="subject or query (query by default)")
    parser.add_argument("--blast_delimiter", default=BLAST_DEFAULT_DELIMITER,
                        help=f"delimiter of BLAST file (default: \"{BLAST_DEFAULT_DELIMITER}\")")
    parser.add_argument("--blast_format", default=BLAST_DEFAULT_FORMAT,
                        help=f"BLAST file format (default: \"{BLAST_DEFAULT_FORMAT}\")")
    return parser.parse_args()


def load_intervals_from_blast(path, format, delimiter, target='query'):
    id_key = 'qseqid' if target == 'query' else 'sseqid'
    start_key = 'qstart' if target == 'query' else 'sstart'
    end_key = 'qend' if target == 'query' else 'send'

    blast_dict = {}
    with open(path) as f:
        for line in f:
            entry = parse_blast_entry(line, format, delimiter)
            seqid = entry[id_key]
            interval = I.closed(int(entry[start_key]), int(entry[end_key]))

            if seqid not in blast_dict:
                blast_dict[seqid] = interval
            else:
                blast_dict[seqid] = blast_dict[seqid] | interval

    return blast_dict


def parse_blast_entry(entry, format, delimiter):
    return dict(zip(format.split(' '), entry.split(delimiter)))


def calc_intervals_sum_len(intervals):
    len = 0
    for i in intervals:
        len += i.upper - i.lower + 1
    return len


if __name__ == '__main__':
    options = parse_options()

    cnt = int(subprocess.check_output('grep -c ">" %s' % options.fasta, shell=True))
    pb = progressbar.ProgressBar(max_value=cnt)

    intervals = load_intervals_from_blast(options.hits, options.blast_format, options.blast_delimiter,
                                          target=options.target)
    fasta_sequences = SeqIO.parse(open(options.fasta), 'fasta')

    with open(options.output, 'w') as out_f:
        for i, record in enumerate(fasta_sequences):
            hits_len = 0

            if record.id in intervals:
                hits_len = calc_intervals_sum_len(intervals[record.id])
            else:
                continue

            pct = float(hits_len) / float(len(record)) * 100
            threshold = float(options.threshold)
            satisfies = False

            if options.condition == 'b':
                satisfies = pct > threshold
            elif options.condition == 'be':
                satisfies = pct >= threshold
            elif options.condition == 'e':
                satisfies = pct == threshold
            elif options.condition == 'l':
                satisfies = pct < threshold
            elif options.condition == 'le':
                satisfies = pct <= threshold

            if satisfies:
                SeqIO.write(record, out_f, "fasta")

            pb.update(i + 1)
