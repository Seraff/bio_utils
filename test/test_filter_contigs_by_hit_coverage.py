#!/usr/bin/env python3
# Author: Serafim Nenarokov, 2018

import unittest, subprocess, os, sys
from pathlib import Path

ROOT_PATH = str(Path(os.path.dirname(os.path.realpath(__file__))).parent)
sys.path.append(ROOT_PATH)

class DataManagerTest(unittest.TestCase):
    def setUp(self):
        self.tmp_fasta_path = 'tmp.fasta'
        self.tmp_csv_path = 'tmp.csv'
        self.tmp_out_path = 'filtered.fasta'

        self.command = f"../filter_contigs_by_hit_coverage.py "
        self.command += f"-f {self.tmp_fasta_path} "
        self.command += f"-b {self.tmp_csv_path}"

    def tearDown(self):
        for p in [self.tmp_csv_path, self.tmp_fasta_path, self.tmp_out_path]:
            subprocess.call(f"rm {p}", shell=True)

    def append_to_command(self, what):
        self.command += f" {what}"

    def call_command(self):
        subprocess.call(self.command, shell=True)

    def output_file_content(self):
        with open(self.tmp_out_path) as f:
            return f.read().strip()

    def add_contig(self, header, content):
        with open(self.tmp_fasta_path, 'a+') as fasta_f:
            fasta_f.write(f">{header}\n{content}\n")

    def add_hit(self, qseqid, start, finish):
        with open(self.tmp_csv_path, 'a+') as csv_f:
            csv_f.write(f"{qseqid},0,bla,0,0,0.001,80,0,0,0,{start},{finish},0,0\n")

    def test_b(self):
        self.add_contig('one', 'GTAGTACCGA')
        self.add_contig('two', 'GTCGTACGGA')

        self.add_hit("one", 1, 9)
        self.add_hit("two", 1, 8)

        threshold = 80

        self.append_to_command(f"--threshold {threshold}")
        self.append_to_command("--condition b")
        self.call_command()

        self.assertEqual(self.output_file_content(), ">one\nGTAGTACCGA")

    def test_l(self):
        self.add_contig('one', 'GTAGTACCGA')
        self.add_contig('two', 'GTCGTACGGA')

        self.add_hit("one", 1, 9)
        self.add_hit("two", 1, 5)

        threshold = 80

        self.append_to_command(f"--threshold {threshold}")
        self.append_to_command("--condition l")
        self.call_command()

        self.assertEqual(self.output_file_content(), ">two\nGTCGTACGGA")

    def test_e(self):
        self.add_contig('one', 'GTAGTACCGA')
        self.add_contig('two', 'GTCGTACGGA')

        self.add_hit("one", 1, 9)
        self.add_hit("two", 1, 5)

        threshold = 50

        self.append_to_command(f"--threshold {threshold}")
        self.append_to_command("--condition e")
        self.call_command()

        self.assertEqual(self.output_file_content(), ">two\nGTCGTACGGA")

    def test_be(self):
        self.add_contig('one', 'GTAGTACCGA')
        self.add_contig('two', 'GTCGTACGGA')
        self.add_contig('three', 'GTCGTACGGA')

        self.add_hit("one", 1, 9)
        self.add_hit("two", 1, 5)
        self.add_hit("three", 1, 6)

        threshold = 60

        self.append_to_command(f"--threshold {threshold}")
        self.append_to_command("--condition be")
        self.call_command()

        self.assertEqual(self.output_file_content(), ">one\nGTAGTACCGA\n>three\nGTCGTACGGA")

    def test_le(self):
        self.add_contig('one', 'GTAGTACCGA')
        self.add_contig('two', 'GTCGTACGGA')
        self.add_contig('three', 'GTCGTACGGA')

        self.add_hit("one", 1, 9)
        self.add_hit("two", 1, 5)
        self.add_hit("three", 1, 2)

        threshold = 50

        self.append_to_command(f"--threshold {threshold}")
        self.append_to_command("--condition le")
        self.call_command()

        self.assertEqual(self.output_file_content(), ">two\nGTCGTACGGA\n>three\nGTCGTACGGA")

if __name__ == '__main__':
    unittest.main()
