import os
import sys
import errno
import glob
import unittest
import subprocess
import contextlib
import collections

import facs
from facs.utils import helpers, galaxy

class DrassBasicTest(unittest.TestCase):
    """Build and query some simple bloom filters.
    """
    def setUp(self):
        self.progs = os.path.join(os.path.dirname(__file__), "bin")
        self.data_dir  = os.path.join(os.path.dirname(__file__), "data")
        self.reference = os.path.join(os.path.dirname(__file__), "data", "reference")
        self.bloom_dir = os.path.join(os.path.dirname(__file__), "data", "bloom")
        self.custom_dir = os.path.join(os.path.dirname(__file__), "data", "custom")
        self.synthetic_fastq = os.path.join(os.path.dirname(__file__), "data", "synthetic_fastq")

        self.fastq_nreads = [1, 8, 200]

        helpers._mkdir_p(self.data_dir)
        helpers._mkdir_p(self.reference)
        helpers._mkdir_p(self.bloom_dir)
        helpers._mkdir_p(self.custom_dir)
        helpers._mkdir_p(self.synthetic_fastq)

	# Check if 2bit decompressor is available
	if not os.path.exists(os.path.join(self.progs, "twoBitToFa")):
		_download_twoBitToFa_bin(self.progs)
		sys.path.append(self.progs)

        # Downloads reference genome(s)
        galaxy.rsync_genomes(self.reference, ["phix"])

    def test_1_build_ref(self):
        """ Build bloom filters out of the reference genomes directory.
        """
        for ref in os.listdir(self.reference):
            facs.build(os.path.join(self.reference, ref),
		os.path.join(self.bloom_dir, os.path.splitext(ref)[0]+".bloom"))

    def test_2_query(self):
        """ Generate dummy fastq files and query them against reference bloom filters.
        """
        for nreads in self.fastq_nreads:
            test_fname = "test%s.fastq" % nreads
            helpers._generate_dummy_fastq(os.path.join(self.synthetic_fastq, test_fname), nreads)

	for ref in self.reference:
		facs.query(os.path.join(self.synthetic_fastq, test_fname),
			   os.path.join(self.bloom_dir, os.path.splitext(ref)[0]+".bloom"))


    def test_3_query_custom(self):
        """ Query against the uncompressed FastQ files files manually deposited in data/custom folder.
        """
        for sample in glob.glob(os.path.join(self.custom_dir, "*.fastq")):
    	    print "\nQuerying against uncompressed sample %s" % sample
		for ref in self.reference:
			facs.query(os.path.join(self.synthetic_fastq, test_fname),
				   os.path.join(self.bloom_dir, os.path.splitext(ref)[0]+".bloom"))


    def test_4_query_custom_small_compressed(self):
	""" Query gzip compressed fastq files (less than 20MB).
	"""
        for sample in glob.glob(os.path.join(self.custom_dir, "*.fastq.gz")): 
    	    print "\nQuerying against compressed sample %s" % sample
            if os.path.getsize(os.path.join(self.custom_dir, sample)) < 20*1024*1204:
		facs.query(os.path.join(self.custom_dir, sample),os.path.join(self.bloom_dir, "U00096.2.bloom"))
