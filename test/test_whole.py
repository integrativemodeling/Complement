#!/usr/bin/env python

import unittest
import os
import sys
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))

class Tests(unittest.TestCase):
    def test_mmcif(self):
        """Test generation of mmCIF output"""
        os.chdir(os.path.join(TOPDIR, 'c3-template'))
        if os.path.exists("complement.cif"):
            os.unlink("complement.cif")
        p = subprocess.check_call(
                ["python", "modeling.py", "--mmcif", "--dry-run"])
        # Check size of output file
        with open("complement.cif") as fh:
            wcl = len(fh.readlines())
        self.assertEqual(wcl, 23474)


if __name__ == '__main__':
    unittest.main()
