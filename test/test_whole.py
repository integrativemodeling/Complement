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
        # Potentially override methods that need network access
        env = os.environ.copy()
        env['PYTHONPATH'] = os.path.join(TOPDIR, 'test', 'mock') \
                            + ':' + env.get('PYTHONPATH', '')
        p = subprocess.check_call(
                ["python", "modeling.py", "--mmcif", "--dry-run"], env=env)
        # Check size of output file
        with open("complement.cif") as fh:
            wcl = len(fh.readlines())
        self.assertTrue(wcl >= 23470)


if __name__ == '__main__':
    unittest.main()
