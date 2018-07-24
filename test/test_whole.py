#!/usr/bin/env python

import unittest
import os
import sys
import subprocess
import ihm.reader

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
        # Check output file
        self._check_mmcif_file('complement.cif')

    def _check_mmcif_file(self, fname):
        with open(fname) as fh:
            s, = ihm.reader.read(fh)
        self.assertEqual(len(s.citations), 1)
        self.assertEqual(s.citations[0].doi, '10.1074/mcp.M115.056473')
        self.assertEqual(len(s.software), 2)
        self.assertEqual(len(s.orphan_starting_models), 2)
        # Should be 3 states
        self.assertEqual(len(s.state_groups), 1)
        state1, state2, state3 = s.state_groups[0]
        self.assertEqual(state1.name,
                         'Human complement component C3')
        self.assertEqual(state2.name,
                         'Activated form C3b of human complement')
        self.assertEqual(state3.name,
                         'Hydrolyzed human complement, iC3 or C3N')
        # Check # of models in each state
        self.assertEqual(sum(len(x) for x in state1), 1)
        self.assertEqual(sum(len(x) for x in state2), 1)
        self.assertEqual(sum(len(x) for x in state3), 2)
        # Check # of spheres and atoms in each model
        models = [g[0] for g in state1+state2+state3]
        self.assertEqual([len(m._spheres) for m in models],
                         [1637, 1560, 1637, 1637])
        self.assertEqual([len(m._atoms) for m in models], [0, 0, 0, 0])
         # Should be 4 ensembles (clusters)
        self.assertEqual([e.num_models for e in s.ensembles],
                         [200, 200, 89, 111])
        # Check localization densities
        self.assertEqual([len(e.densities) for e in s.ensembles],
                         [11, 10, 11, 11])
        self.assertEqual([len(e.sequence) for e in s.entities], [645, 992])
        self.assertEqual([a.details for a in s.asym_units], ['beta', 'alpha'])
        # Just one restraint - crosslinks
        xl, = s.restraints
        self.assertEqual(len(xl.experimental_cross_links), 115)
        self.assertEqual(len(xl.cross_links), 115)
        self.assertEqual(xl.dataset.location.path,
                      'integrativemodeling-Complement-a6a1494/data/'
                      'QCLMS_iC3-Domain-Architecture_Rappsilber_TableS2-1.csv')
        self.assertEqual(xl.dataset.parents[0].location.access_code,
                         'PXD003486')
        # No psi/sigma available
        self.assertEqual(sum(len(x.fits) for x in xl.cross_links), 0)


if __name__ == '__main__':
    unittest.main()
