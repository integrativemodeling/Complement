[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1285940.svg)](https://doi.org/10.5281/zenodo.1285940)

# Complement

This repository contains the modeling files and the analysis related to the
article ["Structure of Complement C3(H2O) Revealed By Quantitative
Cross-Linking/Mass Spectrometry And Modeling"](https://www.ncbi.nlm.nih.gov/pubmed/27250206)
by Zhuo A. Chen et al. in Molecular Cell Proteomics 2016. The directory
structure is the following:


```
c3-template
c3b-template
ic3-template
c3-analysis
c3b-analysis
ic3-analysis
data
```

The directories are organized by system, thereby `c3`, `c3b` and `ic3` correspond to the three different states of the complement, as discussed in the article.

`template` directories contain the
[IMP](https://integrativemodeling.org)
`modeling.py` script, which is run simply by

```
python modeling.py
```

on a single core and

```
mpirun -np 16 python modeling.py
```

on multiple cores (eg 16 in this case). The run produces files using the
[PMI](https://github.com/salilab/pmi) PMI high-level interface.
Refer to the [IMP tutorial](https://integrativemodeling.org/nightly/doc/manual/rnapolii_stalk.html)
and [Nup84](https://salilab.org/nup84) for fuller descriptions of the files.

`analysis` directories contain the scripts and the results of the analysis
on the actual production runs (which are available as `traj-*.tar.xz`
[at Zenodo](https://doi.org/10.5281/zenodo.1285940) - note that two independent
runs were carried out for each state).

`clustering.py` is the first script that needs to be run. It produces directories with the corresponding structural cluster data. The output directory is `kmeans_weight_0_500_1`.

`color_model.py` assigns coded colors to a structure to finalize the image. The output file is `colored.rmf3`.

`make_native_aligned.py` aligns the cluster structures against a given X-ray structure.

`plot_cross_links.py` displays the box plot for the crosslinks. The output file is `distances.pdf`.

`rmsd.py` computes the root mean square distance of the cluster structures from the cluster center. The output file is `rmsd.out`.

`rmsf_precision.py` computes the root mean squared fluctuation of and the domain-wise precision of a cluster.

`show_localization.py` is a [Chimera](http://www.cgl.ucsf.edu/chimera/)
 session script to display the localization densities with the right threshold.

`xl_matrix.py` produces the contact map of the crosslinks; the output file is `XL_table.pdf`.

## Information

_Author(s)_: Riccardo Pellarin

_Date_: August 2016

_License_: [CC-BY-SA-4.0](https://creativecommons.org/licenses/by-sa/4.0/legalcode).
This work is freely available under the terms of the Creative Commons
Attribution-ShareAlike 4.0 International License.

_Last known good IMP version_: [![build info](https://integrativemodeling.org/systems/20/badge.svg?branch=master)](https://integrativemodeling.org/systems/) [![build info](https://integrativemodeling.org/systems/20/badge.svg?branch=develop)](https://integrativemodeling.org/systems/)

_Publications_:
 - Chen ZA, Pellarin R, Fischer L, Sali A, Nilges M, Barlow PN, Rappsilber J.,
   [Structure of Complement C3(H2O) Revealed By Quantitative Cross-Linking/Mass Spectrometry And Modeling](https://www.ncbi.nlm.nih.gov/pubmed/27250206), Mol Cell Proteomics, 2016, 10.1074/mcp.M115.056473.
