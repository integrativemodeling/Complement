# Complement

This repository contains the modeling files and the analysis related to the article: "Protein conformational changes that drives alternative complement activation illuminated by quantitative cross-linking/mass spectrometry" by Zhuo A. Chen et al. in Molecular Cell Proteomics 2016. The directory structure is the following:


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

`template` directories contain the `modeling.py` script, which is run simply by 

```
$IMP-BUILD-DIRECTORY/setup_environment.py python modeling.py
```

on a single core and

```
mpirun -np 16 $IMP-BUILD-DIRECTORY/setup_environment.py python modeling.py
```

on multiple cores (eg 16 in this case). The run produces files using the `IMP.pmi` higher interface. Refer to `IMP.tutorial` and`Nup84` examples for a better descriptions of the files.

`analysis` directories contain the scripts and the results of the analysis on the actual production runs (which are not released here).

`clustering.py` first script that needs to be run. It produces the directories with the corresponding structural clusters data. The output directory is `kmeans_weight_0_500_1`

`color_model.py` assign coded colors to a structure to finalise the image. The output file is `colored.rmf3`

`make_native_aligned.py` align the cluster structures against a given X-ray structure

`plot_cross_links.py` display the box plot for the crosslinks. The output file is `distances.pdf`

`rmsd.py` compute the root mean square distance of the cluster structures from the cluster center. The output file is `rmsd.out`.

`rmsf_precision.py` compute the root mean squared fluctuation of and the domain-wise precision of a cluster

`show_localization.py` chimera session script to display the localization densities with the right threshold

`xl_matrix.py` produces the contact map of the crosslinks, the output file is `XL_table.pdf`

