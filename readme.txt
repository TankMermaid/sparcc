
********************************
General Notes:
********************************
- SparCC is mainly used for calculating correlations in compositional data.

- Additionally, it contains a script for calculating the distance between samples using the JSD metric, its square-root, and many other distance measures.

- SparCC has some other useful capabilities (e.g. calculating sample diversity, rank abundance plots, etc'). However, there are no convenience wrappers for these capabilities, and they are not as well tested. If you known your python (especially, object oriented python), and want to use these, take a peek at the MatrixDictionary and SurveyMatrix classes in the 'lib' subdirectory.

- This module is still under development, and may not work as advertised.  

- Questions, comments, complaints and praise should be send to yonatanf@mit.edu



********************************
Usage Notes:
********************************
- Scripts in the root SparCC directory can be called from the terminal command-line either by explicitly calling python (as is done in the usage examples below), or simply as an executable. The latter will require having execution permission for these file (e.g. chmod +x SparCC.py).

- Help for any one for the scripts in the root SparCC directory is available by typing 'python [script_name] - h' in the command line. e.g.: :: 

   python SparCC.py -h .

- SparCC is implemented in pure python and requires a working version of python (=>2.3, tested with 2.6.6) and numpy (tested with versions 1.4.0 and 1.6.0).

- The optional distance calculations also require the scipy.cluster.hierarchy module from scipy (test with version 0.9.0).
       

********************************
Usage example:
********************************
- The following lists the commands required for analyzing the included 'fake' dataset using the SparCC package, and generating all the files present in the subfolders of the example folder.

- The fake dataset contains simulated abundances of 50 otus in 200 samples, drawn at random from a multinomial log-normal distribution. The true basis correlations used to generate the data are listed in 'true_basis_cor.txt' in the example folder.

- Note that otu 0 is very dominant, and thus, using Pearson or Spearman correlations, appears to be negatively correlated with most other OTUs, though it is in fact not negatively correlated with any OTU.

- In the following, commands are enclosed between quotes "". These lines can be copied (don't copy the quotes) to the command-line and should run. 

---------------------------------
Correlation Calculation:
---------------------------------
First, we'll quantify the correlation between all OTUs, using SparCC, Pearson, and Spearman correlations:

::

   python SparCC.py example/fake_data.txt -i 5 --cor_file=example/basis_corr/cor_sparcc.out
   python SparCC.py example/fake_data.txt -i 5 --cor_file=example/basis_corr/cor_pearson.out -a pearson
   python SparCC.py example/fake_data.txt -i 5 --cor_file=example/basis_corr/cor_spearman.out -a spearman


---------------------------------
Pseudo p-value Calculation:
---------------------------------
Calculating pseudo p-values is done via a bootstrap procedure.
First make shuffled (w. replacement) datasets:
::

   python MakeBootstraps.py example/fake_data.txt -n 5 -o example/pvals/boot

This will generate 5 shuffled datasets, which is clearly not enough to get meaningful p-values, and is used here for convenience.
A more appropriate number of shuffles should be at least a 100, which is the default value. 

Next, you'll have to run SparCC on each of the shuffled data sets. 
Make sure to use the exact same parameters which you used when running SparCC on the real data, name all the output files consistently, numbered sequentially, and with a '.txt' extension.
::

   python SparCC.py example/pvals/boot_0.txt -i 5 --cor_file=example/pvals/sim_cor_0.txt
   python SparCC.py example/pvals/boot_1.txt -i 5 --cor_file=example/pvals/sim_cor_1.txt
   python SparCC.py example/pvals/boot_2.txt -i 5 --cor_file=example/pvals/sim_cor_2.txt
   python SparCC.py example/pvals/boot_3.txt -i 5 --cor_file=example/pvals/sim_cor_3.txt
   python SparCC.py example/pvals/boot_4.txt -i 5 --cor_file=example/pvals/sim_cor_4.txt

Above I'm simply called SparCC 5 separate times. However, it is much more efficient and convenient to write a small script that automates this, and submits these runs as separate jobs to a cluster (if one is available to you. Otherwise, this may take a while to run on a local machine...).

Now that we have all the correlations computed from the shuffled datasets, we're ready to get the pseudo p-values.
Remember to make sure all the correlation files are in the same folder, are numbered sequentially, and have a '.txt' extension.
The following will compute both one and two sided p-values.
::

   python PseudoPvals.py example/basis_corr/cor_sparcc.out example/pvals/sim_cor 5 -o example/pvals/pvals_one_sided.txt -t 'one_sided'
   python PseudoPvals.py example/basis_corr/cor_sparcc.out example/pvals/sim_cor 5 -o example/pvals/pvals_two_sided.txt -t 'two_sided'


---------------------------------
Sample distances:
---------------------------------
Another common task is to compute all pairwise distances between samples. This is useful when clustering the samples (e.g. using UPGMA, k-means, etc'), and when performing dimension reductions (e.g. using metric or non-metric Multi Dimensional Scaling = PCoA).
The code below calculates the distance matrix between all samples of the 'fake' data set using Euclidean distance, and the square-root of the Jensen-Shannon Divergence. 
::

   python SampleDist.py example/fake_data.txt -o example/sample_dist/sample_dist_JSsqrt.out
   python SampleDist.py example/fake_data.txt -m euclidean -o example/sample_dist/sample_dist_euclidean.out

