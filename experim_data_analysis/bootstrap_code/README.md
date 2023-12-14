These codes have been written by Francesco Mambretti (2021-2023).
They are meant to analyze experimental FASTQ files from the SEDES experiment.
This version performs bootstrap, i.e. replicated analyses of subsamples coming from the same dataset.

------------------------------------- REQUIREMENTS -------------------------------------

- `python3` with `numpy, matplotlib, itertools, more_itertools, biopython, pandas, difflib`
- `C++` installed and a `C++` compiler supporting (at least) `C++ - 2011`
- `bash`

------------------------------------- input_params.py -------------------------------------

1) first, modify `input_params.py` setting:
- `key1`: "oligo1", "oligo2", "negative" or "seriesN" -> decide which dataset (others can be added)
- `key2`: "R1", "R2", "R1R2" -> different reading directions (should give similar but not identical results, due to experimental imperfections)
- `key3`: "fw", "rev" or "all" -> select only forward/reverse/all sequences
- `key\_filter`: True/False -> whether to apply (`True`) or not (`False`) a special criterion to filter data. Default criterion here is to exclude PCR by-products.
- `key\_no\_cut`: True/False -> whether to print either full sequences (`True`) or cut sequences, deleting the primer bases (`False`)
- `use_stop` -> decide whether to really do it (`True/False`) --> for bootstrapping this MUST be set True
- `subset_steps`-> size of each subsample (this MUST also be set)

Optionally, other parameters can be modified:
- colors of RSA histograms
- True/False for creating (or not) abundance histograms for unique strands
- minimum threshold quality of the reads
- `l` -> resource length (defaults to 20 bases)
- `n` -> number of top-n strands for the related analysis of dominant individuals
- `random_seq=50` -> number of random nucleotides, by default; not used, currently, apart from `L` ordinary definition
- `cap_size=25` -> size of fixed sequences at the two ends; not really used
- `extra_end=1` -> sometimes there is an extra base, old code versions needed it, currently it is ignored
- `L=random_seq+cap_size+extra_end` -> max length, with cap and last one - length of predators -> can be edited (e.g. for seriesN)
- `lower_bound=L-6` -> discard strand with less than `lower_bound` bases -> can be edited

Another editable parameter is `results_folder`, which can be changed in case one needs to save some results separately.

Note that `input_params.py.old` and `input_params.py.new` contain two examples of specific settings to analyze cycle 9 of experiments (old) and cycles 9bis and 10 (new), as required during the review process of the paper. 

------------------------------------- compilation -------------------------------------

2) `make all` to generate C++ executables (C++-17 is used, but C++-11 compatibility should be enough)

------------------------------------- read\_fastq.py -------------------------------------

3) execute `bash launch_many_bootstrap.sh` which starts `n_replicas` instances of the command `python3 read_fastq.py`, each in a separated subfolder. This, in turn, processes the FASTQ files and generates text files and plots with the outcomes of the performed analyses.
`read_fastq.py` calls itself:
- `find_MCO_serial.x` (executable of the corresponing `C++` code for Maximum Consecutive Overlap calculation between strands - see https://www.mdpi.com/1099-4300/24/4/458 for its definition and related discussions). 
- `find_equal_pair.x` detects the number of consecutive identical bases between two strands passed by command line. Based on the same routines of `find_MCO_serial`, simplified version, used to detect aliens
- `module_functions.py`: process FASTQ files, filter sequences, sort them by abundance, reverse and complement strands, track the abundance of the top-`n` most abundant ones across cycles and compute their cross-MCO matrix
- `main_plot.py`: generate text files and plots for RSA histograms, Shannon entropy associated to them, evolution of top-`n` strands, the fraction of total population covered by top-`n` individuals and the 2D histogram of (MCO,MCO_2nd). It calls `module_plots.py`.
