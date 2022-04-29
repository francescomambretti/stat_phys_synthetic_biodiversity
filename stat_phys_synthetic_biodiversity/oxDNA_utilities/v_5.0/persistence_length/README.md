******************************

Tools to compute and plot the correlation function for tangent vectors along a DNA chain (also with caps), and compute the persistence length.

******************************

`persistence_length.py` is the main file. Requires 6 input arguments: `trajectoryFile steps blocks start_nucl end_nucl tot_nucl`

The `trajectoryFile` is split in blocks which are analyzed separately and the combined for global average (data blocking technique).

Only nucleotides between `start_nucl` and `end_nucl` are included in the computation; this allows, for example but not only, to remove caps & fixed sequences from the computation. The total number of nucleotides `tot_nucl` (of the system) is exploited to extract from the whole trajectory only the coordinates of the relevant nucleotides.

`blocks` sets the number of blocks for data blocking.
