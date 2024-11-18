# Copyright Gleb Bezgin 2024
Unzip the package; then execute line-by-line the following code in MATLAB:

demo_createMcaProfiles.txt

Note that line 1 should be modified by specifying your directory, and line ~8 should be modified in such way that it adds SurfStat to the path (the only external dependency here), such as:

addpath(genpath(/path/to/surfstat))

The code will perform the main MCA steps and output 3 figures: the original map, the subparcellation, and the profiles. Data used are grey matter density from a sample of subjects from openly available NKI-RS dataset; after this paper's acceptance, an example with tau data across Braak stages will be provided, as described in the paper.

This code has been tested on MATLAB R2018a.
