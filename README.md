This is a repo to accompany the manuscript
Skelly et al. (2020) Cell Stem Cell 

To run the example code, first clone this repository and
`cd` to the `figshare_data/` directory.
Download the figshare
data from https://doi.org/10.6084/m9.figshare.12233570 
and unzip in the `figshare_data/` directory.

Next `cd` to the `../code/` directory, pull
the singularity image:
```bash
singularity pull library://daskelly/remote-builds/tidyqtl2-r4.0.0:1.0.0
```
(You may need to log in to singularity using a token
with `singularity remote login` first).

Finally, execute the code
```bash
bash run_code_singularity.sh
```
