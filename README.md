# CMS-BWT
Tool for computing the BWT of a set of highly similar strings using compressed matching statistics.

## Installation

```sh
git clone https://github.com/fmasillo/CMS-BWT.git
cd CMS-BWT
git submodule update --init --recursive
mkdir build
cd build
cmake ..
make
```

## Usage

Usage: 
```sh
./cms_bwt [options] <input filename>
<input filename> is the name of the file containing paths to the reference sequence (in the first line) and to the collection file (in the second line).
  Options: 
        -p      read only a prefix of the file expressed in number of characters, def. whole file
        -b      size for the additional memory buffer in GB, def. 2 
        -r      outputs the run-length encoded BWT, def. false 
        -o      basename for the output files, def. <input filename>
```

Command example for running the implementation and outputting to my_output.rl_bwt the run-length encoded BWT of the first 100000000 characters using 1GB of extra space:
```sh
./cms_bwt -p 100000000 -b 1 -r -o my_output file_example.txt
```

file_example.txt content should look like this:
```
/data/reference.fa
/data/collection.fa

```

## ESA 2023 version of the code

If you are looking for the ESA 2023 version of the code, please check the branch `ESA_version`.

## Citation

Conference paper:

Francesco Masillo. Matching Statistics speed up BWT construction. In Proc. of the 31st Annual European Symposium on Algorithms (ESA 2023), volume 274 of LIPIcs, pages 83:1-83:15. Schloss Dagstuhl - Leibniz-Zentrum für Informatik, 2023: https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ESA.2023.83

```
@inproceedings{Masillo23,
  author       = {Francesco Masillo},
  editor       = {Inge Li G{\o}rtz and
                  Martin Farach{-}Colton and
                  Simon J. Puglisi and
                  Grzegorz Herman},
  title        = {Matching Statistics Speed up {BWT} Construction},
  booktitle    = {Proc. of the 31st Annual European Symposium on Algorithms, {ESA} 2023, September
                  4-6, 2023, Amsterdam, The Netherlands},
  series       = {LIPIcs},
  volume       = {274},
  pages        = {83:1--83:15},
  publisher    = {Schloss Dagstuhl - Leibniz-Zentrum f{\"{u}}r Informatik},
  year         = {2023},
  doi          = {10.4230/LIPICS.ESA.2023.83},
}
```


