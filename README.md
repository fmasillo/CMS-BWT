# CMS-BWT
Tool for computing the BWT of a set of highly similar strings using compressed matching statistics.

## Installation

```sh
git clone https://github.com/fmasillo/CMS-BWT.git
git submodule update --init --recursive
cd CMS-BWT
make
```

## Usage

Usage: ./cms_bwt [options] <input filename>
<input filename> is the name of the file containing paths to the reference sequence (in the first line) and to the collection file (in the second line).
  Options: 
        -p      read only a prefix of the file expressed in number of characters, def. whole file
        -b      size for the additional memory buffer in GB, def. 2 
        -r      outputs the run-length encoded BWT, def. false 
        -m      memory saving implementation, def. false 
        -o      basename for the output files, def. <input filename>

Command example:
```sh
./cms_bwt -p 100000000 -b 1 -r -m -o my_output file_example.txt
```

file_example.txt content should look like this:
```
/data/reference.fa
/data/collection.fa
```


