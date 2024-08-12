# cvmmlst


```
                                  __     __
  ______   ______ ___  ____ ___  / /____/ /_
 / ___/ | / / __ `__ \/ __ `__ \/ / ___/ __/
/ /__ | |/ / / / / / / / / / / / (__  ) /_
\___/ |___/_/ /_/ /_/_/ /_/ /_/_/____/\__/


```

cvmmlst is a bacteria mlst analysis tool that could run on Windows, Linux and MAC os. Some of the code ideas in cvmmlst draw on Torsten Seemanns excellent [mlst](https://github.com/tseemann/mlst) tool.






## Installation
### Using pip
pip3 install cvmmlst

### Using conda
comming soon...

## Dependency
- BLAST+ >2.7.0

**you should add BLAST in your PATH**


## Blast installation
### Windows


Following this tutorial:
[Add blast into your windows PATH](http://82.157.185.121:22300/shares/BevQrP0j8EXn76p7CwfheA)

### Linux/Mac
The easyest way to install blast is:

```
conda install -c bioconda blast
```



## Introduction

### Initialize reference database

After finish installation, you should first initialize the reference database using following command
```
cvmmlst init
```


### Usage
```
usage: cvmmlst -i <genome assemble directory> -o <output_directory>

Author: Qingpo Cui(SZQ Lab, China Agricultural University)

options:
  -h, --help            show this help message and exit
  -i I                  <input_path>: the PATH to the directory of assembled genome files. Could not use with -f
  -f F                  <input_file>: the PATH of assembled genome file. Could not use with -i
  -o O                  <output_directory>: output PATH
  -scheme SCHEME        <mlst scheme want to use>, cvmmlst show_schemes command could output all available schems
  -minid MINID          <minimum threshold of identity>, default=90
  -mincov MINCOV        <minimum threshold of coverage>, default=60
  -t T                  <number of threads>: default=8
  -v, --version         Display version

cvmmlst subcommand:
  {init,show_schemes,add_scheme}
    init                <initialize the reference database>
    show_schemes        <show the list of all available schemes>
    add_scheme          <add custome scheme, use cvmmlst add_scheme -h for help>
```

### Show available schemes

```
cvmmlst show_schemes
```



### Add custome scheme
```
usage: cvmmlst -i <genome assemble directory> -o <output_directory>

Author: Qingpo Cui(SZQ Lab, China Agricultural University) add_scheme
       [-h] [-name NAME] [-path PATH]

optional arguments:
  -h, --help  show this help message and exit
  -name NAME  <the custome scheme name>
  -path PATH  <the path of scheme>
```

-name: str -> the scheme name you want to use with -scheme options

-path: str -> the path of the directory that contains the fasta files of locus in schemes and the profile file

#### Example
```
cvmmlst add_scheme -name my_scheme -path path_to_my_scheme
```



### Output

you will get a text file and a summray file in csv format in the output directory.

The text file like
|dat | bglA | cat |ldh |abcZ | dapE | lhkA | ST | Scheme | FILE|
|---|---|---|---|---|---|---|---|---|---|
|3 |1 |4| 39 | 12 | 14 | 4 |87 | listeria_2 | 665|

The content in csv summary file like
|dat | bglA | cat |ldh |abcZ | dapE | lhkA | ST | Scheme | FILE|
|---|---|---|---|---|---|---|---|---|---|
|3 |1 |4| 39 | 12 | 14 | 4 |87 | listeria_2 | sample01|
|2 |4 |4 |1 |4 |3 |5 |3 |listeria_2 | sample02|
|6 |6| 8 |37 | 7 |8 |1 |121| listeria_2 | sample03|
|3 |1 |4| 39 | 12 | 14 | 4 |87 | listeria_2 | sample04|
|2 |4 |4 |1 |4 |3 |5 |3 |listeria_2 | sample05|
|6 |6| 8 |37 | 7 |8 |1 |121| listeria_2 | sample06|




## Update logs
|Date|Content|
|---|---|
|2024-08-12|Add three subcommand (init, show_schems, add_scheme)|
