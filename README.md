# cvmmlst


```
                                  __     __
  ______   ______ ___  ____ ___  / /____/ /_
 / ___/ | / / __ `__ \/ __ `__ \/ / ___/ __/
/ /__ | |/ / / / / / / / / / / / (__  ) /_
\___/ |___/_/ /_/ /_/_/ /_/ /_/_/____/\__/


```

cvmmlst is a bacteria mlst analysis tool that could run on Windows, Linux and MAC os. Some of the code ideas in cvmmlst draw on Torsten Seemanns excellent [mlst](https://github.com/tseemann/mlst) tool.






## 1. Installation
```
pip3 install cvmmlst
```




## 2. Dependency
- BLAST+ >2.7.0

**you should add BLAST in your PATH**


## 3. Blast installation
### 3.1 Windows


Following this tutorial:
[Add blast into your windows PATH](http://82.157.185.121:22300/shares/BevQrP0j8EXn76p7CwfheA)

### 3.2 Linux/Mac
The easyest way to install blast is:

```
conda install -c bioconda blast
```



## 4. Introduction

### 4.1 Initialize reference database

After finish installation, you should first initialize the reference database using following command
```
cvmmlst init
```


### 4.2 Usage
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

### 4.3 Show available schemes

```
cvmmlst show_schemes
```



### 4.4 Add custome scheme
```
usage: cvmmlst -i <genome assemble directory> -o <output_directory>

Author: Qingpo Cui(SZQ Lab, China Agricultural University) add_scheme
       [-h] [-name NAME] [-path PATH]

optional arguments:
  -h, --help  show this help message and exit
  -name NAME  <the custome scheme name>
  -path PATH  <the path to the files of custome scheme>
```

-name: str -> the scheme name you want to use with -scheme options

-path: str -> the path of the directory that contains the fasta files of locus in schemes and the profile file

#### Example
```
cvmmlst add_scheme -name my_scheme -path PATH_TO_my_scheme
```

The structure of scheme directory should looks like:
```
own_scheme
├── locus1.fasta
├── locus2.fasta
├── locus3.fasta
├── locus4.fasta
├── locus5.fasta
├── locus6.fasta
├── locus7.fasta
└── own_scheme.txt
```

The fasta file of corresponding locus is a multifasta file.

The multifasta file looks like:
```
>locus1_1
ATGATAGGTGAAGATATACAAAGAGTATTAG
>locus1_2
ATGATAGGTGAAGATATACAAAGAGTATTAG
>locus1_3
ATGATAGGTGAAGATATACAAAGAGTATTAG
>locus1_4
ATGATAGGCGAAGATATACAAAGAGTATTAG
>alocus1_5
ATGATAGGCGAAGATATACAAAGAGTATTAG
>locus1_6
ATGATAGGTGAAGATATACAAAGAGTATTAG
```

The own_scheme.txt is a tab-delimited text file.

The profile looks like:


|ST|locus1|locus2|locus3|locus4|locus5|locus6|locus7|clonal_complex|
|---|---|---|---|---|---|---|---|---|
|1|2|1|54|3|4|1|5|ST-21 complex|
|2|4|7|51|4|1|7|1|ST-45 complex|
|3|3|2|5|10|11|11|6|ST-49 complex|
|4|10|11|16|7|10|5|7|ST-403 complex|
|5|7|2|5|2|10|3|6|ST-353 complex|
|6|63|34|27|33|45|5|7||
|7|8|10|2|2|14|12|6|ST-354 complex|


### 4.5 Output

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




## 5. Update logs
|Date|Content|
|---|---|
|2024-08-12|Add three subcommand (init, show_schems, add_scheme)|
