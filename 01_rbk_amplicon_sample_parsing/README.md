# RBK sample ID parsing and summary statistics

## Aims
* Create scripts to parse sample barcodes, IDIS, gene and expected amplicon length from text file
* Link sample information from text file with fastq files
* Write scripts in Python3
* Create unit tests for each function

**USAGE**
```bash
## requires python3 environment
./sample_parsing.py -f testdata/sample-id.txt

## unittest:
pytest -xv
```

## Update 2023-01-09 -- DONE:
* Summary stat table for all BCs added
  
## TO DO (2023-01-09):
* Parallelization
* Check if codes runs for other users

## Update 2023-01-05 -- DONE:
* Create directory for each BC listed in `sample-id.txt` file
* Put concatinated raw fastq files in BC directories
* Make read lengths histogram and summary statistics table in BC directories
* Implement unit testing

## TO DO (2023-01-05):
* Parallelization
* Add a summary stat table for all BCs
  
