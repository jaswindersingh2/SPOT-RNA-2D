# SPOT-RNA-2D
Predicting RNA distance-based contact maps by integrated deep learning on physics-inferred base-pairing and evolutionary-derived coupling

## Installation of SPOT-RNA-2D dependencies

1. `git clone https://github.com/jaswindersingh2/SPOT-RNA-2D.git && cd SPOT-RNA-2D`
2. `virtualenv -p python3.6 venv && source ./venv/bin/activate`
3. `pip install -r requirements.txt`

## Usage

### Run SPOT-RNA-2D for single RNA:

4. `./run.py --rna_id 6ol3_C --input_feats input_features/ --outputs outputs/`

### Run SPOT-RNA-2D for multiple RNAs:

5. `./run.py --list_rna_ids datasets/TS1_ids --input_feats input_features/ --outputs outputs/`

### Run SPOT-RNA-2D-Single:

6. `./run.py --list_rna_ids datasets/TS3_ids --input_feats input_features/ --single_seq 1 --outputs outputs/`


For more options:

6. `./run.py --help`

```
usage: run.py [-h] [--rna_id] [--list_rna_ids] [--input_feats] [--single_seq]
              [--outputs] [--gpu] [--cpu]

optional arguments:
  -h, --help       show this help message and exit
  --rna_id         name of RNA sequence file; default = 6p2h_A
  --list_rna_ids   file consists of list name of RNA sequence files for batch
                   prediction; default = datasets/TS1_ids
  --input_feats    Path to directory consists of input features files with the
                   same name as specified with --rna_id flag; default =
                   inputs_features/
  --single_seq     set equal to 1 for SPOT-RNA-2D-Single; default = 0 for
                   SPOT-RNA-2D
  --outputs        Path to directory to save output files; default = outputs/
  --gpu            To run on GPU, specifiy GPU number. If only one GPU in
                   computer specifiy 0; default = -1 (no GPU)
  --cpu            Specify number of cpu threads that program can use; default
                   = 16
```

