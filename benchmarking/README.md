# SPOT-RNA-2D
Predicting RNA distance-based contact maps by integrated deep learning on physics-inferred base-pairing and evolutionary-derived coupling


1. `conda activate venv_spotrna_2d`
2. `python ss_and_non_ss_accuracy.py`

Old top L contacts 

```
 	 		 TS1 	      TS2 	    TS3


 	 	     ss   non-ss   ss   non-ss   ss   non-ss
GREMLIN             0.194  0.164  0.195  0.145  0.241  0.277 
plmDCA              0.211  0.173  0.228  0.177  0.268  0.299 
mfDCA               0.218  0.174  0.239  0.183  0.263  0.296 
PLMC                0.220  0.194  0.228  0.183  0.263  0.299 
RNAfold             0.243  0.173  0.322  0.235  0.342  0.333 
LinearPartition     0.246  0.167  0.321  0.226  0.348  0.350 
SPOT-RNA-2D-Single  0.369  0.455  0.503  0.516  0.509  0.656 
SPOT-RNA-2D         0.558  0.570  0.559  0.553  0.514  0.649
```

3. `python ss_and_non_ss_accuracy.py`

top L contacts 

```
 	 		 TS1 	      TS2 	    TS3


 	 	     ss   non-ss   ss   non-ss   ss   non-ss
RNAfold             0.201  0.177  0.247  0.236  0.268  0.334 
LinearPartition     0.205  0.173  0.254  0.233  0.270  0.349 
SPOT-RNA-2D-Single  0.218  0.454  0.269  0.529  0.284  0.646 
SPOT-RNA-2D         0.254  0.590  0.284  0.573  0.281  0.643
```

4. `python ss_and_non_ss_accuracy.py`

top L/2 contacts 

```
 	 		 TS1 	      TS2 	    TS3


 	 	     ss   non-ss   ss   non-ss   ss   non-ss
RNAfold             0.336  0.146  0.412  0.208  0.445  0.294 
LinearPartition     0.336  0.142  0.412  0.186  0.448  0.315 
SPOT-RNA-2D-Single  0.343  0.427  0.425  0.464  0.474  0.628 
SPOT-RNA-2D         0.421  0.476  0.448  0.481  0.468  0.610
```

5. `python ss_and_non_ss_accuracy.py`

top L/5 contacts 

```
 	 		 TS1 	      TS2 	    TS3


 	 	     ss   non-ss   ss   non-ss   ss   non-ss
RNAfold             0.544  0.086  0.651  0.130  0.756  0.197 
LinearPartition     0.551  0.080  0.642  0.096  0.759  0.196 
SPOT-RNA-2D-Single  0.519  0.249  0.609  0.253  0.696  0.347 
SPOT-RNA-2D         0.649  0.260  0.642  0.263  0.669  0.380
```

6. `python ss_and_non_ss_accuracy.py`

top L/10 contacts 

```
 	 		 TS1 	      TS2 	    TS3


 	 	     ss   non-ss   ss   non-ss   ss   non-ss
RNAfold             0.560  0.047  0.525  0.074  0.617  0.120 
LinearPartition     0.576  0.039  0.519  0.055  0.617  0.102 
SPOT-RNA-2D-Single  0.506  0.133  0.509  0.137  0.495  0.188 
SPOT-RNA-2D         0.575  0.139  0.494  0.147  0.567  0.210
```


