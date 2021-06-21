# SPOT-RNA-2D
SPOT-RNA-2D: Predicting RNA distance contact-map by deep ensembled learning on evolutionary coupling

1. `git clone https://github.com/jaswindersingh2/SPOT-RNA-2D.git && cd SPOT-RNA-2D`
2. `virtualenv -p python3.6 venv && source ./venv/bin/activate`
3. `pip install -r requirements.txt`
4. `python ss_and_non_ss_accuracy.py`


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


5. `python ss_accuracy_only.py`


```
 		          TS1   TS2   TS3
RNAfold          0.575 0.687 0.838
LinearPartition  0.610 0.674 0.826
```


After removing SPOT-RNA training data overlap with test sets (TS1, TS2)


6. `python ss_accuracy_only.py`

```
		 		  TS1   TS2   TS3
SPOT-RNA         0.684 0.744 0.823
SPOT-RNA2        0.703 0.736 0.828
RNAfold          0.575 0.687 0.838
LinearPartition  0.610 0.674 0.826
```

