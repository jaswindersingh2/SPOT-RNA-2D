#!/usr/bin/env python

import tensorflow as tf
import numpy as np
import pandas as pd
import os
import argparse, tqdm
from utils.utils import get_data, prob_to_secondary_structure
#from utils.FastaMLtoSL import FastaMLtoSL
from pathlib import Path
import time
import subprocess

start = time.time()
from argparse import RawTextHelpFormatter

base_path = os.path.dirname(os.path.realpath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument('--rna_id', default='', type=str, help='name of RNA sequence file; default = ''6p2h_A''\n', metavar='')
parser.add_argument('--list_rna_ids', default='datasets/TS1_ids', type=str, help='file consists of list name of RNA sequence files for batch prediction; default = ''datasets/TS1_ids''\n', metavar='')
parser.add_argument('--input_feats', default='input_features/', type=str, help='Path to directory consists of input features files with the same name as specified with --rna_id flag; default = ''inputs_features/''\n', metavar='')
parser.add_argument('--single_seq', default=0, type=int, help='set equal to 1 for SPOT-RNA-2D-Single; default = ''0 for SPOT-RNA-2D''\n', metavar='')
parser.add_argument('--outputs',default='outputs/', type=str, help='Path to directory to save output files; default = ''outputs/\n', metavar='')
parser.add_argument('--gpu', default=-1, type=int, help='To run on GPU, specifiy GPU number. If only one GPU in computer specifiy 0; default = -1 (no GPU)\n', metavar='')
parser.add_argument('--cpu',default=16, type=int, help='Specify number of cpu threads that program can use; default = 16\n', metavar='')
args = parser.parse_args()


if args.single_seq == 1:
    feat_mean = [0.24416392, 0.19836862, 0.30642843, 0.24948534, 0.24416392, 0.19836862, 0.30642843, 0.24948534, 0.004500812852233642]
    feat_std  = [0.43031312, 0.39954973, 0.46168336, 0.43343267, 0.43031312, 0.39954973, 0.46168336, 0.43343267, 0.058931698637260006]
else:
    feat_mean = [0.24416392, 0.19836862, 0.30642843, 0.24948534, 2.58825501, 2.74345347, 2.53821291, 2.84523057, 0.24416392, 0.19836862, 0.30642843, 0.24948534, 2.58825501, 2.74345347, 2.53821291, 2.84523057, 0.004500812852233642, 6.693295123250723e-05]
    feat_std  = [0.43031312, 0.39954973, 0.46168336, 0.43343267, 1.99199797, 1.99193953, 2.17897853, 2.18739682, 0.43031312, 0.39954973, 0.46168336, 0.43343267, 1.99199797, 1.99193953, 2.17897853, 2.18739682, 0.058931698637260006, 0.1624832931753003]

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

base_path = os.path.dirname(os.path.realpath(__file__))

if args.rna_id:
    list_rna_ids = [args.rna_id]
else:
    with open(args.list_rna_ids) as file:
        list_rna_ids = [line.strip() for line in file.read().splitlines() if line.strip()]

print('\nChecking for inputs features files in path ' + args.input_feats + '\n')

for rna_id in list_rna_ids:

	print('\n', rna_id+':')	

	if Path(args.input_feats + '/' + rna_id).is_file():
		print('RNA sequence  \u2713')
		with open(args.input_feats + '/' + rna_id) as f:
			temp_1 = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0]).values
		seq_ref = ''.join([j.upper() for j in temp_1[0, 0]])
	else: raise ValueError('RNA sequence file does not exists in path ' + args.input_feats + '/' + rna_id + '\n')

#	if args.single_seq == 1:
	if Path(args.input_feats + '/' + rna_id + '_dp.ps').is_file(): print('base pair probability \u2713')
	else:
		print('base pair probability \u2717')
		print('base pair probability file does not exists in path ' + args.input_feats + '/' + rna_id + '_dp.ps')
		print('Generating RNAfold base-pair probability')
		process1 = subprocess.Popen(["RNAfold", "-p", "-i", args.input_feats + '/' + rna_id], stdout=subprocess.PIPE); time.sleep(0.5)
		process2 = subprocess.Popen(["mv", rna_id + "_dp.ps", args.input_feats + '/'], stdout=subprocess.PIPE)
		process2 = subprocess.Popen(["mv", rna_id + "_ss.ps", args.input_feats + '/'], stdout=subprocess.PIPE)
		if Path(args.input_feats + '/' + rna_id + '_dp.ps').is_file(): print('base pair probability \u2713')

	if args.single_seq != 1:
		if Path(args.input_feats + '/' + rna_id + '.pssm').is_file(): print('PSSM features \u2713')
		else: raise ValueError('PSSM file does not exists in path ' + args.input_feats + '/' + rna_id + '.pssm\n')

		if Path(args.input_feats + '/' + rna_id + '.dca').is_file(): print('DCA features \u2713')
		else: raise ValueError('DCA file does not exists in path '  + args.input_feats + '/' + rna_id + '.dca\n')
	time.sleep(0.1)
print()


os.environ["CUDA_VISIBLE_DEVICES"]= str(args.gpu)
os.environ['KMP_WARNINGS'] = 'off'

if args.single_seq==1:
	NUM_MODELS = [4]
else:
	NUM_MODELS = [0, 1, 2, 3]

outputs = {}
mask = {}
sequences = {}
def sigmoid(x):
    return 1/(1+np.exp(-np.array(x, dtype=np.float128)))

for MODEL in NUM_MODELS:

    if args.gpu==-1:
            config = tf.ConfigProto(intra_op_parallelism_threads=args.cpu, inter_op_parallelism_threads=args.cpu)
    else:
	    config = tf.compat.v1.ConfigProto()
	    config.allow_soft_placement=True
	    config.log_device_placement=False
        
    if args.single_seq == 1: print('\nPredicting for SPOT-RNA-2D-Single model')
    else: print('\nPredicting for SPOT-RNA-2D model '+str(MODEL))

    with tf.compat.v1.Session(config=config) as sess:
        saver = tf.compat.v1.train.import_meta_graph(os.path.join(base_path, 'checkpoints', 'model_' + str(MODEL) + '.meta'))
        saver.restore(sess, os.path.join(base_path, 'checkpoints', 'model_' + str(MODEL)))
        graph = tf.compat.v1.get_default_graph()
#        for op in graph.get_operations():
#            print(op.name)

        tmp_out = graph.get_tensor_by_name('output_FC/fully_connected/BiasAdd:0')

        for rna_id in tqdm.tqdm(list_rna_ids):
            #print(rna_id)
            seq, seq_len,feature,zero_mask,label_mask = get_data(args, rna_id)
            out = sess.run([tmp_out], feed_dict={'input_feature:0':(feature-feat_mean)/feat_std, 'zero_mask:0':zero_mask, 'label_mask:0':label_mask, 'seq_len:0':seq_len, 'dropout:0':1})

            if MODEL == 0 or MODEL == 4:
                outputs[rna_id] = [sigmoid(out[0])]
                mask[rna_id] = label_mask
                sequences[rna_id] = seq
            else:
                outputs[rna_id].append(sigmoid(out[0]))

    tf.compat.v1.reset_default_graph()

print('\nSaving Output to directory ' + args.outputs)
for rna_ids in tqdm.tqdm(list_rna_ids):

    ensemble_outputs = np.mean(outputs[rna_ids],0)
    prob_to_secondary_structure(ensemble_outputs, mask[rna_ids], sequences[rna_ids], rna_ids, args, base_path)

print('\nFinished!')
end = time.time()
print('\nProcesssing Time {} seconds'.format(end - start))
