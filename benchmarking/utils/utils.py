import pandas as pd
import numpy as np
import os

def get_pairs(labels):
    pairs = []
    stack = []
    stack_2 = []
    stack_3 = []
    stack_4 = []
    stack_5 = []
    stack_6 = []
    stack_7 = []
    for i, I in enumerate(labels):
        if I == '(':
            stack.append(i)
            #print(I, i)
        elif I == ')':
            pairs.append(sorted([stack[-1], i]))
            #print(I, i)
            del stack[-1]
        elif I == '<':
            stack_2.append(i)
        elif I == '>':
            pairs.append(sorted([stack_2[-1], i]))
            del stack_2[-1]
        elif I == '{':
            stack_3.append(i)
        elif I == '}':
            pairs.append(sorted([stack_3[-1], i]))
            del stack_3[-1]
        elif I == '[':
            stack_4.append(i)
        elif I == ']':
            pairs.append(sorted([stack_4[-1], i]))
            del stack_4[-1]
        elif I == 'A':
            stack_5.append(i)
        elif I == 'a':
            pairs.append(sorted([stack_5[-1], i]))
            del stack_5[-1]
        elif I == 'B':
            stack_6.append(i)
        elif I == 'b':
            pairs.append(sorted([stack_6[-1], i]))
            del stack_6[-1]
        elif I == 'C':
            stack_7.append(i)
        elif I == 'c':
            pairs.append(sorted([stack_7[-1], i]))
            del stack_7[-1]
        elif I == '.' or I == '-':
            continue
        else:
            print(I)
    return pairs



########--------------------- parse GREMLIN output ---------------------------- #########################
def GREMLIN_prob(base_path, id, seq):
		
############ read gremlin output ##############################
	try:
		with open(base_path + '/predictions/GREMLIN/' + id + '.dca') as f:
			temp4 = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0], usecols=[0,1,2]).values

		dca = np.zeros((len(seq), len(seq)))
		for k in temp4:
			if abs(int(k[0]) - int(k[1])) < 4:
			    dca[int(k[0]), int(k[1])] = 1*k[2]
			else:
			    dca[int(k[0]), int(k[1])] = k[2]
	except:
		#print("dca missing", id)
		dca = np.zeros((len(seq), len(seq)))

	return dca

########--------------------- parse plmc output output ---------------------------- #########################
def plmc_prob(base_path, id, seq):
		
	with open(base_path + '/predictions/PLMC/' + id + '.dca') as f:
		temp4 = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, usecols=[0,2,5]).values

	dca = np.zeros((len(seq), len(seq)))
	for k in temp4:
		if abs(int(k[0]) - int(k[1])) < 4:
		    dca[int(k[0]-1), int(k[1]-1)] = 1*k[2]
		else:
		    dca[int(k[0]-1), int(k[1]-1)] = k[2]

	return dca

########--------------------- parse mfdca output output ---------------------------- #########################
def mfdca_prob(base_path, id, seq):
		
	with open(base_path + '/predictions/mfDCA/' + 'MFDCA_apc_fn_scores_' + id + '.txt') as f:
		temp4 = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0,1,2,3,4,5,6,7,8,9,10], usecols=[0,1,2]).values

	dca = np.zeros((len(seq), len(seq)))
	for k in temp4:
		if abs(int(k[0]) - int(k[1])) < 4:
			dca[int(k[0]-1), int(k[1]-1)] = 1*k[2]
		else:
			dca[int(k[0]-1), int(k[1]-1)] = k[2]

	return dca

########--------------------- parse plmdca output output ---------------------------- #########################
def plmdca_prob(base_path, id, seq):
		
	with open(base_path + '/predictions/plmDCA/PLMDCA_apc_fn_scores_' + id + '.txt') as f:
		temp4 = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0,1,2,3,4,5,6,7,8,9,10,11], usecols=[0,1,2]).values
	dca = np.zeros((len(seq), len(seq)))
	for k in temp4:
		if abs(int(k[0]) - int(k[1])) < 4:
			dca[int(k[0]-1), int(k[1]-1)] = 1*k[2]
		else:
			dca[int(k[0]-1), int(k[1]-1)] = k[2]

	return dca

######## --------------------- parse base-pair probability RNAfold output ---------------------------- #########################
def RNAfold_bp_prob(base_path, id, seq):

    with open(base_path + '/predictions/RNAfold/' + str(id) + '.prob') as f:
        temp = pd.read_csv(f, comment='#', header=None).values
    output_pred = np.zeros((len(seq), len(seq)))
    for i in temp[:,0]:
        a = i.split(' ')
        if abs(int(a[0]) - int(a[1])) < 4:
            output_pred[int(a[0]) - 1, int(a[1]) - 1] = 1*float(a[2])**2
        else:
            output_pred[int(a[0]) - 1, int(a[1]) - 1] = float(a[2])**2

    return output_pred

######## --------------------- parse base-pair probability SPOT-RNA output ---------------------------- #########################
def spotrna(id, seq):

    output_pred = np.loadtxt(base_path + '/predictions/SPOT-RNA/' + str(id) + '.prob')

    return output_pred

######## --------------------- parse base-pair probability SPOT-RNA2 output ---------------------------- #########################
def spotrna2(base_path, id, seq):

    output_pred = np.loadtxt(base_path + '/predictions/SPOT-RNA2/' + str(id) + '.prob')

    return output_pred

######## --------------------- parse LinearPartition output ---------------------------- #########################
def LinearPartition(base_path, id, seq):

	with open(base_path + '/predictions/LinearPartition/' + id + '.prob', 'r') as f:
		prob = pd.read_csv(f, delimiter=None, delim_whitespace=True, header=None).values
	y_pred =  np.zeros((len(seq), len(seq)))
	for i in prob:
		y_pred[int(i[0])-1, int(i[1])-1] = i[2]

	return y_pred


######## --------------------- parse base-pair probability SPOT-RNA-2D-Single output ---------------------------- #########################
def spotrna_2d_single(base_path, id, seq):

    output_pred = np.loadtxt(base_path + '/predictions/SPOT-RNA-2D-Single/' + str(id) + '.prob')

    return output_pred


######## --------------------- parse base-pair probability SPOT-RNA-2D output ---------------------------- #########################
def spotrna_2d(base_path, id, seq):
    #print(base_path)
    output_pred = np.loadtxt(base_path + '/predictions/SPOT-RNA-2D/' + str(id) + '.prob')

    return output_pred


######## --------------------- parse RNAfold output ---------------------------- #########################
def RNAfold_bps(base_path, id, seq):

    with open(base_path + '/predictions/RNAfold/' + str(id) + '.dbn') as f:
        temp = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, usecols=[0], skiprows=[0,3,4,5]).values
    seq_pred = [i for i in temp[0,0]]
    labels = [I for i,I in enumerate(temp[1, 0])]
    assert len(seq) == len(seq_pred) == len(labels)
    pred_pairs = get_pairs(labels)

    return pred_pairs

######## --------------------- parse LinearPartition output ---------------------------- #########################
def LinearPartition_bps(base_path, id, seq):

	with open(base_path + '/predictions/LinearPartition/' + id + '.prob', 'r') as f:
		prob = pd.read_csv(f, delimiter=None, delim_whitespace=True, header=None).values
	y_pred =  np.zeros((len(seq), len(seq)))
	for i in prob:
		y_pred[int(i[0])-1, int(i[1])-1] = i[2]

	tri_inds = np.triu_indices(y_pred.shape[0], k=1)

	out_pred = y_pred[tri_inds]
	outputs = out_pred[:, None]
	seq_pairs = [[tri_inds[0][j], tri_inds[1][j], ''.join([seq[tri_inds[0][j]], seq[tri_inds[1][j]]])] for j in
		         range(tri_inds[0].shape[0])]

	outputs_T = np.greater_equal(outputs, 0.198)
	pred_pairs = [i for I, i in enumerate(seq_pairs) if outputs_T[I]]
	pred_pairs = [i[:2] for i in pred_pairs]

	return pred_pairs


############ load base-pair prob form SPOT-RNA ##############################
def spot_rna_bps(base_path, id, seq):

    y_pred = np.loadtxt(base_path + '/predictions/SPOT-RNA/' + str(id) + '.prob')
    tri_inds = np.triu_indices(y_pred.shape[0], k=1)

    out_pred = y_pred[tri_inds]
    outputs = out_pred[:, None]
    seq_pairs = [[tri_inds[0][j], tri_inds[1][j], ''.join([seq[tri_inds[0][j]], seq[tri_inds[1][j]]])] for j in
                 range(tri_inds[0].shape[0])]

    outputs_T = np.greater_equal(outputs, 0.335)
    pred_pairs = [i for I, i in enumerate(seq_pairs) if outputs_T[I]]
    pred_pairs = [i[:2] for i in pred_pairs]

    return pred_pairs


############ load base-pair prob form SPOT-RNA ##############################
def spot_rna2_bps(base_path, id, seq):

    y_pred = np.loadtxt(base_path + '/predictions/SPOT-RNA2/' + str(id) + '.prob')
    tri_inds = np.triu_indices(y_pred.shape[0], k=1)

    out_pred = y_pred[tri_inds]
    outputs = out_pred[:, None]
    seq_pairs = [[tri_inds[0][j], tri_inds[1][j], ''.join([seq[tri_inds[0][j]], seq[tri_inds[1][j]]])] for j in
                 range(tri_inds[0].shape[0])]

    outputs_T = np.greater_equal(outputs, 0.795)
    pred_pairs = [i for I, i in enumerate(seq_pairs) if outputs_T[I]]
    pred_pairs = [i[:2] for i in pred_pairs]

    return pred_pairs

