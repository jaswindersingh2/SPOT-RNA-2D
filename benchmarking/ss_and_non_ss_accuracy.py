import pandas as pd
import numpy as np
import os
from utils.utils import *
from pathlib import Path

base_path = os.path.dirname(os.path.realpath(__file__))
base_path = str(Path(base_path).parent.absolute())

with open(base_path + '/datasets/TS1_ids') as f:
    TS1_ids = f.read().splitlines()
with open(base_path + '/datasets/TS2_ids') as f:
    TS2_ids = f.read().splitlines()
with open(base_path + '/datasets/TS3_ids') as f:
    TS3_ids = f.read().splitlines()

n = 1  # top L/n contacts

#predictors = ['GREMLIN', 'plmDCA', 'mfDCA', 'PLMC', 'RNAfold', 'LinearPartition', 'SPOT-RNA-2D-Single', 'SPOT-RNA-2D']
predictors = ['RNAfold', 'LinearPartition', 'SPOT-RNA-2D-Single', 'SPOT-RNA-2D']
sets = [TS1_ids, TS2_ids, TS3_ids]


print('\n \t \t\t TS1 \t      TS2 \t    TS3\n' )
print('\n \t \t     ss   non-ss   ss   non-ss   ss   non-ss' )

for number, pred in enumerate(predictors[0:]):

	for num, ids in enumerate(sets[0:]):

		# initialize empty count and lists to count no. of RNAs and save f1 scores
		count = 0;		save_base_pairs_f1 = [];	save_non_base_pairs_f1 = []

		for K,k in enumerate(ids[0:]):

			# read RNA sequence
			with open(base_path + '/datasets/all_sequences/' + str(k)) as f:
				temp_1 = pd.read_csv(f, comment='#', delim_whitespace=True, header=None, skiprows=[0]).values
			seq = [j.upper() for j in temp_1[0, 0]]

			# read list of native secondary structure base pairs (canonical, non-canonical, pseudoknots, multiplets, lone-pairs)
			native_ss_base_pairs = np.loadtxt(base_path + '/datasets/ss_base_pairs_labels/' + k + '.bps')
			native_ss_base_pairs = [ [int(i[0]), int(i[1])] for i in native_ss_base_pairs]  # convert float base-pairs [1.0, 67.0] to int base-pair [1, 67]
			native_ss_base_pairs = [i for i in native_ss_base_pairs if abs(i[0]-i[1]) >= 24]  # extract long-range base-pairs only

			# read native distance-based native contacts LxL matrix
			native_dist_contacts = np.loadtxt(base_path + '/datasets/distance_contact_labels/' + k + '.true')
			native_dist_contacts[native_dist_contacts==-1] = 0
			tri_inds = np.where(native_dist_contacts==1)
			native_dist_contacts = [[i,j] for i,j in zip(tri_inds[0], tri_inds[1]) if abs(i-j) >= 24] # extract long-range base-pairs only

			native_non_ss_pairs = [i for i in native_dist_contacts if i not in native_ss_base_pairs]

			# read predicted contacts LxL matrix from different predictors
			if pred == 'RNAfold':
				pred_contacts_matrix = RNAfold_bp_prob(base_path, k, seq)
			elif pred == 'LinearPartition':
				pred_contacts_matrix = LinearPartition(base_path, k, seq)
			elif pred == 'CentroidFold':
				pred_contacts_matrix = CentroidFold_mfe_bp_prob(base_path, k, seq)
			elif pred == 'SPOT-RNA':
				pred_contacts_matrix = spotrna(base_path, k, seq)
			elif pred == 'SPOT-RNA2':
				pred_contacts_matrix = spotrna2(base_path, k, seq)
			elif pred == 'GREMLIN':
				pred_contacts_matrix = GREMLIN_prob(base_path, k, seq)
			elif pred == 'PLMC':
				pred_contacts_matrix = plmc_prob(base_path, k, seq)
			elif pred == 'mfDCA':
				pred_contacts_matrix = mfdca_prob(base_path, k, seq)
			elif pred == 'plmDCA':
				pred_contacts_matrix = plmdca_prob(base_path, k, seq)
			elif pred == 'SPOT-RNA-2D-Single':
				pred_contacts_matrix = spotrna_2d_single(base_path, k, seq)
			elif pred == 'SPOT-RNA-2D':
				pred_contacts_matrix = spotrna_2d(base_path, k, seq)
				
			# extract upper triangular top L long-range contacts predicted proability LxL matrix
#			tri_inds = np.triu_indices(pred_contacts_matrix.shape[0], k=1)
#			all_pred_contacts = np.array([[i,j,pred_contacts_matrix[i,j]] for i,j in zip(tri_inds[0], tri_inds[1]) if abs(i-j) >= 24])

			tri_inds = np.triu_indices(pred_contacts_matrix.shape[0], k=24)
			all_pred_contacts = np.array([[i,j,pred_contacts_matrix[i,j]] for i,j in zip(tri_inds[0], tri_inds[1]) if [i, j] not in native_non_ss_pairs])

			all_pred_contacts_sorted = all_pred_contacts[all_pred_contacts[:,2].argsort()[::-1]]
			pred_contacts = [[int(i[0]), int(i[1])] for i in all_pred_contacts_sorted[0:int(len(seq)/n)]]


########################## seconadary structure base pair perfornace measure ##################
			true_pairs = native_ss_base_pairs

			pred_correctly = [i for i in pred_contacts if i in true_pairs]; #print(len(pred_correctly))
			pred_wrongly = [i for i in pred_contacts if i not in true_pairs]; #print(len(pred_wrongly))
#			pred_wrongly = [i for i in pred_contacts if i not in native_dist_contacts]; #print(len(pred_wrongly))

			tp = len(pred_correctly)
			fp = len(pred_wrongly)
			fn = len([i for i in true_pairs if i not in pred_contacts])


			try:
				pre = tp / (tp + fp)
				sen = tp / (tp + fn)
				f1 = 2*((pre*sen)/(pre + sen))
			except:
				f1 = 0

			save_base_pairs_f1.append(f1)


########################## non-base pair perfornace measure ##################

			true_pairs = native_non_ss_pairs

			# read predicted contacts LxL matrix from different predictors
			if pred == 'RNAfold':
				pred_contacts_matrix = RNAfold_bp_prob(base_path, k, seq)
			elif pred == 'LinearPartition':
				pred_contacts_matrix = LinearPartition(base_path, k, seq)
			elif pred == 'CentroidFold':
				pred_contacts_matrix = CentroidFold_mfe_bp_prob(base_path, k, seq)
			elif pred == 'SPOT-RNA':
				pred_contacts_matrix = spotrna(base_path, k, seq)
			elif pred == 'SPOT-RNA2':
				pred_contacts_matrix = spotrna2(base_path, k, seq)
			elif pred == 'GREMLIN':
				pred_contacts_matrix = GREMLIN_prob(base_path, k, seq)
			elif pred == 'PLMC':
				pred_contacts_matrix = plmc_prob(base_path, k, seq)
			elif pred == 'mfDCA':
				pred_contacts_matrix = mfdca_prob(base_path, k, seq)
			elif pred == 'plmDCA':
				pred_contacts_matrix = plmdca_prob(base_path, k, seq)
			elif pred == 'SPOT-RNA-2D-Single':
				pred_contacts_matrix = spotrna_2d_single(base_path, k, seq)
			elif pred == 'SPOT-RNA-2D':
				pred_contacts_matrix = spotrna_2d(base_path, k, seq)


			tri_inds = np.triu_indices(pred_contacts_matrix.shape[0], k=24)
			all_pred_contacts = np.array([[i,j,pred_contacts_matrix[i,j]] for i,j in zip(tri_inds[0], tri_inds[1]) if [i, j] not in native_ss_base_pairs])

			all_pred_contacts_sorted = all_pred_contacts[all_pred_contacts[:,2].argsort()[::-1]]
			pred_contacts = [[int(i[0]), int(i[1])] for i in all_pred_contacts_sorted[0:int(len(seq)/n)]]

			pred_correctly = [i for i in pred_contacts if i in true_pairs]; #print(len(pred_correctly))
			pred_wrongly = [i for i in pred_contacts if i not in true_pairs]; #print(len(pred_wrongly))
#			pred_wrongly = [i for i in pred_wrongly if i not in native_ss_base_pairs]; #print(len(pred_wrongly))

			tp = len(pred_correctly)
			fp = len(pred_wrongly)
			fn = len([i for i in true_pairs if i not in pred_contacts])

			try:
				pre = tp / (tp + fp)
				sen = tp / (tp + fn)
				f1 = 2*((pre*sen)/(pre + sen))
			except:
				f1 = 0

			save_non_base_pairs_f1.append(f1)

			count += 1

		if num==0:
			print(pred + " "*(len('SPOT-RNA-2D-Single')-len(pred)),end=' ')
		print(' {0:.3f}'.format(np.mean(save_base_pairs_f1)), end=' ')
		print(' {0:.3f}'.format(np.mean(save_non_base_pairs_f1)), end=' ')
		#print(' {}'.format(count), end='\t')
	print()

print()


