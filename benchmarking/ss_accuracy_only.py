import pandas as pd
import numpy as np
import pickle as pkl
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

spotrna_TS1 = ['4r4v_A', '4wfa_Y', '3suh_X', '4l81_A', '3ktw_C', '3adb_C', '5u3g_B', '4qg3_B', '3amt_B', '3d2v_A', '1c0a_B', '1qf6_B', '1gax_C', '2csx_C', '2dr2_B', '2zni_C', '2oiu_P', '5kpy_A', '5tpy_A', '4pqv_A', '1un6_E', '1mzp_B', '1u63_B', '1i6u_C', '4pcj_A']
spotrna_TS2 = ['3pdr_A', '3oww_A', '4mgn_A', '3q3z_A','2qus_A', '3slm_A', '6u8k_A', '3r4f_A', '2qwy_A', '4c7o_E', '4rmo_B', '4pmi_A', '5bjp_E', '5vof_A']

#TS1_ids = [i for i in TS1_ids if i not in spotrna_TS1]
#TS2_ids = [i for i in TS2_ids if i not in spotrna_TS2]


#predictors = ['SPOT-RNA', 'SPOT-RNA2', 'RNAfold', 'LinearPartition', 'SPOT-RNA-2D-Single', 'SPOT-RNA-2D']
predictors = ['RNAfold', 'LinearPartition', 'SPOT-RNA-2D-Single', 'SPOT-RNA-2D']
sets = [TS1_ids, TS2_ids, TS3_ids]
#sets = [['6p2h_A']]
print('\n \t\t     TS1   TS2   TS3' )

for number, pred in enumerate(predictors):

	for num, ids in enumerate(sets[0:]):

		#print()

		count = 0; native_ss_count = 0; native_non_ss_count = 0
		f1_list = []

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
#			print(len(native_dist_contacts), len(native_non_ss_pairs), len(native_ss_base_pairs))

			if pred == 'RNAfold':
				pred_pairs = RNAfold_bps(base_path, k, seq)
			elif pred == 'LinearPartition':
				pred_pairs = LinearPartition_bps(base_path, k, seq)
			elif pred == 'SPOT-RNA':
				pred_pairs = spot_rna_bps(base_path, k, seq)
			elif pred == 'SPOT-RNA2':
				pred_pairs = spot_rna2_bps(base_path, k, seq)
			elif pred == 'SPOT-RNA-2D-Single':
				pred_contacts_matrix = spotrna_2d_single(base_path, k, seq)
				# extract upper triangular top L long-range contacts predicted proability LxL matrix
				tri_inds = np.triu_indices(pred_contacts_matrix.shape[0], k=24)
				all_pred_contacts = np.array([[i,j,pred_contacts_matrix[i,j]] for i,j in zip(tri_inds[0], tri_inds[1]) if [i, j] not in native_non_ss_pairs])
#				print(len(tri_inds[0]), len(all_pred_contacts))
				all_pred_contacts_sorted = all_pred_contacts[all_pred_contacts[:,2].argsort()[::-1]]
				pred_pairs = [[int(i[0]), int(i[1])] for i in all_pred_contacts_sorted[0:int(len(seq)/n)]]
			elif pred == 'SPOT-RNA-2D':
				pred_contacts_matrix = spotrna_2d(base_path, k, seq)
				# extract upper triangular top L long-range contacts predicted proability LxL matrix
				tri_inds = np.triu_indices(pred_contacts_matrix.shape[0], k=24)
				all_pred_contacts = np.array([[i,j,pred_contacts_matrix[i,j]] for i,j in zip(tri_inds[0], tri_inds[1]) if [i, j] not in native_non_ss_pairs])
#				print(len(tri_inds[0]), len(all_pred_contacts))
				all_pred_contacts_sorted = all_pred_contacts[all_pred_contacts[:,2].argsort()[::-1]]
				pred_pairs = [[int(i[0]), int(i[1])] for i in all_pred_contacts_sorted[0:int(len(seq)/n)]]

			pred_contacts = pred_pairs
			pred_contacts = [i for i in pred_contacts if abs(i[0]-i[1]) >= 24]  # extract long-range base-pairs only

			true_pairs = native_ss_base_pairs


			pred_correctly = [i for i in pred_contacts if i in true_pairs]; #print(len(pred_correctly))
			pred_wrongly = [i for i in pred_contacts if i not in true_pairs]; #print(len(pred_wrongly))

			tp = len(pred_correctly)
			fp = len(pred_wrongly)
			fn = len([i for i in true_pairs if i not in pred_contacts])

			try:
				pre = tp / (tp + fp)
				sen = tp / (tp + fn)
				f1 = 2*((pre*sen)/(pre + sen))
			except:
				f1 = 0

			f1_list.append(f1)

			count += 1

		if num==0:
			print(pred + " "*(len('SPOT-RNA-2D-Single')-len(pred)),end=' ')
		print(' {0:.3f}'.format(np.mean(f1_list)), end='')

	print()

#print(count)
print()

