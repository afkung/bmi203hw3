"""
Andrew Kung
BMI203 HW #3
Answering questions using the Smith-Waterman algorithm
Last modified: 2/23/18
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import .algs

def run_stuff():

	"""
	# 1-1 Testing Gap Penalties
	scoring_matrix, index_dict = loadMatrix("BLOSUM50")

	score_array = [[0]*5 for x in range(20)]
	print("FP rates:")
	print("	1	2	3	4	5")
	for open_penalty in range(1,21): # iterating across possible penalties
		print(open_penalty, end = '	')
		for extend_penalty in range(1,6):
			pos, neg = scoreMatrix(scoring_matrix, index_dict, open_penalty, extend_penalty)
			score = findFP(0.7, pos, neg)
			print(score,end = '	')
			score_array[open_penalty-1][extend_penalty-1] = score
		print()
	# heat map
	plt.imshow(score_array, cmap = 'hot')
	plt.show()
	"""


	"""
	# 1-2 Testing Matrices
	matrices = ["BLOSUM50","BLOSUM62","MATIO","PAM100","PAM250"]

	open_penalty = 6 #from 1-1
	extend_penalty = 3 #from 1-1

	for matrix in matrices: # iterating across different matrices
		scoring_matrix, index_dict = loadMatrix(matrix)
		pos, neg = scoreMatrix(scoring_matrix, index_dict, open_penalty, extend_penalty)
		roc_values = []
		for FP_threshold in np.arange(0.0,1.0,0.02): # step size of 1/50 since there are 50 pos/neg pairs
			roc_values.append(1-findFP(FP_threshold, pos, neg))
		rocCurve(np.arange(0.0,1.0,0.02),roc_values, matrix)
	"""

	"""
	# 1-3 Changing Scores
	open_penalty = 6 # penalties from Part 1
	extend_penalty = 3
	scoring_matrix, index_dict = loadMatrix("BLOSUM50")
	pos, neg = scoreMatrix(scoring_matrix, index_dict, open_penalty, extend_penalty)
	roc_values = []
	for FP_threshold in np.arange(0.0,1.0,0.02): # step size of 1/50 since there are 50 pos/neg pairs
		roc_values.append(1-findFP(FP_threshold, pos, neg))
	rocCurve(np.arange(0.0,1.0,0.02),roc_values, "BLOSUM50 Norm")
	"""

	# 2-1 Optimization Functions




	# 2-2 Optimizing Best scoring matrix using above functions
	open_penalty = 6 # penalties from Part 1
	extend_penalty = 3
	scoring_matrix, index_dict = loadMatrix("BLOSUM50")
	optimized_matrix = optimizeMatrixFast(scoring_matrix, index_dict, open_penalty, extend_penalty)
	print(index_dict)
	for row in optimized_matrix:
		for item in row:
			print(item, end = '	')
		print()
	print("Original Score: " + str(findFitness(scoring_matrix, index_dict, open_penalty, extend_penalty)))
	print("Optimized Score: " + str(findFitness(optimized_matrix, index_dict, open_penalty, extend_penalty)))

	# plotting ROC curves
	pos_orig, neg_orig = scoreMatrix(scoring_matrix, index_dict, open_penalty, extend_penalty)
	pos_opt, neg_opt = scoreMatrix(optimized_matrix, index_dict, open_penalty, extend_penalty)
	roc_values_orig = []
	roc_values_opt = []
	for FP_threshold in np.arange(0.0,1.0,0.02): # step size of 1/50 since there are 50 pos/neg pairs
		roc_values_orig.append(1-findFP(FP_threshold, pos_orig, neg_orig))
		roc_values_opt.append(1-findFP(FP_threshold, pos_opt, neg_opt))
	rocCurve(np.arange(0.0,1.0,0.02),roc_values_orig, "Original BLOSUM50")
	rocCurve(np.arange(0.0,1.0,0.02),roc_values_opt, "Optimized BLOSUM50")



	"""
	# 2-3 Optimizing MATIO using above functions
	open_penalty = 6 # penalties from Part 1
	extend_penalty = 3
	scoring_matrix, index_dict = loadMatrix("MATIO")
	optimized_matrix = optimizeMatrix(scoring_matrix, index_dict, open_penalty, extend_penalty)
	print("Original Score: " + str(findFitness(scoring_matrix, index_dict, open_penalty, extend_penalty)))
	print("Optimized Score: " + str(findFitness(optimized_matrix, index_dict, open_penalty, extend_penalty)))

	# plotting ROC curves
	roc_values = []
	opt_roc_values = []
	for FP_threshold in np.arange(0.0,1.0,0.02): # since there are 50 pos/neg pairs
		roc_values.append(findFP(FP_threshold, scoring_matrix))
		opt_roc_values.append(findFP(FP_threshold, optimized_matrix))
	rocCurve(np.arange(0.0,1.0,0.02),roc_values, "MATIO Classic")
	rocCurve(np.arange(0.0,1.0,0.02),opt_roc_values, "New MATIO")

"""