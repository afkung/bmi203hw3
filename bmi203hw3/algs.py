"""
Andrew Kung
BMI203 HW #3
Algorithms for implementing the Smith-Waterman algorithm
Last modified: 2/23/18
"""

# importing useful libraries
import numpy as np
import matplotlib.pyplot as plt
import copy
import random


# loadMatrix: function to load scoring matrix
# input: file name of matrix to use (must be in directory)
# ouput: scoring matrix, dictionary of indices for amino acids
def loadMatrix(file_name):
	f = open(file_name, 'r')
	no_title = True
	titles = []
	scoring_matrix = []
	for line in f:
		if line[0] != '#':
			if no_title:
				no_title = False
				titles = line.strip().split()
			else:
				scoring_matrix.append(line.strip().split())
	index_dict = {}
	for index in range(len(titles)):
		index_dict[titles[index]] = index
	return scoring_matrix, index_dict

# loadSequence: function to load sequence
# input: file name of sequence (must be in directory)
# output: sequence as string
def loadSequence(file_path):
	f = open(file_path, 'r')
	sequence = ""
	for line in f:
		if line[0] != '>':
			sequence += str(line).strip()
	f.close()
	return sequence

# fillTable: function to fill out alignment table using Smith-Waterman algorithm
# input: seq1 and seq2 as strings, scoring matrix as list of lists and index dictionary, gap penalties
# output: filled out table as list of lists
def fillTable(seq1, seq2, scoring_matrix, index_dict, open_penalty, extend_penalty):
	table = [[0]*len(seq2) for x in range(len(seq1))]
	for index1 in range(len(seq1)):
		for index2 in range(len(seq2)):
			if index1 == 0 or index2 == 0:
				table[index1][index2] = 0
			else:
				if index2 > 1 and table[index1][index2-1] == (table[index1][index2-2] - open_penalty): # if gap from left is already open, extend. Otherwise open
					left = table[index1][index2-1] - extend_penalty
				else:
					left = table[index1][index2-1] - open_penalty						
				if index1 > 1 and table[index1-1][index2] == (table[index1-2][index2] - open_penalty): # if gap from above is already open, extend. Otherwise open
					up = table[index1-1][index2] - extend_penalty
				else:
					up = table[index1-1][index2] - open_penalty
				diagonal = table[index1-1][index2-1] + float(scoring_matrix[index_dict[seq1[index1].upper()]][index_dict[seq2[index2].upper()]]) # match, with score from scoring matrix
				table[index1][index2] = max(left, up, diagonal, 0)
	return table

# traceBack: function to trace back optimal score given filled out table
# input: table with values from S-W implementation
# output: alignment score as float
def traceBack(table):
	max_score = 0.0
	for row in table:
		if max(row) > max_score:
			max_score = max(row)
	return max_score

# align: operator function that aligns two sequences using functions above
# input: seq1 and seq2 as strings, scoring matrix as list of lists and associated index dictionary, gap penalties
# output: alignment score as float
def align(seq1, seq2, scoring_matrix, index_dict, open_penalty, extend_penalty):
	return traceBack(fillTable(seq1, seq2, scoring_matrix, index_dict, open_penalty, extend_penalty))

# function to score the known positive and negative pairs using a given scoring matrix
# Input: Scoring matrix, associated dictionary of indices, gap penalties
# Output: list of positive and negative scores as floats
def scoreMatrix(scoring_matrix, index_dict, open_penalty, extend_penalty):
	pos_scores = []
	f = open("Pospairs.txt",'r')
	for item in f:
		path1 = item.split()[0]
		seq1 = loadSequence(path1)
		path2 = item.split()[1]
		seq2 = loadSequence(path2)
		pos_scores.append(align(seq1, seq2, scoring_matrix, index_dict, open_penalty, extend_penalty))
# uncomment following line for normalized score
#		pos_scores.append(align(seq1, seq2, scoring_matrix, index_dict, open_penalty, extend_penalty)/min(len(seq1),len(seq2)))
	f.close()
	# finding scores for known negatives
	neg_scores = []
	g = open("Negpairs.txt",'r')
	for item in g:
		path1 = item.split()[0]
		seq1 = loadSequence(path1)
		path2 = item.split()[1]
		seq2 = loadSequence(path2)
		neg_scores.append(align(seq1, seq2, scoring_matrix, index_dict, open_penalty, extend_penalty))
# uncomment following line for normalized score
#		neg_scores.append(align(seq1, seq2, scoring_matrix, index_dict, open_penalty, extend_penalty)/min(len(seq1),len(seq2)))
	g.close()
	return pos_scores, neg_scores
	
# function to find the false positive rate based on known positive and negative pairs
# Input: true positive threshold as float, scores for known +/- pairs as lists of floats
# Output: false positive rate from known negative pairs
def findFP(TP_threshold, pos_scores, neg_scores):
	pos_scores.sort()
	score_cutoff = pos_scores[int(TP_threshold * len(pos_scores))-1]
	# counting how many negative pairs scored above the threshold
	FP_count = 0
	for item in neg_scores:
		if item > score_cutoff:
			FP_count += 1
	return FP_count / len(neg_scores)

# function to find the true positive rate based on known positive and negative pairs
# Input: false positive threshold as float, scores for known +/- pairs as lists of floats
# Output: true positive rate from known positive pairs
def findTP(FP_threshold, pos_scores, neg_scores):
	neg_scores.sort()
	score_cutoff = neg_scores[int((1-FP_threshold) * len(neg_scores))-1]
	# counting how many positive pairs scored above the threshold
	TP_count = 0
	for item in pos_scores:
		if item > score_cutoff:
			TP_count += 1
	return TP_count / len(pos_scores)

# function to plot a ROC curve, with false positive values on the x-axis and false negative values on the y-axis
# Input: list of FP values, list of TP values, name of matrix as string
# Output: ROC plot, displayed in interface as saved as .png
def rocCurve(FP_values, TP_values, name):
	area = 0.0
	for item in TP_values:
		area += item * (1 / len(FP_values))  # area under curve estimate using rectangle method
	fig, ax = plt.subplots()
	ax.set(xlabel = "False Positives", ylabel = "True Positives", title = "ROC " + name + " Area = " + str(area))
	ax.grid()
	plt.ylim(0,1)
	plt.xlim(0,1)
	ax.plot(FP_values, TP_values)
	fig.savefig(name+"ROC.png")
	
# findFitness: function to find fitness score (defined as sum of TP rates for FP = 0.0,0.1,0.2,0.3) given a scoring matrix
# input: scoring matrix, index dictionary to test, gap penalties
# output: fitness score as float
def findFitness(scoring_matrix, index_dict, open_penalty, extend_penalty):
	pos, neg = scoreMatrix(scoring_matrix, index_dict, open_penalty, extend_penalty)
	fitness_score = findTP(0.0, pos, neg) + findTP(0.1, pos, neg) + findTP(0.2, pos, neg) + findTP(0.3, pos, neg)
	return fitness_score

# optimizeMatrix: function to optimize a given scoring matrix using one iteration of greedy pseudo-gradient descent
# input: scoring matrix, associated index dict, gap penalties
# output: optimized matrix
def optimizeMatrix(scoring_matrix, index_dict, open_penalty, extend_penalty):
	test_matrix = copy.deepcopy(scoring_matrix)
	step_sizes = [0.1,-0.1] # trying steps in both directions
	for step_size in step_sizes:
		for index1 in range(len(scoring_matrix)):
			for index2 in range(index1,len(scoring_matrix)):
				while findFitness(test_matrix, index_dict, open_penalty, extend_penalty) >= findFitness(scoring_matrix, index_dict, open_penalty, extend_penalty): # try stepping in a direction, if better then keep going
					test_matrix[index1][index2] = float(test_matrix[index1][index2]) + step_size
					test_matrix[index2][index1] = float(test_matrix[index2][index1]) + step_size
					print("Add 1")
				scoring_matrix = copy.deepcopy(test_matrix)
				print("Index 1:"+str(index1)+" Index 2:"+str(index2))
	return scoring_matrix, index_dict

# optimizeMatrixFast: function to optimize given scoring matrix, similar to optimizeMatrix but instead of iterating across all values, randomly selects 100 to test
# input: scoring matrix, associated index dict, gap penalties
# output: optimized matrix
def optimizeMatrixFast(scoring_matrix, index_dict, open_penalty, extend_penalty):
	test_matrix = copy.deepcopy(scoring_matrix)
	step_size = 0.1
	iterations = 100
	for iter in range(iterations):
		step_size = step_size * (random.randint(0,1)*2-1) # either +1 or -1, step in either direction
		index1 = random.randint(0,len(scoring_matrix)-1)
		index2 = random.randint(0,len(scoring_matrix)-1)
		test_matrix[index1][index2] = float(test_matrix[index1][index2]) + step_size
		test_matrix[index2][index1] = float(test_matrix[index2][index1]) + step_size
		test = findFitness(test_matrix, index_dict, open_penalty, extend_penalty)
		scoring = findFitness(scoring_matrix, index_dict, open_penalty, extend_penalty)
		print(str(test) + '  ' + str(scoring))
#		if findFitness(test_matrix, index_dict, open_penalty, extend_penalty) > findFitness(scoring_matrix, index_dict, open_penalty, extend_penalty): # try stepping in a direction, if better then save
		if test > scoring:
			scoring_matrix = copy.deepcopy(test_matrix)
			print("CHANGE")
		print(iter)
	return scoring_matrix, index_dict