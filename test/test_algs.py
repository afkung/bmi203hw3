import numpy as np
from example import algs

def test_bubblesort():
	x = np.array([1,2,4,0,1])
	x,assign,cond = algs.bubblesort(x))
	assert x = [0,1,1,2,4]

def test_quicksort():
	x = np.array([1,2,4,0,1])
	x,assign,cond = algs.quicksort(x,0,len(x)-1,0,0)
	assert x = [0,1,1,2,4]