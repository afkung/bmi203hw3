from bmi203hw3 import algs
	
def test_load_matrix():
	scoring_matrix, index_dict = algs.loadMatrix("PAM100")
	assert index_dict == {'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C': 4, 'Q': 5, 'E': 6, 'G': 7, 'H': 8, 'I': 9, 'L': 10, 'K': 11, 'M': 12, 'F': 13, 'P': 14, 'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19, 'B': 20, 'Z': 21, 'X': 22, '*': 23}
	assert scoring_matrix[0] == ['4', '-3', '-1', '-1', '-3', '-2', '0', '1', '-3', '-2', '-3', '-3', '-2', '-5', '1', '1', '1', '-7', '-4', '0', '-1', '-1', '-1', '-9']

def test_load_sequence():
	seq1 = algs.loadSequence("sequences/prot-0004.fa")
	assert seq1 == "SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM"
	seq2 = algs.loadSequence("sequences/prot-0008.fa")
	assert seq2 == "ANKTRELCMKSLEHAKVDTSNEARQDGIDLYKHMFENYPPLRKYFKSREEYTAEDVQNDPFFAKQGQKILLACHVLCATYDDRETFNAYTRELLDRHARDHVHMPPEVWTDFWKLFEEYLGKKTTLDEPTKQAWHEIGREFAKEINK"
	
def test_align():
	scoring_matrix, index_dict = algs.loadMatrix("PAM100")
	seq1 = "BAA"
	seq2 = "AAA"
	assert algs.traceBack(algs.fillTable(seq1, seq2, scoring_matrix, index_dict, 0, 0)) == 2 * float(scoring_matrix[index_dict['A']][index_dict['A']])
"""
Expected alignment: BAA or BAA
                    AAA or -AAA
Expected score: (match A) * 2
"""
