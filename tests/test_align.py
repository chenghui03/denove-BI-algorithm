import setwd

import unittest

from modules.seq import MySeq
from modules.align_utils import Table, AlignResult, TraceBack, SubstMatrix
from modules.align import PairwiseAlign, ProgressiveAlign

sm = SubstMatrix.create_submat(1, -1, "ATCG")
align = PairwiseAlign(sm, -1)

seq2_case1 = list(map(MySeq, ["AT","T"]))
seq2_case2 = list(map(MySeq, ["ATCGATCG","CGAT"]))

# traceback
seq2_case3 = list(map(MySeq, ["TATA","TA"]))
seq2_case4 = list(map(MySeq, ["GCGCGCGC","GC"]))
seq2_case5 = list(map(MySeq, ["ATCG","ATCG"]))

# mutiple alignment
seq3_case1 = list(map(MySeq, ["ATC", "T", "TC"]))
seq3_case2 = list(map(MySeq, ["ATCG","ATCG","ATCG"]))

# 
seq4_case1 = [AlignResult(["ATCG", "A_CG","ATCG"], "DNA"),
              MySeq("TCG", "DNA")]

seq4_case2 = [MySeq("TCAAGAG", "DNA"),
              MySeq("TCGAA", "DNA"),
              MySeq("TATTCG", "DNA"),
              MySeq("TCAGTT", "DNA"),
              MySeq("TCAGTAAT", "DNA")]

seq4_case3 = [MySeq("ATCG", "DNA"),
              MySeq("ATCG", "DNA")]

class test_PairwiseAlign(unittest.TestCase):

        
    def test_after_align(self):
        with self.assertRaises(ValueError):
            align.get_align()

    def test_get_match_score_x(self):
        seqs = [MySeq("ATTC"), MySeq("ATG")]
        self.assertEqual(PairwiseAlign.get_match_score_x(seqs, [3,2], "11", -1, sm), 1)
        self.assertEqual(PairwiseAlign.get_match_score_x(seqs, [2,2], "11", -1, sm), 1)
        self.assertEqual(PairwiseAlign.get_match_score_x(seqs, [4,3], "10", -1, sm), 0)
        self.assertEqual(PairwiseAlign.get_match_score_x(seqs, [4,3], "01", -1, sm), 0)
    
    def test_init_table(self):
        _align = PairwiseAlign(sm, -1, False)
        _align.seqs = seq2_case1
        _align._init_table()

        self.assertEqual(_align._S.to_list(), [[0, -1, -2],[-1, 0, 0]])
        self.assertEqual(_align._T.to_list(), [[[],["10"],["10"]],[["01"], [], []]])
        
        #TODO 3seq
    
    def test_align(self):
        # needleman_Wunsch
        align = PairwiseAlign(sm, -1, is_local=False) # AT T
        align.align(seq2_case1)
        self.assertEqual(align._S.to_list(),
                         [[0, -1, -2],
                          [-1, -1, 0]])
        
        self.assertEqual(align._T.to_list(),
                         [[[], ['10'], ['10']],
                          [['01'], ['11'], ['11']]])

        # smith_Waterman(local align)
        align = PairwiseAlign(sm, -1, is_local=True)

        align.align(seq2_case1)
        
        self.assertEqual(align._S.to_list(), [[0,0,0],
                                              [0,0,1]])
        self.assertEqual(align._T.to_list(), [[[],[],[]],
                                              [[],[],["11"]]])
        self.assertEqual(align.max_pos_list, [[2,1],])

    def test_get_align(self):
        # needleman_Wunsch
        align = PairwiseAlign(sm, -1, is_local=False)
        
        ## no tie
        align.align(seq2_case1)
        res = align.get_align()
        self.assertEqual(res[0][0], "AT")
        self.assertEqual(res[0][1], "_T")

        align.align(seq2_case2)
        res = align.get_align()
        self.assertEqual(res[0][0], "ATCGATCG")
        self.assertEqual(res[0][1], "__CGAT__")
        
        ## tie
        align.align(seq2_case3)
        res = align.get_align()
        res = [res[i][1] for i in range(len(res))]
        expected = ["TA__", "__TA", "T__A",] 
        self.assertEqual(sorted(res), sorted(expected))

        align.align(seq2_case4) # GCGCGCGC GC
        res = align.get_align()
        res = [res[i][1] for i in range(len(res))]
        
        def __func(i,j):
            res = []
            for x in range(8):
                if x==i:
                    res.append("G")
                elif x==j:
                    res.append("C")
                else:
                    res.append("_")
            return "".join(res)
        
        expected = []
        for i in range(0, 8, 2):
            for j in range(i+1, 8, 2):
                expected.append(__func(i,j))

        self.assertEqual(sorted(res), sorted(expected))

        # smith_Waterman(local align)
        align = PairwiseAlign(sm, -1, is_local=True)

        # no tie
        align.align(seq2_case1) #AT T
        
        res = align.get_align()
        self.assertEqual(res[0][0],"T")
        self.assertEqual(res[0][1],"T")
    
    def test_evaluate_align(self):
        align = PairwiseAlign(sm, -1, is_local=True)
        align.align(seq2_case5)
        self.assertEqual(align.evaluate_align(), {"score": 4, "identity": 1})

class test_ProgrossiveAlign(unittest.TestCase):
    def test_after_align(self):
        align = ProgressiveAlign(sm, -1)
        with self.assertRaises(ValueError):
            align.get_align()

    def test_align(self):
        align = ProgressiveAlign(sm, -1)
        align.align(seq4_case1)
        self.assertEqual(align.get_align().num_seqs(), 4)

        align.align(seq4_case2)
        self.assertEqual(align.get_align().num_seqs(), 5)

    def test_evaluate_align(self):
        align = ProgressiveAlign(sm, -1)
        align.align(seq4_case3)
        self.assertEqual(align.evaluate_align(), 4)

if __name__ == "__main__":
     unittest.main()