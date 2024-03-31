import setwd

import unittest

from modules.seq import MySeq
from modules.align_utils import Table, AlignResult, TraceBack, SubstMatrix

class test_SubstMatrix(unittest.TestCase):
    sm = SubstMatrix.create_submat(1, -1, "ATCG")

    def test_create_sm(self):

        self.assertEqual(test_SubstMatrix.sm["CC"], 1)
        self.assertEqual(test_SubstMatrix.sm["CG"], -1)
        with self.assertRaises(ValueError):
            test_SubstMatrix.sm["X"]


class test_Align(unittest.TestCase):
    def test_consensus(self):
        case1 = AlignResult(["ATCG","AC_G","AT_C"], "DNA")
        res = case1.consensus()
        self.assertEqual(res, "ATCG")


        case2 = AlignResult(["AC_G","AC_G"], "DNA")
        with self.assertRaises(ValueError):
            case2.consensus()
            

class test_Table(unittest.TestCase):
    # seq1 xxx
    # seq2 xx
    case1 = Table([[0, -3, -2],
                   [-1, 1, 2]])
    
    # seq1,2,3 = xxx, xx, xx
    # seq_len_s 3,2,2
    # shape: seq_len_s[::-1]
    case2 = [[[0, 1, 2], [3, 4, 5]],
             [[0, -1, -2], [-3, -4, -5]]]
    case2 = Table(case2)

    def test_get_item(self):
        self.assertEqual(test_Table.case1[(2, 0)], -2)
        self.assertEqual(test_Table.case1[(2, 1)], 2)
        self.assertEqual(test_Table.case1[(0, 1)], -1)
        self.assertEqual(test_Table.case1[(1, 0)], -3)

        self.assertEqual(test_Table.case2[(2, 1, 1)], -5)
        self.assertEqual(test_Table.case2[(2, 0, 1)], -2)


class test_Recover_record(unittest.TestCase):
    case = list(map(MySeq,
                    ["TT","T"]))
    trace_table = Table([[[], ['10'], ['10']],
                         [['01'], ['11'], ['10', '11']]]) 
    seq2case1 = TraceBack(case, trace_table)

    case = list(map(MySeq,
                    ["ATAT","AT"]))
    trace_table = Table([[[], ['10'], ['10'], ['10'], ['10']],
                         [['01'], ['11'], ['10'], ['10', '11'], ['10']],
                         [['01'], ['01'], ['11'], ['10'], ['10', '11']]]) 
    seq2case2 = TraceBack(case, trace_table)

    def test_is_done(self):
        res = test_Recover_record.seq2case1.is_done()
        self.assertEqual(res, False)
        
    def test_step(self):
        # seq2case1.T
        # [[], ['10'], ['10']]
        # [['01'], ['11'], ['10','11']]
        res = test_Recover_record.seq2case1.step()

        _res = [res[0].get_pos(), res[1].get_pos()]
        expected = [[2,1], [2,1]]
        self.assertEqual(sorted(_res), sorted(expected))

        _res = [res[0].prior_distance, res[1].prior_distance]
        expected = [['10'], ['11']]
        self.assertEqual(sorted(_res), sorted(expected))

        _res = res[1].step()[0].get_pos()
        expected = [[1,0], [1,1]]
        self.assertIn(_res, expected)

        # seq2case1
        res = test_Recover_record.seq2case2.step()

        _res = [res[0].prior_distance, res[1].prior_distance]
        expected = [['11'], ['10']]
        self.assertEqual(sorted(_res), sorted(expected))

        _res = [res[0].step()[0].prior_distance, res[1].step()[0].prior_distance]
        expected = ['10'] #['11', '10']
        self.assertIn(expected, _res)

    def test_recover_align(self):
        res = test_Recover_record.seq2case1.recover_align()
        res = [res[0][1], res[1][1]]
        expected = ["_T","T_"]
        self.assertEqual(sorted(res), sorted(expected))

        res = test_Recover_record.seq2case2.recover_align()
        res = [res[i][1] for i in range(3)]
        expected = ["AT__", "__AT", "A__T"]
        self.assertEqual(sorted(res), sorted(expected))

if __name__=="__main__":
    unittest.main()