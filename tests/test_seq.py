import setwd
import unittest
from modules.seq import MySeq

class test_MySeq(unittest.TestCase):
    s1 = MySeq("ATCG")

    def test_validate(self):
        self.assertEqual(test_MySeq.s1.validate(), True)

        seq = MySeq("XXXX","DNA")
        with self.assertRaises(ValueError):
            seq.validate()

    def test_trans(self):

        # transcribe
        self.assertEqual(test_MySeq.s1.transcribe().seq, "AUCG")

        # retrotranscribe
        seq = MySeq("AUGUAA","RNA")
        self.assertEqual(seq.retro_transcribe().seq, "ATGTAA")
        
        # translate
        seq = MySeq("ATGATGATG","DNA")
        self.assertEqual(seq.translate().seq, "MMM")

        seq = MySeq("ATGTAA","DNA")
        self.assertEqual(seq.translate().seq, "M_")

        seq = MySeq("AUGUAA","RNA")
        self.assertEqual(seq.translate().seq, "M_")

    def test_repr(self):
        repr_res = []
        seq = MySeq("ATCGCGCG", "DNA")
        repr_res = [seq.get_biotype(), seq.alphabet(), repr(seq)]
        for i in repr_res:
            print(i)
            
if __name__=="__main__":
    unittest.main()