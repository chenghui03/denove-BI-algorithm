import setwd

from abc import ABC,abstractmethod
from itertools import product
import functools 

from modules.seq import MySeq
from modules.align_utils import Table, AlignResult, TraceBack, SubstMatrix

class Align(ABC):
    def __init__(self, substMatrix: SubstMatrix) -> None:
        self.sm = substMatrix
        self.seqs = []
    
    @abstractmethod
    def after_align(func):
        pass
    
    @abstractmethod
    def align(self, seqs: list[MySeq]) -> None:
        pass
    
    @abstractmethod
    def get_align(self):
        pass
    
    @abstractmethod
    def evaluate_align(self):
        pass


class PairwiseAlign(Align):
    """smith and waterman"""
    def __init__(self, substMatrix: SubstMatrix, gap: int, is_local: bool=False) -> None:
        self.is_local = is_local
        self.gap = gap
        self._S = None
        self._T = None
        super().__init__(substMatrix)

    def _set_ST(self, seq1, seq2):
        #TODO support 3
        S = [[ 0 for i in range(len(seq1) + 1)] for j in range(len(seq2) + 1)]
        T = [[ [] for i in range(len(seq1) + 1)] for j in range(len(seq2) + 1)]
        self._S, self._T = Table(S), Table(T)

    def _init_table(self):
        if None: 
            raise ValueError("not same seq type") #TODO left

        def onehot(pos, vec_len) -> str:
            """apply to 2: what about 3"""
            return ''.join(['1' if i == pos else '0' for i in range(vec_len)])

        seq1, seq2 = self.seqs
        self._set_ST(seq1, seq2)
        if self.is_local:
            return None

        # seqs = [AT, T]
        # len_seq_s = [2, 1]
        # i_seq = [0, 1]
        # for x [1,2  1]

        n_seqs = len(self.seqs)

        for i_seq, seq in enumerate(self.seqs):
            for x in range(1, len(seq) + 1):
                # index: [0, 0, x, 0, ...]

                index = [0 if i == i_seq else x for i in range(n_seqs)][::-1]

                temp = self.gap * x
                self._S[index] = temp
                self._T[index] = [onehot(i_seq, n_seqs),]
    
    @staticmethod
    def get_match_score_x(seqs: list[MySeq],
                    pos: tuple,
                    not_gap_flags: str,
                    gap: int,
                    sm: SubstMatrix
                    ) -> int:
        # TODO
        """
        get sum of match score of x nums of chars

        pos:
        seqs = [ATC, T, TC]
        for compare: T T T
        pos =  [2, 1, 1]
        warning : pos is (row, col) not (row, col)
        """
        def get_score_2(c1, c2, gap, sm):
            if any(c == "-" for c in [c1, c2]):
                return gap
            else:
                return sm[c1, c2]

        not_gap = []
        for flag, i, seq in zip(not_gap_flags, pos, seqs):
            if flag == "1":
                temp = seq[i-1]
                not_gap.append(temp)

        not_gap = [seq[i-1] for flag, i, seq in zip(not_gap_flags, pos, seqs) if flag == "1"] 
        # pre 0111
        # pos 3,4,5,2
        # 4,5,2

        match_score = 0
        for i in range(len(not_gap)):
            for j in range(i+1, len(not_gap)):
                match_score += get_score_2(not_gap[i], not_gap[j], gap, sm)

        return match_score

    @staticmethod
    def get_gap_score_x(prior_distance: str, gap: int) -> int:
        """get sum of gap score of x nums of chars"""
        gap_n = 0
        for i in range(len(prior_distance)):
            for j in range(i+1, len(prior_distance)):
                if any(x == "0" for x in [prior_distance[i], 
                                          prior_distance[j]]):
                    gap_n +=1
        
        return gap_n * gap
    
    @staticmethod
    def after_align(func):
        @functools.wraps(func)
        def wrapper(self, *arg, **kwargs):
            if any(x == None for x in [self._S, self._T]):
                raise ValueError(f"method '{func.__name__}' should only be execuated after align")

            return func(self, *arg, **kwargs)
        return wrapper
    
    def get_max_nbr(self, score_list: list) -> tuple:
        max_score = max(score_list)

        prior_distance = []

        if self.is_local:
            max_score = max(max_score, 0)
            if max_score == 0:
                return 0, []
        
        for i, score in enumerate(score_list):
            if score == max_score:
                nbr_dist = format(i + 1, "b").zfill(len(self.seqs))
                prior_distance.append(nbr_dist)

        return max_score, prior_distance
    
    def align(self, seqs: list[MySeq]) -> None:
        self.seqs = seqs
        self._init_table()

        if self.is_local:
            self.max_S_score = 0
            self.max_pos_list = []

        s_list = [None for i in range(2**len(seqs) - 1)] 

        pos_list = product(*(range(1, len(x) + 1) for x in seqs))

        for i, pos in enumerate(pos_list):
            for prior_i in range(1, 2**len(seqs)):
                prior_distance = format(prior_i, "b").zfill(len(seqs))
                prior_index = [i-int(j) for i,j in zip(pos, prior_distance)]

                prior_score = self._S[prior_index]
                match_score = self.get_match_score_x(self.seqs, pos, prior_distance, self.gap, self.sm) 
                gap_score = self.get_gap_score_x(prior_distance, self.gap)

                s_list[prior_i - 1] = prior_score + match_score + gap_score

            self._S[pos], self._T[pos] = self.get_max_nbr(s_list) # is_local 

            if self.is_local:
                if self._S[pos] > self.max_S_score:
                    self.max_S_score = self._S[pos]
                    self.max_pos_list = [list(pos)]
                elif self._S[pos] == self.max_S_score:
                    self.max_pos_list.append(list(pos))
                else:
                    pass

    @after_align
    def get_align(self):
        return TraceBack(self.seqs, self._T).recover_align()
    
    @after_align
    def evaluate_align(self):
        def get_last_index(seqs):
            return [len(seq) for seq in self.seqs]

        def identity(seqs: list[MySeq]):
            sm = SubstMatrix.create_submat(1, 0, seqs[0].alphabet())
            align = PairwiseAlign(sm, 0, is_local=False)
            align.align(seqs)

            matches = align._S[get_last_index(seqs)]
            identity_score = matches / max(map(len, seqs))
            return identity_score

        score = self._S[get_last_index(self.seqs)]
        identity_s = identity(self.seqs)
        
        return {"score":score, "identity": identity_s}
        