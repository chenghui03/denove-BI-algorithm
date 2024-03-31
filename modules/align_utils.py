import setwd

from collections import Counter
import copy

from modules.seq import MySeq
from modules.constants import * 


class SubstMatrix:
    def __init__(self, alphabet: list, sm: dict) -> None:
        self.alphabet = alphabet
        self.sm = sm

    def __repr__(self) -> str:
        return f"Alphabet: {self.alphabet}. sm: ..."

    def __score_pair(self,c1,c2) -> int:
        if all(c in self.alphabet for c in (c1, c2)):
            return self.sm[c1+c2]
        else:
            raise ValueError('char not in alphabet')

    def __getitem__(self,ft) -> int:
        from_word,to_word = ft
        return self.__score_pair(from_word,to_word)

    @classmethod
    def from_dict(cls, d: dict):
        sm = d
        alphabet = list({key[0] for key in d.keys()})
        return cls(alphabet, sm)
    
    @classmethod
    def create_submat(cls, match: int, mismatch: int, alphabet: str):
        alphabet = list(alphabet)
        sm = {}
        for c1 in alphabet:
            for c2 in alphabet:
                sm[c1+c2] = match if c1 == c2 else mismatch

        return cls(alphabet, sm)
    
    @classmethod
    def from_file(self):
        raise NotImplementedError


class AlignResult:
    def __init__(self,
                 seqs: list[str],
                 seq_type: str="DNA") -> None:
        '''
        seq_type: str["DNA"|"RNA"|"protein"]
        '''
        self.seqs = seqs
        self.seq_type = seq_type

    def __len__(self):
        return len(self.seqs[0])
    
    def __getitem__(self, args):
        if isinstance(args,int):
            return self.seqs[args]
        elif isinstance(args,tuple) and len(args) == 2:
            return self.seqs[args[0]][args[1]]
        else:
            raise TypeError(f"takes more than 2 positional argument")

    def __repr__(self) -> str:
        return "\n".join(self.seqs) + "\n"
    
    def __str__(self) -> str:
        return "\n".join(self.seqs) + "\n"

    def num_seqs(self):
        return len(self.seqs)

    def column(self,index):
        return [seq[index] for seq in self.seqs]
    
    def __most_frequent(lst):
        
        counter = Counter((char for char in lst if char != "_"))
        if counter:
            return counter.most_common(1)[0][0]
        else:
            raise ValueError("only gap in one position")
    
    def consensus(self) -> str:
        '''
        the most frequent char in each column, ignoring gaps
        '''
        return "".join((AlignResult.__most_frequent(self.column(i))
                        for i in range(len(self))))


class Table:
    """
    class for S,T table
    index start from 0
    """
    def __init__(self, data: list):
        self.data = data
    
    def __getitem__(self, index: tuple):
        """seq1pos, seq2pos, ..."""
        index_str = "".join([f"[{i}]" for i in index[::-1]])
        return eval(f"self.data{index_str}")
    
    def __setitem__(self, index: tuple, value: int):
        index_str = "".join([f"[{i}]" for i in index[::-1]])
        exec(f"self.data{index_str} = {value}")
    
    def __repr__(self) -> str:
        # ↖(U+2196)
        # ← (U+2190)
        # ↑ (U+2191)
        rows = [repr(row) for row in self.data]

        return "\n".join(rows) + "\n"
    
    def to_list(self):
        return self.data

class TraceBack:
    def __init__(self,
                    seqs: list[MySeq],
                    trace_table: Table,
                    pos: list = None,
                    res_seqs: list = None,
                    prior_distance: list[str] = None) -> None:
        self.seqs = seqs
        self.T = trace_table
        self.seq_type = self.seqs[0].seq_type # TODO raise error

        if pos:
            self.pos = pos
        else:
            self.pos = [len(seq) for seq in seqs]

        if prior_distance:
            self.prior_distance = prior_distance
        else:
            self.prior_distance = self.T[self.pos]

        if res_seqs:
            self.res_seqs = res_seqs
        else:
            self.res_seqs = ["" for i in range(len(self.seqs))]

    def __repr__(self) -> str:
        return f"""pos: {self.pos} prior_dist:{self.prior_distance}""" 
    
    def is_done(self) -> bool:
        if self.prior_distance:
            return False
        return True
        # if any(i > 0 for i in self.pos):
        #     return False
        # return True
    
    def get_pos(self) -> tuple:
        return self.pos
    
    def get_align(self) -> AlignResult:
        return AlignResult(self.res_seqs, seq_type=None) #TODO

    def step(self):
        """retrace or split when tie encounted"""

        if len(self.prior_distance) > 1:
            task_list = []
            for pd_i in self.prior_distance:
                task_list.append(TraceBack(self.seqs, self.T,
                                          copy.deepcopy(self.pos), # warning copy 
                                          copy.deepcopy(self.res_seqs), #TODO 
                                          [pd_i,]))

            return task_list

        for i,distance in enumerate(self.prior_distance[0]):
            res_seq = self.res_seqs[i]
            seq = self.seqs[i]

            if distance == "0":
                self.res_seqs[i] = "_" + res_seq
                pass # pos[i] stay

            elif distance == "1":
                self.res_seqs[i] = seq[self.pos[i] - 1] + res_seq

                self.pos[i] -= 1 
            else:
                raise ValueError("distance must be 0 or 1")

        self.prior_distance = self.T[self.pos]

        return [TraceBack(self.seqs, self.T, self.pos, self.res_seqs, self.prior_distance),]
    
    def recover_align(self) -> list[AlignResult]:
        if None:
            raise NotImplementedError

        recover_tasks = [TraceBack(self.seqs, self.T),]

        i = 0
        while i < len(recover_tasks):
            current = recover_tasks[i]
            if not current.is_done():
                recover_tasks.extend(recover_tasks[i].step())
                recover_tasks.pop(i)
            else:
                i += 1

        return list(map(lambda x: x.get_align(), recover_tasks))