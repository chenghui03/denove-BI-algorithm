""" Chapter 4"""

import setwd

import re
import functools

from modules.constants import * 


class MySeq:
    def __init__(self,seq:str,seq_type:str="DNA") -> None:
        self.seq = seq.upper()
        self.seq_type = seq_type

    def __len__(self) -> int:
        return len(self.seq)
    
    def __getitem__(self,index) -> str:
        return self.seq[index]
    
    def __getslice__(self,i,j) -> str:
        return self.seq[i:j]
    
    def __str__(self) -> str:
        return self.seq
    
    def __repr__(self) -> str:
        return f"Sequence: {self.seq} type: {self.seq_type}"
    
    def get_biotype(self) -> str:
        return self.seq_type
    
    def get_seq(self) -> str:
        return self.seq
    
    def alphabet(self) -> str|None:
        return ALPHABET.get(self.seq_type,None)
    
    def validate(self) -> bool:
        '''raise error when can't pass validation'''
        match_iter = re.finditer(
            ERROR_ALPHABET_PATTERN[self.seq_type],self.seq)
        
        for error in match_iter:
            raise ValueError(
                "{} at {} index is unsupported in alphabet".format(
                    error.group()[0],error.span()[0]))
        return True
    
    @staticmethod
    def check_type(supported_type: list[str]):
        def decorator(func):
            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                self = args[0]
                if self.seq_type not in supported_type:
                    raise ValueError(f"The sequence type '{self.seq_type}' is not supported. Expected {','.join(supported_type)} sequence.")
                return func(*args, **kwargs)
            return wrapper
        return decorator
        
    @check_type(["DNA"])
    def transcribe(self) -> "MySeq":
        """without reverse"""
        return MySeq(self.seq.replace("T", "U"), "RNA")
    
    @check_type(["DNA"])
    def reverse_comp(self) -> "MySeq":
        comp_table = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd',
                                    'TGCAtgcaYRKMyrkmBVDHbvdh')
        seq = self.seq.translate(comp_table)
        return MySeq(seq[::-1],"DNA")

    def __translate_codon(self, codon) -> str:
        return CODEN.get(codon, 'X')
    
    @check_type(["RNA"])
    def retro_transcribe(self) -> "MySeq":
        """without reverse"""
        return MySeq(self.seq.replace("U", "T"), "DNA")

    @check_type(["DNA", "RNA"])
    def translate(self, ini_pos:int=0) -> "MySeq":
        if self.seq_type == "DNA" :
            seq = ''.join((self.__translate_codon(self.seq[pos:pos+3]) \
                            for pos in range(ini_pos, len(self.seq)-2, 3)))
            return MySeq(seq,"protein")
        
        elif self.seq_type == "RNA" :
            return self.retro_transcribe().translate()
        
