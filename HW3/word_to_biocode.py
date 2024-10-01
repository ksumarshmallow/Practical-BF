import string
from dataclasses import dataclass
from Bio.Seq import Seq

@dataclass
class WordToBiocode:
    lang: str = 'eng'
    nucl: tuple = ('A', 'T', 'G', 'C')

    def __post_init__(self):
        if self.lang not in ('eng', 'rus'):
            raise ValueError('lang must be either "eng" or "rus"')

        if self.lang == 'eng':
            self.alphabet = string.ascii_lowercase
        elif self.lang == 'rus':
            self.alphabet = 'абвгдеёжзийклмнопрстуфхцчшщъыьэюя'
        
        self.alphabet_to_binary = {letter: bin(i)[2:].zfill(6) for i, letter in enumerate(self.alphabet)}

        self.fourth_to_binary = {str(idx): binary for idx, binary in enumerate(['00', '01', '10', '11'])}
        self.binary_to_fourth = {binary: idx for idx, binary in self.fourth_to_binary.items()}

        self.fourth_to_nucl = {str(idx): nucl for idx, nucl in enumerate(self.nucl)}
        self.nucl_to_fourth = {nucl: idx for idx, nucl in self.fourth_to_nucl.items()}

    def _convert_to_binary(self, word) -> list:
        binary_word = [self.alphabet_to_binary[letter] for letter in word]
        return binary_word
    
    def _convert_to_fourth(self, binary_word) -> list:
        fourth_word = []
        for i in range(len(binary_word)):
            letter_binary = binary_word[i]
            letter_binary_pairs = [letter_binary[j:j+2] for j in range(0, len(letter_binary), 2)]
            letter_fourth = ''.join([self.binary_to_fourth[j] for j in letter_binary_pairs])
            fourth_word.append(letter_fourth)
        return fourth_word

    def _convert_to_seq(self, fourth_word) -> str:
        seq = ''.join([self.fourth_to_nucl[fourth] for fourth in ''.join(fourth_word)])
        return seq

    def get_frequency(self, word: str, genome_length: float = 3.2 * 10**9):
        prob_seq = 1 / 4**(3 * len(word))
        return genome_length * prob_seq

    def get_probability(self, word: str, genome_length: float = 3.2 * 10**9):
        prob_seq = 1 / 4**(3 * len(word))
        return 1 - (1 - prob_seq)**(genome_length)

    def convert_to_nucleotide(self, word):
        word = word.lower()

        if set(word) - set(self.alphabet):
            raise ValueError('word contains invalid characters')

        binary_word = self._convert_to_binary(word)
        fourth_word = self._convert_to_fourth(binary_word)
        seq = self._convert_to_seq(fourth_word)
        return seq

    def convert_to_protein(self, word):
        sequence = Seq(self.convert_to_nucleotide(word))    
        protein_sequence = "".join(sequence.translate())
        return protein_sequence