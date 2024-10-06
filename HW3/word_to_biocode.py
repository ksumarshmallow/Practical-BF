import string
from dataclasses import dataclass
from Bio.Seq import Seq
from Bio.Data import CodonTable

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
    
    def get_frequency(self, word: str, genome_length: float = 3.2 * 10**9):
        prob_seq = 1 / 4**(3 * len(word))
        return genome_length * prob_seq

    def get_probability(self, word: str, genome_length: float = 3.2 * 10**9):
        prob_seq = 1 / 4**(3 * len(word))
        return 1 - (1 - prob_seq)**(genome_length)

    def _convert_to_binary(self, word) -> list:
        binary_word = [self.alphabet_to_binary[letter] for letter in word]
        return binary_word

    def _convert_from_binary(self, binary_word: list) -> str:
        decoded_letters = []
        for i in range(0, len(binary_word), 6):
            binary_letter = binary_word[i:i+6]
            letter_index = int(binary_letter, 2)
            if letter_index >= len(self.alphabet):
                decoded_letter = ''
            else:
                decoded_letter = self.alphabet[letter_index]
            decoded_letters.append(decoded_letter)
        return ''.join(decoded_letters)
    
    def _convert_to_fourth(self, binary_word) -> list:
        fourth_word = []
        for i in range(len(binary_word)):
            letter_binary = binary_word[i]
            letter_binary_pairs = [letter_binary[j:j+2] for j in range(0, len(letter_binary), 2)]
            letter_fourth = ''.join([self.binary_to_fourth[j] for j in letter_binary_pairs])
            fourth_word.append(letter_fourth)
        return fourth_word

    def _convert_from_fourth(self, fourth_word: list) -> list:
        binary_word = ''.join([self.fourth_to_binary[fourth] for fourth in fourth_word])
        return binary_word
    
    def _convert_to_seq(self, fourth_word) -> str:
        seq = ''.join([self.fourth_to_nucl[fourth] for fourth in ''.join(fourth_word)])
        return seq
    
    def _convert_from_seq(self, seq: str) -> list:
        fourth_seq = ''.join([self.nucl_to_fourth[nucl] for nucl in seq])
        return fourth_seq

    def reverse_translation(self, protein_seq: str) -> str:
        standard_table = CodonTable.unambiguous_dna_by_id[1]

        nucleotide_seq = []
        for aa in protein_seq:
            if aa in standard_table.back_table:
                codon = standard_table.back_table[aa][0]  # первый кодон для АК
            else:
                codon = 'N'
            nucleotide_seq.append(codon)
        nucleotide_seq = ''.join(nucleotide_seq)
        return nucleotide_seq

    def encode(self, word: str, type_seq='DNA'):
        if type_seq not in ('DNA', 'PROTEIN'):
            raise ValueError('type_seq must be either "DNA" or "PROTEIN"')

        word = word.lower()
        if set(word) - set(self.alphabet):
            raise ValueError('word contains invalid characters')

        binary_word = self._convert_to_binary(word)
        fourth_word = self._convert_to_fourth(binary_word)
        seq = self._convert_to_seq(fourth_word)

        if type_seq == 'PROTEIN':
            seq = Seq(seq)    
            seq = "".join(seq.translate())
        
        return seq

    def decode(self, seq: str, type_seq='DNA'):
        if type_seq not in ('DNA', 'PROTEIN'):
            raise ValueError('type_seq must be either "DNA" or "PROTEIN"')

        seq = seq.upper()
        if len(seq) % 3 != 0:
            seq = seq[:-(len(seq) % 3 + 1)]

        if type_seq == 'PROTEIN':
            seq = self.reverse_translation(protein_seq=seq)
        
        fourth_seq = self._convert_from_seq(seq)
        binary_word = self._convert_from_fourth(fourth_seq)
        decoded_word = self._convert_from_binary(binary_word)

        return decoded_word