import sys


class Base:
    minimal_orf_lenght = 120
    translation_dict = {'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'UUA': 'L', 'UUG': 'L',
                        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'CGU': 'R', 'CGC': 'R',
                        'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 'AAA': 'K', 'AAG': 'K',
                        'AAU': 'N', 'AAC': 'N', 'AUG': 'M', 'GAU': 'D', 'GAC': 'D', 'UUU': 'F',
                        'UUC': 'F', 'UGU': 'C', 'UGC': 'C', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P',
                        'CCG': 'P', 'CAA': 'Q', 'CAG': 'Q', 'UCU': 'S', 'UCC': 'S', 'UCA': 'S',
                        'UCG': 'S', 'AGU': 'S', 'AGC': 'S', 'GAA': 'E', 'GAG': 'E', 'ACU': 'T',
                        'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G',
                        'GGG': 'G', 'UGG': 'W', 'CAU': 'H', 'CAC': 'H', 'UAU': 'Y', 'UAC': 'Y',
                        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V',
                        'GUG': 'V', 'UAG': '*', 'UGA': '*', 'UAA': '*'}
    weights = {}
    water_weight = 18

    def __init__(self, s):
        self.s = s.upper()

    def __str__(self):
        return self.s

    def as_dna(self):
        raise NotImplementedError()

    def as_rna(self):
        raise NotImplementedError()

    def as_protein(self):
        raise NotImplementedError()

    def weight(self):
        return sum([self.weights[m] - self.water_weight for m in str(self)]) + self.water_weight

    @staticmethod
    def check_sequence(s):
        raise NotImplementedError()

    @staticmethod
    def find_longest_orf(s):
        raise NotImplementedError()


class DNA(Base):
    weights = {'A': 331.2218, 'T': 322.2085, 'G': 347.2212, 'C': 307.1971}

    def as_dna(self):
        return self


class RNA(Base):
    back_transcription_dict = {'U': 'A', 'A': 'T', 'G': 'C', 'C': 'G'}
    start_codon = 'AUG'
    stop_codons = ['UAG', 'UGA', 'UAA']
    weights = {'A': 347.2212, 'G': 363.2206, 'C': 323.1965, 'U': 324.1813}

    def as_dna(self):
        return CodingDNA(''.join([self.back_transcription_dict[letter] for letter in self.s]))

    def as_rna(self):
        return self

    def as_protein(self):
        s = str(self.as_rna())
        protein = str()
        for n in range(0, len(s), 3):
            protein += self.translation_dict[s[n:n + 3]]
        return Protein(protein)

    @staticmethod
    def check_sequence(s):
        if 'U' in s:
            return True
        else:
            return False

    @staticmethod
    def find_longest_orf(s):
        indexes = None
        max_len = RNA.minimal_orf_lenght
        for i in range(0, len(s) - 3):
            if s[i:i + 3] == RNA.start_codon:
                for j in range(i + 3, len(s), 3):
                    if s[j:j + 3] in RNA.stop_codons:
                        if len(s[i:j + 3]) > max_len:
                            max_len = len(s[i:j + 3])
                            indexes = [i, j + 3]
                        break
        if indexes:
            return s[indexes[0]:indexes[1]]
        else:
            return None


class CodingDNA(DNA):
    transcription_dict = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}
    start_codon = 'TAC'
    stop_codons = ['ATC', 'ACT', 'ATT']

    def as_rna(self):
        return RNA(''.join([self.transcription_dict[letter] for letter in self.s]))

    def as_protein(self):
        s = str(self.as_rna())
        protein = ''
        for n in range(0, len(s), 3):
            protein += self.translation_dict[s[n:n + 3]]
        return Protein(protein)

    @staticmethod
    def check_sequence(s):
        if 'U' not in s:
            return True
        else:
            return False

    @staticmethod
    def find_longest_orf(s):
        indexes = None
        max_len = CodingDNA.minimal_orf_lenght
        for i in range(0, len(s) - 3):
            if s[i:i + 3] == CodingDNA.start_codon:
                for j in range(i + 3, len(s), 3):
                    if s[j:j + 3] in CodingDNA.stop_codons:
                        if len(s[i:j + 3]) > max_len:
                            max_len = len(s[i:j + 3])
                            indexes = [i, j + 3]
                        break
        if indexes:
            return s[indexes[0]:indexes[1]]
        else:
            return None


class NoncodingDNA(DNA):
    transcription_dict = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}
    start_codon = 'ATG'
    stop_codons = ['TAG', 'TGA', 'TAA']

    def as_rna(self):
        return RNA(self.s.replace('T', 'U'))

    def as_protein(self):
        s = str(self.as_rna())
        protein = ''
        for n in range(0, len(s), 3):
            protein += self.translation_dict[s[n:n + 3]]
        return Protein(protein)

    @staticmethod
    def check_sequence(s):
        if 'U' not in s:
            return True
        else:
            return False

    @staticmethod
    def find_longest_orf(s):
        indexes = None
        max_len = NoncodingDNA.minimal_orf_lenght
        for i in range(0, len(s) - 3):
            if s[i:i + 3] == NoncodingDNA.start_codon:
                for j in range(i + 3, len(s), 3):
                    if s[j:j + 3] in NoncodingDNA.stop_codons:
                        if len(s[i:j + 3]) > max_len:
                            max_len = len(s[i:j + 3])
                            indexes = [i, j + 3]
                        break
        if indexes:
            return s[indexes[0]:indexes[1]]
        else:
            return None


class Protein(Base):
    weights = {'A': 89.0932, 'L': 131.1729, 'R': 174.201, 'K': 146.1876, 'N': 132.1179, 'M': 149.2113,
               'D': 133.1027, 'F': 165.1891, 'C': 121.1582, 'P': 115.1305, 'Q': 146.1445, 'S': 105.0926,
               'E': 147.1293, 'T': 119.1192, 'G': 75.0666, 'W': 204.2252, 'H': 155.1546, 'Y': 181.1885,
               'I': 131.1729, 'V': 117.1463, '*': 0}

    def as_protein(self):
        return self


def create_orf_obj(sequence):
    possible = (NoncodingDNA, CodingDNA, RNA)
    for t in possible:
        ORF = t.find_longest_orf(sequence)
        if t.check_sequence(sequence) and ORF:
            return t(ORF)


data = sys.stdin.read()
longest_ORFs = [create_orf_obj(''.join(seq.split('\n')[1:])) for seq in data.split('>')[1:]]

for orf in longest_ORFs:
    dna = orf.as_dna()
    rna = orf.as_rna()
    protein = orf.as_protein()
    print('DNA Sequence: {:>5}  Weight: {:>5}\n'.format(str(dna), dna.weight()))
    print('RNA Sequence: {:>5}  Weight: {:>5}\n'.format(str(rna), rna.weight()))
    print('Protein Sequence: {:>5}  Weight: {:>5}\n'.format(str(protein), protein.weight()))

# with open('C:/Users/Andrey/PycharmProjects/lecture1/lecture2/lec2_hw_sample1.fsa', 'r') as inp:
