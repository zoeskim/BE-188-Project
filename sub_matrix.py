""" Adapted from code by Adam Weiner: https://github.com/adamcweiner/flu-vax/blob/master/fluv/sub_matrix.py
"""

""" Code for substitution matrices is inspired by https://github.com/benhid/pyMSA/blob/master/pymsa/core/substitution_matrix.py """

class SubstitutionMatrix:
    """ Class representing a substitution matrix, such as PAM250, Blosum62, etc. """

    def __init__(self, gap_penalty: int, gap_character: str):
        self.gap_penalty = gap_penalty
        self.gap_character = gap_character

    def get_distance(self, char1, char2) -> int:
        """ Returns the distance between two characters.
        :param char1: First character.
        :param char2: Second character.
        :return: the distance value """

        if char1 is self.gap_character and char2 is self.gap_character:
            distance = 1
        elif char1 is self.gap_character or char2 is self.gap_character:
            distance = self.gap_penalty
        else:
            matrix = self.get_distance_matrix()
            try:
                distance = matrix[(char1, char2)] if (char1, char2) in matrix else matrix[(char2, char1)]
            except KeyError:
                print("The pair ({0},{1}) couldn't be found in the substitution matrix.".format(char1, char2))
                raise

        return distance

    def get_distance_matrix(self) -> None:
        pass

class PAM250(SubstitutionMatrix):
    """ Class implementing the PAM250 substitution matrix
    Reference: https://en.wikipedia.org/wiki/Point_accepted_mutation"""
    def __init__(self, gap_penalty=-8, gap_character: str = '-'):
        super(PAM250, self).__init__(gap_penalty, gap_character)
        self.distance_matrix = \
            {('W', 'F'): 0, ('L', 'R'): -3, ('S', 'P'): 1, ('V', 'T'): 0, ('Q', 'Q'): 4, ('N', 'A'): 0, ('Z', 'Y'): -4,
             ('W', 'R'): 2, ('Q', 'A'): 0, ('S', 'D'): 0, ('H', 'H'): 6, ('S', 'H'): -1, ('H', 'D'): 1, ('L', 'N'): -3,
             ('W', 'A'): -6, ('Y', 'M'): -2, ('G', 'R'): -3, ('Y', 'I'): -1, ('Y', 'E'): -4, ('B', 'Y'): -3,
             ('Y', 'A'): -3,
             ('V', 'D'): -2, ('B', 'S'): 0, ('Y', 'Y'): 10, ('G', 'N'): 0, ('E', 'C'): -5, ('Y', 'Q'): -4,
             ('Z', 'Z'): 3,
             ('V', 'A'): 0, ('C', 'C'): 12, ('M', 'R'): 0, ('V', 'E'): -2, ('T', 'N'): 0, ('P', 'P'): 6, ('V', 'I'): 4,
             ('V', 'S'): -1, ('Z', 'P'): 0, ('V', 'M'): 2, ('T', 'F'): -3, ('V', 'Q'): -2, ('K', 'K'): 5,
             ('P', 'D'): -1,
             ('I', 'H'): -2, ('I', 'D'): -2, ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2,
             ('P', 'H'): 0,
             ('F', 'Q'): -5, ('Z', 'G'): 0, ('X', 'L'): -1, ('T', 'M'): -1, ('Z', 'C'): -5, ('X', 'H'): -1,
             ('D', 'R'): -1,
             ('B', 'W'): -5, ('X', 'D'): -1, ('Z', 'K'): 0, ('F', 'A'): -3, ('Z', 'W'): -6, ('F', 'E'): -5,
             ('D', 'N'): 2,
             ('B', 'K'): 1, ('X', 'X'): -1, ('F', 'I'): 1, ('B', 'G'): 0, ('X', 'T'): 0, ('F', 'M'): 0, ('B', 'C'): -4,
             ('Z', 'I'): -2, ('Z', 'V'): -2, ('S', 'S'): 2, ('L', 'Q'): -2, ('W', 'E'): -7, ('Q', 'R'): 1,
             ('N', 'N'): 2,
             ('W', 'M'): -4, ('Q', 'C'): -5, ('W', 'I'): -5, ('S', 'C'): 0, ('L', 'A'): -2, ('S', 'G'): 1,
             ('L', 'E'): -3,
             ('W', 'Q'): -5, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 1, ('N', 'R'): 0, ('H', 'C'): -3,
             ('Y', 'N'): -2,
             ('G', 'Q'): -1, ('Y', 'F'): 7, ('C', 'A'): -2, ('V', 'L'): 2, ('G', 'E'): 0, ('G', 'A'): 1, ('K', 'R'): 3,
             ('E', 'D'): 3, ('Y', 'R'): -4, ('M', 'Q'): -1, ('T', 'I'): 0, ('C', 'D'): -5, ('V', 'F'): -1,
             ('T', 'A'): 1,
             ('T', 'P'): 0, ('B', 'P'): -1, ('T', 'E'): 0, ('V', 'N'): -2, ('P', 'G'): 0, ('M', 'A'): -1, ('K', 'H'): 0,
             ('V', 'R'): -2, ('P', 'C'): -3, ('M', 'E'): -2, ('K', 'L'): -3, ('V', 'V'): 4, ('M', 'I'): 2,
             ('T', 'Q'): -1,
             ('I', 'G'): -3, ('P', 'K'): -1, ('M', 'M'): 6, ('K', 'D'): 0, ('I', 'C'): -2, ('Z', 'D'): 3,
             ('F', 'R'): -4,
             ('X', 'K'): -1, ('Q', 'D'): 2, ('X', 'G'): -1, ('Z', 'L'): -3, ('X', 'C'): -3, ('Z', 'H'): 2,
             ('B', 'L'): -3,
             ('B', 'H'): 1, ('F', 'F'): 9, ('X', 'W'): -4, ('B', 'D'): 3, ('D', 'A'): 0, ('S', 'L'): -3, ('X', 'S'): 0,
             ('F', 'N'): -3, ('S', 'R'): 0, ('W', 'D'): -7, ('V', 'Y'): -2, ('W', 'L'): -2, ('H', 'R'): 2,
             ('W', 'H'): -3,
             ('H', 'N'): 2, ('W', 'T'): -5, ('T', 'T'): 3, ('S', 'F'): -3, ('W', 'P'): -6, ('L', 'D'): -4,
             ('B', 'I'): -2,
             ('L', 'H'): -2, ('S', 'N'): 1, ('B', 'T'): 0, ('L', 'L'): 6, ('Y', 'K'): -4, ('E', 'Q'): 2, ('Y', 'G'): -5,
             ('Z', 'S'): 0, ('Y', 'C'): 0, ('G', 'D'): 1, ('B', 'V'): -2, ('E', 'A'): 0, ('Y', 'W'): 0, ('E', 'E'): 4,
             ('Y', 'S'): -3, ('C', 'N'): -4, ('V', 'C'): -2, ('T', 'H'): -1, ('P', 'R'): 0, ('V', 'G'): -1,
             ('T', 'L'): -2,
             ('V', 'K'): -2, ('K', 'Q'): 1, ('R', 'A'): -2, ('I', 'R'): -2, ('T', 'D'): 0, ('P', 'F'): -5,
             ('I', 'N'): -2,
             ('K', 'I'): -2, ('M', 'D'): -3, ('V', 'W'): -6, ('W', 'W'): 17, ('M', 'H'): -2, ('P', 'N'): 0,
             ('K', 'A'): -1,
             ('M', 'L'): 4, ('K', 'E'): 0, ('Z', 'E'): 3, ('X', 'N'): 0, ('Z', 'A'): 0, ('Z', 'M'): -2, ('X', 'F'): -2,
             ('K', 'C'): -5, ('B', 'Q'): 1, ('X', 'B'): -1, ('B', 'M'): -2, ('F', 'C'): -4, ('Z', 'Q'): 3,
             ('X', 'Z'): -1,
             ('F', 'G'): -5, ('B', 'E'): 3, ('X', 'V'): -1, ('F', 'K'): -5, ('B', 'A'): 0, ('X', 'R'): -1,
             ('D', 'D'): 4,
             ('W', 'G'): -7, ('Z', 'F'): -5, ('S', 'Q'): -1, ('W', 'C'): -8, ('W', 'K'): -3, ('H', 'Q'): 3,
             ('L', 'C'): -6,
             ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -4, ('W', 'S'): -2, ('S', 'E'): 0, ('H', 'E'): 1,
             ('S', 'I'): -1,
             ('H', 'A'): -1, ('S', 'M'): -2, ('Y', 'L'): -1, ('Y', 'H'): 0, ('Y', 'D'): -4, ('E', 'R'): -1,
             ('X', 'P'): -1,
             ('G', 'G'): 5, ('G', 'C'): -3, ('E', 'N'): 1, ('Y', 'T'): -3, ('Y', 'P'): -5, ('T', 'K'): 0, ('A', 'A'): 2,
             ('P', 'Q'): 0, ('T', 'C'): -2, ('V', 'H'): -2, ('T', 'G'): 0, ('I', 'Q'): -2, ('Z', 'T'): -1,
             ('C', 'R'): -4,
             ('V', 'P'): -1, ('P', 'E'): -1, ('M', 'C'): -5, ('K', 'N'): 1, ('I', 'I'): 5, ('P', 'A'): 1,
             ('M', 'G'): -3,
             ('T', 'S'): 1, ('I', 'E'): -2, ('P', 'M'): -2, ('M', 'K'): 0, ('I', 'A'): -1, ('P', 'I'): -2,
             ('R', 'R'): 6,
             ('X', 'M'): -1, ('L', 'I'): 2, ('X', 'I'): -1, ('Z', 'B'): 2, ('X', 'E'): -1, ('Z', 'N'): 1, ('X', 'A'): 0,
             ('B', 'R'): -1, ('B', 'N'): 2, ('F', 'D'): -6, ('X', 'Y'): -2, ('Z', 'R'): 0, ('F', 'H'): -2,
             ('B', 'F'): -4,
             ('F', 'L'): 2, ('X', 'Q'): -1, ('B', 'B'): 3}

    def get_distance_matrix(self) -> dict:
        return self.distance_matrix
