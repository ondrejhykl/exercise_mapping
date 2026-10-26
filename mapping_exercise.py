def brute_force(text, p):
    indexes = []
    for i in range(len(text) - len(p) + 1):
        pat = text[i:i+len(p)]
        if pat == p:
            indexes.append(i)
    return indexes


def suffix_array(text):
    suffixes = [(text, 0)]
    for i in range(1, len(text)):
        suffixes.append((text[i:], i))
    suffixes = sorted(suffixes)
    sa = [index for _, index in suffixes]
    return sa


def find_block(text, sa, p):
    # First row whose suffix prefix is >= p
    lo = 0
    hi = len(sa)
    while lo < hi:
        mid = (lo + hi) // 2
        current = text[sa[mid]:sa[mid] + len(p)]

        if current < p:
            lo = mid + 1
        else:
            hi = mid
    start = lo

    lo = start
    hi = len(sa)
    while lo < hi:
        mid = (lo + hi) // 2
        current = text[sa[mid]:sa[mid] + len(p)]

        if current <= p:
            lo = mid + 1
        else:
            hi = mid
    end = lo

    return (start, end)


def bwt(text, sa):
    bwt = ""
    for i in sa:
        bwt = bwt + text[i - 1]
    return bwt


def first_occurrence(letters):
    result = {}
    position = 0

    for letter in sorted(set(letters)):
        result[letter] = position
        position += letters.count(letter)

    return result


def inverse_bwt(bwt):
    text = ""
    row = 0
    fo = first_occurrence(bwt)
    for i in range(len(bwt) - 1):
        k = bwt[:row].count(bwt[row]) + 1
        text = text + bwt[row]
        row = fo[bwt[row]] + k - 1
    return text[::-1] + "$"


class FMIndex:
    def __init__(self, text, s):
        self.bwt = bwt(text, suffix_array(text))
        self.first_occurrence = first_occurrence(self.bwt)
        self.counts = {}
        for X in self.bwt:
            if X not in self.counts:
                self.counts[X] = []
        for i in range(len(text)+1):
            for X in self.counts.keys():
                self.counts[X].append(self.bwt[:i].count(X))
        self.sampled = {}
        for row, pos in enumerate(suffix_array(text)):
            if pos % s == 0:
                self.sampled[row] = pos
        self.s = s

    def lf(self, row):
        X = self.bwt[row]
        return self.first_occurrence[X] + self.counts[X][row]

    def count(self, p):
        top = 0
        bottom = len(self.bwt)
        for letter in p[::-1]:
            ct = self.counts[letter][top]
            cb = self.counts[letter][bottom]
            top = self.first_occurrence[letter] + ct
            bottom = self.first_occurrence[letter] + cb

        return top, bottom


if __name__ == "__main__":
    sequence = "AAGGCGGCTTCA$"
    pattern = "GGC"
    print(brute_force(sequence, pattern))
    sa = suffix_array(sequence)
    print(sa)
    print(find_block(sequence, sa, pattern))
    transform = bwt(sequence, sa)
    print(transform)
    fo = first_occurrence(transform)
    print(fo)
    inv_bwt = inverse_bwt(transform)
    print(inv_bwt)
    fmi = FMIndex("CAT$", 2)
    fmi2 = FMIndex(sequence, 2)
    print(fmi.sampled, fmi.bwt, fmi.counts["T"])
    print(fmi.lf(1))
    print(fmi.count("AT"))
    print(fmi2.count("GGC"))
