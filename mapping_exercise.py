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


if __name__ == "__main__":
    sequence = "AAGGCGGCTTCA$"
    pattern = "GGC"
    print(brute_force(sequence, pattern))
    sa = suffix_array(sequence)
    print(sa)
    print(find_block(sequence, sa, pattern))
    bwt = bwt(sequence, sa)
    print(bwt)
    fo = first_occurrence(bwt)
    print(fo)
    inv_bwt = inverse_bwt(bwt)
    print(inv_bwt)


