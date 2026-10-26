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


if __name__ == "__main__":
    sequence = "AAGGCGGCTTCA$"
    pattern = "GGC"
    print(brute_force(sequence, pattern))
    sa = suffix_array(sequence)
    print(sa)
    print(find_block(sequence, sa, pattern))

