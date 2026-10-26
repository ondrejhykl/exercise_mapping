from Bio import SeqIO


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

    def locate(self, p):
        locations = []
        top, bottom = self.count(p)
        for start_row in range(top, bottom):
            row = start_row
            step = 0
            while row not in self.sampled:
                row = self.lf(row)
                step += 1

            position = self.sampled[row] + step
            locations.append(position)

        return locations


def reverse_complement(sequence):
    complement = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }

    return "".join(complement[base] for base in sequence)[::-1]


def hamming(a, b):
    mismatches = 0
    for i in range(len(a)):
        if a[i] != b[i]:
            mismatches += 1
    return mismatches


def seeds(read, n):
    result = []

    for i in range(n):
        start = i * len(read) // n
        end = (i + 1) * len(read) // n
        piece = read[start:end]
        result.append((start, piece))

    return result


def map_read(read, fm, ref, d=2, n_seeds=None):
    if n_seeds is None:
        n_seeds = d + 1

    hits = set()
    total_candidates = 0

    for strand, sequence in [
        ("+", read),
        ("-", reverse_complement(read))
    ]:
        candidates = set()

        for offset, seed in seeds(sequence, n_seeds):
            for pos in fm.locate(seed):
                start = pos - offset

                # Read must fit completely inside the reference
                if start >= 0 and start + len(sequence) <= len(ref):
                    candidates.add(start)

        total_candidates += len(candidates)

        for start in candidates:
            mismatches = hamming(
                ref[start:start + len(sequence)],
                sequence
            )

            if mismatches <= d:
                hits.add((strand, start, mismatches))

    return sorted(hits), total_candidates


def brute_map(read, ref, d=2):
    hits = []

    for strand, sequence in [
        ("+", read),
        ("-", reverse_complement(read))
    ]:
        for start in range(len(ref) - len(sequence) + 1):
            mismatches = hamming(
                ref[start:start + len(sequence)],
                sequence
            )

            if mismatches <= d:
                hits.append((strand, start, mismatches))

    return sorted(hits)


if __name__ == "__main__":
    """
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
    print(fmi2.locate("GGC"))
    """
    ref_record = SeqIO.read("reference.fasta", "fasta")
    ref = str(ref_record.seq)

    reads = [
        (record.id, str(record.seq))
        for record in SeqIO.parse("reads.fasta", "fasta")
    ]

    fm = FMIndex(ref + "$", 32)
    d = 2

    results = []

    for name, read in reads:
        hits, candidates = map_read(read, fm, ref, d)
        results.append((name, read, hits, candidates))

    mapped = sum(bool(hits) for _, _, hits, _ in results)
    unmapped = len(results) - mapped
    multiple = sum(len(hits) > 1 for _, _, hits, _ in results)

    print("Total reads:", len(reads))
    print("Mapped:", mapped)
    print("Unmapped:", unmapped)
    print("More than one hit:", multiple)

    # Print hits of the first five reads
    print("\nHits of the first five reads:")

    for name, read, hits, _ in results[:5]:
        if hits:
            for strand, position, mismatches in hits:
                print(name, strand, position, mismatches)
        else:
            print(name, "unmapped")

    # Self-check 1: validate every hit
    valid_hits = True

    for name, read, hits, _ in results:
        for strand, position, mismatches in hits:
            sequence = read if strand == "+" else reverse_complement(read)

            observed = hamming(
                sequence,
                ref[position:position + len(read)]
            )

            if not (0 <= position <= len(ref) - len(read)):
                valid_hits = False
            if observed != mismatches or mismatches > d:
                valid_hits = False

    print("\nSelf-check 1 - Hit positions and mismatches:",
              "OK" if valid_hits else "ERROR")

    # Self-check 2: compare 30 random reads with brute force
    import random

    random.seed(42)
    sample = random.sample(results, min(30, len(results)))
    matching_reads = 0

    for name, read, hits, _ in sample:
        expected = brute_map(read, ref, d)
        if hits == expected:
            matching_reads += 1

    print("Self-check 2 - FM-index vs brute force:",
          f"{matching_reads}/{len(sample)} reads match")

    # Self-check 3: check all unmapped reads
    unmapped_reads = [
        (name, read)
        for name, read, hits, _ in results
        if not hits
    ]

    false_unmapped = 0

    for name, read in unmapped_reads:
        if brute_map(read, ref, d):
            false_unmapped += 1

    print("Self-check 3 - Unmapped reads with brute-force hits:",
          false_unmapped)

    # Self-check 4: reference ends
    first_hits, _ = map_read(ref[:100], fm, ref, d=0)
    first_ok = ("+", 0, 0) in first_hits

    last_read = reverse_complement(ref[-100:])
    last_hits, _ = map_read(last_read, fm, ref, d=0)
    last_ok = ("-", len(ref) - 100, 0) in last_hits

    print("Self-check 4 - Reference ends:")
    print("  First 100 bp:", "OK" if first_ok else "ERROR")
    print("  Reverse complement of last 100 bp:",
          "OK" if last_ok else "ERROR")
