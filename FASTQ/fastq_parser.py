import sys

PHRED_CONST = 33
WINDOW_SIZE = 5
QUALITY_BORDER = 30
MIN_LENGTH = 60

def phred_score(char):
    """
    Convert FASTQ quality character to Phred score
    """
    return ord(char) - PHRED_CONST


def read_fastq(filename):
    """
    Read FASTQ file and return list of (sequence, quality)
    """
    reads = []

    with open(filename) as f:
        while True:
            header = f.readline()
            if not header:
                break

            sequence = f.readline().strip()
            f.readline()
            quality = f.readline().strip()

            if not sequence or not quality:
                break

            reads.append((sequence, quality))

    return reads


def gc_content(reads):
    """
    Calculates GC %
    """
    gc = 0
    total = 0

    for sequence, _ in reads:
        gc += sequence.count("G")
        gc += sequence.count("C")
        total += len(sequence)

    if total == 0:
        return 0

    return round(100 * gc / total, 2)


def phred_position_average(reads, pos):
    """
    Average Phred score for a specific position
    """
    scores = []

    for _, quality in reads:
        if len(quality) >= pos:
            scores.append(phred_score(quality[pos - 1]))

    if not scores:
        return 0

    return round(sum(scores) / len(scores))


def sliding_window_trim(sequence, quality):
    """
    Sliding window trimming
    """
    length = len(sequence)

    for i in range(length - WINDOW_SIZE + 1):
        window = quality[i:i + WINDOW_SIZE]

        avg_quality = sum(phred_score(c) for c in window) / WINDOW_SIZE

        if avg_quality < QUALITY_BORDER:
            return sequence[:i], quality[:i]

    return sequence, quality


def quality_trimming(reads):
    """
    Apply sliding window trimming
    """
    trimmed_reads = []
    lengths = []
    removed_reads = 0

    for sequence, quality in reads:

        new_seq, new_qual = sliding_window_trim(sequence, quality)

        if len(new_seq) == 0:
            removed_reads += 1
            continue

        trimmed_reads.append((new_seq, new_qual))
        lengths.append(len(new_seq))

    return trimmed_reads, lengths, removed_reads


def length_filter(reads):
    """
    Filter reads shorter than MIN_LENGTH
    """
    result = []

    for sequence, quality in reads:
        if len(sequence) >= MIN_LENGTH:
            result.append((sequence, quality))

    return result


def length_stats(lengths):
    """
    Return min, average and max length
    """
    if not lengths:
        return 0, 0, 0

    min_val = min(lengths)
    max_val = max(lengths)
    aver = round(sum(lengths) / len(lengths))

    return min_val, aver, max_val


def main():
    file_name = sys.argv[1]
    reads = read_fastq(file_name)

    lengths_of_seqs = [len(seq) for seq, _ in reads]

    total_reads = len(reads)
    min_len, avg_len, max_len = length_stats(lengths_of_seqs)

    print("Общее число прочтений -", total_reads)
    print("Минимальная длина -", min_len)
    print("Средняя длина -", avg_len)
    print("Максимальная длина -", max_len)

    gc = gc_content(reads)
    print("GC-состав -", gc)

    phred_avg = phred_position_average(reads, 10)
    print("Среднее качество по шкале Phred -", phred_avg)

    trimmed_reads, trimmed_lengths, removed_reads = quality_trimming(reads)
    print("Сколько прочтений подверглось триммингу -", removed_reads)
    
    min_trim, avg_trim, max_trim = length_stats(trimmed_lengths)
    print("Минимальная длина -", min_trim)
    print("Средняя длина -", avg_trim)
    print("Максимальная длина -", max_trim)


    final_reads = length_filter(trimmed_reads)
    print("Оставшееся число прочтений -", len(final_reads))


if __name__ == "__main__":
    main()
