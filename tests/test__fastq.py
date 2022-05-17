from fastqdedup._fastq import average_error_rate


def test_average_error_rate():
    # (0.1 + 0.001) / 2 == 0.0505
    assert average_error_rate(chr(10) + chr(30), phred_offset=0) == 0.0505