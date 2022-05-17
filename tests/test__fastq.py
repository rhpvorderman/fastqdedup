from fastqdedup._fastq import average_error_rate

import pytest


def test_average_error_rate():
    # (0.1 + 0.001) / 2 == 0.0505
    assert average_error_rate(chr(10) + chr(30), phred_offset=0) == 0.0505


def test_average_error_rate_with_default_offset():
    assert average_error_rate(chr(43) + chr(63)) == 0.0505


@pytest.mark.parametrize("i", list(range(33)) + [127])
def test_average_error_rate_out_of_range(i):
    with pytest.raises(ValueError) as error:
        average_error_rate(chr(i))
    error.match(f"{chr(i)} outside of valid phred range")


def test_average_error_rate_non_ascii():
    with pytest.raises(ValueError) as error:
        average_error_rate(chr(128))
    error.match("phred_scores must be ASCII encoded")

