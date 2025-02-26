import math


def round_up(x):
    if math.log(x) >= 1:
        num_digits = len(str(x))
        num_zeros = num_digits-1
        multiply_by =10**num_zeros
        return int(math.ceil(x/multiply_by)*multiply_by)

    else:
        raise Exception