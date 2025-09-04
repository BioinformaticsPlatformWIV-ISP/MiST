from collections.abc import Sequence


def unique_preserve_order(seq: Sequence) -> list:
    """
    Only keeps the unique elements from the input list while preserving the order.
    :param seq: Input sequence
    :return: List of unique elements
    """
    seen = set()
    result = []
    for item in seq:
        if item in seen:
            continue
        seen.add(item)
        result.append(item)
    return result
