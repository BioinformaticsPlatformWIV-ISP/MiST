from mist.app import model

# Real allele ids are always non-negative, so these can't collide with an actual allele.
CODE_ABSENT = 0
CODE_WILDCARD = -1
CODE_UNREADABLE = -2  # missing/non-numeric allele value -> never matches


def encode_allele(raw: str) -> int:
    """
    Encodes an allele string into a compact int code.
    :param raw: Allele value
    :return: Encoded allele code
    """
    if raw == model.ALLELE_WILDCARD:
        return CODE_WILDCARD
    try:
        return int(raw)
    except ValueError:
        return CODE_UNREADABLE


def decode_allele(code: int) -> str:
    """
    Decodes an int allele code back into its original string representation.
    :param code: Encoded allele
    :return: Allele string
    """
    if code == CODE_WILDCARD:
        return model.ALLELE_WILDCARD
    if code == CODE_UNREADABLE:
        return 'n/a'
    return str(code)


def candidate_codes(res: model.QueryResult | None) -> set[int]:
    """
    Determines which encoded allele codes a detected result could match at a single locus. Candidates that
    aren't plain integers (e.g. a novel-allele marker like "12*") can never match a stored profile allele, so
    they're dropped.
    :param res: Detected result for the locus (None if nothing was detected)
    :return: Candidate allele codes
    """
    if res is None:
        return {CODE_ABSENT}
    candidates = [res.allele_str] if len(res.allele_results) == 1 else res.allele_str.split('__')
    codes = set()
    for candidate in candidates:
        try:
            codes.add(int(candidate))
        except ValueError:
            continue
    return codes
