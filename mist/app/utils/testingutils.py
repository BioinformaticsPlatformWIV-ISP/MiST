import tempfile


def get_temp_dir() -> tempfile.TemporaryDirectory:
    """
    Returns a temporary directory for testing.
    :return: Temporary directory
    """
    return tempfile.TemporaryDirectory()
