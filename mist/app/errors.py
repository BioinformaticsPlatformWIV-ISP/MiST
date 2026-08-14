class DependencyError(RuntimeError):
    """
    Error that is raised when a dependency is missing.
    """

    pass

class LinCodeError(Exception):
    """
    Raised when a LIN code cannot be extracted for the given input.
    """

    pass
