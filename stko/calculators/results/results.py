"""
Results
=======

#. :class:`.Results`

Base results class for extracting molecule properties.

"""

from abc import ABC


class Results(ABC):
    """
    An abstract base class for results.
    ...

    """

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'<{self.__class__.__name__} at {id(self)}>'
