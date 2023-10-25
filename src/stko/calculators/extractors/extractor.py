"""
Extractor
=========

#. :class:`.Extractor`

Base results class for extracting molecule properties.

"""

from abc import ABC


class Extractor(ABC):
    """
    An abstract base class for extractors.

    """

    def __init__(self, output_file):
        """
        Initializes :class:`Extractor`

        Parameters
        ----------
        output_file : :class:`str`
            Output file to extract properties from.

        """

        self.output_file = output_file
        # Explictly set encoding to UTF-8 because default encoding on
        # Windows will fail to read the file otherwise.
        with open(self.output_file, 'r', encoding='UTF-8') as f:
            self.output_lines = f.readlines()

        self._extract_values()

    def _extract_values(self):
        """
        Extract values from `path`.

        Parameters
        ----------
        path : :class:`str`
            Path to file to extract values from.

        """

        raise NotImplementedError()

    def _properties_dict(self):
        """
        Dictionary of defined properties and extraction strings.

        """

        raise NotImplementedError()

    def _check_line(self, line, option):
        """
        Checks a line for a string based on option.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to check.

        option : :class:`str`
            Define which property and string being checked for.
            They are defined in :meth:`_properties_dict`.

        Returns
        -------
        :class:`bool`
            Returns ``True`` if the desired string is present.

        """

        options = self._properties_dict()

        if options[option] in line:
            return True
