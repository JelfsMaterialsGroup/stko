def check_line(
    self,
    line: str,
    option: str,
    options_dict: dict[str, str],
) -> bool:
    """
    Checks a line for a string based on option.

    Parameters:

        line:
            Line of output file to check.

        option:
            Define which property and string being checked for.
            They are defined in :meth:`_properties_dict` of the
            Extractor.

        options_dict:
            The :meth:`_properties_dict` of the Extractor.

    Returns:

        Returns ``True`` if the desired string is present.

    """

    options = self._properties_dict()

    if options[option] in line:
        return True

    return False
