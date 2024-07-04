def check_line(line: str, option: str, options_dict: dict[str, str]) -> bool:
    """Checks a line for a string based on option.

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
    return options_dict[option] in line
