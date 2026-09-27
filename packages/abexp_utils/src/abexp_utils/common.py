def url_encode(raw: str):
    """
    Returns a url-encoded version of a raw string
    """
    from urllib.parse import quote

    return quote(raw.strip())


def url_decode(raw: str):
    """
    Returns a url-decoded version of a raw string
    """
    from urllib.parse import quote, unquote

    return unquote(raw.strip())
