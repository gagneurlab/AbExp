import sys


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# from snakemk_util import recursive_format
def recursive_format(data, params, fail_on_unknown=False):
    """
    format a (nested) dictionary of strings with a set of params
    """
    if isinstance(data, str):
        try:
            return data.format_map(params)
        except ValueError as e:
            eprint(f"Failed to format '{data}' with params '{params}'!")
            raise e
    elif isinstance(data, dict):
        return {k: recursive_format(v, params) for k, v in data.items()}
    elif isinstance(data, list):
        return [recursive_format(v, params) for v in data]
    else:
        if fail_on_unknown:
            raise ValueError("Handling of data type not implemented: %s" % type(data))
        else:
            return data


class SafeDict(dict):
    def __missing__(self, key):
        return '{' + key + '}'
