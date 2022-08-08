from scipy.io import loadmat
import numpy as np


def load_model(file: str, squeeze_me: bool = True):
    """Load discdata in .mat format into structured array.
    Args:
        file (str): Path to .mat data file.
        squeeze_me (bool, optional): Remove length 1 dimension in array.
    Returns:
        MD_dict (dict): .mat data file with type dictionary.

    Sources: https://numpy.org/doc/stable/user/basics.rec.html

    """

    matdata = loadmat(
        file, squeeze_me=squeeze_me
    )  # raw data, squeeze_me removes length 1 dimensions

    MD = matdata["MD"]
    MD_dict = strucArr2dict(MD, keys=get_keys(MD))

    return MD_dict


def get_keys(strucArr):
    """Get names of dtypes in structured array.

    Args:
        strucArr (array): Structured array with named dtype.

    Returns:
        (tuple): Names of dtype in strucArr.

    Example:
    >>> x = np.array([('Rex', 9, 81.0), ('Fido', 3, 27.0)], dtype=[('name', 'U10'), ('age', 'i4'), ('weight', 'f4')])
    >>> get_keys(x)
    ('name', 'age', 'weight')
    """
    return strucArr.dtype.names


def strucArr2dict(strucArr, keys: list):
    """Convert structured array to dictionary.

    Args:
        strucArr (array): Structured array containing magnetodisc data.
        keys (list): List of field names for the magnetodisc structured array.

    Returns:
        dataDict (dict): Dictionary containing data of the specified fields.
    """
    dataDict = {}
    for key in keys:
        subkeys = strucArr[key].item().dtype.names

        if subkeys:
            dataDict.update({key: dict.fromkeys(subkeys)})
            for subkey in subkeys:
                dataDict[key][subkey] = strucArr[key].item()[subkey].item()
        else:
            dataDict.update({key: strucArr[key].item()})

    return dataDict
