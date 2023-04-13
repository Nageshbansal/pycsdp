"""API defining checkpoint."""
import os
import sys
_ROOT = os.path.abspath(os.path.dirname(__file__))


def get_pca(file):
    """Return the path of data file from data directory."""
    return os.path.normpath(os.path.join(_ROOT, 'pca', file))


if sys.version_info[:2] >= (3, 8):
    from importlib.metadata import PackageNotFoundError, version  # pragma: no cover
else:
    from importlib_metadata import PackageNotFoundError, version  # pragma: no cover

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = version(dist_name)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError
