import shlex
import subprocess
import sys

from vtam.utils.PathManager import PathManager


def pip_install_vtam_for_tests():

    """This function is used in the tests when the vtam command is run"""

    cmd = '{} -m pip install . -q --upgrade'.format(sys.executable)
    if sys.platform.startswith("win"):
        args = cmd
    else:
        args = shlex.split(cmd)
    subprocess.run(args=args, check=True, cwd=PathManager.instance().get_project_path())

def tqdm_hook(t):
    """Wraps tqdm instance.
    Don't forget to close() or __exit__()
    the tqdm instance once you're done with it (easiest using `with` syntax).
    Example
    -------
    >>> with tqdm(...) as t:
    ...     reporthook = my_hook(t)
    ...     urllib.urlretrieve(..., reporthook=reporthook)
    """
    last_b = [0]

    def update_to(b=1, bsize=1, tsize=None):
        """
        b  : int, optional
            Number of blocks transferred so far [default: 1].
        bsize  : int, optional
            Size of each block (in tqdm units) [default: 1].
        tsize  : int, optional
            Total size (in tqdm units). If [default: None] remains unchanged.
        """
        if tsize is not None:
            t.total = tsize
        t.update((b - last_b[0]) * bsize)
        last_b[0] = b

    return update_to

