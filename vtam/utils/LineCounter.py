import gzip 
import bz2
from functools import partial

class LineCounter():

    def __init__(self, filename):
        self.filename = filename

    def sequence_counter(self):

        def _make_gen(reader):
            b = reader( 1024 * 1024 )
            while b:
                yield b
                b = reader( 1024 * 1024 )

        def rawgencount(filename):

            if filename.endswith(".gz"):      
                _open = partial(gzip.open) 
            elif filename.endswith(".bz2"):
                _open = partial(bz2.open)
            else:
                _open = open

            f = open(filename, 'rb')
            f_gen = _make_gen(f.raw.read)
            return sum( buf.count(b'>') for buf in f_gen )
        
        return rawgencount(self.filename)
