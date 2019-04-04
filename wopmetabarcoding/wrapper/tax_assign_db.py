'''
A simple one-file solution for those times when you want to check if one organism is
a descended from another, but don't need a full phylogenetic tree manipulation library.

The library is just a single file that depends only upon the standard library.
You can easily embed it in another library by copying this script.
'''
from unittest import TestCase

import sqlite3
import tarfile
import os
import re
import shutil
from tempfile import mkdtemp
from os import remove
import urllib.request


def _strip_tab(name):
    return name.replace('\t', '')

#: The default location to download taxonomy information from
SOURCE_URL = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
#: The separator used to tokenize the in-database lineage string
SEP_TOKEN = "zzz"


def _from_ftp(url=SOURCE_URL):
    try:
        tempdir = mkdtemp()
        if url is not None:
            with open(os.path.join(tempdir, "taxdump.tar.gz"), 'wb') as outhandle:
                outhandle.writelines(urllib.request.urlopen(url))
        archive_path = os.path.join(tempdir, "taxdump.tar.gz")
        with tarfile.open(archive_path, 'r:gz') as tarchive:
            tarchive.extract("names.dmp", tempdir)
            tax2name = {}

            for line in open(os.path.join(tempdir, "names.dmp")):
                tax_id, name, unique_name, name_class, _ = map(_strip_tab, line.split("|"))
                if name_class == "scientific name":
                    tax2name[tax_id] = unique_name if name == "" else name
            remove(os.path.join(tempdir, "names.dmp"))
            tarchive.extract("nodes.dmp", tempdir)

            for line in open(os.path.join(tempdir, "nodes.dmp")):
                parts = map(_strip_tab, line.split("|"))
                tax_id, parent_tax_id, rank = parts[:3]
                name = tax2name[tax_id]
                yield tax_id, name, parent_tax_id, rank
            remove(os.path.join(tempdir, "nodes.dmp"))
        remove(archive_path)
    finally:
        shutil.rmtree(tempdir)


class Taxonomy(object):
    """Operate on taxonomic hierarchies downloaded from the NCBI Taxonomy database
    using a compact SQLite database.

    Parameters
    ----------
    store_path: str
        Path to the sqlite database containing the hierarchies


    Attributes
    ----------
    connection: sqlite3.Connection
        The underlying connection to the sqlite database
    """
    @classmethod
    def from_source(cls, store_path='taxonomy.db', url=SOURCE_URL):
        """Construct a new :class:`Taxonomy` instance and associated database file
        from source data downloaded from NCBI's FTP servers.

        If `url` is :const:`None`, then it will look for the source information in the
        current directory at the name "taxdump.tar.gz".


        Parameters
        ----------
        store_path: str
            Path to construct the database at. Defaults to "taxonomy.db" in the current
            directory
        url: str
            The URL to download the taxonomy information from. Defaults to :data:`SOURCE_URL`

        Returns
        -------
        :class:`Taxonomy`
        """
        store = cls(store_path)
        store._init_schema()
        store.executemany('INSERT INTO taxonomy VALUES (?,?,?,?,"");', _from_ftp(url))
        store._construct_lineage()
        store._init_index()
        store.commit()
        return store

    def __init__(self, store_path):
        self.store_path = store_path
        self.connection = sqlite3.connect(store_path)
        try:
            base_lineage = self.execute("SELECT lineage FROM taxonomy where taxa_id = 1;").next()
            self.sep = re.split(r'\d', base_lineage)[0]
        except:
            self.sep = 'zzz'

    def _init_schema(self):
        self.execute('DROP TABLE IF EXISTS taxonomy')
        self.execute('''CREATE TABLE taxonomy (taxa_id INTEGER PRIMARY KEY,
                                               taxa_name VARCHAR(50),
                                               parent_taxa INTEGER,
                                               rank VARCHAR(20),
                                               lineage VARCHAR(200));''')
        self.commit()

    def _init_index(self):
        self.execute('''CREATE INDEX IF NOT EXISTS taxname ON taxonomy(taxa_name);''')
        self.execute('''CREATE INDEX IF NOT EXISTS parent_id ON taxonomy(parent_taxa);''')
        self.execute('''CREATE INDEX IF NOT EXISTS lineage ON taxonomy(lineage);''')
        self.commit()

    def _construct_lineage(self):
        tax2lineage = {}
        for row in self.execute("SELECT * FROM taxonomy"):
            lineage = self.lineage(row[0])
            tax2lineage[row[0]] = self.sep + self.sep.join(map(str, lineage)) + self.sep

        self.executemany("UPDATE taxonomy SET lineage = ?2 WHERE taxa_id = ?1;", tax2lineage.iteritems())

    def execute(self, stmt, args=""):
        """Execute raw SQL against the underlying database.

        See :meth:`sqlite3.Connection.execute`
        """
        return self.connection.execute(stmt, args)

    def executemany(self, stmt, args=""):
        return self.connection.executemany(stmt, args)

    def close(self):
        """Close the underlying database connection.

        See :meth:`sqlite3.Connection.close`
        """
        self.connection.close()

    def commit(self):
        """Save pending changes to the underlying database

        See :meth:`sqlite3.Connection.commit`
        """
        self.connection.commit()

    def name_to_tid(self, name):
        """Translates a scientific name `name` string into its equivalent taxonomic id number

        Parameters
        ----------
        name: str
            A scientific name like "Homo sapiens"

        Returns
        -------
        tid: int
        """
        result = self.execute(
            "SELECT taxa_id FROM taxonomy WHERE taxa_name = ?", (name,)).fetchone()
        if result is not None:
            result = result[0]
        return result

    def tid_to_name(self, tid):
        """Translates a taxonomic id number `tid` into its equivalent scientific name

        Parameters
        ----------
        tid: int
            A taxonomic id number like 9606
        Returns
        -------
        name: str
            A scientific name like "Homo sapiens"
        """

        tid = (tid,)
        result = self.execute(
            "SELECT taxa_name FROM taxonomy WHERE taxa_id = ?", tid).fetchone()
        if result is not None:
            result = result[0]
        return result

    def tid_to_rank(self, tid):
        tid = (tid,)
        result = self.execute(
            "SELECT rank FROM taxonomy WHERE taxa_id = ?", tid).fetchone()
        if result is not None:
            result = result[0]
        return result

    def is_parent(self, child_tid, parent_tid):
        """Test if `parent_tid` is a parent taxa of `child_tid`

        Parameters
        ----------
        child_tid: int
        parent_tid: int

        Returns
        -------
        bool
        """
        parent_tid = "%{}{}{}%".format(self.sep, parent_tid, self.sep)
        child_tid = child_tid
        try:
            self.execute('''SELECT taxa_id FROM taxonomy WHERE taxa_id = ?
                                                         AND lineage LIKE ?''',
                         (child_tid, parent_tid)).next()
            return True
        except StopIteration:
            return False

    def lineage(db, tid):
        """Construct the taxonomic "path" from `tid` to the root of the
        phylogenetic hierarchy

        Parameters
        ----------
        tid: int

        Returns
        -------
        list of ints
        """
        path = []
        path.append(tid)
        while tid != 1:
            tid = db.parent(tid)
            if tid is None:
                break
            path.append(tid)
        return path[::-1]

    def children(self, tid, deep=False):
        """Retrieve all child taxonomic id numbers of `tid`. If `deep` is `True`, retrieve all descendants

        Parameters
        ----------
        tid: int
        deep: bool
            Retrieve all descendants, not just direct children

        Returns
        -------
        list of ints
        """
        tid = (tid,)
        children = []
        if deep:
            tid = tid[0]
            tid_str = ("%{}{}{}%".format(self.sep, tid, self.sep),)
            for row in self.execute("SELECT taxa_id FROM taxonomy WHERE lineage LIKE ?", tid_str):
                if row[0] != tid:
                    children.append(row[0])
        else:
            for row in self.execute("SELECT taxa_id FROM taxonomy WHERE parent_taxa = ?", tid):
                children.append(row[0])
        return children

    def parent(self, tid):
        """Extract the taxonomic id number of the parent of `tid`

        Parameters
        ----------
        tid: int

        Returns
        -------
        int
        """
        tid = (tid,)
        parent = self.execute("SELECT parent_taxa FROM taxonomy WHERE taxa_id = ?", tid).fetchone()
        return parent if parent is None else parent[0]

    def siblings(self, tid):
        """Extract the taxonomic id numbers of the siblings (same parent) of `tid`

        Parameters
        ----------
        tid: int

        Returns
        -------
        list of ints
        """
        parent = self.parent(tid)
        siblings = self.children(parent)
        siblings.remove(tid)
        return siblings

    def relatives(self, tid, degree=1):
        """Retrieve relatives of `tid` out to `degree` steps removed

        Parameters
        ----------
        tid: int
        degree: int

        Returns
        -------
        list of ints
        """
        root = self.parent(tid)
        for i in range(degree):
            root = self.parent(root)

        relatives = []
        current_layer = [root]
        next_layer = []
        for i in range(degree * 2):
            for entry in current_layer:
                next_layer.extend(row for row in self.children(entry))
            relatives.extend(current_layer)
            current_layer = next_layer
            next_layer = []
        return relatives

    def nearest_common_ancestor(self, a, b):
        alineage = self.lineage(a)
        blineage = self.lineage(b)
        for i, tida in enumerate(alineage[::-1]):
            for j, tidb in enumerate(blineage[::-1]):
                if tida == tidb:
                    return i + j, tida



