# -*- coding: utf-8 -*-

import os
import pandas
from unittest import TestCase

from sqlalchemy import create_engine

from vtam.models.Variant import Variant
from vtam.models.VariantReadCount import VariantReadCount
from vtam.utils.KnownOccurrences import KnownOccurrences
from vtam.utils.PathManager import PathManager

from vtam.models.Run import Run
from vtam.models.Marker import Marker
from vtam.models.Biosample import Biosample
from wopmars.Base import Base
from sqlalchemy.orm import sessionmaker


class TestKnownOccurrences(TestCase):

    @classmethod
    def setUp(self):

        # ################################################################################################################
        # #
        # # Creates memory database
        # #
        # ################################################################################################################
        #
        # self.engine = create_engine('sqlite://')  # memory db
        #
        # Session = sessionmaker(bind=self.engine)
        # self.session = Session()
        #
        # Base.metadata.create_all(self.engine)
        #
        # readinfo_tsv = os.tsv_path.join(PathManager.get_test_path(), 'test_files', 'readinfo.tsv')
        # df = pandas.read_csv(readinfo_tsv, sep="\t", header=0)
        #
        # for run_name in df.Run.unique().tolist():
        #     run = Run(name=run_name)
        #     self.session.add(run)
        #
        # for marker_name in df.Marker.unique().tolist():
        #     marker = Marker(name=marker_name)
        #     self.session.add(marker)
        #
        # for biosample_name in df.Biosample.unique().tolist():
        #     biosample = Biosample(name=biosample_name)
        #     self.session.add(biosample)
        #
        # self.session.commit()
        #
        # self.tsv_path = os.tsv_path.join(PathManager.get_test_path(), 'test_files', 'known_occurrences.tsv')
        # self.known_occurrences_df = pandas.read_csv(self.tsv_path, sep="\t", header=0)
        #
        # self.readinfo_tsv = os.tsv_path.join(PathManager.get_test_path(), 'test_files', 'readinfo.tsv')

        known_occurrences_tsv_path = os.path.join(PathManager.get_package_path(), 'doc/data/dryad.f40v5/known_occurrences.tsv')
        known_occurrences_df = KnownOccurrences.read_tsv_into_df(known_occurrences_tsv_path)
        import pdb; pdb.set_trace()
        #
        readinfo_tsv = os.path.join(PathManager.get_test_path(), 'test_files', 'readinfo.tsv')
        df = pandas.read_csv(readinfo_tsv, sep="\t", header=0)

    def test_variant_in_db(self):

        variant = Variant(sequence="TCTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT")
        self.session.add(variant)
        self.session.commit()

        variant_read_count = VariantReadCount(run_id=1, marker_id=1, biosample_id=1, replicate=1, variant_id=1, read_count=10)
        self.session.add(variant_read_count)
        self.session.commit()

        # Does not raises error
        known_occurrences_inst = KnownOccurrences(self.known_occurrences_df, readinfo_tsv=self.readinfo_tsv, engine=self.engine)

    def test_variant_not_in_db(self):

        variant = Variant(sequence="TAGC")
        self.session.add(variant)
        self.session.commit()

        with self.assertRaises(SystemExit):
            # Does raises error
            known_occurrences_inst = KnownOccurrences(self.known_occurrences_df, readinfo_tsv=self.readinfo_tsv, engine=self.engine)

    def test_run_marker_biosample_variant_not_in_db(self):

        variant = Variant(sequence="TCTATATTTCATTTTTGGTGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGAGCCGAGTTGGGTAACCCGGGTTCATTAATTGGGGACGATCAAATTTATAACGTAATCGTAACTGCTCATGCCTTTATTATGATTTTTTTTATAGTGATACCTATTATAATT")
        self.session.add(variant)

        variant_read_count = VariantReadCount(run_id=1, marker_id=1, biosample_id=2, replicate=1, variant_id=1, read_count=10)
        self.session.add(variant_read_count)
        self.session.commit()

        with self.assertRaises(SystemExit):
            # Does raises error
            known_occurrences_inst = KnownOccurrences(self.known_occurrences_df, readinfo_tsv=self.readinfo_tsv, engine=self.engine)
