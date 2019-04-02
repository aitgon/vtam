# -*- coding: utf-8 -*-

from unittest import TestCase

from bin.create_taxonomy_db import create_parser, f_create_taxonomy_db

class TestCreateTaxonomyDBSqlite(TestCase):

    def setUp(self):
        self.parser = create_parser()

    def test_parser(self):
        args = self.parser.parse_args(['-o', 'output.sqlite'])
        result = f_create_taxonomy_db(args.output)
        self.assertTrue(False)

