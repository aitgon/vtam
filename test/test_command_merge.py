# -*- coding: utf-8 -*-
import io
from unittest import TestCase

import yaml

from vtam.CommandMerge import VSearchMergeRunner

import tempfile


class TestCommandMerge(TestCase):

    def setUp(self):
        pass

    @classmethod
    def setUpClass(cls):
        pass

    def test_params(self):

        params_str = """fastq_ascii: 10"""
        params_yml = tempfile.NamedTemporaryFile()
        # Open the file for writing.
        with open(params_yml.name, 'w') as fout:
            fout.write(params_str)  # where `stuff` is, y'know... stuff to write (a string)

        vsearch_merge_runner = VSearchMergeRunner('fastq_fw_abspath', 'fastq_rv_abspath', 'fasta_abspath',
                                                  params_yml=params_yml.name)
        vsearch_merge_runner.load_parameters()
        self.assertTrue(vsearch_merge_runner.parameters['fastq_ascii'], 10)

    def test_wrong_params(self):

        wrong_params_str = """lfn_variant_replicate_threshold: 0.001"""
        wrong_params_yml = tempfile.NamedTemporaryFile()
        # Open the file for writing.
        with open(wrong_params_yml.name, 'w') as fout:
            fout.write(wrong_params_str)  # where `stuff` is, y'know... stuff to write (a string)

        with self.assertRaises(SystemExit) as se:
            vsearch_merge_runner = VSearchMergeRunner('fastq_fw_abspath', 'fastq_rv_abspath', 'fasta_abspath', params_yml=wrong_params_yml.name)
            vsearch_merge_runner.load_parameters()
        self.assertEqual(se.exception.code, 1)
