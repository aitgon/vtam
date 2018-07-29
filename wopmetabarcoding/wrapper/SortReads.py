import inspect
import tempfile

import os
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
from wopmetabarcoding.utils.PathFinder import PathFinder

from wopmetabarcoding.utils.VSearch import VSearch1
from wopmetabarcoding.utils.logger import logger
from wopmetabarcoding.wrapper.SortReadsUtilities import \
    create_primer_tag_fasta_for_vsearch, discard_tag_primer_alignment_with_low_sequence_quality,  trim_reads, \
    convert_trimmed_tsv_to_fasta, annotate_reads, gather_files, count_reads, insert_variant

import errno
# import subprocess

from wopmetabarcoding.utils.constants import tempdir


class SortReads(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.SortReads"
    }
    # Input
    # Input table
    __input_table_sample_information = "SampleInformation"
    __input_table_file = "File"
    __input_table_marker = "Marker"
    # Output
    # Output file
    __output_sortreads_samplecount_csv = "SortReads_sample_counts"
    # __output_file_vsearch_output_tsv = "vsearch_output_tsv"
    # __output_file_primer_tag_fasta = "primer_tag_fasta"
    # __output_file_checked_vsearch_output_tsv = "checked_vsearch_output_tsv"
    # Output table
    # __output_table_obifasta = 'ObiFasta'
    __output_table_variant = 'Variant'

    def specify_input_table(self):
        return [
            SortReads.__input_table_sample_information,
            SortReads.__input_table_file,
            SortReads.__input_table_marker
        ]

    def specify_output_file(self):
        return [
            SortReads.__output_sortreads_samplecount_csv
        ]

    def specify_output_table(self):
        return [
            SortReads.__output_table_variant,
        ]

    def specify_params(self):
        """

        :return:
        """
        return {
            "sort_reads_output_dir": 'str',
            "min_id": "float",
            "minseqlength": "int",
            "overhang": "int"

        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        con = engine.connect()
        #
        # Input tables models
        sample_information_model = self.input_table(SortReads.__input_table_sample_information)
        file_model = self.input_table(SortReads.__input_table_file)
        marker_model = self.input_table(SortReads.__input_table_marker)
        # Temp variables:
        annoted_tsv_list = []
        run_list = {}
        primer_tag_fasta = os.path.join(tempdir, 'primer_tag.fasta')
        checked_vsearch_output_tsv = os.path.join(tempdir, 'checked_vsearch_output.tsv')

        # Output file models
        sortreads_samplecount = self.output_file(SortReads.__output_sortreads_samplecount_csv)

        # Output tables models
        variant_model = self.output_table(SortReads.__output_table_variant)
        #
        # Parameters
        sort_reads_output_dir = self.option("sort_reads_output_dir")
        try:
            os.makedirs(sort_reads_output_dir)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        #  TODO: Later we will see in case files are already trimmed
        # for file_obj in session.query(file_model).all():
        #     if file_obj.trimmed_status is "1":
        #         raise Exception('Read trimming cannot be skipped')
        #         # read_counter(session, file_obj.name, readcount_model)
        #         # fasta_writer(session, readcount_model, merged_fasta)
        # else:
        for file_obj_i,file_obj in enumerate(session.query(file_model).order_by('name').all()):
            merged_fasta = file_obj.name
            logger.debug(
                "file: {}; line: {}; FASTA {} {}".format(__file__, inspect.currentframe().f_lineno, file_obj_i+1, merged_fasta))
            file_id = file_obj.id
            sample_information_obj = session.query(sample_information_model).filter(sample_information_model.file_id==file_id).all()
            PathFinder.mkdir_p(os.path.join(tempdir, "SortReads", os.path.basename(merged_fasta)))
            # 
            ############################################
            # First (reverse) trim
            ############################################
            is_forward_strand = True
            #
            ############################################
            # Vsearch --db primer_tag_fasta --usearch_global merged_fasta
            ############################################
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}".format(__file__, inspect.currentframe().f_lineno, file_obj_i+1, merged_fasta, is_forward_strand))
            # Logger.instance().info("Creating a fasta query file to align on the merged reads fasta for forward trimming.")
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; FASTA for forward trimming: {}".format(__file__, inspect.currentframe().f_lineno, file_obj_i+1, merged_fasta, is_forward_strand, primer_tag_fasta))
            create_primer_tag_fasta_for_vsearch(sample_information_obj, is_forward_strand, primer_tag_fasta)
            # Logger.instance().info("Processing Vsearch for forward trimming.")
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; VSearch forward trimming".format(__file__, inspect.currentframe().f_lineno, file_obj_i+1, merged_fasta, is_forward_strand))
            #
            # self.vsearch_subprocess(merged_fasta, is_forward_strand, primer_tag_fasta, vsearch_output_tsv)
            vsearch_output_tsv = os.path.join(tempdir, "SortReads", os.path.basename(merged_fasta), "vsearch_output_tsv")

            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; vsearch_output_tsv".format(__file__, inspect.currentframe().f_lineno, file_obj_i+1, merged_fasta, is_forward_strand))
            vsearch_params = {'db': primer_tag_fasta,
                              'usearch_global': merged_fasta,
                              'id': str(self.option("min_id")),
                              'maxhits': 1,
                              'maxrejects': 0,
                              'maxaccepts': 0,
                              'minseqlength': str(self.option("minseqlength")),
                              'userfields': "query+target+tl+qilo+qihi+tilo+tihi+qrow",
                              'userout': vsearch_output_tsv,
                              }
            vsearch1 = VSearch1(**vsearch_params)
            vsearch1.run()
            del vsearch1
            #
            ############################################
            # discard_tag_primer_alignment_with_low_sequence_quality
            ############################################
            # Logger.instance().info("Eliminating non SRS conforms reads for forward trimming.")
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; Eliminating non SRS conforms reads for forward trimming".format(__file__, inspect.currentframe().f_lineno, file_obj_i+1, merged_fasta, is_forward_strand))
            discard_tag_primer_alignment_with_low_sequence_quality(vsearch_output_tsv, checked_vsearch_output_tsv, self.option("overhang"))
            #
            ############################################
            # Trim reads and write to sqlite
            ############################################
            # Logger.instance().info("Trimming reads for forward trimming.")
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; Trimming reads for forward trimming".format(__file__, inspect.currentframe().f_lineno, file_obj_i+1, merged_fasta, is_forward_strand))
            # trimmed_tsv = os.path.join(tempdir, os.path.basename(merged_fasta), 'forward_trimmed.tsv')
            trimmed_tsv = os.path.join(tempdir, "SortReads", os.path.basename(merged_fasta), 'forward_trimmed.tsv')
            temp_db_sqlite = os.path.join(tempdir, "SortReads", os.path.basename(merged_fasta), 'forward_trimmed.sqlite')
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; Trimming reads for forward trimming: {}".format(__file__, inspect.currentframe().f_lineno, file_obj_i+1, merged_fasta, is_forward_strand, temp_db_sqlite))
            trim_reads(checked_vsearch_output_tsv, merged_fasta, trimmed_tsv, temp_db_sqlite)
            #
            ############################################
            # convert_trimmed_tsv_to_fasta
            ############################################
            # trimmed_fasta = trimmed_tsv.replace('.tsv', '.fasta')
            Logger.instance().info("Writing fasta file for forward trimming.")
            trimmed_fasta = os.path.join(tempdir, "SortReads", os.path.basename(merged_fasta), 'forward_trimmed.fasta')
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; Writing fasta file for trimming.: {}".format(__file__, inspect.currentframe().f_lineno, file_obj_i+1, merged_fasta, is_forward_strand, trimmed_fasta))
            convert_trimmed_tsv_to_fasta(trimmed_tsv, trimmed_fasta)
            #
            ############################################
            #
            # Second (reverse) trim
            #
            ############################################
            #
            is_forward_strand = False
            #
            ############################################
            # Vsearch --db primer_tag_fasta --usearch_global merged_fasta
            ############################################
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}".format(__file__, inspect.currentframe().f_lineno,
                                                                     file_obj_i+1, merged_fasta, is_forward_strand))
            # Logger.instance().info("Creating a fasta query file to align on the merged reads fasta for forward trimming.")
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; FASTA for forward trimming: {}".format(__file__,
                                                                                                  inspect.currentframe().f_lineno,
                                                                                                     file_obj_i+1,
                                                                                                     merged_fasta,
                                                                                                  is_forward_strand,
                                                                                                  primer_tag_fasta))
            create_primer_tag_fasta_for_vsearch(sample_information_obj, is_forward_strand, primer_tag_fasta)
            # Logger.instance().info("Processing Vsearch for forward trimming.")
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; VSearch forward trimming".format(__file__,
                                                                                            inspect.currentframe().f_lineno,
                                                                                               file_obj_i+1, merged_fasta,
                                                                                            is_forward_strand))
            #
            # self.vsearch_subprocess(merged_fasta, is_forward_strand, primer_tag_fasta, vsearch_output_tsv)
            vsearch_output_tsv = os.path.join(tempdir, "SortReads", os.path.basename(merged_fasta),
                                              "vsearch_output_tsv")

            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; vsearch_output_tsv".format(__file__,
                                                                                      inspect.currentframe().f_lineno,
                                                                                         file_obj_i+1, merged_fasta, is_forward_strand))
            vsearch_params = {'db': primer_tag_fasta,
                              'usearch_global': trimmed_fasta,
                              'id': str(self.option("min_id")),
                              'maxhits': 1,
                              'maxrejects': 0,
                              'maxaccepts': 0,
                              'minseqlength': str(self.option("minseqlength")),
                              'userfields': "query+target+tl+qilo+qihi+tilo+tihi+qrow",
                              'userout': vsearch_output_tsv,
                              }
            vsearch1 = VSearch1(**vsearch_params)
            vsearch1.run()
            del vsearch1
            #
            ############################################
            # discard_tag_primer_alignment_with_low_sequence_quality
            ############################################
            # Logger.instance().info("Eliminating non SRS conforms reads for forward trimming.")
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; Eliminating non SRS conforms reads for forward trimming".format(
                    __file__, inspect.currentframe().f_lineno, file_obj_i+1, merged_fasta, is_forward_strand))
            discard_tag_primer_alignment_with_low_sequence_quality(vsearch_output_tsv, checked_vsearch_output_tsv,
                                                                   self.option("overhang"))
            #
            ############################################
            # Trim reads and write to sqlite
            ############################################
            # Logger.instance().info("Trimming reads for forward trimming.")
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; Trimming reads for reverse trimming".format(__file__,
                                                                                                       inspect.currentframe().f_lineno,
                                                                                                          file_obj_i+1,
                                                                                                          merged_fasta,
                                                                                                       is_forward_strand))
            # trimmed_tsv = os.path.join(tempdir, os.path.basename(merged_fasta), 'forward_trimmed.tsv')
            trimmed_tsv = os.path.join(tempdir, "SortReads", os.path.basename(merged_fasta), 'reverse_trimmed.tsv')
            temp_db_sqlite = os.path.join(tempdir, "SortReads", os.path.basename(merged_fasta),
                                          'reverse_trimmed.sqlite')
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; Trimming reads for reverse trimming: {}".format(__file__,
                                                                                                           inspect.currentframe().f_lineno,
                                                                                                              file_obj_i+1,
                                                                                                              merged_fasta,
                                                                                                           is_forward_strand,
                                                                                                           temp_db_sqlite))
            # trim_reads(checked_vsearch_output_tsv, trimmed_fasta, trimmed_tsv, tempdir)
            trim_reads(checked_vsearch_output_tsv, trimmed_fasta, trimmed_tsv, temp_db_sqlite)
            #
            # trimmed_fasta = trimmed_tsv.replace('.tsv', '.fasta')
            # for file_obj in session.query(file_model).all():
            Logger.instance().info("Annotating reads with Sample Information.")
            # annotated_reads_tsv = (merged_fasta).replace(".fasta", "_annotated_reads.tsv")
            annotated_reads_tsv = os.path.join(tempdir, os.path.basename(merged_fasta).replace('.fasta', '_annotated_reads.tsv'))
            annoted_tsv_list.append(annotated_reads_tsv)
            run_list[annotated_reads_tsv] = file_obj.run_name
            logger.debug(
                "file: {}; line: {}; trimmed_tsv {}".format(__file__, inspect.currentframe().f_lineno, trimmed_tsv))
            annotate_reads(session, sample_information_model, trimmed_tsv,
                           file_id=file_id, annotated_reads_tsv=annotated_reads_tsv)
        #
        # # For each marker_id, concatenate its files with the annotated reads and count unique reads per marker_id
        with open(sortreads_samplecount, 'w') as fout_sortread_samplecount:
            for marker_obj in session.query(marker_model).all():
                marker_name = marker_obj.name
                logger.debug(
                    "file: {}; line: {}; marker_name {}".format(__file__, inspect.currentframe().f_lineno,
                                                                marker_name))
                marker_id = marker_obj.id
                gathered_marker_file = os.path.join(tempdir, marker_name + "_file.tsv")
                sample_count_tsv = os.path.join(sort_reads_output_dir, marker_name + "_sample_count.tsv")
                fout_sortread_samplecount.write(marker_name + "\t" + str(marker_id) + "\t" + sample_count_tsv + "\n")
                count_reads_marker = gathered_marker_file.replace(".tsv", ".sqlite")
                Logger.instance().info("Gathering all files from annotated files the same marker_id into one.")
                gather_files(marker_name, gathered_marker_file, annoted_tsv_list, run_list)
                gathered_marker_file = os.path.join(tempdir, marker_name + "_file_prerun.tsv")
                # session.commit()
                # read_count_per_marker_sqlite = conn_name = marker_obj.name + ".sqlite"
                Logger.instance().info("Counting reads for each marker_id.")
                count_reads(gathered_marker_file, count_reads_marker, marker_name, sample_count_tsv)
                # session.commit()
                Logger.instance().info("Inserting variant in the Variant table of the database.")
                insert_variant(con, count_reads_marker, variant_model)
                session.commit()
