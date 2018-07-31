import inspect

import os

import sqlalchemy
from sqlalchemy.sql import select

import pandas
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
from wopmetabarcoding.utils.PathFinder import PathFinder

from wopmetabarcoding.utils.VSearch import VSearch1
from wopmetabarcoding.utils.logger import logger
from wopmetabarcoding.wrapper.SortReadsUtilities import \
    create_primer_tag_fasta_for_vsearch, discard_tag_primer_alignment_with_low_sequence_quality,  trim_reads, \
    convert_trimmed_tsv_to_fasta, annotate_reads, gather_files, count_reads, insert_variant

import errno

from wopmetabarcoding.utils.constants import tempdir


class SortReads(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.SortReads"
    }
    # Input
    # Input table
    __input_table_sample_information = "SampleInformation"
    __input_table_fasta = "Fasta"
    __input_table_marker = "Marker"
    # Output
    # Output file
    # __output_read_count_tsv = "read_count"
    # __output_file_vsearch_output_tsv = "vsearch_output_tsv"
    # __output_file_primer_tag_fasta = "primer_tag_fasta"
    # __output_file_checked_vsearch_output_tsv = "checked_vsearch_output_tsv"
    # Output table
    # __output_table_obifasta = 'ObiFasta'
    __output_table_variant = 'Variant'
    __output_table_variant_read_count = 'VariantReadCount'

    def specify_input_table(self):
        return [
            SortReads.__input_table_sample_information,
            SortReads.__input_table_fasta,
            SortReads.__input_table_marker
        ]

    def specify_output_file(self):
        return [
            # SortReads.__output_read_count_tsv
        ]

    def specify_output_table(self):
        return [
            SortReads.__output_table_variant,
            SortReads.__output_table_variant_read_count,
        ]

    def specify_params(self):
        """

        :return:
        """
        return {
            # "sort_reads_output_dir": 'str',
            "min_id": "float",
            "minseqlength": "int",
            "overhang": "int"

        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        conn = engine.connect()
        #
        # Input tables models
        sample_information_model = self.input_table(SortReads.__input_table_sample_information)
        fasta_model = self.input_table(SortReads.__input_table_fasta)
        marker_model = self.input_table(SortReads.__input_table_marker)
        # Temp variables:
        tsv_file_list_with_read_annotations = [] # List of TSV files with read annotations
        run_list = {}
        primer_tag_fasta = os.path.join(tempdir, 'primer_tag.fasta')
        checked_vsearch_output_tsv = os.path.join(tempdir, 'checked_vsearch_output.tsv')

        # Output file models
        # read_count_tsv = self.output_file(SortReads.__output_read_count_tsv)

        # Output tables models
        variant_model = self.output_table(SortReads.__output_table_variant)
        variant_read_count_model = self.output_table(SortReads.__output_table_variant_read_count)
        #
        # Parameters
        sort_reads_output_dir = self.option("sort_reads_output_dir")
        # try:
        #     os.makedirs(sort_reads_output_dir)
        # except OSError as exception:
        #     if exception.errno != errno.EEXIST:
        #         raise
        #  TODO: Later we will see in case files are already trimmed
        # for file_obj in session.query(file_model).all():
        #     if file_obj.trimmed_status is "1":
        #         raise Exception('Read trimming cannot be skipped')
        #         # read_counter(session, file_obj.name, readcount_model)
        #         # fasta_writer(session, readcount_model, merged_fasta)
        # else:
        marker2fasta2readannotationtsv_dict = {} # Dict of dicts, where for each marker, there are fasta and readannotationtsv
        for fasta_obj in session.query(fasta_model).order_by('name').all():
            fasta_id = fasta_obj.id
            fasta_name = fasta_obj.name
            # Get marker of this fasta file
            marker_id = session.query(sample_information_model).filter(
                sample_information_model.fasta_id == fasta_id).first().marker_id
            marker_name = session.query(marker_model).filter(
                marker_model.id == marker_id).first().name
            logger.debug(
                "file: {}; line: {}; FASTA {} {}".format(__file__, inspect.currentframe().f_lineno, fasta_id, fasta_name))
            # file_id = fasta_obj.id
            sample_information_obj = session.query(sample_information_model).filter(sample_information_model.fasta_id==fasta_id).all()
            PathFinder.mkdir_p(os.path.join(tempdir, "SortReads", os.path.basename(fasta_name)))
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
                "file: {}; line: {}; FASTA {} {}; forward {}".format(__file__, inspect.currentframe().f_lineno, fasta_id, fasta_name, is_forward_strand))
            # Logger.instance().info("Creating a fasta query file to align on the merged reads fasta for forward trimming.")
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; FASTA for forward trimming: {}".format(__file__, inspect.currentframe().f_lineno, fasta_id, fasta_name, is_forward_strand, primer_tag_fasta))
            create_primer_tag_fasta_for_vsearch(sample_information_obj, is_forward_strand, primer_tag_fasta)
            # Logger.instance().info("Processing Vsearch for forward trimming.")
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; VSearch forward trimming".format(__file__, inspect.currentframe().f_lineno, fasta_id, fasta_name, is_forward_strand))
            #
            # self.vsearch_subprocess(merged_fasta, is_forward_strand, primer_tag_fasta, vsearch_output_tsv)
            vsearch_output_tsv = os.path.join(tempdir, "SortReads", os.path.basename(fasta_name), "vsearch_output_tsv")

            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; vsearch_output_tsv".format(__file__, inspect.currentframe().f_lineno, fasta_id, fasta_name, is_forward_strand))
            vsearch_params = {'db': primer_tag_fasta,
                              'usearch_global': fasta_name,
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
                "file: {}; line: {}; FASTA {} {}; forward {}; Eliminating non SRS conforms reads for forward trimming".format(__file__, inspect.currentframe().f_lineno, fasta_id, fasta_name, is_forward_strand))
            discard_tag_primer_alignment_with_low_sequence_quality(vsearch_output_tsv, checked_vsearch_output_tsv, self.option("overhang"))
            #
            ############################################
            # Trim reads and write to sqlite
            ############################################
            # Logger.instance().info("Trimming reads for forward trimming.")
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; Trimming reads for forward trimming".format(__file__, inspect.currentframe().f_lineno, fasta_id, fasta_name, is_forward_strand))
            # trimmed_tsv = os.path.join(tempdir, os.path.basename(merged_fasta), 'forward_trimmed.tsv')
            trimmed_tsv = os.path.join(tempdir, "SortReads", os.path.basename(fasta_name), 'forward_trimmed.tsv')
            temp_db_sqlite = os.path.join(tempdir, "SortReads", os.path.basename(fasta_name), 'forward_trimmed.sqlite')
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; Trimming reads for forward trimming: {}".format(__file__, inspect.currentframe().f_lineno, fasta_id, fasta_name, is_forward_strand, temp_db_sqlite))
            trim_reads(checked_vsearch_output_tsv, fasta_name, trimmed_tsv, temp_db_sqlite)
            #
            ############################################
            # convert_trimmed_tsv_to_fasta
            ############################################
            # trimmed_fasta = trimmed_tsv.replace('.tsv', '.fasta')
            Logger.instance().info("Writing fasta file for forward trimming.")
            trimmed_fasta = os.path.join(tempdir, "SortReads", os.path.basename(fasta_name), 'forward_trimmed.fasta')
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; Writing fasta file for trimming.: {}".format(__file__, inspect.currentframe().f_lineno, fasta_id, fasta_name, is_forward_strand, trimmed_fasta))
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
                                                                     fasta_id, fasta_name, is_forward_strand))
            # Logger.instance().info("Creating a fasta query file to align on the merged reads fasta for forward trimming.")
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; FASTA for forward trimming: {}".format(__file__,
                                                                                                  inspect.currentframe().f_lineno,
                                                                                                     fasta_id,
                                                                                                     fasta_name,
                                                                                                  is_forward_strand,
                                                                                                  primer_tag_fasta))
            create_primer_tag_fasta_for_vsearch(sample_information_obj, is_forward_strand, primer_tag_fasta)
            # Logger.instance().info("Processing Vsearch for forward trimming.")
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; VSearch forward trimming".format(__file__,
                                                                                            inspect.currentframe().f_lineno,
                                                                                               fasta_id, fasta_name,
                                                                                            is_forward_strand))
            #
            # self.vsearch_subprocess(merged_fasta, is_forward_strand, primer_tag_fasta, vsearch_output_tsv)
            vsearch_output_tsv = os.path.join(tempdir, "SortReads", os.path.basename(fasta_name),
                                              "vsearch_output_tsv")

            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; vsearch_output_tsv".format(__file__,
                                                                                      inspect.currentframe().f_lineno,
                                                                                         fasta_id, fasta_name, is_forward_strand))
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
                    __file__, inspect.currentframe().f_lineno, fasta_id, fasta_name, is_forward_strand))
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
                                                                                                          fasta_id,
                                                                                                          fasta_name,
                                                                                                       is_forward_strand))
            # trimmed_tsv = os.path.join(tempdir, os.path.basename(merged_fasta), 'forward_trimmed.tsv')
            trimmed_tsv = os.path.join(tempdir, "SortReads", os.path.basename(fasta_name), 'reverse_trimmed.tsv')
            temp_db_sqlite = os.path.join(tempdir, "SortReads", os.path.basename(fasta_name),
                                          'reverse_trimmed.sqlite')
            logger.debug(
                "file: {}; line: {}; FASTA {} {}; forward {}; Trimming reads for reverse trimming: {}".format(__file__,
                                                                                                           inspect.currentframe().f_lineno,
                                                                                                              fasta_id,
                                                                                                              fasta_name,
                                                                                                           is_forward_strand,
                                                                                                           temp_db_sqlite))
            # trim_reads(checked_vsearch_output_tsv, trimmed_fasta, trimmed_tsv, tempdir)
            trim_reads(checked_vsearch_output_tsv, trimmed_fasta, trimmed_tsv, temp_db_sqlite)
            #
            # trimmed_fasta = trimmed_tsv.replace('.tsv', '.fasta')
            # for file_obj in session.query(file_model).all():
            Logger.instance().info("Annotating reads with Sample Information.")
            # read_annotation_tsv = (merged_fasta).replace(".fasta", "_annotated_reads.tsv")
            # read_annotation_tsv = os.path.join(tempdir, os.path.basename(merged_fasta).replace('.fasta', '_annotated_reads.tsv'))
            # One TSV file with read annotation per merged FASTA Fasta
            read_annotation_tsv = os.path.join(tempdir, "SortReads", os.path.basename(fasta_name), 'read_annotation.tsv')
            tsv_file_list_with_read_annotations.append(read_annotation_tsv)
            run_list[read_annotation_tsv] = fasta_obj.run_name
            logger.debug(
                "file: {}; line: {}; trimmed_tsv {}".format(__file__, inspect.currentframe().f_lineno, trimmed_tsv))
            ################################################################
            # Annotated reads
            ################################################################
            marker2fasta2readannotationtsv_dict[marker_name] = {}
            marker2fasta2readannotationtsv_dict[marker_name][fasta_name] = read_annotation_tsv
            annotate_reads(session, sample_information_model, trimmed_tsv,
                           fasta_id=fasta_id, out_tsv=read_annotation_tsv)
            read_annotation_df = pandas.read_csv(read_annotation_tsv, sep='\t',
                                 header=None,
                                 names=['read_id', 'marker_id', 'run', 'tag_forward', 'tag_reverse', 'biosample_id',
                                        'replicate_id', 'variant_sequence'])
            fasta_variant_count_df = read_annotation_df.groupby(['run', 'biosample_id', 'replicate_id', 'variant_sequence']).size().reset_index(name='count')
            logger.debug(
                "file: {}; line: {};  Insert variants: marker {} fasta {}".format(__file__, inspect.currentframe().f_lineno, marker_name, fasta_name))
            for row in fasta_variant_count_df.itertuples():
                biosample_id = row[2]
                replicate_id = row[3]
                variant_sequence = row[4]
                read_count = row[5]
                # stmt = variant_model.__table__.insert().values(marker_id=marker_id, biosample_id=biosample_id, replicate_id=replicate_id,
                #                                                sequence=variant_sequence,
                #                                                read_count=read_count)
                # conn.execute(stmt)
                # variant_table = variant_model.__table__
                # variant_id = None
                try:
                    stmt_ins_var = variant_model.__table__.insert().values(sequence=variant_sequence)
                    stmt_result_var = conn.execute(stmt_ins_var)
                    variant_id = stmt_result_var.inserted_primary_key[0]
                except sqlalchemy.exc.IntegrityError:
                    stmt_select_var = select([variant_model.__table__.c.id]).where(variant_model.__table__.c.sequence==variant_sequence)
                    variant_id = conn.execute(stmt_select_var).first()[0]
                try:
                    #
                    stmt_ins_read_count = variant_read_count_model.__table__.insert().values(variant_id=variant_id, marker_id=marker_id, biosample_id=biosample_id, replicate_id=replicate_id, read_count=read_count)
                    conn.execute(stmt_ins_read_count)
                except sqlalchemy.exc.IntegrityError:
                    stmt_upd = variant_read_count_model.__table__.update()\
                        .where(variant_read_count_model.__table__.c.variant_id==variant_id)\
                        .where(variant_read_count_model.__table__.c.marker_id==marker_id)\
                        .where(variant_read_count_model.__table__.c.biosample_id==biosample_id)\
                        .where(variant_read_count_model.__table__.c.replicate_id==replicate_id)\
                        .values(read_count=read_count)
                    # import pdb; pdb.set_trace()
                    conn.execute(stmt_upd)
            logger.debug("file: {}; line: {};  Insert variants: Finished".format(__file__, inspect.currentframe().f_lineno))
            # variant_data = conn.execute("SELECT id, marker_id, seq, count FROM count_read")
            # check_read = conn.execute('SELECT EXISTS (SELECT id FROM count_read WHERE seq=?)', (line[7],))
            # from sqlalchemy.sql import select
            # s = select([users])
            # for row in conn.execute(s):
            #     print(row)
        # for marker_obj in session.query(marker_model).all():
        #     marker_name = marker_obj.name
        #     marker_read_count_tsv = os.path.join(sort_reads_output_dir, "read_count_{}.tsv".format(marker_name))
        #     for marker_name_fasta in marker2fasta2readannotationtsv_dict[marker_name]:
        #
        # For each marker_id, concatenate its files with the annotated reads and count unique reads per marker_id
        # # instead of function gather_files
        # for marker_obj in session.query(marker_model).all():
        #     marker_name = marker_obj.name
        #     PathFinder.mkdir_p(os.path.join(tempdir, "SortReads", marker_name))
        #     trimmed_marker_tsv = os.path.join(tempdir, "SortReads", marker_name, 'trimmed_marker.tsv')
        #     with open(trimmed_marker_tsv, 'a') as fout:
        #         for marker_fasta in marker2fasta2readannotationtsv_dict[marker_name]:
        #             read_annotation_tsv = marker2fasta2readannotationtsv_dict[marker_name][marker_fasta]
        #             with open(read_annotation_tsv, 'r') as marker_fasta_fin:
        #                 fout.write(marker_fasta_fin.read())
        #     read_count_marker_tsv = os.path.join(tempdir, "SortReads", marker_name, 'read_count.tsv')
        #     read_count_marker_sqlite = os.path.join(tempdir, "SortReads", marker_name, 'read_count.sqlite')
        #     out_tsv = os.path.join(tempdir, "SortReads", marker_name, 'out.sqlite')
        #     count_reads(read_count_marker_tsv, read_count_marker_sqlite, marker_name, out_tsv)
        #     # import pdb; pdb.set_trace()
        # # import pdb; pdb.set_trace()
        # with open(read_count_tsv, 'w') as fout_sortread_samplecount:
        #     for marker_obj in session.query(marker_model).all():
        #         marker_name = marker_obj.name
        #         logger.debug(
        #             "file: {}; line: {}; marker_name {}".format(__file__, inspect.currentframe().f_lineno,
        #                                                         marker_name))
        #         marker_id = marker_obj.id
        #         read_count_marker_tsv = os.path.join(tempdir, marker_name, "read_count.tsv")
        #         # import pdb; pdb.set_trace()
        #         # sample_count_tsv = os.path.join(sort_reads_output_dir, marker_name + "_sample_count.tsv")
        #         sample_count_tsv = "sample_count_tsv"
        #         fout_sortread_samplecount.write(marker_name + "\t" + str(marker_id) + "\t" + sample_count_tsv + "\n")
        #         count_reads_marker = read_count_marker_tsv.replace(".tsv", ".sqlite")
        #         Logger.instance().info("Gathering all files from annotated files the same marker_id into one.")
        #         ################################################################
        #         # gather_files
        #         ################################################################
        #         # gather_files(marker_name, read_count_marker_tsv, tsv_file_list_with_read_annotations, run_list)
        #         import pdb; pdb.set_trace()
        #         read_count_marker_tsv = os.path.join(tempdir, marker_name + "_file_prerun.tsv")
        #         # session.commit()
        #         # read_count_per_marker_sqlite = conn_name = marker_obj.name + ".sqlite"
        #         Logger.instance().info("Counting reads for each marker_id.")
        #         ################################################################
        #         # count_reads
        #         ################################################################
        #         count_reads(read_count_marker_tsv, count_reads_marker, marker_name, sample_count_tsv)
        #         # session.commit()
        #         Logger.instance().info("Inserting variant in the Variant table of the database.")
        #         insert_variant(conn, count_reads_marker, variant_model)
        #         session.commit()
