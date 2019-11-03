import inspect

import os
import shutil

import pandas
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper

from vtam.utils.FastaInformation import FastaInformation
from vtam.utils.OptionManager import OptionManager
from vtam.utils.PathManager import PathManager
from vtam.utils.SortReadsRunner import SortReadsRunner

class SortReads(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.SortReads"
    }
    # Input
    # Input table
    __input_table_sample_information = "SampleInformation"
    __input_table_fasta = "Fasta"
    __input_table_run = "Run"
    __input_table_marker = "Marker"
    __input_table_biosample = "Biosample"
    __input_table_replicate = "Replicate"
    #Input file
    __input_file_fastainfo = "fastainfo"
    # Output
    # Output file
    # Output table
    __output_file_sort_reads = 'sortreads'

    def specify_input_table(self):
        return [
            SortReads.__input_table_sample_information,
            SortReads.__input_table_fasta,
            SortReads.__input_table_run,
            SortReads.__input_table_marker,
            SortReads.__input_table_biosample,
            SortReads.__input_table_replicate,

        ]

    def specify_input_file(self):
        return[
            SortReads.__input_file_fastainfo,
        ]

    def specify_output_file(self):
        return [
            SortReads.__output_file_sort_reads,
        ]

    def specify_params(self):
        """

        :return:
        """
        return {
            "fasta_dir": "str",
            "min_id": "float",
            "minseqlength": "int",
            "overhang": "int",
            "log_verbosity": "int",
            "log_file": "str",
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        OptionManager.instance()['log_verbosity'] = int(self.option("log_verbosity"))
        if not self.option("log_verbosity") is None:
            OptionManager.instance()['log_file'] = str(self.option("log_file"))
        this_step_tmp_dir = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__))
        PathManager.mkdir_p(this_step_tmp_dir)

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ##########################################################
        #
        # Input tables
        sample_information_model = self.input_table(SortReads.__input_table_sample_information)
        fasta_model = self.input_table(SortReads.__input_table_fasta)
        marker_model = self.input_table(SortReads.__input_table_marker)
        run_model = self.input_table(SortReads.__input_table_run)
        biosample_model = self.input_table(SortReads.__input_table_biosample)
        replicate_model = self.input_table(SortReads.__input_table_replicate)
        #
        # Input files
        input_file_fastainfo = self.input_file(SortReads.__input_file_fastainfo)
        #
        # Output files
        # sortreads with columns: read_name, fasta_id, run_id, marker_id, biosample_id, replicate_id
        sort_reads_tsv = self.output_file(SortReads.__output_file_sort_reads)
        #
        # Options
        fasta_dir = str(self.option("fasta_dir"))
        min_id = str(self.option("min_id"))
        minseqlength = str(self.option("minseqlength"))

        ##########################################################
        #
        # 1. Read fastainfo to get run_id, marker_id, biosample_id, replicate_id for current analysis
        #
        ##########################################################

        fasta_information_obj = FastaInformation(input_file_fastainfo, engine, run_model, marker_model, biosample_model, replicate_model)
        fasta_information_record_list = fasta_information_obj.get_fasta_information_record_list(extended_information=True)
        fasta_information_df = pandas.DataFrame(data=fasta_information_record_list)

        ############################################
        #
        # For each fasta_path file output in the DB (Table Fasta)
        #
        # 1. Trimming (Forward and reverse): Remove primer and tag sequence from each read sequence (Each sequence in Fasta)
        # 2. Store read count of each variant in table 'VariantReadCount'
        # 3. Eliminate singleton: Variants found one time throughout all biosample-replicates
        #
        ############################################

        sort_reads_tsv_list = []

        for fasta_file_name in sorted(fasta_information_df.fasta_file_name.unique().tolist()):
            fasta_information_per_fasta_df = fasta_information_df.loc[fasta_information_df.fasta_file_name==fasta_file_name]
            fasta_file_id = session.query(fasta_model.id).filter(fasta_model.name == fasta_file_name).one()[0]
            fasta_path = os.path.join(fasta_dir, fasta_file_name)

            alignement_parameters = {'min_id': min_id, 'minseqlength': minseqlength}

            sort_reads_tsv_i = os.path.join(this_step_tmp_dir, "sort_reads_fasta_id_{}.tsv".format(fasta_file_id))

            # Create SortReadsRunner
            sort_reads_runner = SortReadsRunner(fasta_path=fasta_path, fasta_id=fasta_file_id, alignement_parameters=alignement_parameters,
                                                    fasta_information_df=fasta_information_per_fasta_df, sort_reads_tsv=sort_reads_tsv_i)

            sort_reads_runner.run()
            sort_reads_tsv_list.append(sort_reads_tsv_i)

        ########################################################
        #
        # Concatenate sortreads files of different fasta_path files
        #
        ########################################################
        with open(sort_reads_tsv, 'w') as fout: # header
            fout.write('read_id\trun_id\tfasta_id\tmarker_id\tbiosample_id\treplicate_id\tread_sequence\n')
        with open(sort_reads_tsv, 'ab') as wfd:
            for f in sort_reads_tsv_list:
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd)

        # for fasta_obj in session.query(fasta_model).order_by('name').all():
        #     fasta_id = fasta_obj.id
        #     fasta_name = fasta_obj.name
        #
        #     alignement_parameters = {'min_id': min_id, 'minseqlength': minseqlength}
        #     import pdb; pdb.set_trace()
        #
        #     # Create SortReadsRunner
        #     sort_reads_runner = SortReadsRunner(fasta_path=fasta_name, fasta_id=fasta_id, alignement_parameters=alignement_parameters,
        #                                             fasta_information_df=fasta_information_df, sort_reads_tsv=sort_reads_tsv)







        #     PathManager.mkdir_p(os.path.join(this_step_tmp_dir, os.path.basename(fasta_name)))
        #     # Get marker of this fasta_path file
        #     marker_id = session.query(sample_information_model).filter(
        #         sample_information_model.fasta_id == fasta_id).first().marker_id
        #     marker_name = session.query(marker_model).filter(
        #         marker_model.id == marker_id).first().name
        #     Logger.instance().debug(
        #         "file: {}; line: {}; FASTA {} {}".format(__file__, inspect.currentframe().f_lineno, fasta_id, fasta_name))
        #     # file_id = fasta_obj.id
        #     sample_information_obj = session.query(sample_information_model).filter(sample_information_model.fasta_id==fasta_id).all()
        #     ############################################
        #     #
        #     # First (forward) trim: create_primer_tag_fasta_for_vsearch
        #     #
        #     # 1. Create fasta_path file for primer-tag sequences
        #     # 2. Run vsearch with 'db' parameter: primer_tag_fasta and 'usearch_global' parameter: fasta_path with the reads
        #     ############################################
        #     is_forward_strand = True
        #     #
        #     ############################################
        #     # Vsearch --db primer_tag_fasta --usearch_global merged_fasta
        #     # Vsearch --db primer_tag_fasta --usearch_global merged_fasta
        #     ############################################
        #     Logger.instance().debug("file: {}; line: {}; FASTA {} {}; forward {}".format(__file__,
        #                                      inspect.currentframe().f_lineno, fasta_id, fasta_name, is_forward_strand))
        #     Logger.instance().debug(
        #         "file: {}; line: {}; FASTA {} {}; forward {}; FASTA for forward trimming: {}".format(__file__,
        #                      inspect.currentframe().f_lineno, fasta_id, fasta_name, is_forward_strand, primer_tag_fasta))
        #     #
        #     # This create the primer + tag fasta_path file
        #     create_primer_tag_fasta_for_vsearch(sample_information_obj, is_forward_strand, primer_tag_fasta)
        #     Logger.instance().debug(
        #         "file: {}; line: {}; FASTA {} {}; forward {}; VSearch forward trimming".format(__file__,
        #                                     inspect.currentframe().f_lineno, fasta_id, fasta_name, is_forward_strand))
        #     vsearch_output_tsv = os.path.join(this_step_tmp_dir, os.path.basename(fasta_name), "vsearch_output_fwd.tsv")
        #     Logger.instance().debug(
        #         "file: {}; line: {}; FASTA {} {}; forward {}; vsearch_output_tsv: {}".format(__file__,
        #                 inspect.currentframe().f_lineno, fasta_id, fasta_name, is_forward_strand, vsearch_output_tsv))
        #     #
        #     ############################################
        #     # Run vsearch (Trim)
        #     # 
        #     # 1. Define vsearch parameters
        #     # 2. Run vsearch: output written to 'vsearch_output_tsv'
        #     ############################################
        #     # vsearch_params = {'db': primer_tag_fasta,
        #     #                   'usearch_global': fasta_name,
        #     #                   'id': str(self.option("min_id")),
        #     #                   'maxhits': 1,
        #     #                   'maxrejects': 0,
        #     #                   'maxaccepts': 0,
        #     #                   'minseqlength': str(self.option("minseqlength")),
        #     #                   'userfields': "query+target+tl+qilo+qihi+tilo+tihi+qrow",
        #     #                   'userout': vsearch_output_tsv,
        #     #                   }
        #     # vsearch1 = VSearchUsearchGlobal(**vsearch_params)
        #     # vsearch1.run()
        #     # del vsearch1
        #
        #     #
        #     # Create object and run vsearch
        #     vsearch_parameters = {'--db': primer_tag_fasta,
        #                                           '--usearch_global': fasta_name,
        #                                           '--id': str(self.option("min_id")),
        #                                           '--maxhits': 1,
        #                                           '--maxrejects': 0,
        #                                           '--maxaccepts': 0,
        #                                           '--minseqlength': str(self.option("minseqlength")),
        #                                           '--userfields': "query+target+tl+qilo+qihi+tilo+tihi+qrow",
        #                                           '--userout': vsearch_output_tsv,
        #                           }
        #     vsearch_cluster = VSearch(parameters=vsearch_parameters)
        #     vsearch_cluster.run()
        #
        #
        #     #
        #     ############################################
        #     # discard_tag_primer_alignment_with_low_sequence_quality
        #     ############################################
        #     Logger.instance().debug(
        #         "file: {}; line: {}; FASTA {} {}; forward {}; Eliminating non SRS conforms reads for forward trimming"
        #             .format(__file__, inspect.currentframe().f_lineno, fasta_id, fasta_name, is_forward_strand))
        #     discard_tag_primer_alignment_with_low_sequence_quality(vsearch_output_tsv, checked_vsearch_output_tsv, self.option("overhang"))
        #     #
        #     ############################################
        #     # Trim reads and write to sqlite
        #     ############################################
        #     Logger.instance().debug(
        #         "file: {}; line: {}; FASTA {} {}; forward {}; Trimming reads for forward trimming".format(__file__,
        #                                   inspect.currentframe().f_lineno, fasta_id, fasta_name, is_forward_strand))
        #     trimmed_tsv = os.path.join(this_step_tmp_dir, os.path.basename(fasta_name), 'trimmed_fwd.tsv')
        #     temp_db_sqlite = os.path.join(this_step_tmp_dir, os.path.basename(fasta_name), 'trimmed_fwd.sqlite')
        #     Logger.instance().debug(
        #         "file: {}; line: {}; FASTA {} {}; forward {}; Trimming reads for forward trimming: {}"
        #     .format(__file__, inspect.currentframe().f_lineno, fasta_id, fasta_name, is_forward_strand, temp_db_sqlite))
        #     trim_reads(checked_vsearch_output_tsv, fasta_name, trimmed_tsv, temp_db_sqlite)
        #     #
        #     ############################################
        #     # convert_trimmed_tsv_to_fasta
        #     ############################################
        #     Logger.instance().info("Writing fasta_path file for forward trimming.")
        #     trimmed_fasta = os.path.join(this_step_tmp_dir, os.path.basename(fasta_name), 'trimmed_fwd.fasta_path')
        #     Logger.instance().debug(
        #         "file: {}; line: {}; FASTA {} {}; forward {}; Writing fasta_path file for trimming.: {}".format(__file__,
        #                        inspect.currentframe().f_lineno, fasta_id, fasta_name, is_forward_strand, trimmed_fasta))
        #     convert_trimmed_tsv_to_fasta(trimmed_tsv, trimmed_fasta)
        #     #
        #     ############################################
        #     #
        #     # Second (reverse) trim
        #     #
        #     ############################################
        #     #
        #     is_forward_strand = False
        #     #
        #     ############################################
        #     # Vsearch --db primer_tag_fasta --usearch_global merged_fasta
        #     ############################################
        #     Logger.instance().debug(
        #         "file: {}; line: {}; FASTA {} {}; forward {}".format(__file__, inspect.currentframe().f_lineno,
        #                                                              fasta_id, fasta_name, is_forward_strand))
        #     Logger.instance().info("Creating a fasta_path query file to align on the merged reads fasta_path for forward trimming.")
        #     Logger.instance().debug(
        #         "file: {}; line: {}; FASTA {} {}; forward {}; FASTA for forward trimming: {}".format(__file__,
        #                   inspect.currentframe().f_lineno, fasta_id, fasta_name, is_forward_strand, primer_tag_fasta))
        #     create_primer_tag_fasta_for_vsearch(sample_information_obj, is_forward_strand, primer_tag_fasta)
        #     # Logger.instance().info("Processing Vsearch for forward trimming.")
        #     Logger.instance().debug(
        #         "file: {}; line: {}; FASTA {} {}; forward {}; VSearch forward trimming".format(__file__,
        #                                                                                     inspect.currentframe().f_lineno,
        #                                                                                        fasta_id, fasta_name,
        #                                                                                     is_forward_strand))
        #     #
        #     # self.vsearch_subprocess(merged_fasta, is_forward_strand, primer_tag_fasta, vsearch_output_tsv)
        #     vsearch_output_tsv = os.path.join(this_step_tmp_dir, os.path.basename(fasta_name), "vsearch_output_rev.tsv")
        #
        #     Logger.instance().debug(
        #         "file: {}; line: {}; FASTA {} {}; forward {}; vsearch_output_tsv".format(__file__,
        #                                                                               inspect.currentframe().f_lineno,
        #                                                                                  fasta_id, fasta_name, is_forward_strand))
        #     # vsearch_params = {'db': primer_tag_fasta,
        #     #                   'usearch_global': trimmed_fasta,
        #     #                   'id': str(self.option("min_id")),
        #     #                   'maxhits': 1,
        #     #                   'maxrejects': 0,
        #     #                   'maxaccepts': 0,
        #     #                   'minseqlength': str(self.option("minseqlength")),
        #     #                   'userfields': "query+target+tl+qilo+qihi+tilo+tihi+qrow",
        #     #                   'userout': vsearch_output_tsv,
        #     #                   }
        #     # vsearch1 = VSearchUsearchGlobal(**vsearch_params)
        #     # vsearch1.run()
        #     # del vsearch1
        #
        #     #
        #     # Create object and run vsearch
        #     vsearch_parameters = {'--db': primer_tag_fasta,
        #                                           '--usearch_global': trimmed_fasta,
        #                                           '--id': str(self.option("min_id")),
        #                                           '--maxhits': 1,
        #                                           '--maxrejects': 0,
        #                                           '--maxaccepts': 0,
        #                                           '--minseqlength': str(self.option("minseqlength")),
        #                                           '--userfields': "query+target+tl+qilo+qihi+tilo+tihi+qrow",
        #                                           '--userout': vsearch_output_tsv,
        #                           }
        #     vsearch_cluster = VSearch(parameters=vsearch_parameters)
        #     vsearch_cluster.run()
        #     #
        #     ############################################
        #     # discard_tag_primer_alignment_with_low_sequence_quality
        #     ############################################
        #     Logger.instance().debug(
        #         "file: {}; line: {}; FASTA {} {}; forward {}; Eliminating non SRS conforms reads for forward trimming".format(
        #             __file__, inspect.currentframe().f_lineno, fasta_id, fasta_name, is_forward_strand))
        #     discard_tag_primer_alignment_with_low_sequence_quality(vsearch_output_tsv, checked_vsearch_output_tsv,
        #                                                            self.option("overhang"))
        #     #
        #     ############################################
        #     # Trim reads in reverse strand and write to sqlite
        #     ############################################
        #     Logger.instance().debug(
        #         "file: {}; line: {}; FASTA {} {}; forward {}; Trimming reads for reverse trimming".format(__file__,
        #                                                                                                inspect.currentframe().f_lineno,
        #                                                                                                   fasta_id,
        #                                                                                                   fasta_name,
        #                                                                                                is_forward_strand))
        #     trimmed_tsv = os.path.join(this_step_tmp_dir, os.path.basename(fasta_name), 'trimmed_rev.tsv')
        #     temp_db_sqlite = os.path.join(this_step_tmp_dir, os.path.basename(fasta_name), 'trimmed_rev.sqlite')
        #     Logger.instance().debug(
        #         "file: {}; line: {}; FASTA {} {}; forward {}; Trimming reads for reverse trimming: {}".format(__file__,
        #                                                                                                    inspect.currentframe().f_lineno,
        #                                                                                                       fasta_id,
        #                                                                                                       fasta_name,
        #                                                                                                    is_forward_strand,
        #                                                                                                    temp_db_sqlite))
        #     trim_reads(checked_vsearch_output_tsv, trimmed_fasta, trimmed_tsv, temp_db_sqlite)
        #     #
        #     Logger.instance().info("Annotating reads with Sample Information.")
        #     # One TSV file with read annotation per merged FASTA Fasta
        #     read_annotation_tsv = os.path.join(this_step_tmp_dir, os.path.basename(fasta_name), 'read_annotation.tsv')
        #     tsv_file_list_with_read_annotations.append(read_annotation_tsv)
        #     # run_list[read_annotation_tsv] = fasta_obj.run_id
        #     Logger.instance().debug(
        #         "file: {}; line: {}; trimmed_tsv {}".format(__file__, inspect.currentframe().f_lineno, trimmed_tsv))
        #     ################################################################
        #     # Count number of times a given variant (Sequence unique) is observed in fasta_path file
        #     # In other words, store in 'VariantReadCount' for each 'variant_id' -> 'read count'
        #     ################################################################
        #     fasta_sort_reads_tsv = os.path.join(this_step_tmp_dir, os.path.basename(fasta_name), 'sortreads.tsv')
        #     fasta_sort_reads_tsv_list.append(fasta_sort_reads_tsv)
        #     annotate_reads(session, sample_information_model, trimmed_tsv, fasta_id=fasta_id, out_tsv=fasta_sort_reads_tsv)
        #
        # ########################################################
        # #
        # # Concatenate sortreads files of different fasta_path files
        # #
        # ########################################################
        # with open(sort_reads_tsv, 'wb') as wfd:
        #     for f in fasta_sort_reads_tsv_list:
        #         with open(f, 'rb') as fd:
        #             shutil.copyfileobj(fd, wfd)
