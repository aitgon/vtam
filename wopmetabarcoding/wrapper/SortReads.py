import tempfile

import os
from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger

from wopmetabarcoding.utils.VSearch import VSearch1
from wopmetabarcoding.wrapper.SortReadsUtilities import \
    create_primer_tag_fasta_for_vsearch, check_criteria_in_vsearch_output,  trim_reads, \
    convert_trimmed_tsv_to_fasta, annotate_reads, gather_files, count_reads, insert_variant

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
            "min_id": "float",
            "minseqlength": "int",
            "overhang": "int"

        }

    # def vsearch_subprocess(self, file, fasta_db, vsearch_output_tsv):
    #     """
    #
    #     :param file:
    #     :param fasta_db:
    #     :param vsearch_output_tsv:
    #     :return:
    #     """
    #     # subprocess.call(
    #     #     "vsearch --usearch_global " + file + " --db " + fasta_db + " --id " + str(self.option("min_id")) +
    #     #     " --maxhits 1 --maxrejects 0 --maxaccepts 0 --minseqlength " + str(self.option("minseqlength")) +
    #     #     " --userout " + vsearch_output_tsv + " --userfields query+target+tl+qilo+qihi+tilo+tihi+qrow", shell=True
    #     # )
    #     vsearch_params = {'db': fasta_db,
    #                       'usearch_global': file,
    #                       'id0': str(self.option("min_id")),
    #                       'maxhits': 1,
    #                       'maxrejects': 0,
    #                       'maxaccepts': 0,
    #                       'minseqlength': str(self.option("minseqlength")),
    #                       'userfields': "query+target+tl+qilo+qihi+tilo+tihi+qrow",
    #                       'userout': vsearch_output_tsv,
    #                       }
    #     vsearch1 = VSearch1(**vsearch_params)
    #     vsearch1.run()

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        conn = engine.connect()
        #
        # Input tables models
        sample_information_model = self.input_table(SortReads.__input_table_sample_information)
        file_model = self.input_table(SortReads.__input_table_file)
        marker_model = self.input_table(SortReads.__input_table_marker)

        # Temp variables:
        annoted_tsv_list = []
        run_list = {}
        vsearch_output_tsv = os.path.join(tempdir, 'vsearch_output.tsv')
        primer_tag_fasta = os.path.join(tempdir, 'primer_tag.fasta')
        checked_vsearch_output_tsv = os.path.join(tempdir, 'checked_vsearch_output.tsv')

        # Output file models
        sortreads_samplecount = self.output_file(SortReads.__output_sortreads_samplecount_csv)

        # Output tables models
        variant_model = self.output_table(SortReads.__output_table_variant)

        #  TODO: Later we will see in case files are already trimmed
        # for file_obj in session.query(file_model).all():
        #     if file_obj.trimmed_status is "1":
        #         raise Exception('Read trimming cannot be skipped')
        #         # read_counter(session, file_obj.name, readcount_model)
        #         # fasta_writer(session, readcount_model, file_obj.name)
        # else:
        for file_obj in session.query(file_model).all():
            file_id = file_obj.id
            sample_information_obj = session.query(sample_information_model).filter(sample_information_model.file_id==file_id).all()
            # 
            # First: Forward trim
            is_forward_strand = True
            Logger.instance().info("Creating a fasta query file to align on the merged reads fasta for forward trimming.")
            create_primer_tag_fasta_for_vsearch(sample_information_obj, is_forward_strand, primer_tag_fasta)
            Logger.instance().info("Processing Vsearch for forward trimming.")
            #
            # self.vsearch_subprocess(file_obj.name, primer_tag_fasta, vsearch_output_tsv)
            vsearch_params = {'db': primer_tag_fasta,
                              'usearch_global': file_obj.name,
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
            Logger.instance().info("Eliminating non SRS conforms reads for forward trimming.")
            check_criteria_in_vsearch_output(vsearch_output_tsv, checked_vsearch_output_tsv, self.option("overhang"))
            Logger.instance().info("Trimming reads for forward trimming.")
            trimmed_tsv = os.path.join(tempdir, os.path.basename(file_obj.name).replace('.fasta', '_forward_trimmed.tsv'))
            trim_reads(checked_vsearch_output_tsv, file_obj.name, trimmed_tsv, is_forward_strand, tempdir)
            trimmed_fasta = trimmed_tsv.replace('.tsv', '.fasta')
            Logger.instance().info("Writing fasta file for forward trimming.")
            convert_trimmed_tsv_to_fasta(trimmed_tsv, trimmed_fasta)
            # 
            # Second: Reverse trim
            is_forward_strand = False
            Logger.instance().info("Creating a fasta query file to align on the merged reads fasta for reverse trimming.")
            # create_primer_tag_fasta_for_vsearch(session, sample_information_model,primer_tag_fasta, file_obj.name,is_forward_strand)
            create_primer_tag_fasta_for_vsearch(sample_information_obj, is_forward_strand, primer_tag_fasta)
            Logger.instance().info("Processing Vsearch for reverse trimming.")
            #
            # self.vsearch_subprocess(trimmed_fasta, primer_tag_fasta, vsearch_output_tsv)
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
            Logger.instance().info("Eliminating non SRS conforms reads for reverse trimming.")
            check_criteria_in_vsearch_output(vsearch_output_tsv, checked_vsearch_output_tsv, self.option("overhang"))
            Logger.instance().info("Trimming reads for reverse trimming.")
            trimmed_tsv = trimmed_fasta.replace('_forward_trimmed.fasta', '_reverse_trimmed.tsv')
            trim_reads(checked_vsearch_output_tsv, trimmed_fasta, trimmed_tsv, is_forward_strand, tempdir)
            trimmed_fasta = trimmed_tsv.replace('.tsv', '.fasta')
            # for file_obj in session.query(file_model).all():
            Logger.instance().info("Annotating reads with Sample Information.")
            # annotated_reads_tsv = (file_obj.name).replace(".fasta", "_annotated_reads.tsv")
            annotated_reads_tsv = os.path.join(tempdir, os.path.basename(file_obj.name).replace('.fasta', '_annotated_reads.tsv'))
            annoted_tsv_list.append(annotated_reads_tsv)
            run_list[annotated_reads_tsv] = file_obj.run_name
            annotate_reads(session, sample_information_model, trimmed_tsv,
                           file_id=file_id, annotated_reads_tsv=annotated_reads_tsv)
        #
        # # For each marker_id, concatenate its files with the annotated reads and count unique reads per marker_id
        with open(sortreads_samplecount, 'w') as fout_sortread_samplecount:
            for marker_obj in session.query(marker_model).all():
                marker_name = marker_obj.name
                marker_id = marker_obj.id
                gathered_marker_file = os.path.join(tempdir, marker_name + "_file.tsv")
                sample_count_tsv = os.path.join(tempdir, marker_name + "_sample_count.tsv")
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
                insert_variant(session, count_reads_marker, variant_model)
                session.commit()
