from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
from wopmetabarcoding.wrapper.VsearchSortReads_functions import \
    create_primer_tag_fasta_for_vsearch, check_criteria_in_vsearch_output, read_counter, fasta_writer, trim_reads, \
    convert_trimmed_tsv_to_fasta, annotate_reads, gather_files, count_reads, insert_variant

import subprocess



class VsearchSortReads(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.VsearchSortReads"
    }
    # Input
    # Input table
    __input_table_fileinformation = "SampleInformation"
    __input_table_file = "File"
    __input_table_marker = "Marker"
    # Output
    # Output file
    # __output_file_vsearch_output_tsv = "vsearch_output_tsv"
    # __output_file_primer_tag_fasta = "primer_tag_fasta"
    # __output_file_checked_vsearch_output_tsv = "checked_vsearch_output_tsv"
    # Output table
    __output_table_readcount = "ReadCount"
    # __output_table_obifasta = 'ObiFasta'
    __output_table_variant = 'Variant'

    def specify_input_table(self):
        return [
            VsearchSortReads.__input_table_fileinformation,
            VsearchSortReads.__input_table_file,
            VsearchSortReads.__input_table_marker
        ]

    # def specify_output_file(self):
    #     return [
    #         VsearchSortReads.__output_file_vsearch_output_tsv,
    #         VsearchSortReads.__output_file_primer_tag_fasta,
    #         VsearchSortReads.__output_file_checked_vsearch_output_tsv
    #     ]

    def specify_output_table(self):
        return [
            VsearchSortReads.__output_table_readcount,
            VsearchSortReads.__output_table_variant
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

    def vsearch_subprocess(self, file, fasta_db, vsearch_output_tsv):
        """

        :param file:
        :param fasta_db:
        :param vsearch_output_tsv:
        :return:
        """
        subprocess.call(
            "vsearch --usearch_global " + file + " --db " + fasta_db + " --id " + str(self.option("min_id")) +
            " --maxhits 1 --maxrejects 0 --maxaccepts 0 --minseqlength " + str(self.option("minseqlength")) +
            " --userout " + vsearch_output_tsv + " --userfields query+target+tl+qilo+qihi+tilo+tihi+qrow", shell=True
        )


    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        conn = engine.connect()

        # Input tables models
        file_information_model = self.input_table(VsearchSortReads.__input_table_fileinformation)
        file_model = self.input_table(VsearchSortReads.__input_table_file)
        marker_model = self.input_table(VsearchSortReads.__input_table_marker)

        # # Output file
        # vsearch_output_tsv = self.output_file(VsearchSortReads.__output_file_vsearch_output_tsv)
        # primer_tag_fasta = self.output_file(VsearchSortReads.__output_file_primer_tag_fasta)
        # checked_vsearch_output_tsv = self.output_file(VsearchSortReads.__output_file_checked_vsearch_output_tsv)

        # Temp variables:
        annoted_tsv_list = []
        run_list = {}
        vsearch_output_tsv = 'data/output/vsearch_output.tsv'
        primer_tag_fasta = 'data/output/primer_tag.fasta'
        checked_vsearch_output_tsv = 'data/output/checked_vsearch_output.tsv'

        # Output tables models
        readcount_model = self.output_table(VsearchSortReads.__output_table_readcount)
        variant_model = self.output_table(VsearchSortReads.__output_table_variant)

        #  TODO: Later we will see in case files are already trimmed
        # for file_obj in session.query(file_model).all():
        #     if file_obj.trimmed_status is "1":
        #         raise Exception('Read trimming cannot be skipped')
        #         # read_counter(session, file_obj.name, readcount_model)
        #         # fasta_writer(session, readcount_model, file_obj.name)
        # else:
        for file_obj in session.query(file_model).all():
            # 
            # First: Forward trim
            is_forward_strand = True
            Logger.instance().info("Creating a fasta query file to align on the merged reads fasta for forward trimming.")
            create_primer_tag_fasta_for_vsearch(session, file_information_model, primer_tag_fasta, file_obj.name,is_forward_strand)
            Logger.instance().info("Processing Vsearch for forward trimming.")
            self.vsearch_subprocess(file_obj.name, primer_tag_fasta, vsearch_output_tsv)
            Logger.instance().info("Eliminating non SRS conforms reads for forward trimming.")
            check_criteria_in_vsearch_output(vsearch_output_tsv, checked_vsearch_output_tsv, self.option("overhang"))
            Logger.instance().info("Trimming reads for forward trimming.")
            trimmed_tsv = (file_obj.name).replace('.fasta', '_forward_trimmed.tsv')
            trim_reads(checked_vsearch_output_tsv, file_obj.name, trimmed_tsv, is_forward_strand)
            trimmed_fasta = trimmed_tsv.replace('.tsv', '.fasta')
            Logger.instance().info("Writing fasta file for forward trimming.")
            convert_trimmed_tsv_to_fasta(trimmed_tsv, trimmed_fasta)
            # 
            # Second: Reverse trim
            is_forward_strand = False
            Logger.instance().info("Creating a fasta query file to align on the merged reads fasta for reverse trimming.")
            create_primer_tag_fasta_for_vsearch(session, file_information_model,primer_tag_fasta, file_obj.name,is_forward_strand)
            Logger.instance().info("Processing Vsearch for reverse trimming.")
            self.vsearch_subprocess(trimmed_fasta, primer_tag_fasta, vsearch_output_tsv)
            Logger.instance().info("Eliminating non SRS conforms reads for reverse trimming.")
            check_criteria_in_vsearch_output(vsearch_output_tsv, checked_vsearch_output_tsv, self.option("overhang"))
            Logger.instance().info("Trimming reads for reverse trimming.")
            trimmed_tsv = trimmed_fasta.replace('_forward_trimmed.fasta', '_reverse_trimmed.tsv')
            trim_reads(checked_vsearch_output_tsv, trimmed_fasta, trimmed_tsv, is_forward_strand)
            trimmed_fasta = trimmed_tsv.replace('.tsv', '.fasta')
            # for file_obj in session.query(file_model).all():
            Logger.instance().info("Annotating reads with Sample Information.")
            annotated_reads_tsv = (file_obj.name).replace(".fasta", "_annotated_reads.tsv")
            annoted_tsv_list.append(annotated_reads_tsv)
            run_list[annotated_reads_tsv] = file_obj.run_name
            annotate_reads(session, file_information_model, trimmed_tsv,
                           merged_fasta_file_name=file_obj.name, annotated_reads_tsv=annotated_reads_tsv)
        #
        # # For each marker, concatenate its files with the annotated reads and count unique reads per marker
        for marker_obj in session.query(marker_model).all():
            marker_name = marker_obj.marker_name
            gathered_marker_file = marker_name + "_file.tsv"
            count_reads_marker = gathered_marker_file.replace(".tsv", ".sqlite")
            Logger.instance().info("Gathering all files from annotated files the same marker into one.")
            gather_files(marker_name, gathered_marker_file, annoted_tsv_list)
            # session.commit()
            # read_count_per_marker_sqlite = conn_name = marker_obj.marker_name + ".sqlite"
            Logger.instance().info("Counting reads for each marker.")
            count_reads(gathered_marker_file, count_reads_marker)
            # session.commit()
            Logger.instance().info("Inserting variant in the Variant table of the database.")
            insert_variant(session, count_reads_marker, variant_model)
            # session.commit()
