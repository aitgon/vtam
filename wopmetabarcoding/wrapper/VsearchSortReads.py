from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmetabarcoding.wrapper.VsearchSortReads_functions import \
    create_fastadb, dereplicate, read_counter, fasta_writer, insert_read, trim_reads, reverse_complement, \
    create_fasta, attribute_combination, singleton_removing, variant_table

import subprocess


class VsearchSortReads(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.VsearchSortReads"
    }
    __input_table_fileinformation = "FileInformation"
    __input_table_file = "File"
    __output_file_csvsearch = "csvsearch"
    __output_file_fastadb = "fastadb"
    __output_table_readcount = "ReadCount"
    __output_table_outputcsv = "OutputCsv"
    __output_table_reversetrimmed = "ReverseTrimmed"
    __output_file_reversefastadb = "reversefastadb"
    __output_file_fowardtrimmed = "fowardtrimmed"
    __output_file_outputforward = 'outputforward'
    __output_table_obifasta = 'ObiFasta'
    __output_table_variant = 'Variant'

    def specify_input_table(self):
        return [
            VsearchSortReads.__input_table_fileinformation,
            VsearchSortReads.__input_table_file
        ]

    def specify_output_file(self):
        return [
            VsearchSortReads.__output_file_csvsearch,
            VsearchSortReads.__output_file_fastadb,
            VsearchSortReads.__output_file_reversefastadb,
            VsearchSortReads.__output_file_fowardtrimmed,
            VsearchSortReads.__output_file_outputforward
        ]

    def specify_output_table(self):
        return [
            VsearchSortReads.__output_table_readcount,
            VsearchSortReads.__output_table_outputcsv,
            VsearchSortReads.__output_table_reversetrimmed,
            VsearchSortReads.__output_table_obifasta,
            VsearchSortReads.__output_table_variant
        ]

    def specify_params(self):
        """

        :return:
        """
        return {
            "min_id": "float",
            "minseqlength": "int"
        }

    def vsearch_sr(self, file, fasta_db, output_csv):
        """

        :param file:
        :param fasta_db:
        :param output_csv:
        :return:
        """
        subprocess.call(
            "vsearch --usearch_global " + file + " --db " + fasta_db + " --id " + str(self.option("min_id")) +
            " --maxhits 1 --maxrejects 0 --maxaccepts 0 --minseqlength " + str(self.option("minseqlength")) +
            " --userout " + output_csv + " --userfields query+target+tl+qilo+qihi+tilo+tihi+qrow", shell=True
        )

    def algo_trim(self, session, output_csv, csv_file, database, fasta_file):
        """

        :param session:
        :param output_csv:
        :param csv_file:
        :param database:
        :param fasta_file:
        :return:
        """
        self.vsearch_sr(fasta_file, database, csv_file)
        print('ok')
        dereplicate(session, output_csv, csv_file)
        print('ok1')
        insert_read(session, output_csv, fasta_file)
        print('ok2')
        trim_reads(session, output_csv, fasta_file)

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        conn = engine.connect()

        # Output file
        csv_output = self.output_file(VsearchSortReads.__output_file_csvsearch)
        fasta_db = self.output_file(VsearchSortReads.__output_file_fastadb)
        reverse_fasta_db = self.output_file(VsearchSortReads.__output_file_reversefastadb)
        foward_trimmed = self.output_file(VsearchSortReads.__output_file_fowardtrimmed)
        outputforward = self.output_file(VsearchSortReads.__output_file_outputforward)

        # Table model
        file_information_model = self.input_table(VsearchSortReads.__input_table_fileinformation)
        file_model = self.input_table(VsearchSortReads.__input_table_file)
        readcount_model = self.output_table(VsearchSortReads.__output_table_readcount)
        outputcsv_model = self.output_table(VsearchSortReads.__output_table_outputcsv)
        reversetrimmed_model = self.output_table(VsearchSortReads.__output_table_reversetrimmed)
        obifasta_model = self.output_table(VsearchSortReads.__output_table_obifasta)
        variant_model = self.output_table(VsearchSortReads.__output_table_variant)

        for element in session.query(file_model).all():
            if element.dereplicate_status is "1":
                read_counter(session, element.file_name, readcount_model)
                fasta_writer(session, readcount_model, element.file_name)
            else:
                create_fastadb(session, file_information_model, fasta_db, reverse_fasta_db)
        #         self.algo_trim(session, outputcsv_model, csv_output, fasta_db, element.file_name)
        # reverse_complement(session, outputcsv_model)
        # create_fasta(session, outputcsv_model, foward_trimmed)
        # self.algo_trim(session, reversetrimmed_model, csv_output, reverse_fasta_db, foward_trimmed)
        # reverse_complement(session, reversetrimmed_model)
        # create_fasta(session, reversetrimmed_model, outputforward)
        # attribute_combination(session, file_information_model, outputcsv_model,  obifasta_model, outputforward)
        # singleton_removing(session, obifasta_model)
        # variant_table(session, obifasta_model, variant_model)
        # session.commit()


