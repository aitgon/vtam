import os
import pathlib
import shutil

import pandas
from wopmars.models.ToolWrapper import ToolWrapper

from vtam.utils.FastaInformation import FastaInformation
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
            "min_id": "float",
            "minseqlength": "int",
            "overhang": "int",
        }

    def run(self):
        session = self.session
        engine = session._session().get_bind()
        fasta_dir = str(os.getenv('VTAM_FASTA_DIR'))

        this_temp_dir = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__))
        pathlib.Path(this_temp_dir).mkdir(exist_ok=True)

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
        #
        # Input files
        input_file_fastainfo = self.input_file(SortReads.__input_file_fastainfo)
        #
        # Output files
        # sortreads with columns: read_name, fasta_id, run_id, marker_id, biosample_id, replicate
        sort_reads_tsv = self.output_file(SortReads.__output_file_sort_reads)
        #
        # Options
        min_id = str(self.option("min_id"))
        minseqlength = str(self.option("minseqlength"))
        overhang = int(self.option("overhang"))

        ##########################################################
        #
        # 1. Read fastainfo to get run_id, marker_id, biosample_id, replicate for current analysis
        #
        ##########################################################

        fasta_information_obj = FastaInformation(input_file_fastainfo, engine, run_model, marker_model, biosample_model)
        sample_information_df = fasta_information_obj.get_sample_information_df(add_tag_primer_fasta=True)
        # sample_information_df = fasta_information_obj.get_fasta_information_record_list(
        #     tag_primer_fasta_information=True)
        sample_information_df = pandas.DataFrame(data=sample_information_df)

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

        for fasta_file_name in sorted(sample_information_df.fasta_file_name.unique().tolist()):
            fasta_information_per_fasta_df = sample_information_df.loc[sample_information_df.fasta_file_name==fasta_file_name]
            fasta_file_id = session.query(fasta_model.id).filter(fasta_model.name == fasta_file_name).one()[0]
            fasta_path = os.path.join(fasta_dir, fasta_file_name)

            alignement_parameters = {'min_id': min_id, 'minseqlength': minseqlength, 'overhang': overhang}

            sort_reads_tsv_i = os.path.join(this_temp_dir, "sort_reads_fasta_id_{}.tsv".format(fasta_file_id))

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
            fout.write('read_id\trun_id\tfasta_id\tmarker_id\tbiosample_id\treplicate\tread_sequence\n')
        with open(sort_reads_tsv, 'ab') as wfd:
            for f in sort_reads_tsv_list:
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
