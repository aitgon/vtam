from wopmars.framework.database.tables.ToolWrapper import ToolWrapper

from vtam.utils.utilities import get_or_create
from vtam.utils.OptionManager import OptionManager


import os

class SampleInformation(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.SampleInformation"
    }
    __input_file_csv = "fastainfo"
    #
    __output_table_biosample = "Biosample"
    __output_table_fasta = "Fasta"
    __output_table_marker = "Marker"
    __output_table_primerpair = "PrimerPair"
    __output_table_replicate = "Replicate"
    __output_table_run = "Run"
    __output_table_sample_information = "SampleInformation"
    __output_table_tagpair = "TagPair"

    def specify_input_file(self):
        return [SampleInformation.__input_file_csv]

    def specify_output_table(self):
        return [
            SampleInformation.__output_table_biosample,
            SampleInformation.__output_table_fasta,
            SampleInformation.__output_table_marker,
            SampleInformation.__output_table_primerpair,
            SampleInformation.__output_table_replicate,
            SampleInformation.__output_table_run,
            SampleInformation.__output_table_sample_information,
            SampleInformation.__output_table_tagpair,
        ]

    def specify_params(self):
        return{
            "fasta_dir": "str",
            "log_verbosity": "int",
            "log_file": "str",
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        OptionManager.instance()['log_verbosity'] = int(self.option("log_verbosity"))
        if not self.option("log_verbosity") is None:
            OptionManager.instance()['log_file'] = str(self.option("log_file"))

        ##########################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        ##########################################################

        # input file paths
        csv_path = self.input_file(SampleInformation.__input_file_csv)
        #
        # Models
        biosample_model = self.output_table(SampleInformation.__output_table_biosample)
        fasta_model = self.output_table(SampleInformation.__output_table_fasta)
        marker_model = self.output_table(SampleInformation.__output_table_marker)
        primerpair_model = self.output_table(SampleInformation.__output_table_primerpair)
        replicate_model = self.output_table(SampleInformation.__output_table_replicate)
        run_model = self.output_table(SampleInformation.__output_table_run)
        sampleinformation_model = self.output_table(SampleInformation.__output_table_sample_information)
        tagpair_model = self.output_table(SampleInformation.__output_table_tagpair)

        with open(csv_path, 'r') as fin:
            next(fin) # skip header
            for line in fin:
                tag_forward = line.strip().split('\t')[0]
                primer_forward = line.strip().split('\t')[1]
                tag_reverse = line.strip().split('\t')[2]
                primer_reverse = line.strip().split('\t')[3]
                marker_name = line.strip().split('\t')[4]
                biosample_name = line.strip().split('\t')[5]
                replicate_name = line.strip().split('\t')[6]
                if len(line.strip().split('\t')) == 9: # No fastq file columns
                    file_name = os.path.join(self.option("fasta_dir"), line.strip().split('\t')[8])
                elif len(line.strip().split('\t')) == 11: # Yes fastq file columns
                    file_name = os.path.join(self.option("fasta_dir"), line.strip().split('\t')[10])
                else:
                    raise IndexError("{} error. Verify nb of columns in this input file: {}".format(self.__class__.__name__, os.path.join(os.getcwd(), csv_path)))
                if not os.path.isfile(os.path.join(os.getcwd(), file_name)):
                    raise FileNotFoundError("{} error. Verify this file output: {}".format(self.__class__.__name__, os.path.join(os.getcwd(), file_name)))
                run_name = line.strip().split('\t')[7]
                #
                # Insert run
                run_obj = {'name': run_name}
                run_instance = get_or_create(session, run_model, **run_obj)
                run_id = run_instance.id
                #
                # Insert marker_id
                marker_obj = {'name': marker_name}
                marker_instance = get_or_create(session, marker_model, **marker_obj)
                marker_id = marker_instance.id
                #
                # Insert Biosamples
                biosample_obj = {'name': biosample_name}
                biosample_instance = get_or_create(session, biosample_model, **biosample_obj)
                biosample_id = biosample_instance.id
                #
                # Insert replicate
                replicate_obj = {'name': replicate_name}
                replicate_instance = get_or_create(session, replicate_model, **replicate_obj)
                replicate_id = replicate_instance.id
                #
                # Insert file output
                is_trimmed = False # Default
                fasta_obj = {'name': file_name, 'is_trimmed': is_trimmed}
                fasta_instance = get_or_create(session, fasta_model, **fasta_obj)
                fasta_id = fasta_instance.id

                # Insert sample_information
                sample_information_obj = {'tag_forward': tag_forward, 'primer_forward': primer_forward, 'tag_reverse': tag_reverse, 'primer_reverse': primer_reverse, 'biosample_id': biosample_id, 'replicate_id': replicate_id, 'run_id': run_id}
                sample_information_obj['fasta_id'] = fasta_id
                sample_information_obj['marker_id'] = marker_id
                get_or_create(session, sampleinformation_model, **sample_information_obj)

                # Primer pair
                primerpair_obj = {'primer_forward': primer_forward, 'primer_reverse': primer_reverse}
                get_or_create(session, primerpair_model, **primerpair_obj)

                # Tag pair
                tagpair_obj = {'tag_forward': tag_forward, 'tag_reverse': tag_reverse}
                get_or_create(session, tagpair_model, **tagpair_obj)

