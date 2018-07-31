from wopmars.framework.database.tables.ToolWrapper import ToolWrapper

from wopmetabarcoding.utils.utilities import get_or_create

import os

class FastaInformation(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.FastaInformation"
    }
    __input_file_csv = "sample2fasta"
    #
    __output_table_fasta = "Fasta"
    __output_table_marker = "Marker"
    __output_table_primerpair = "PrimerPair"
    __output_table_biosample = "Biosample"
    __output_table_tagpair = "TagPair"
    __output_table_replicate = "Replicate"
    __output_table_sample_information = "SampleInformation"

    def specify_input_file(self):
        return [FastaInformation.__input_file_csv]

    def specify_output_table(self):
        return [
            FastaInformation.__output_table_sample_information,
            FastaInformation.__output_table_fasta,
            FastaInformation.__output_table_marker,
            FastaInformation.__output_table_primerpair,
            FastaInformation.__output_table_biosample,
            FastaInformation.__output_table_tagpair,
            FastaInformation.__output_table_replicate
        ]

    def specify_params(self):
        return{
            "fasta_dir": "str"
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        # conn = engine.connect()
        #
        # input file paths
        csv_path = self.input_file(FastaInformation.__input_file_csv)
        #
        # Input models
        marker_model = self.output_table(FastaInformation.__output_table_marker)
        fasta_model = self.output_table(FastaInformation.__output_table_fasta)
        sampleinformation_model = self.output_table(FastaInformation.__output_table_sample_information)
        replicate_model = self.output_table(FastaInformation.__output_table_replicate)
        biosample_model = self.output_table(FastaInformation.__output_table_biosample)
        primerpair_model = self.output_table(FastaInformation.__output_table_primerpair)
        tagpair_model = self.output_table(FastaInformation.__output_table_tagpair)

        with open(csv_path, 'r') as fin:
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
                    # TODO: Add it to a exception logger
                    raise FileNotFoundError("{} error. Verify this file path: {}".format(self.__class__.__name__, os.path.join(os.getcwd(), file_name)))
                run_name = line.strip().split('\t')[7]
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
                biosample_name = biosample_instance.name
                #
                # Insert replicate
                replicate_obj = {'biosample_id': biosample_id, 'marker_id': marker_id, 'file_name': file_name, 'name': replicate_name}
                replicate_instance = get_or_create(session, replicate_model, **replicate_obj)
                replicate_id = replicate_instance.id
                #
                # Insert file path
                is_trimmed = False # Default
                fasta_obj = {'name': file_name, 'run_name': run_name, 'is_trimmed': is_trimmed}
                fasta_instance = get_or_create(session, fasta_model, **fasta_obj)
                fasta_id = fasta_instance.id

                # Insert sample_information
                sample_information_obj = {'tag_forward': tag_forward, 'primer_forward': primer_forward, 'tag_reverse': tag_reverse, 'primer_reverse': primer_reverse, 'biosample_id': biosample_id, 'replicate_id': replicate_id, 'run_name': run_name}
                sample_information_obj['fasta_id'] = fasta_id
                sample_information_obj['marker_id'] = marker_id
                get_or_create(session, sampleinformation_model, **sample_information_obj)

                # Primer pair
                primerpair_obj = {'primer_forward': primer_forward, 'primer_reverse': primer_reverse}
                get_or_create(session, primerpair_model, **primerpair_obj)

                # Tag pair
                tagpair_obj = {'tag_forward': tag_forward, 'tag_reverse': tag_reverse}
                get_or_create(session, tagpair_model, **tagpair_obj)

