from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger

from wopmetabarcoding.model import FileAG
from wopmetabarcoding.utils.utilities import insert_table, get_or_create
from wopmetabarcoding.wrapper.SampleInformationUtilities import \
    insert_marker, insert_fileinformation, insert_file, insert_tagpair, insert_primer, insert_sample, insert_replicate,\
    insert_replicate_marker

# import models
import wopmetabarcoding.model.Marker
import wopmetabarcoding.model.PrimerPair


class SampleInformation(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.SampleInformation"
    }
    __input_file_csv = "csv"
    #
    __output_table_file = "File"
    __output_table_fileag = "FileAG"
    __output_table_sampleinformation = "SampleInformation"
    __output_table_marker = "Marker"
    __output_table_primerpair = "PrimerPair"
    __output_table_biosample = "Biosample"
    __output_table_tagpair = "TagPair"
    __output_table_replicate = "Replicate"
    __output_table_replicatemarker = "Replicatemarker"

    def specify_input_file(self):
        return [SampleInformation.__input_file_csv]

    def specify_output_table(self):
        return [
            SampleInformation.__output_table_file,
            SampleInformation.__output_table_fileag,
            SampleInformation.__output_table_sampleinformation,
            SampleInformation.__output_table_marker,
            SampleInformation.__output_table_primerpair,
            SampleInformation.__output_table_biosample,
            SampleInformation.__output_table_tagpair,
            SampleInformation.__output_table_replicate,
            SampleInformation.__output_table_replicatemarker
        ]

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        # file paths
        csv_path = self.input_file(SampleInformation.__input_file_csv)
        # Models
        marker_model = self.output_table(SampleInformation.__output_table_marker)
        # primerpair_model = self.output_table(SampleInformation.__output_table_primerpair)
        # tagpair_model = self.output_table(SampleInformation.__output_table_tagpair)
        file_model = self.output_table(SampleInformation.__output_table_file)
        # biosample_model = self.output_table(SampleInformation.__output_table_biosample)
        sampleinformation_model = self.output_table(SampleInformation.__output_table_sampleinformation)
        replicate_model = self.output_table(SampleInformation.__output_table_replicate)
        replicatemarker_model = self.output_table(SampleInformation.__output_table_replicatemarker)
        fileag_model = self.output_table(SampleInformation.__output_table_fileag)
        with open(csv_path, 'r') as fin:
            next(fin)  # skip header
            for line in fin:
                file_name = line.strip().split(',')[7]
                run_name = line.strip().split(',')[8]
                # AG
                is_trimmed = False # Default
                fileag_obj = {'name': file_name, 'run_name': run_name, 'is_trimmed': is_trimmed}
                instance = get_or_create(session, fileag_model, **fileag_obj)
                file_id = instance.id

                sample_information_ag_obj_keys = ['tag_forward', 'primer_forward', 'tag_reverse', 'primer_reverse', 'marker_name', 'sample_name', 'replicate_name', 'file_name', 'run_name']
                sample_information_ag_obj_vals = line.strip().split(',')
                sample_information_ag_obj = dict(zip(sample_information_ag_obj_keys, sample_information_ag_obj_vals))
                sample_information_ag_obj['file_id'] = file_id
                get_or_create(session, sampleinformation_model, **sample_information_ag_obj)
                #AG
                # AG
                # # Marker table
                insert_marker(session, marker_model, line)
                # # Primer table
                # insert_primer(session, primerpair_model,  line)
                # # TagPair table
                # insert_tagpair(session, tagpair_model, line)
                # # File table
                insert_file(session, file_model, line)
                # # Biosample table
                # insert_sample(session, biosample_model, line)
                # # File_information
                # insert_fileinformation(session, sampleinformation_model, line)
                # # Replicate table
                # insert_replicate(session, replicate_model, line)
                # # Replicatemarker table
                # insert_replicate_marker(session, replicatemarker_model, line)
                # session.commit()

