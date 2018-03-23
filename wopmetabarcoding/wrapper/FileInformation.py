from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
from wopmetabarcoding.wrapper.FileInformation_functions import \
    insert_marker, insert_fileinformation, insert_file, insert_tag, insert_primer, insert_sample

# import models
import wopmetabarcoding.model.Marker
import wopmetabarcoding.model.PrimerPair


class FileInformation(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.SampleInformation"
    }
    __input_file_csv = "csv"
    #
    __output_table_file = "File"
    __output_table_sampleinformation = "SampleInformation"
    __output_table_marker = "Marker"
    __output_table_primerpair = "PrimerPair"
    __output_table_biosample = "Biosample"
    __output_table_tag = "Tag"

    def specify_input_file(self):
        return [FileInformation.__input_file_csv]

    def specify_output_table(self):
        return [
            FileInformation.__output_table_file,
            FileInformation.__output_table_sampleinformation,
            FileInformation.__output_table_marker,
            FileInformation.__output_table_primerpair,
            FileInformation.__output_table_biosample,
            FileInformation.__output_table_tag,
        ]

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        conn = engine.connect()
        # file paths
        csv_path = self.input_file(FileInformation.__input_file_csv)
        # Models
        marker_model = self.output_table(FileInformation.__output_table_marker)
        primerpair_model = self.output_table(FileInformation.__output_table_primerpair)
        tag_model = self.output_table(FileInformation.__output_table_tag)
        file_model = self.output_table(FileInformation.__output_table_file)
        biosample_model = self.output_table(FileInformation.__output_table_biosample)
        sampleinformation_model = self.output_table(FileInformation.__output_table_sampleinformation)
        try:
            with open(csv_path, 'r') as fin:
                next(fin)  # skip header
                for line in fin:
                    # Marker table
                    insert_marker(session, marker_model, line)
                    # Primer table
                    insert_primer(session, primerpair_model,  line)
                    # Tag table
                    insert_tag(session, tag_model, line)
                    # File table
                    insert_file(session, file_model, line)
                    # Biosample table
                    insert_sample(session, biosample_model, line)
                    # File_information
                    insert_fileinformation(session, sampleinformation_model, line)
                    session.commit()
        except FileNotFoundError:
            context = \
                'While searching for the csv file or the fastq directory, any file or directory are found at the pointed directory'
            Logger.instance().exception(context + ' Please check if the csv file or the fastq directory are present')

