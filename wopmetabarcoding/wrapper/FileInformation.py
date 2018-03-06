from wopmars.framework.database.tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
from wopmetabarcoding.wrapper.functions import insert_table

# import models
import wopmetabarcoding.model.Marker
import wopmetabarcoding.model.PrimerPair


class FileInformation(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.FileInformation"
    }
    __input_file_csv = "csv"
    __output_table_marker = "Marker"
    __output_table_primerpair = "PrimerPair"
    __output_table_tag = "Tag"
    __output_table_file = "File"
    __output_table_sample = "Sample"
    __output_table_fileinformation = "FileInformation"

    def specify_input_file(self):
        return [FileInformation.__input_file_csv]

    def specify_output_table(self):
        return [
            FileInformation.__output_table_marker, FileInformation.__output_table_primerpair,
            FileInformation.__output_table_tag, FileInformation.__output_table_file,
            FileInformation.__output_table_fileinformation, FileInformation.__output_table_sample,
        ]

    def insert_marker(self, session, model, line):
        marker_name = line.split(',')[4]
        obj_marker = {'marker_name': marker_name}
        insert_table(session, model, obj_marker)

    def insert_primer(self, session, model, line):
        primer_forward = line.split(',')[1]
        primer_reverse = line.split(',')[3]
        obj_primer = {'primer_forward': primer_forward, 'primer_reverse': primer_reverse}
        insert_table(session, model, obj_primer)

    def insert_tag(self, session, model, line):
        tag_forward = line.split(',')[0]
        tag_reverse = line.split(',')[2]
        obj_tag = {'tag_forward': tag_forward, 'tag_reverse': tag_reverse}
        insert_table(session, model, obj_tag)

    def insert_file(self, session, model, line):
        file_name = line.split(',')[7]
        run_name = line.split(',')[8].strip()
        if line.split(',')[0] == "" or line.split(',')[2] == "" or line.split(',')[3] == "" or line.split(',')[4] == "":
            dereplicate = True
        else:
            dereplicate = False
        obj_file = {'file_name': file_name, 'run_name': run_name, 'dereplicate_status': dereplicate}
        insert_table(session, model, obj_file)

    def insert_sample(self, session, model, line):
        sample_name = line.split(',')[5]
        obj_sample = {'sample_name': sample_name}
        insert_table(session, model, obj_sample)

    def insert_fileinformation(self, session, model, line):
        marker_name = line.split(',')[4]
        primer_forward = line.split(',')[1]
        primer_reverse = line.split(',')[3]
        tag_forward = line.split(',')[0]
        tag_reverse = line.split(',')[2]
        file_name = line.split(',')[7]
        run_name = line.split(',')[8].strip()
        sample_name = line.split(',')[5]
        obj_fileinformation = {
            'marker_name': marker_name, 'tag_forward': tag_forward, 'tag_reverse': tag_reverse,
            'primer_forward': primer_forward, 'primer_reverse': primer_reverse, 'file_name': file_name,
            'run_name': run_name, 'sample_name': sample_name,
        }
        insert_table(session, model, obj_fileinformation)

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
        sample_model = self.output_table(FileInformation.__output_table_sample)
        fileinformation_model = self.output_table(FileInformation.__output_table_fileinformation)
        try:
            with open(csv_path, 'r') as fin:
                next(fin)  # skip header
                for line in fin:
                    # Marker table
                    self.insert_marker(session, marker_model, line)
                    # Primer table
                    self.insert_primer(session, primerpair_model,  line)
                    # Tag table
                    self.insert_tag(session, tag_model, line)
                    # File table
                    self.insert_file(session, file_model, line)
                    # Sample table
                    self.insert_sample(session, sample_model, line)
                    # File_information
                    self.insert_fileinformation(session, fileinformation_model, line)
                    session.commit()

        except FileNotFoundError:
            context = \
                'While searching for the csv file or the fastq directory, any file or directory are found at the pointed directory'
            Logger.instance().exception(context + ' Please check if the csv file or the fastq directory are present')

