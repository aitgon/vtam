from wopmars.models.ToolWrapper import ToolWrapper
from vtam.utils.FileSampleInformation import FileSampleInformation


class SampleInformation(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.SampleInformation"
    }
    __input_file_csv = "sortedinfo"
    #
    __output_table_sample = "Sample"
    __output_table_sorted_read_file = "SortedReadFile"
    __output_table_marker = "Marker"
    __output_table_primerpair = "PrimerPair"
    __output_table_run = "Run"
    __output_table_sample_information = "SampleInformation"
    __output_table_tagpair = "TagPair"

    def specify_input_file(self):
        return [SampleInformation.__input_file_csv]

    def specify_output_table(self):
        return [
            SampleInformation.__output_table_sample,
            SampleInformation.__output_table_sorted_read_file,
            SampleInformation.__output_table_marker,
            SampleInformation.__output_table_run,
            SampleInformation.__output_table_sample_information,
        ]

    def specify_params(self):
        return{
            "fasta_dir": "str",
        }

    def run(self):
        session = self.session

        #######################################################################
        #
        # Wrapper inputs, outputs and parameters
        #
        #######################################################################

        # input file paths
        csv_path = self.input_file(SampleInformation.__input_file_csv)

        FileSampleInformation(csv_path).to_sqlite(session=session)

        #######################################################################
        #
        # Touch output tables, to update modification date
        #
        #######################################################################

        for output_table_i in self.specify_output_table():
            declarative_meta_i = self.output_table(output_table_i)
            obj = session.query(declarative_meta_i).order_by(
                declarative_meta_i.id.desc()).first()
            session.query(declarative_meta_i).filter_by(
                id=obj.id).update({'id': obj.id})
            session.commit()

    @staticmethod
    def get_or_create(session, model, **kwargs):
        instance = session.query(model).filter_by(**kwargs).first()
        if instance:  # get
            return instance
        else:  # create
            instance = model(**kwargs)
            session.add(instance)
            session.commit()
            return instance
