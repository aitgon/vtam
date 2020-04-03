from wopmars.models.ToolWrapper import ToolWrapper
import pandas


class SampleInformation(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.SampleInformation"
    }
    __input_file_csv = "readinfo"
    #
    __output_table_biosample = "Biosample"
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
            SampleInformation.__output_table_biosample,
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
        engine = session._session().get_bind()

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
        fasta_model = self.output_table(SampleInformation.__output_table_sorted_read_file)
        marker_model = self.output_table(SampleInformation.__output_table_marker)
        run_model = self.output_table(SampleInformation.__output_table_run)
        sampleinformation_model = self.output_table(SampleInformation.__output_table_sample_information)

        sorted_read_info_df = pandas.read_csv(csv_path, sep="\t", header=0)
        sorted_read_info_df.columns = sorted_read_info_df.columns.str.lower()
        for row in sorted_read_info_df.itertuples():

            run_name = row.run
            marker_name = row.marker
            biosample_name = row.biosample
            replicate = row.replicate
            sorted_read_file = row.fastasorted
            #
            # Insert run
            run_obj = {'name': run_name}
            run_instance = SampleInformation.get_or_create(session, run_model, **run_obj)
            run_id = run_instance.id
            #
            # Insert marker_id
            marker_obj = {'name': marker_name}
            marker_instance = SampleInformation.get_or_create(session, marker_model, **marker_obj)
            marker_id = marker_instance.id
            #
            # Insert Biosamples
            biosample_obj = {'name': biosample_name}
            biosample_instance = SampleInformation.get_or_create(session, biosample_model, **biosample_obj)
            biosample_id = biosample_instance.id
            #
            # Insert file output
            fasta_obj = {'name': sorted_read_file, 'run_id': run_id}
            fasta_instance = SampleInformation.get_or_create(session, fasta_model, **fasta_obj)
            fasta_id = fasta_instance.id

            # Insert sample_information
            sample_information_obj = {'biosample_id': biosample_id, 'replicate': replicate, 'run_id': run_id}
            sample_information_obj['sortedreadfile_id'] = fasta_id
            sample_information_obj['marker_id'] = marker_id
            SampleInformation.get_or_create(session, sampleinformation_model, **sample_information_obj)

    @staticmethod
    def get_or_create(session, model, **kwargs):
        instance = session.query(model).filter_by(**kwargs).first()
        if instance:
            return instance
        else:
            instance = model(**kwargs)
            session.add(instance)
            session.commit()
            return instance
