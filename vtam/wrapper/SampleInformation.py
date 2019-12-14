from wopmars.models.ToolWrapper import ToolWrapper
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
            SampleInformation.__output_table_run,
            SampleInformation.__output_table_sample_information,
            SampleInformation.__output_table_tagpair,
        ]

    def specify_params(self):
        return{
            "fasta_dir": "str",
        }

    def run(self):
        session = self.session
        engine = session._session().get_bind()

        fasta_dir = str(os.getenv('VTAM_FASTA_DIR'))

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
                # replicate_name = line.strip().split('\t')[6]
                replicate = line.strip().split('\t')[6]
                if len(line.strip().split('\t')) == 9: # No fastq file columns
                    fasta_file_name = line.strip().split('\t')[8]
                    # fasta_file_full_path = os.path.join(self.option("fasta_dir"), line.strip().split('\t')[8])
                    fasta_file_full_path = os.path.join(fasta_dir, line.strip().split('\t')[8])
                elif len(line.strip().split('\t')) == 11: # Yes fastq file columns
                    fasta_file_name = line.strip().split('\t')[10]

                    # fasta_file_full_path = os.path.join(self.option("fasta_dir"), line.strip().split('\t')[10])
                    fasta_file_full_path = os.path.join(fasta_dir, line.strip().split('\t')[10])
                else:
                    raise IndexError("{} error. Verify nb of columns in this input file: {}".format(self.__class__.__name__, os.path.join(os.getcwd(), csv_path)))
                if not os.path.isfile(os.path.join(os.getcwd(), fasta_file_full_path)):
                    raise FileNotFoundError("{} error. Verify this file output: {}".format(self.__class__.__name__, os.path.join(os.getcwd(), fasta_file_full_path)))
                run_name = line.strip().split('\t')[7]
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
                is_trimmed = False # Default
                fasta_obj = {'name': fasta_file_name, 'is_trimmed': is_trimmed, 'run_id': run_id}
                fasta_instance = SampleInformation.get_or_create(session, fasta_model, **fasta_obj)
                fasta_id = fasta_instance.id

                # Insert sample_information
                sample_information_obj = {'tag_forward': tag_forward, 'primer_forward': primer_forward, 'tag_reverse':
                    tag_reverse, 'primer_reverse': primer_reverse, 'biosample_id': biosample_id, 'replicate': replicate,
                                          'run_id': run_id}
                sample_information_obj['fasta_id'] = fasta_id
                sample_information_obj['marker_id'] = marker_id
                SampleInformation.get_or_create(session, sampleinformation_model, **sample_information_obj)

                # Primer pair
                primerpair_obj = {'primer_forward': primer_forward, 'primer_reverse': primer_reverse}
                SampleInformation.get_or_create(session, primerpair_model, **primerpair_obj)

                # Tag pair
                tagpair_obj = {'tag_forward': tag_forward, 'tag_reverse': tag_reverse}
                SampleInformation.get_or_create(session, tagpair_model, **tagpair_obj)

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
