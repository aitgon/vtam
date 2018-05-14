from wopmars.framework.database. tables.ToolWrapper import ToolWrapper
from sqlalchemy import select
from wopmetabarcoding.utils.constants import tempdir
from wopmetabarcoding.wrapper.FilterUtilities import Variant2Sample2Replicate2Count
from wopmetabarcoding.utils.VSearch import VSearch1
import os, pickle


class Taxassign(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.Taxassign"
    }
    __input_table_marker = "Marker"
    __input_table_variant = "Variant"
    __input_file_taxassign_db = "taxassign_db_fasta"
    __filtered_dataframe_path = "filtered_dataframe_path"

    def specify_input_table(self):
        return [
            Taxassign.__input_table_marker,
            Taxassign.__input_table_variant
        ]

    def specify_input_file(self):
        return [
            Taxassign.__input_file_taxassign_db,
            Taxassign.__filtered_dataframe_path

        ]

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        # Input table models
        marker_model = self.input_table(Taxassign.__input_table_marker)
        variant_model = self.input_table(Taxassign.__input_table_variant)
        # Input files
        taxassign_db_fasta = self.input_file(Taxassign.__input_file_taxassign_db)
        filter_output = self.input_file(Taxassign.__filtered_dataframe_path)
        #
        with open(filter_output, 'r') as fin:
            for line in fin:
                line = line.strip().split('\t')
                marker_name = line[0]
                dataframe_path = line[1]
                filtered_variants_fasta = line[2]
                variant2sample2replicate2count_df = pickle.load(open(dataframe_path, 'rb'))
                output_tsv = filtered_variants_fasta.replace('.fasta', '.tsv')
                # vsearch_usearch_global_args = {'db': taxassign_db_fasta,
                #                                'usearch_global': filtered_variants,
                #                                'id': 0,
                #                                'maxrejects': 0,
                #                                'maxaccepts': 0,
                #                                'userout': output_tsv,
                #                                'userfields': "--userfields query+target+id",
                #                                }
                # vsearch_1 = VSearch1(**vsearch_usearch_global_args)
                # vsearch_1.run()
                os.system(
                    "vsearch --usearch_global " + filtered_variants_fasta +" --db "+ taxassign_db_fasta +
                    " --maxaccept 0 --maxreject 0  --userout " + output_tsv + " --userfields query+target+id --id "
                    + str(0.8)
                )



