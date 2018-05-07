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

    def specify_input_table(self):
        return [
            Taxassign.__input_table_marker,
            Taxassign.__input_table_variant
        ]

    def specify_input_file(self):
        return [
            Taxassign.__input_file_taxassign_db
        ]

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        # Input table models
        marker_model = self.input_table(Taxassign.__input_table_marker)
        variant_model = self.input_table(Taxassign.__input_table_variant)
        # Input files
        taxassign_db_fasta = self.input_file(Taxassign.__input_file_taxassign_db)
        marker_select = select([marker_model.id, marker_model.name])
        marker_obj = engine.execute(marker_select)
        for marker in marker_obj:
            filtered_variants = os.path.join(tempdir, (marker.name + "_filtered_variants.fasta"))
            filtered_variant_dataframe = os.path.join(tempdir, (marker.name + "_filtered_dataframe.pkl"))
            variant2sample2replicate2count_df = pickle.load(open(filtered_variant_dataframe, 'rb'))
            variants = list(set(variant2sample2replicate2count_df['sequence'].tolist()))
            Variant2Sample2Replicate2Count.filter_fasta(variants, filtered_variants, False)
            output_tsv = filtered_variants.replace('.fasta', '.tsv')
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
                "vsearch --usearch_global " + filtered_variants +" --db "+ taxassign_db_fasta +
                " --maxaccept 0 --maxreject 0  --userout " + output_tsv + " --userfields query+target+id --id "
                + str(0.1)
            )



