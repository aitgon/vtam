import inspect
import os
import pandas
import pathlib

from sqlalchemy import create_engine, select

from vtam.models.TaxAssign import TaxAssign as tax_assign_declarative
from vtam.models.Variant import Variant as variant_declarative
from vtam.utils.Logger import Logger
from vtam.utils.PathManager import PathManager
from vtam.utils.TaxAssignRunner import TaxAssignRunner
from vtam.utils.VariantDFutils import VariantDFutils


class CommandTaxAssign(object):
    """Class for the Pool Marker wrapper"""

    @classmethod
    def main(cls, db, mode, variants_tsv, taxonomy_tsv, blasdb_dir_path, blastdbname_str, ltg_rule_threshold, include_prop,
             min_number_of_taxa, num_threads):

        this_temp_dir = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__))
        pathlib.Path(this_temp_dir).mkdir(exist_ok=True)

        ################################################################################################################
        #
        # Load db and tables as classes and delete taxassign in reset mode
        #
        ################################################################################################################

        engine = create_engine('sqlite:///{}'.format(db), echo=False)

        variants_declarative_table = variant_declarative.__table__
        tax_assign_declarative_table = tax_assign_declarative.__table__

        if mode == 'reset':
            with engine.connect() as conn:
                conn.execute(tax_assign_declarative_table.delete())

        ################################################################################################################
        #
        # Use variants that are not already already assigned in TaxAssign
        #
        ################################################################################################################

        variants_df = pandas.read_csv(variants_tsv, sep="\t", header=0)
        variant_sequence_list = variants_df.ix[:,-1].tolist()  # last column

        Logger.instance().debug(
            "file: {}; line: {}; Get variants and sequences that passed the filters".format(__file__,
                                                                                            inspect.currentframe().f_lineno,
                                                                                            'TaxAssign'))

        ################################################################################################################
        #
        # Select if not already taxassigned
        #
        ################################################################################################################

        tax_assign_declarative_table.create(bind=engine, checkfirst=True)
        stmt_variant_tax_assign = select([tax_assign_declarative_table.c.variant_id])\
            .where(tax_assign_declarative_table.c.ltg_tax_id.isnot(None))

        # These are the variants that are already in taxassign and do not need recalculate
        variant_tax_assign_list = []
        with engine.connect() as conn:
            variant_tax_assign_list = [variant_id[0] for variant_id in conn.execute(stmt_variant_tax_assign).fetchall()]

        variants_declarative_table.create(bind=engine, checkfirst=True)

        stmt_variant = select([variants_declarative_table.c.id, variants_declarative_table.c.sequence]) \
            .where(variants_declarative_table.c.sequence.in_(variant_sequence_list)) \
            .where(variants_declarative_table.c.id.notin_(variant_tax_assign_list)) \
            .distinct()\
            .order_by("id")

        # Select to DataFrame the variants that will be tried to be assigned to taxa
        variant_list = []
        with engine.connect() as conn:
            for row in conn.execute(stmt_variant).fetchall():
                variant_list.append({'id': row.id, 'sequence': row.sequence})
        variant_df = pandas.DataFrame.from_records(variant_list, index='id')

        ################################################################################################################
        #
        # 2 Create FASTA file with Variants
        #
        ################################################################################################################

        Logger.instance().debug(
            "file: {}; line: {}; Create Fasta from Variants".format(__file__, inspect.currentframe().f_lineno))
        variant_fasta = os.path.join(this_temp_dir, 'variant.fasta')
        variant_df_utils = VariantDFutils(variant_df)
        variant_df_utils.to_fasta(variant_fasta)

        ################################################################################################################
        #
        # Run TaxAssignRunner
        #
        ################################################################################################################

        taxonomy_df = pandas.read_csv(taxonomy_tsv, sep="\t", header=0,
                                      dtype={'tax_id': 'int', 'parent_tax_id': 'int', 'old_tax_id': 'float'})

        tax_assign_runner = TaxAssignRunner(variant_df=variant_df, taxonomy_df=taxonomy_df,
                                            blast_db_dir=blasdb_dir_path, blast_db_name=blastdbname_str,
                                            ltg_rule_threshold=ltg_rule_threshold,
                                            include_prop=include_prop, min_number_of_taxa=min_number_of_taxa,
                                            num_threads=num_threads)
        ltg_df = tax_assign_runner.ltg_df
        ltg_df['blastdb'] = blastdbname_str

        ################################################################################################################
        #
        # 6. Insert data
        #
        ################################################################################################################

        Logger.instance().debug("file: {}; line: {}; Insert variant_id, ltg_tax_id, ltg_rank to DB".format(__file__,
                                                                                                           inspect.currentframe().f_lineno))
        with engine.connect() as conn:
            conn.execute(tax_assign_declarative_table.insert(), ltg_df.to_dict('records'))
