import inspect
import os
import sqlalchemy

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
    def main(cls, db, mode, variants_tsv, variant_taxa_tsv, taxonomy_tsv, blasdb_dir_path, blastdbname_str, ltg_rule_threshold, include_prop,
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
        variants_declarative_table.create(bind=engine, checkfirst=True)
        tax_assign_declarative_table = tax_assign_declarative.__table__
        tax_assign_declarative_table.create(bind=engine, checkfirst=True)

        if mode == 'reset':
            with engine.connect() as conn:
                conn.execute(tax_assign_declarative_table.delete())

        ################################################################################################################
        #
        # Use variants that are not already already assigned in TaxAssign
        #
        ################################################################################################################

        variant_input_df = pandas.read_csv(variants_tsv, sep="\t", header=0)
        variant_sequence_list = variant_input_df.ix[:, -1].tolist()  # get variant sequences as list

        # Add variant to DB if not already there
        for variant_sequence in variant_sequence_list:
            with engine.connect() as conn:
                row_variant = conn.execute(sqlalchemy.select([variants_declarative_table.c.id])
                                          .where(variants_declarative_table.c.sequence == variant_sequence)).first()
                if row_variant is None:  # variant_sequence IS NOT in the database, so INSERT it
                    conn.execute(variants_declarative_table.insert().values(sequence=variant_sequence))

        ################################################################################################################
        #
        # Select if not already taxassigned
        #
        ################################################################################################################

        stmt_variant_tax_assign = select([tax_assign_declarative_table.c.variant_id])\
            .where(tax_assign_declarative_table.c.ltg_tax_id.isnot(None))

        # These are the variants that are already in taxassign and do not need recalculate
        variant_tax_assign_list = []
        with engine.connect() as conn:
            variant_tax_assign_list = [variant_id[0] for variant_id in conn.execute(stmt_variant_tax_assign).fetchall()]

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

        ################################################################################################################
        #
        # Merge variant_df and ltg_df and write to DB
        #
        ################################################################################################################

        # if not (ltg_df is None) and ltg_df.shape[0] > 0:

        ltg_df = variant_df.merge(ltg_df, left_index=True, right_on='variant_id', how='outer')

        ltg_df['blast_db'] = blastdbname_str

        ############################################################################################################
        #
        # Insert or update data into DB
        #
        ############################################################################################################

        Logger.instance().debug("file: {}; line: {}; Insert variant_id, ltg_tax_id, ltg_rank to DB".format(__file__,
                                                                                   inspect.currentframe().f_lineno))

        for ltg_row in ltg_df.itertuples():
            variant_id = ltg_row.variant_id
            with engine.connect() as conn:
                select_row = conn.execute(sqlalchemy.select([tax_assign_declarative])
                                          .where(tax_assign_declarative_table.c.variant_id == variant_id)).first()
                if select_row is None:  # variant_id IS NOT in the database, so INSERT it
                    conn.execute(tax_assign_declarative_table.insert(), dict(ltg_row._asdict()))
                else:  # variant_sequence IS in the database, so update row
                    tax_assign_declarative_table.update()\
                        .where(tax_assign_declarative_table.c.variant_id == variant_id).values()

        ############################################################################################################
        #
        # Update data of variant input file
        #
        ############################################################################################################

        variant_output_df = variant_input_df.copy()
        # Add ltg columns to variant_df if it do not exist
        for ltg_df_col in ['ltg_tax_id', 'ltg_tax_name', 'ltg_rank', 'identity', 'blast_db']:
            if not (ltg_df_col in variant_output_df.columns):
                variant_output_df[ltg_df_col] = None
        # Move sequence column to end
        variant_df_columns = variant_output_df.columns.tolist()
        variant_df_columns.append(variant_df_columns.pop(variant_df_columns.index('sequence')))
        variant_output_df = variant_output_df[variant_df_columns]

        for variant_row in variant_output_df.itertuples():
            variant_id = variant_row.variant_id
            with engine.connect() as conn:
                select_row = conn.execute(sqlalchemy.select([tax_assign_declarative.ltg_tax_id,
                                                             tax_assign_declarative.ltg_tax_name,
                                                             tax_assign_declarative.ltg_rank,
                                                             tax_assign_declarative.identity,
                                                             tax_assign_declarative.blast_db,
                                                             ])
                                          .where(tax_assign_declarative_table.c.variant_id == variant_id)).first()
            tax_assign_dict = dict(zip(['ltg_tax_id', 'ltg_tax_name', 'ltg_rank', 'identity', 'blast_db'],
                                        select_row))
            for k in tax_assign_dict:
                variant_output_df.loc[variant_output_df.variant_id == variant_id, k] = tax_assign_dict[k]

        variant_output_df.to_csv(variant_taxa_tsv, sep='\t', index=False, header=True)

