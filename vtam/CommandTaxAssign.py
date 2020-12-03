import inspect
import multiprocessing
import os
import pandas
import pathlib
import sqlalchemy

from vtam.models.TaxAssign import TaxAssign
from vtam.models.Variant import Variant
from vtam.utils.Logger import Logger
from vtam.utils.PathManager import PathManager
from vtam.utils.RunnerTaxAssign import RunnerTaxAssign
from vtam.utils.TaxLineage import TaxLineage
from vtam.utils.Taxonomy import Taxonomy


class CommandTaxAssign(object):
    """Class for the Pool Marker wrapper"""

    @classmethod
    def main(cls, db, mode, asvtable_tsv, output, taxonomy_tsv, blastdb_dir_path, blastdbname_str,
             num_threads=multiprocessing.cpu_count(), params=None):
        """

        Parameters
        ----------
        db: str
            Path to SQLITE database with Variant and Taxassign tables
        mode
        asvtable_tsv
        output
        taxonomy_tsv
        blastdb_dir_path
        blastdbname_str
        num_threads
        params

        Returns
        -------

        """

        this_temp_dir = os.path.join(
            PathManager.instance().get_tempdir(),
            os.path.basename(__file__))
        pathlib.Path(this_temp_dir).mkdir(exist_ok=True)

        #######################################################################
        #
        # Parameters
        #
        #######################################################################

        # params_dic = constants.get_params_default_dic()
        # params_dic = FileParams(params).get_params_dic()

        # ltg_rule_threshold = params_dic['ltg_rule_threshold']
        # include_prop = params_dic['include_prop']
        # min_number_of_taxa = params_dic['min_number_of_taxa']
        # qcov_hsp_perc = params_dic['qcov_hsp_perc']

        #######################################################################
        #
        # Load db and tables as classes and delete taxassign in reset mode
        #
        #######################################################################

        engine = sqlalchemy.create_engine('sqlite:///{}'.format(db), echo=False)

        variant_declarative_table = Variant.__table__
        variant_declarative_table.create(bind=engine, checkfirst=True)
        tax_assign_declarative_table = TaxAssign.__table__
        tax_assign_declarative_table.create(bind=engine, checkfirst=True)

        if mode == 'reset':
            with engine.connect() as conn:
                conn.execute(tax_assign_declarative_table.delete())

        #######################################################################
        #
        # Use variants that are not already already assigned in TaxAssign
        #
        #######################################################################

        variant_input_df = pandas.read_csv(asvtable_tsv, sep="\t", header=0)
        # get list of variant sequences
        variant_sequence_list = variant_input_df.sequence.tolist()

        # Add variant to DB if not already there
        for variant_sequence in variant_sequence_list:
            with engine.connect() as conn:
                row_variant = conn.execute(sqlalchemy.select([variant_declarative_table.c.id]) .where(
                    variant_declarative_table.c.sequence == variant_sequence)).first()
                if row_variant is None:  # variant_sequence IS NOT in the database, so INSERT it
                    conn.execute(
                        variant_declarative_table.insert().values(
                            sequence=variant_sequence))

        #######################################################################
        #
        # Get already tax-assigned variants with all informations including sequence
        #
        #######################################################################

        stmt_variant_tax_assign = sqlalchemy.select([
            tax_assign_declarative_table.c.variant_id,
            tax_assign_declarative_table.c.identity,
            tax_assign_declarative_table.c.ltg_rank,
            tax_assign_declarative_table.c.ltg_tax_id,
            tax_assign_declarative_table.c.ltg_tax_name,
            tax_assign_declarative_table.c.blast_db,
            variant_declarative_table.c.sequence,
        ])\
            .where(tax_assign_declarative_table.c.ltg_tax_id.isnot(None))\
            .where(tax_assign_declarative_table.c.variant_id == variant_declarative_table.c.id)\
            .where(variant_declarative_table.c.sequence.in_(variant_sequence_list))\
            .distinct()

        # These are the variants that are already in taxassign and do not need
        # recalculate
        ltg_from_db_list = []
        with engine.connect() as conn:
            for row in conn.execute(stmt_variant_tax_assign).fetchall():
                ltg_from_db_list.append(dict(zip(row.keys(), row.values())))
        """(Pdb) pandas.DataFrame.from_records(ltg_from_db_list)
   identity ltg_rank  ltg_tax_id              ltg_tax_name                                           sequence  variant_id
0       100  species     2028017  Orthocladiinae sp. BAP34  AGCATGATCTGGAATAGTAGGTACTTCCCTTAGTATCTTAATTCGA...         325
1        99  species     2028029   Rheocricotopus sp. DH90  GGCTTGATCCGGAATAGTAGGAACTTCTTTAAGAATTCTAATTCGA...        1203
2       100  species     1592914            Caenis pusilla  GGCTTGATCCGGAATGCTGGGCACCTCTCTAAGCCTTCTAATTCGT...        1443
3       100  species     2028029   Rheocricotopus sp. DH90  TGCTTGATCAGGAATAGTAGGAACTTCTTTAAGAATTCTAATTCGA...        2298
4        90   family        7149              Chironomidae  TGCTTGATCAGGGATAGTGGGAACTTCTTTAAGAATTCTTATTCGA...        2498
5       100  species      189839            Baetis rhodani  TGCTTGGGCAGGTATGGTAGGTACCTCATTAAGACTTTTAATTCGA...        2610"""
        ltg_db_df = pandas.DataFrame.from_records(ltg_from_db_list)
        ltg_db_df = ltg_db_df.reindex(
            sorted(
                ltg_db_df.columns),
            axis=1)  # sort columns

        #######################################################################
        #
        # Get list of variants (id and sequence) that need blast for tax assignation
        #
        #######################################################################

        stmt_variant = sqlalchemy.select([variant_declarative_table.c.id, variant_declarative_table.c.sequence]) \
            .where(variant_declarative_table.c.sequence.in_(variant_sequence_list)) \

        if ltg_db_df.shape[0] > 0:
            stmt_variant = stmt_variant.where(
                variant_declarative_table.c.id.notin_(
                    ltg_db_df.variant_id.tolist()))
        stmt_variant = stmt_variant.distinct().order_by("id")

        variant_not_tax_assigned = []
        with engine.connect() as conn:
            for row in conn.execute(stmt_variant).fetchall():
                variant_not_tax_assigned.append(
                    dict(zip(row.keys(), row.values())))

        #######################################################################
        #
        # Run RunnerTaxAssign for variant_not_tax_assigned
        #
        #######################################################################

        blast_variant_df = pandas.DataFrame()
        ltg_blast_df = pandas.DataFrame()

        if len(variant_not_tax_assigned) > 0:  # Run blast for variants that need tax assignation

            blast_variant_df = pandas.DataFrame.from_records(variant_not_tax_assigned, index='id')
            taxonomy = Taxonomy(tsv=taxonomy_tsv)
            sequence_list = blast_variant_df.sequence.tolist()
            tax_assign_runner = RunnerTaxAssign(
                sequence_list=sequence_list,
                taxonomy=taxonomy,
                blast_db_dir=blastdb_dir_path,
                blast_db_name=blastdbname_str,
                num_threads=num_threads,
                params = None)
            ltg_blast_df = tax_assign_runner.ltg_df

            ######################################################
            # Uncomment to debug because blast is slow
            # pandas.to_pickle(ltg_df, "ltg_df.pkl")
            # ltg_df = pandas.read_pickle("ltg_df.pkl")
            # import pdb; pdb.set_trace()
            ######################################################

            ltg_blast_df.rename({'variant_id': 'sequence'},
                                inplace=True, axis=1)

            ltg_blast_df = blast_variant_df.merge(
                ltg_blast_df, on='sequence', how='outer')

            ltg_blast_df['blast_db'] = blastdbname_str

            ltg_blast_df = ltg_blast_df.reindex(
                sorted(ltg_blast_df.columns), axis=1)  # sort columns
        del(blast_variant_df)  

        #######################################################################
        #
        # Concatenate tax-assigned variants from DB and from Blast
        # Merge variant_df and ltg_df and write to DB
        #
        #######################################################################

        if ltg_db_df.shape[0] > 0 and ltg_blast_df.shape[0] > 0:
            ltg_df = pandas.concat([ltg_db_df[["blast_db",
                                               "identity",
                                               "ltg_rank",
                                               "ltg_tax_id",
                                               "ltg_tax_name",
                                               "sequence"]],
                                    ltg_blast_df],
                                   axis=0)
        elif ltg_db_df.shape[0] > 0:
            ltg_df = ltg_db_df.copy()
        elif ltg_blast_df.shape[0] > 0:
            ltg_df = ltg_blast_df.copy()
        del(ltg_blast_df)

        #######################################################################
        #
        # Insert or update variant and taxassign tables
        #
        #######################################################################

        Logger.instance().debug(
            "file: {}; line: {}; Insert variant_id, ltg_tax_id, ltg_rank to DB".format(
                __file__, inspect.currentframe().f_lineno))

        for ltg_row in ltg_df.itertuples():
            variant_sequence = ltg_row.sequence
            with engine.connect() as conn:
                variant_id = conn.execute(sqlalchemy.select([variant_declarative_table.c.id]) .where(
                    variant_declarative_table.c.sequence == variant_sequence)).first()[0]
                select_row = conn.execute(
                    sqlalchemy.select(
                        [TaxAssign]) .where(
                        tax_assign_declarative_table.c.variant_id == variant_id)).first()
                # import pdb; pdb.set_trace()
                if select_row is None:  # variant_id IS NOT in the database, so INSERT it
                    ltg_row_dic = ltg_row._asdict()
                    ltg_row_dic['variant_id'] = variant_id
                    conn.execute(
                        tax_assign_declarative_table.insert(),
                        dict(ltg_row_dic))
                else:  # variant_sequence IS in the database, so update row
                    tax_assign_declarative_table.update() .where(
                        tax_assign_declarative_table.c.variant_id == variant_id).values()

        #######################################################################
        #
        # Update LTGs for variant output file
        #
        #######################################################################

        Logger.instance().debug(
            "file: {}; line: {}; Update LTGs for variant output file".format(
                __file__, inspect.currentframe().f_lineno))

        variant_output_df = variant_input_df.copy()
        del(variant_input_df)
        # Add ltg columns to variant_df if it do not exist
        for ltg_df_col in ['ltg_tax_id', 'ltg_tax_name', 'ltg_rank', 'identity', 'blast_db']:
            if not (ltg_df_col in variant_output_df.columns):
                variant_output_df[ltg_df_col] = None
        # Move sequence column to end
        variant_df_columns = variant_output_df.columns.tolist()
        variant_df_columns.append(variant_df_columns.pop(variant_df_columns.index('sequence')))
        variant_output_df = variant_output_df[variant_df_columns]

        for variant_row in variant_output_df.itertuples():
            # variant_id = variant_row.variant_id
            variant_sequence = variant_row.sequence
            with engine.connect() as conn:
                variant_id = conn.execute(sqlalchemy.select([variant_declarative_table.c.id]) .where(
                    variant_declarative_table.c.sequence == variant_sequence)).first()[0]
                select_row = conn.execute(
                    sqlalchemy.select(
                        [
                            TaxAssign.ltg_tax_id,
                            TaxAssign.ltg_tax_name,
                            TaxAssign.ltg_rank,
                            TaxAssign.identity,
                            TaxAssign.blast_db,
                        ]) .where(
                        tax_assign_declarative_table.c.variant_id == variant_id)).first()
            tax_assign_dict = dict(zip(
                ['ltg_tax_id', 'ltg_tax_name', 'ltg_rank', 'identity', 'blast_db'], select_row))
            for k in tax_assign_dict:
                variant_output_df.loc[variant_output_df.sequence ==
                                      variant_sequence, k] = tax_assign_dict[k]
        # do not move. required because sometimes tax_id is none
        variant_output_df = variant_output_df.astype({'ltg_tax_id': 'object'})

        #######################################################################
        #
        # Update tax lineages for variant output file
        #
        #######################################################################

        Logger.instance().debug(
            "file: {}; line: {}; Update tax lineages for variant output file".format(
                __file__, inspect.currentframe().f_lineno))

        tax_id_list = variant_output_df.ltg_tax_id.unique().tolist()  # unique list of tax ids
        tax_lineage = TaxLineage(taxonomic_tsv_path=taxonomy_tsv)
        tax_lineage_df = tax_lineage.create_lineage_from_tax_id_list(
            tax_id_list=tax_id_list, tax_name=True)

        # Merge
        variant_output_df = variant_output_df.merge(
            tax_lineage_df, left_on='ltg_tax_id', right_on='tax_id', how='left')
        variant_output_df.drop('tax_id', axis=1, inplace=True)

        Logger.instance().debug(
            "file: {}; line: {}; Reorder columns".format(
                __file__, inspect.currentframe().f_lineno))
        # Move sequence column to end
        variant_df_columns = variant_output_df.columns.tolist()
        variant_df_columns.append(variant_df_columns.pop(variant_df_columns.index('sequence')))
        variant_output_df = variant_output_df[variant_df_columns]
        Logger.instance().debug(
            "file: {}; line: {}; Write to TSV".format(
                __file__, inspect.currentframe().f_lineno))
        variant_output_df.to_csv(output, sep='\t', index=False, header=True)

