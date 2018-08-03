import multiprocessing
from itertools import repeat

import resource

import numpy
from wopmars.framework.database. tables.ToolWrapper import ToolWrapper
from wopmetabarcoding.wrapper.TaxAssignUtilities import vsearch_command, sub_fasta_creator, vsearch_output_to_sqlite, \
    get_vsearch_output_for_variant_as_df, f_variant_vsearch_output_to_ltg, convert_fileinfo_to_otu_df, f_taxid2taxname
import pandas,os
from wopmetabarcoding.utils.constants import tempdir
from wopmetabarcoding.utils.logger import logger, LOGGER_LEVEL

import inspect
from wopmetabarcoding.utils.PathFinder import PathFinder

from sqlalchemy import select

def f_variant2taxid(variant_id, tax_assign_sqlite, tax_assign_pars_tsv):
    # marker_name = variant_class[variant_seq][0]
    logger.debug(
        "file: {}; line: {}; variant_id {}...".format(__file__, inspect.currentframe().f_lineno,
                                                                   variant_id))
    # vsearch_per_variant_df_pkl = variant_class[variant_seq][1]

    vsearch_per_variant_df_pkl = os.path.join(tempdir, "TaxAssign",
                                              "vsearch_df_for_variant_id_{}.pkl".format(variant_id))
    tax_id = f_variant_vsearch_output_to_ltg(variant_id, vsearch_per_variant_df_pkl, tax_assign_sqlite, tax_assign_pars_tsv)
    return (variant_id, tax_id)


class TaxAssign(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.TaxAssign"
    }
    # Input file
    __input_file_db_udb = "db_udb"
    __input_file_tax_assign_db_sqlite = "tax_assign_db_sqlite"
    __input_tax_assign_pars_tsv = "tax_assign_pars_tsv"
    # Input table
    __input_table_variant = "Variant"
    __input_table_variant_read_count = "VariantReadCount"
    __input_table_variant_selected = "VariantSelected"
    # Output table
    __output_table_variant_taxa = 'VariantTaxa'

    def specify_input_file(self):
        return [
            TaxAssign.__input_file_db_udb,
            TaxAssign.__input_file_tax_assign_db_sqlite,
            TaxAssign.__input_tax_assign_pars_tsv,
        ]

    def specify_input_table(self):
        return [
            TaxAssign.__input_table_variant,
            TaxAssign.__input_table_variant_read_count,
            TaxAssign.__input_table_variant_selected,
        ]

    def specify_output_table(self):
        return [
            TaxAssign.__output_table_variant_taxa,
        ]

    def specify_params(self):
        return {
            "output_dir_taxassign": "str",
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        #
        # Input file
        db_udb = self.input_file(TaxAssign.__input_file_db_udb)
        tax_assign_pars_tsv = self.input_file(TaxAssign.__input_tax_assign_pars_tsv)
        tax_assign_sqlite = self.input_file(TaxAssign.__input_file_tax_assign_db_sqlite)
        #
        # Input table models
        variant_model = self.input_table(TaxAssign.__input_table_variant)
        variant_read_count_model = self.input_table(TaxAssign.__input_table_variant_read_count)
        variant_selected_model = self.input_table(TaxAssign.__input_table_variant_selected)
        #
        # Output file
        variant_taxa_model = self.output_table(TaxAssign.__output_table_variant_taxa)
        with engine.connect() as conn:
            conn.execute(variant_taxa_model.__table__.delete())
        #
        ########################################
        # Create variant selected df
        ########################################
        variant_model_table = variant_model.__table__
        variant_read_count_model_table = variant_read_count_model.__table__
        variant_selected_model_table = variant_selected_model.__table__
        stmt = select([variant_read_count_model_table.c.marker_id, variant_selected_model_table.c.variant_id, variant_model_table.c.sequence])\
            .distinct()\
            .where(variant_selected_model_table.c.passed==1)\
            .where(variant_model_table.c.id==variant_selected_model_table.c.variant_id)\
            .where(variant_read_count_model_table.c.variant_id==variant_selected_model_table.c.variant_id)\
            .order_by(variant_read_count_model_table.c.marker_id) \
            .order_by(variant_selected_model_table.c.variant_id)
        variant_selected_list = []
        with engine.connect() as conn:
            for row2 in conn.execute(stmt).fetchall():
                variant_selected_list.append(row2)
        variant_selected_df = pandas.DataFrame.from_records(variant_selected_list, columns=['marker_id', 'variant_id', 'variant_sequence'])
        #
        ################################################
        # for subsets of variants, create fasta files for vsearch
        ################################################
        PathFinder.mkdir_p(os.path.join(tempdir, "TaxAssign"))
        fasta_path = os.path.join(tempdir, "TaxAssign", "variants.fasta")
        logger.debug("file: {}; line: {}; Write fasta to {}".format(__file__, inspect.currentframe().f_lineno, fasta_path))
        variant_id_list = []
        for row in variant_selected_df.itertuples(index=True):
            variant_id_list.append(row.variant_id)
            with open(fasta_path, 'a') as fout:
                fout.write(">{}\n{}\n".format(row.variant_id, row.variant_sequence))
            #  every 10 sequences run and analyse variants or at the end
            if ((row.Index+1)%10 == 0) or ((row.Index+1) == variant_selected_df.shape[0]):
                #
                #################################################
                # Run vsearch to find COI targets in the Metazoa database
                #################################################
                vsearch_output_tsv = os.path.join(tempdir, "TaxAssign", "vsearch_variant_number_{}.tsv".format(row.Index))
                logger.debug("file: {}; line: {}; Run vsearch and store here: {}".format(__file__, inspect.currentframe().f_lineno, fasta_path))
                vsearch_command(fasta_path, db_udb, vsearch_output_tsv)
                #
                #################################################
                # Save vsearch output in sqlite DB
                #################################################
                vsearch_output_sqlite = os.path.join(tempdir, "TaxAssign", "vsearch_variant_number_{}.sqlite".format(row.Index))
                logger.debug("file: {}; line: {}; Convert vsearch TSV to SQLITE: {}".format(__file__, inspect.currentframe().f_lineno,
                                                                                               vsearch_output_sqlite))
                vsearch_output_to_sqlite(vsearch_output_tsv, vsearch_output_sqlite)
                #
                #################################################
                # Store pickle DF per variant
                #################################################
                logger.debug(
                    "file: {}; line: {}; Store pickle DF per variant".format(__file__, inspect.currentframe().f_lineno))
                for variant_id in variant_id_list:
                    vsearch_per_variant_df = get_vsearch_output_for_variant_as_df(vsearch_output_sqlite, variant_id)
                    vsearch_per_variant_df_pkl = os.path.join(tempdir, "TaxAssign",
                                                         "vsearch_df_for_variant_id_{}.pkl".format(variant_id))
                    vsearch_per_variant_df.to_pickle(vsearch_per_variant_df_pkl, compression='gzip')
                #
                os.remove(fasta_path)
                variant_id_list = [] # reset variant id list
        #
        # ################################################
        # # Start of Parallel Version: Comment out after debugging
        # ################################################
        # variant_id_list = sorted(variant_selected_df.variant_id.unique().tolist())
        # logger.debug(
        #     "file: {}; line: {}; Parallel. Mem {}".format(__file__, inspect.currentframe().f_lineno, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        # with multiprocessing.Pool(multiprocessing.cpu_count(), maxtasksperchild=1) as p:
        #     variant2marker2taxid_list = p.starmap(f_variant2taxid, zip(variant_id_list, repeat(tax_assign_sqlite),
        #                                                                repeat(tax_assign_pars_tsv)), chunksize=1)
        # # End of Parallel Version: Comment out after debugging
        #
        ################################################
        # Start of Non Parallel Version: Comment out after debugging
        ################################################
        variant2taxid_list = []
        for variant_id in variant_selected_df.variant_id.unique():
            variant_id, tax_id = f_variant2taxid(variant_id, tax_assign_sqlite, tax_assign_pars_tsv)
            variant_id = int(variant_id)
            if numpy.isnan(tax_id):
                tax_id = None
            else:
                tax_id = int(tax_id)
            variant2taxid_list.append({'variant_id': variant_id, 'tax_id': tax_id})
        #
        ################################
        # Insert into db
        ################################
        with engine.connect() as conn:
            conn.execute(variant_taxa_model.__table__.insert(), variant2taxid_list)
        # End of Non Parallel Version: Comment out after debugging
