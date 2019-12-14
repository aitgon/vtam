import errno
import inspect
import os
import pathlib

import pandas
from Bio.Blast.Applications import NcbiblastnCommandline
from sqlalchemy import select
from wopmars.models.ToolWrapper import ToolWrapper

from vtam.utils.VariantDFutils import VariantDFutils
from vtam.utils.Logger import Logger
from vtam.utils.PathManager import PathManager
from vtam.utils.TaxAssignUtilities import f01_taxonomy_sqlite_to_df
from vtam.utils.TaxAssignUtilities import f04_1_tax_id_to_taxonomy_lineage
from vtam.utils.TaxAssignUtilities import f07_blast_result_to_ltg_tax_id


class TaxAssign(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.TaxAssign"}

    # Input file
    __input_file_fastainfo = "fastainfo"
    __input_file_taxonomy = "taxonomy"
    # __input_file_map_taxids = "map_taxids"
    # Input table
    __input_table_marker = "Marker"
    __input_table_run = "Run"
    __input_table_biosample = "Biosample"
    __input_table_filter_codon_stop = "FilterCodonStop"
    __input_table_variant = "Variant"

    # Output table
    __output_table_tax_assign = "TaxAssign"

    def specify_input_file(self):
        return [
            TaxAssign.__input_file_fastainfo,
            TaxAssign.__input_file_taxonomy,
            # TaxAssign.__input_file_map_taxids,
        ]

    def specify_input_table(self):
        return [
            TaxAssign.__input_table_marker,
            TaxAssign.__input_table_run,
            TaxAssign.__input_table_biosample,
            TaxAssign.__input_table_variant,
            TaxAssign.__input_table_filter_codon_stop,
        ]

    def specify_output_table(self):
        return [
            TaxAssign.__output_table_tax_assign,

        ]

    def specify_params(self):
        return {
            "ltg_rule_threshold": "float",  # percentage
            "include_prop": "float",  # percentage
            "min_number_of_taxa": "int",  # count
            "log_verbosity": "int",
            "log_file": "str",
            "blast_db": "str",
            "num_threads": "str",
        }

    def run(self):
        session = self.session
        engine = session._session().get_bind()
        threads = int(os.getenv('VTAM_THREADS'))

        this_temp_dir = os.path.join(PathManager.instance().get_tempdir(), os.path.basename(__file__))
        pathlib.Path(this_temp_dir).mkdir(exist_ok=True)

        #########################################################
        #
        # 1. Wrapper inputs, outputs and parameters
        #
        #########################################################
        Logger.instance().debug(
            "file: {}; line: {}; Wrapper inputs, outputs and parameters.".format(__file__,
                                                                                 inspect.currentframe().f_lineno, ))
        #
        # Input file
        input_file_fastainfo = self.input_file(TaxAssign.__input_file_fastainfo)
        input_file_taxonomy = self.input_file(TaxAssign.__input_file_taxonomy)
        #
        # Input table models
        marker_model = self.input_table(TaxAssign.__input_table_marker)
        run_model = self.input_table(TaxAssign.__input_table_run)
        biosample_model = self.input_table(TaxAssign.__input_table_biosample)
        filter_codon_stop_model = self.input_table(TaxAssign.__input_table_filter_codon_stop)
        variant_model = self.input_table(TaxAssign.__input_table_variant)
        # Output table models
        tax_assign_model = self.output_table(TaxAssign.__output_table_tax_assign)
        #
        # Options
        ltg_rule_threshold = float(self.option("ltg_rule_threshold"))  # percentage
        include_prop = float(self.option("include_prop"))  # percentage
        min_number_of_taxa = int(self.option("min_number_of_taxa"))  # count
        blast_db = str(self.option("blast_db"))  # count
        # num_threads = str(self.option("num_threads"))  # count

        ##########################################################
        #
        # 2. Read fastainfo to get run_id, marker_id, biosample_id, replicate for current analysis
        #
        ##########################################################
        fastainfo_df = pandas.read_csv(input_file_fastainfo, sep="\t", header=0, \
                                       names=['tag_forward', 'primer_forward', 'tag_reverse', 'primer_reverse',
                                              'marker_name', 'biosample_name', \
                                              'replicate', 'run_name', 'fastq_fwd', 'fastq_rev', 'fasta_path'])
        sample_instance_list = []
        for row in fastainfo_df.itertuples():
            marker_name = row.marker_name
            run_name = row.run_name
            biosample_name = row.biosample_name
            replicate = row.replicate
            with engine.connect() as conn:
                # get run_id ###########
                stmt_select_run_id = select([run_model.__table__.c.id]).where(run_model.__table__.c.name == run_name)
                run_id = conn.execute(stmt_select_run_id).first()[0]
                # get marker_id ###########
                stmt_select_marker_id = select([marker_model.__table__.c.id]).where(
                    marker_model.__table__.c.name == marker_name)
                marker_id = conn.execute(stmt_select_marker_id).first()[0]
                # get biosample_id ###########
                stmt_select_biosample_id = select([biosample_model.__table__.c.id]).where(
                    biosample_model.__table__.c.name == biosample_name)
                biosample_id = conn.execute(stmt_select_biosample_id).first()[0]
                # add this sample_instance ###########
                sample_instance_list.append({'run_id': run_id, 'marker_id': marker_id, 'biosample_id': biosample_id,
                                             'replicate': replicate})

        ##########################################################
        #
        # 3a. Select variants from this run/marker/biosample/replicate combination
        # 3b. Delete variants from this run/marker/biosample/replicate combination
        #
        ##########################################################
        #
        # 3a. Select variants from this run/marker/biosample/replicate combination
        #
        variant_id_delete_list = []
        for sample_instance in sample_instance_list:
            stmt_select_variant_id_delete = select([filter_codon_stop_model.__table__.c.variant_id]) \
                .where(filter_codon_stop_model.__table__.c.run_id == sample_instance['run_id']) \
                .where(filter_codon_stop_model.__table__.c.marker_id == sample_instance['marker_id']) \
                .where(filter_codon_stop_model.__table__.c.biosample_id == sample_instance['biosample_id']) \
                .where(filter_codon_stop_model.__table__.c.replicate == sample_instance['replicate'])
            # Select to DataFrame
            with engine.connect() as conn:
                for row in conn.execute(stmt_select_variant_id_delete).fetchall():
                    variant_id = row[0]
                    if not variant_id in variant_id_delete_list:
                        variant_id_delete_list.append(row)
        #
        # 3b. Delete variants from this run/markerbiosample/replicate combination
        #
        variant_instance_list = [{'variant_id': variant_id} for variant_id in variant_id_delete_list]
        with engine.connect() as conn:
            conn.execute(tax_assign_model.__table__.delete(), variant_instance_list)

        #
        ##########################################################
        #
        # Get variants that passed the filter
        #
        ##########################################################
        Logger.instance().debug(
            "file: {}; line: {}; Get variants and sequences that passed the filters".format(__file__,
                                                                                            inspect.currentframe().f_lineno,
                                                                                            'TaxAssign'))

        filter_codon_stop_model_table = filter_codon_stop_model.__table__
        variant_model_table = variant_model.__table__
        stmt_variant = select([filter_codon_stop_model_table.c.variant_id, variant_model_table.c.sequence]) \
            .where(filter_codon_stop_model_table.c.variant_id == variant_model_table.c.id) \
            .where(filter_codon_stop_model_table.c.filter_delete == 0).distinct().order_by("variant_id")
        # Select to DataFrame
        variant_list = []
        with engine.connect() as conn:
            for row in conn.execute(stmt_variant).fetchall():
                variant_list.append({'id': row.variant_id, 'sequence': row.sequence})
        variant_df = pandas.DataFrame.from_records(variant_list, index='id')

        # creation one fasta_path file containing all the variant
        #
        ##########################################################
        #
        # 2 Create FASTA file with Variants
        #
        ##########################################################
        Logger.instance().debug(
            "file: {}; line: {}; Create Fasta from Variants".format(__file__, inspect.currentframe().f_lineno))
        variant_fasta = os.path.join(this_temp_dir, 'variant.fasta')
        variant_df_utils = VariantDFutils(variant_df)
        variant_df_utils.to_fasta(variant_fasta)
        # VariantDFutils.to_fasta(variant_df, variant_fasta)
        #

        ##########################################################
        #
        # 3 Run local blast
        #
        ##########################################################
        Logger.instance().debug(
            "file: {}; line: {}; Running local blast with FASTA input {}".format(__file__,
                                                                                 inspect.currentframe().f_lineno,
                                                                                 variant_fasta))
        #
        # Run and read local blast result
        blast_output_tsv = os.path.join(this_temp_dir, 'blast_output.tsv')
        # blast_output_tsv = "/home/gonzalez/tmp/blast/blast_output.tsv" # uncomment for testing
        # get blast db dir and filename prefix from NHR file
        os.environ['BLASTDB'] = blast_db
        # if os.path.basename(map_taxids_tsv_path) == 'None': # run blast with full NCBI blast db
        #     blastn_cline = NcbiblastnCommandline(query=variant_fasta, db=blast_db_basename, evalue=1e-5,
        #                                          outfmt='"6 qseqid sacc pident evalue qcovhsp staxids"', dust='yes',
        #                                          qcov_hsp_perc=80, num_threads=1, out=blast_output_tsv)
        # else: # run blast with custom blast db
        # map_taxids_tsv_path, coi_blast_db_dir = download_coi_db()
        blastn_cline = NcbiblastnCommandline(query=variant_fasta, db='nt', evalue=1e-5,
                                             outfmt='"6 qseqid sacc pident evalue qcovhsp staxids"', dust='yes',
                                             qcov_hsp_perc=80, num_threads=threads, out=blast_output_tsv)
        Logger.instance().debug(
            "file: {}; line: {}; {}".format(__file__, inspect.currentframe().f_lineno, str(blastn_cline)))
        #
        # Run blast
        stdout, stderr = blastn_cline()

        ##########################################################
        #
        # Process blast reults
        #
        ##########################################################
        # if os.path.basename(map_taxids_tsv_path) == 'None': # Process result from full DB
        Logger.instance().debug(
            "file: {}; line: {}; Reading Blast output from: {}".format(__file__, inspect.currentframe().f_lineno, blast_output_tsv))
        blast_output_df = pandas.read_csv(blast_output_tsv, sep='\t', header=None,
                                          names=['variant_id', 'target_id', 'identity', 'evalue', 'coverage',
                                                 'target_tax_id'])
        # Remove null target tax ids
        blast_output_df = blast_output_df.loc[~blast_output_df.target_tax_id.isnull()]
        # expand multiple target_tax_ids
        blast_output_df.target_tax_id = blast_output_df.target_tax_id.astype('str') # first convert as string
        blast_output_df = (
            pandas.concat([blast_output_df, blast_output_df.target_tax_id.str.split(pat=';', n=1, expand=True)],
                          axis=1))
        # Select first tax_id
        blast_output_df = blast_output_df[['variant_id', 'target_id', 'identity', 'evalue', 'coverage', 0]]
        # rename first tax_id
        blast_output_df = blast_output_df.rename(columns={0: 'target_tax_id'})
        # Convert columns back to int
        blast_output_df.target_tax_id = blast_output_df.target_tax_id.astype('float')
        blast_output_df.target_tax_id = blast_output_df.target_tax_id.astype('int')
        # Blast output extract
        """   variant_id  target_id  identity        evalue  coverage  target_tax_id
0           2  MF7836761    99.429  1.620000e-86       100        1469487
1           2  MF7836761    99.429  1.620000e-86       100         189839
2           2  KY2618191    98.857  7.520000e-85       100         189839
3           2  MF7834791    98.857  7.520000e-85       100         189839
4           2  KU9559321    98.857  7.520000e-85       100         189839
"""
        ##########################################################
        #
        # Read target_tax_id
        # Compute lineages for each unique target_tax_id
        # Create a DF with these columns: tax_id and its lineage in wide format
        # Merge to the blast result
        #
        ##########################################################
        #
        Logger.instance().debug(
            "file: {}; line: {}; Open taxonomy.sqlite DB".format(__file__, inspect.currentframe().f_lineno))
        blast_output_df.target_tax_id = pandas.to_numeric(blast_output_df.target_tax_id)
        # getting the taxonomy_db to df
        # taxonomy_sqlite_path = __download_taxonomy_sqlite()
        taxonomy_sqlite_path = input_file_taxonomy
        taxonomy_db_df = f01_taxonomy_sqlite_to_df(taxonomy_sqlite_path)
        #
        Logger.instance().debug(
            "file: {}; line: {}; Annotate each target_tax_id with its lineage as columns in wide format".format(
                __file__, inspect.currentframe().f_lineno))
        lineage_list = []
        for target_tax_id in blast_output_df.target_tax_id.unique().tolist():
            lineage_list.append(f04_1_tax_id_to_taxonomy_lineage(target_tax_id, taxonomy_db_df))
        tax_id_to_lineage_df = pandas.DataFrame(lineage_list)
        #
        Logger.instance().debug(
            "file: {}; line: {}; Merge blast result including tax_id with their lineages".format(__file__,
                                                                                                 inspect.currentframe().f_lineno))
        # Merge local blast output with tax_id_to_lineage_df
        variantid_identity_lineage_df = blast_output_df.merge(tax_id_to_lineage_df, left_on='target_tax_id',
                                                               right_on='tax_id')
        variantid_identity_lineage_df.drop('tax_id', axis=1, inplace=True)
        variantid_identity_lineage_tsv = os.path.join(this_temp_dir, 'variantid_identity_lineage.tsv')
        variantid_identity_lineage_df.to_csv(variantid_identity_lineage_tsv, sep="\t", header=True)

        ##########################################################
        #
        #  6 test_f05_select_ltg_identity
        #
        ##########################################################
        Logger.instance().debug(
            "file: {}; line: {}; Main loop over variant and identity to"
            "compute the whole set of ltg_tax_id and ltg_rank for each variant_id"
            "to a dataframe".format(__file__, inspect.currentframe().f_lineno))
        #
        # f07_blast_result_to_ltg_tax_id(tax_lineage_df,ltg_rule_threshold=97, include_prop=90, min_number_of_taxa=3):
        # this function return a data frame containing the Ltg rank and Ltg Tax_id for each variant
        #
        ltg_df = f07_blast_result_to_ltg_tax_id(variantid_identity_lineage_df,
                                                ltg_rule_threshold=int(ltg_rule_threshold),
                                                include_prop=int(include_prop), min_number_of_taxa=min_number_of_taxa)
        ##########################################################
        #
        # 6. Insert Filter data
        #
        ##########################################################
        Logger.instance().debug("file: {}; line: {}; Insert variant_id, ltg_tax_id, ltg_rank to DB".format(__file__,
                                                                                                           inspect.currentframe().f_lineno))
        with engine.connect() as conn:
            conn.execute(tax_assign_model.__table__.insert(), ltg_df.to_dict('records'))
