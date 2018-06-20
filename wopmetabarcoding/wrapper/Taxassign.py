from wopmars.framework.database. tables.ToolWrapper import ToolWrapper
from wopmars.utils.Logger import Logger
from wopmetabarcoding.wrapper.TaxassignUtilities import alignment_vsearch, create_phylogenetic_line_df, dataframe2ltgdefinition, rank_hierarchy, seq2tax_db_sqlite_to_df
import pandas
from Bio import SeqIO


class Taxassign(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.Taxassign"
    }
    __input_file_taxassign_db = "taxassign_db_fasta"
    __filtered_dataframe_path = "filtered_dataframe_path"
    __assignlvl2id = "assignlvl2id"
    __tax_assign_db_sqlite = "tax_assign_db_sqlite"
    __default_output = 'default_output'

    def specify_input_file(self):
        return [
            Taxassign.__input_file_taxassign_db,
            Taxassign.__filtered_dataframe_path,
            Taxassign.__assignlvl2id,
            Taxassign.__tax_assign_db_sqlite

        ]

    def specify_output_file(self):
        return [
            Taxassign.__default_output
        ]

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        # Input files
        taxassign_db_fasta = self.input_file(Taxassign.__input_file_taxassign_db)
        filter_output = self.input_file(Taxassign.__filtered_dataframe_path)
        tax_assign_pars_tsv = self.input_file(Taxassign.__assignlvl2id)
        tax_assign_sqlite = self.input_file(Taxassign.__tax_assign_db_sqlite)
        #
        default_output = self.output_file(Taxassign.__default_output)
        with open(filter_output, 'r') as fin:
            for line in fin:
                line = line.strip().split('\t')
                filtered_variants_fasta = line[2]
                output_tsv = filtered_variants_fasta.replace('.fasta', '.tsv')
                # output_tsv = output_tsv.replace(tempdir, '/tmp/tmpe6yiaf0x/')
                for record in SeqIO.parse(filtered_variants_fasta, 'fasta'):
                    sequence_fasta = "data/sequence_fasta.fasta"
                    with open(sequence_fasta, 'w') as fin_sequence_fasta:
                        sequence_id = str(record.description)
                        fin_sequence_fasta.write(">" + sequence_id + "\n")
                        sequence = str(record.description)
                        fin_sequence_fasta.write(sequence + "\n")
                    alignment_vsearch(sequence_fasta, taxassign_db_fasta, output_tsv)
                    # import pdb; pdb.set_trace()
                    vsearch2seq2tax_df = pandas.read_csv(output_tsv, sep="\t", header=None, index_col=1)
                    vsearch2seq2tax_df.columns = ["var_seq", "alignment_identity"]
                    vsearch2seq2tax_df.index.name = 'tax_seq_id'
                    #
                    tax_seq_id_list = vsearch2seq2tax_df.index.tolist()
                    #
                    seq2tax_df = seq2tax_db_sqlite_to_df(tax_assign_sqlite, tax_seq_id_list)
                    #
                    # tax_assign_pars df
                    names = ["identity_threshold", "min_tax_level", "max_tax_resolution", "min_tax_n"]
                    tax_assign_pars_df = pandas.read_csv(tax_assign_pars_tsv, sep="\t", header=None, names=names)
                    #
                    # Merge of the vsearch alignment, the sequence and taxa information
                    vsearch2seq2tax_df = pandas.merge(vsearch2seq2tax_df, seq2tax_df, left_index=True,
                                                      right_on="tax_seq_id")
                    vsearch2seq2tax_df = vsearch2seq2tax_df.assign(
                        rank_id=vsearch2seq2tax_df.rank_name.apply(lambda x: rank_hierarchy.index(x)))
                    #
                    # Loop over each identity threshold
                    for tax_assign_pars_df_row_i, tax_assign_pars_df_row in tax_assign_pars_df.iterrows():
                        # identity_threshold, min_tax_level, max_tax_resolution, min_tax_n = 100, "species", "subspecies", 1
                        # identity_threshold, min_tax_level, max_tax_resolution, min_tax_n = 97, "genus", "species", 1
                        # identity_threshold, min_tax_level, max_tax_resolution, min_tax_n = 95, "family", "species", 3
                        # identity_threshold, min_tax_level, max_tax_resolution, min_tax_n = 90, "order", "family", 3
                        # identity_threshold, min_tax_level, max_tax_resolution, min_tax_n = 85, "order", "order", 3
                        # identity_threshold, min_tax_level, max_tax_resolution, min_tax_n = 80, "class", "order", 5
                        identity_threshold, min_tax_level, max_tax_resolution, min_tax_n = tax_assign_pars_df_row.tolist()
                        Logger.instance().info("Selecting sequences with " + str(identity_threshold) + "% identity.")
                        min_tax_level_id = rank_hierarchy.index(min_tax_level)
                        max_tax_resolution_id = rank_hierarchy.index(max_tax_resolution)
                        #
                        # test identity_threshold
                        vsearch2seq2tax_df_selected = vsearch2seq2tax_df.loc[
                            vsearch2seq2tax_df.alignment_identity >= identity_threshold]
                        if vsearch2seq2tax_df_selected.empty:  #  no lines selected at this alignment identity threshold
                            Logger.instance().info(
                                "Any sequences are selected passing to next identity threshold."
                            )
                            continue  #  next identity threshold
                        #  continue only if selected lines
                        #
                        # test min_tax_level
                        vsearch2seq2tax_df_selected = vsearch2seq2tax_df_selected.loc[
                            vsearch2seq2tax_df_selected.rank_id >= min_tax_level_id]
                        if vsearch2seq2tax_df_selected.empty:  #  no lines selected at this alignment identity threshold
                            Logger.instance().info(
                                "Any sequence with enought detailled taxonomic "
                                "level found, passing to next identity threshold."
                            )
                            continue  #  next identity threshold
                        #  continue only if selected lines
                        #
                        # test min_tax_n
                        if vsearch2seq2tax_df_selected.shape[0] < min_tax_n:
                            Logger.instance().info(
                                "Not enought sequences are selected passing to next identity threshold."
                            )
                            continue  #  next identity threshold
                        # continue only if selected lines
                        tax_seq_id_list = vsearch2seq2tax_df_selected.tax_seq_id.tolist()
                        #
                        # Create lineage df
                        tax_lineage_df = create_phylogenetic_line_df(tax_seq_id_list, tax_assign_sqlite)
                        #
                        #  Search LTG
                        tax_count_perc = dataframe2ltgdefinition(tax_lineage_df)
                        #
                        # test min_tax_n
                        if tax_count_perc.empty:
                            Logger.instance().info(
                                "Any taxonomic level with the given proportion to become LTG."
                            )
                            continue  #  next identity threshold
                        tax_count_perc['rank_index'] = [rank_hierarchy.index(rank_name) for rank_name in
                                                        tax_count_perc.index.tolist()]
                        #
                        #  Criteria: lineage df tax id more detailed than
                        tax_count_perc.loc[tax_count_perc.rank_index >= min_tax_level_id]
                        tax_count_perc_ltg = tax_count_perc.loc[tax_count_perc.rank_index >= min_tax_level_id]
                        if tax_count_perc_ltg.empty:
                            Logger.instance().info(
                                "Nothing survive."
                            )
                            continue
                        ltg_tax_id = tax_count_perc_ltg.tax_id.tolist()[-1]
                        ltg_rank_id = tax_count_perc_ltg.rank_index.tolist()[-1]
                        #
                        if ltg_rank_id > max_tax_resolution_id:  #  go up in lineage of ltg_tax_id up to max_tax_resolution_id
                            ltg_tax_id = tax_lineage_df.loc[
                                tax_lineage_df[rank_hierarchy[ltg_rank_id]] == ltg_tax_id, max_tax_resolution].unique()
                        print(ltg_tax_id)






