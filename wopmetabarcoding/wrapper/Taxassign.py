from wopmars.framework.database. tables.ToolWrapper import ToolWrapper
from sqlalchemy import select
from wopmetabarcoding.utils.constants import tempdir, order, taxonomic_levels
from wopmetabarcoding.wrapper.TaxassignUtilities import alignment_vsearch, most_common
from wopmetabarcoding.wrapper.FilterUtilities import Variant2Sample2Replicate2Count
from wopmetabarcoding.utils.VSearch import VSearch1
import os, pickle, pandas
from Bio import SeqIO


class Taxassign(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "wopmetabarcoding.wrapper.Taxassign"
    }
    __input_table_marker = "Marker"
    __input_table_variant = "Variant"
    __input_file_taxassign_db = "taxassign_db_fasta"
    __filtered_dataframe_path = "filtered_dataframe_path"
    __assignlvl2id = "assignlvl2id"
    __default_output = 'default_output'


    def specify_input_table(self):
        return [
            Taxassign.__input_table_marker,
            Taxassign.__input_table_variant
        ]

    def specify_input_file(self):
        return [
            Taxassign.__input_file_taxassign_db,
            Taxassign.__filtered_dataframe_path,
            Taxassign.__assignlvl2id

        ]

    def specify_output_file(self):
        return [
            Taxassign.__default_output
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
        assignlvl2id = self.input_file(Taxassign.__assignlvl2id)
        df_assignlvl2id = pandas.read_csv(assignlvl2id, sep='\t', names=["id", "min_target_taxlevel", "max_tax_resolution", "min_taxon_n"])
        # df_assignlvl2id.columns = ["id", "min_target_taxlevel", "max_tax_resolution", "min_taxon_n"]
        #
        taxassign_proportion = 0.8
        default_output = self.output_file(Taxassign.__default_output)
        with open(filter_output, 'r') as fin:
            for line in fin:
                line = line.strip().split('\t')
                marker_name = line[0]
                dataframe_path = line[1]
                filtered_variants_fasta = line[2]
                variant2sample2replicate2count_df = pandas.read_csv(dataframe_path, sep='\t')
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
                    import pdb; pdb.set_trace()
                    variants = list(set(variant2sample2replicate2count_df["variant_seq"].tolist()))
                    df = pandas.read_csv(output_tsv, sep='\t', names=['query', 'target', 'id'])
                    print(df)
                    # df_info = create_info_df(taxassign_db_fasta)
                    df2 = df.loc[df['query'] == record.description]
                    for ids in order:
                        df3 = df2.loc[df2['id'] >= ids]
                        targets = list(set(df3['target'].tolist()))
                        with open(taxassign_db_fasta, 'r') as fin:
                            with open(os.path.join(tempdir + 'temp.tsv'), 'w') as fout:
                                for line in fin:
                                    if '>' in line:
                                        linebis = line.split(' ')
                                        seq_id = int(linebis[0].replace(">", ''))
                                        if seq_id in targets:
                                            # 5164832 name=Papestra cristifera tax_id=1485856 rank=species parent_taxid=685411
                                            linefinal = line.strip().split('=')
                                            seq_name = linefinal[0].replace(' name', '')
                                            seq_name = seq_name.replace('>', '')
                                            name = linefinal[1].replace(' tax_id', '')
                                            tax_id =linefinal[2].replace(' rank', '')
                                            rank = linefinal[3].replace(' parent_taxid', '')
                                            parent_id = linefinal[4]
                                            fout.write(seq_name + '\t' + name + "\t" + tax_id + "\t" + rank + "\t" + parent_id + "\n")
                        df4 = pandas.read_csv(os.path.join(tempdir + 'temp.tsv'), sep='\t', names=['sequence_name', 'name', 'tax_id', 'rank', 'parent_id'])
                        # df4.columns = ['sequence_name', 'name', 'tax_id', 'rank', 'parent_id']
                        # rank_list = df4['rank'].tolist()
                        # most_common = max(rank_list, key = rank_list.count)
                        assign_info = df_assignlvl2id.loc[df_assignlvl2id['id'] == ids]
                        min_target_taxlevel = assign_info['min_target_taxlevel'].tolist()
                        min_target_taxlevel = min_target_taxlevel[0]
                        max_tax_resolution = assign_info['max_tax_resolution'].tolist()
                        max_tax_resolution = max_tax_resolution[0]
                        name_taxon = df4['name'].tolist()
                        most_common_name = most_common(name_taxon)
                        print((name_taxon.count(most_common_name)/len(name_taxon) >= 0.80))
                        df5 = df4.loc[df4['rank'] == min_target_taxlevel]

                    del df




