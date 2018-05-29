import pandas
from Bio import SeqIO


def create_info_df(fasta_name):
    columns_name = ['target', 'tax_id', 'name', 'rank', 'parent_id']
    df_info = pandas.DataFrame(columns=columns_name)
    for record in SeqIO.parse(fasta_name, 'fasta'):
        sequence_id = record.description.strip().split('=')
        seq_name = sequence_id[0].replace(' name', '')
        name = sequence_id[1].replace(' tax_id', '')
        tax_id = sequence_id[2].replace(' rank', '')
        rank = sequence_id[3].replace(' parent_taxid', '')
        parent_id = sequence_id[4]
        sequence_info = [seq_name, name, tax_id, rank, parent_id]
        df_info.loc[len(df_info)] = sequence_info
    return df_info

# test1 = open(output_tsv, 'w')
                # test1.close()
                # vsearch_usearch_global_args = {'db': taxassign_db_fasta,
                #                                'usearch_global': filtered_variants_fasta,
                #                                'id': str(0.80),
                #                                'maxrejects': 0,
                #                                'maxaccepts': 0,
                #                                'userout': output_tsv,
                #                                'userfields': "--userfields query+target+id --id",
                #                                }
                # vsearch_1 = VSearch1(**vsearch_usearch_global_args)
                # vsearch_1.run()