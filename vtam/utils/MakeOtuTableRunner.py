import inspect

import pandas
import sqlalchemy

from vtam.utils.Logger import Logger
from vtam.model.Biosample import Biosample
from vtam.utils.FastaInformation import FastaInformation
from vtam.utils.VariantReadCountDF import VariantReadCountDF


class MakeOtuTableRunner(object):

    def __init__(self, engine, input_file_fastainfo, run_model, marker_model, biosample_model, replicate_model, filter_chimera_borderline_model, filter_codon_stop_model, variant_model, tax_assign_model):

        self.engine = engine
        self.input_file_fastainfo = input_file_fastainfo
        self.run_model = run_model
        self.marker_model = marker_model
        self.biosample_model = biosample_model
        self.replicate_model = replicate_model
        self.filter_chimera_borderline_model = filter_chimera_borderline_model
        self.filter_codon_stop_model = filter_codon_stop_model
        self.variant_model = variant_model
        self.tax_assign_model = tax_assign_model

    def run(self):

        fasta_info = FastaInformation(self.input_file_fastainfo, self.engine, self.run_model, self.marker_model, self.biosample_model, self.replicate_model)
        variant_read_count_df = fasta_info.get_variant_read_count_df(self.filter_codon_stop_model)
        variant_df = fasta_info.get_variant_df(variant_read_count_like_model=self.filter_codon_stop_model,
                                               variant_model=self.variant_model)

        biosample_df = fasta_info.get_biosample_df(variant_read_count_like_model=self.filter_codon_stop_model)
        marker_df = fasta_info.get_marker_df(variant_read_count_like_model=self.filter_codon_stop_model)
        run_df = fasta_info.get_run_df(variant_read_count_like_model=self.filter_codon_stop_model)

        # Aggregate replicates
        variant_read_count_obj = VariantReadCountDF(variant_read_count_df)
        N_ij_df = variant_read_count_obj.get_N_ij_df()

        # TOdo need multiindex
        otu_df = N_ij_df.pivot_table(index='variant_id', columns='biosample_id', values='N_ij', fill_value=0)
        import pdb; pdb.set_trace()

        # Add marker id and replace it with name
        otu_df = variant_read_count_df[['marker_id', 'variant_id']].drop_duplicates(inplace=False).merge(otu_df,
                                                                              on='variant_id', validate='one_to_one')
        otu_df = otu_df.merge(marker_df, left_on = 'marker_id', right_index = True)

        # Add run id and replace it with name
        otu_df = otu_df.merge(variant_read_count_df[['run_id', 'variant_id']].drop_duplicates(inplace=False),
                              right_on='variant_id', left_index=True, validate='one_to_one')
        import pdb; pdb.set_trace()
        otu_df = otu_df.merge(run_df, left_on='run_id', right_index = True)

        # Add sequence
        otu_df = otu_df.merge(variant_df, left_on='variant_id', right_index=True, validate='one_to_one')
        otu_df['sequence_length'] = otu_df.sequence.apply(lambda x: len(x))
        import pdb; pdb.set_trace()

        # Chimera
        variant_to_chimera_borderline_df = self.get_chimera_borderline_df()
        otu_df = otu_df.merge(variant_to_chimera_borderline_df, left_on='variant_id', right_index=True, validate='one_to_one')

        # Run name



    def get_chimera_borderline_df(self):
        # #########################################################
        #
        # Get variants that passed the filter
        # Get also chimera borderline information
        #
        # #########################################################
        Logger.instance().debug(
            "file: {}; line: {}; Get variants and sequences that passed the filters".format(__file__, inspect.currentframe().f_lineno,'TaxAssign'))

        filter_codon_stop_model_table = self.filter_codon_stop_model.__table__
        filter_chimera_borderline_model_table = self.filter_chimera_borderline_model.__table__
        variant_model_table = self.variant_model.__table__
        stmt_filter_codon_stop = sqlalchemy.select([
                            filter_codon_stop_model_table.c.variant_id,
                            filter_chimera_borderline_model_table.c.filter_delete, ]) \
            .where(filter_chimera_borderline_model_table.c.variant_id == variant_model_table.c.id) \
            .where(filter_codon_stop_model_table.c.filter_delete == 0).distinct()
        # Select to DataFrame
        variant_to_chimera_borderline_list = []
        with self.engine.connect() as conn:
            for row in conn.execute(stmt_filter_codon_stop).fetchall():
                variant_to_chimera_borderline_list.append(row)
        variant_to_chimera_borderline_df = pandas.DataFrame.from_records(variant_to_chimera_borderline_list, index='variant_id', columns=['variant_id', 'chimera_borderline'])
        return variant_to_chimera_borderline_df
