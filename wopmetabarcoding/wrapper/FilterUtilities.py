from sqlalchemy import select
from wopmetabarcoding.utils.VSearch import VSearch1, Vsearch2, Vsearch3
import pandas, itertools
from Bio import SeqIO
import re
import os
from wopmetabarcoding.utils.constants import tempdir


class Variant2Sample2Replicate2Count():

    def __init__(self, variant2sample2replicate2count_df):
        self.df = variant2sample2replicate2count_df
        if self.df.shape[1] != 5:
            raise Exception('Columns missing in the variant2sample2replicate2count data frame!')
        self.indices_to_drop = []

    def lfn1_per_replicate(self, lfn_per_replicate_threshold):
        """
        Function calculating the Low Frequency Noise per sample replicate
        :param df: dataframe containing the information
        :param lfn_per_replicate_threshold: threshold defined by the user
        :return: List of the index which don't pass the filter
        """
        # Calculating the total of reads by sample replicates
        df2 = self.df.groupby(by=['sample_replicate']).sum()
        # Merge the column with the total reads by sample replicates for calculate the ratio
        df2 = self.df.merge(df2, left_on='sample_replicate', right_index=True)
        df2.columns = ['sequence', 'replicate', 'sample', 'sample_replicate', 'read_count_per_variant_per_sample_replicate', 'read_count_per_sample_replicate']
        # Calculate the ratio
        df2['low_frequence_noice_per_replicate'] = df2.read_count_per_variant_per_sample_replicate / df2.read_count_per_sample_replicate
        # Selecting all the indexes where the ratio is below the ratio
        indices_to_drop= list(df2.loc[df2.low_frequence_noice_per_replicate < lfn_per_replicate_threshold].index)
        self.indices_to_drop.append(indices_to_drop) # select rare reads to discard later
        return self.indices_to_drop

    def lfn2_per_variant(self, lfn_per_variant_threshold):
        """
        Function calculating the Low Frequency Noise per variant
        :param df: dataframe containing the information
        :param lfn_per_variant_threshold: threshold defined by the user
        :return: List of the index which don't pass the filter
        """
        # Calculating the total of reads by variant
        df2 = self.df.groupby(by=['sequence']).sum()
        # Merge the column with the total reads by variant for calculate the ratio
        df2 = self.df.merge(df2, left_on='sequence', right_index=True)
        df2.columns = ['sequence', 'replicate', 'sample', 'sample_replicate', 'read_count_per_variant_per_sample_replicate', 'read_count_per_variant']
        # Calculate the ratio
        df2['low_frequence_noice_per_variant'] = df2.read_count_per_variant_per_sample_replicate / df2.read_count_per_variant
        # Selecting all the indexes where the ratio is below the ratio
        indices_to_drop = list(df2.loc[df2.low_frequence_noice_per_variant < lfn_per_variant_threshold].index) # select rare reads to discard later
        self.indices_to_drop.append(indices_to_drop)
        return self.indices_to_drop

    def lfn2_per_replicate_series(self, lfn_per_replicate_series_threshold):
        """
        Function calculating the Low Frequency Noise per replicate series
        :param df: dataframe containing the information
        :param lfn_per_replicate_series_threshold: threshold defined by the user
        :return: List of the index which don't pass the filter
        """
        # Calculating the total of reads by replicate series
        df2 = self.df.groupby(by=['replicate']).sum()
        # Merge the column with the total reads by replicate series for calculate the ratio
        df2 = self.df.merge(df2, left_on='replicate', right_index=True)
        df2.columns = ['sequence', 'replicate', 'sample', 'sample_replicate',
                       'read_count_per_variant_per_sample_replicate', 'read_count_per_replicate_series']
        df2[
            'low_frequence_noice_per_replicate_series'] = df2.read_count_per_variant_per_sample_replicate / df2.read_count_per_replicate_series
        # Selecting all the indexes where the ratio is below the ratio
        indices_to_drop = list(df2.ix[
                        df2.low_frequence_noice_per_replicate_series < lfn_per_replicate_series_threshold].index)  #  select rare reads to discard later
        self.indices_to_drop.append(indices_to_drop)
        return self.indices_to_drop

    def lfn3_read_count(self, lfn_read_count_threshold):
        """
        Function calculating the Low Frequency Noise per users defined minimal readcount
        :param df: dataframe containing the information
        :param lfn_read_count_threshold: threshold defined by the user
        :return: List of the index which don't pass the filter
        """
        # Selecting all the indexes where the ratio is below the minimal readcount
        indices_to_drop = list(self.df.loc[self.df['count'] <= lfn_read_count_threshold].index)
        self.indices_to_drop.append(indices_to_drop)
        return self.indices_to_drop

    def lfn4_per_variant_with_cutoff(self, cutoff_tsv):
        """
        Function calculating the Low Frequency Noise per variant against cutoff
        :param df: dataframe containing the information
        :param cutoff_tsv: file containing the cutoffs for each variant
        :return: List of the index which don't pass the filter
        """

        cutoff_df = pandas.read_csv(cutoff_tsv, sep="\t", header=0)
        cutoff_df.columns = ['sequence', 'cutoff']
        #
        # Calculating the total of reads by variant
        df2 = self.df.groupby(by=['sequence']).sum()
        # Merge the column with the total reads by variant for calculate the ratio
        df2 = self.df.merge(df2, left_on='sequence', right_index=True)
        df2.columns = ['sequence', 'replicate', 'sample', 'sample_replicate', 'read_count_per_variant_per_sample_replicate', 'read_count_per_variant']
        df2['low_frequence_noice_per_variant'] = df2.read_count_per_variant_per_sample_replicate / df2.read_count_per_variant
        #
        # merge with cutoff
        df2 = df2.merge(cutoff_df, left_on="sequence", right_on="sequence")
        indices_to_drop = list(df2.ix[df2.low_frequence_noice_per_variant < df2.cutoff].index)
        self.indices_to_drop.append(indices_to_drop)
        return self.indices_to_drop

    def lfn4_per_replicate_series_with_cutoff(self, cutoff_tsv):
        """
            Function calculating the Low Frequency Noise per replicate series against cutoff
            :param df: dataframe containing the information
            :param cutoff_tsv: file containing the cutoffs for each variant
            :return: List of the index which don't pass the filter
            """
        cutoff_df = pandas.read_csv(cutoff_tsv, sep="\t", header=0)
        cutoff_df.columns = ['sequence', 'cutoff']
        #
        df2 = self.df.groupby(by=['replicate']).sum()
        df2 = self.df.merge(df2, left_on='replicate', right_index=True)
        df2.columns = ['sequence', 'replicate', 'sample', 'sample_replicate',
                       'read_count_per_variant_per_sample_replicate', 'read_count_per_replicate_series']
        df2[
            'low_frequence_noice_per_replicate_series'] = df2.read_count_per_variant_per_sample_replicate / df2.read_count_per_replicate_series
        #
        # merge with cutoff
        df2 = df2.merge(cutoff_df, left_on="sequence", right_on="sequence")
        indices_to_drop = list(df2.ix[df2.low_frequence_noice_per_replicate_series < df2.cutoff].index)
        self.indices_to_drop.append(indices_to_drop)
        return self.indices_to_drop

    def drop_indices(self):
        indices_to_drop = []
        for lists in self.indices_to_drop:
            indices_to_drop += lists
        self.df.drop(indices_to_drop, inplace=True)

    def min_repln(self, engine, replicate_model, marker_id, min_count):
        """
        Filter watching if a variant is present the minimal number of sample replicates of a sample
        :param engine: engine of the database
        :param replicate_model: model of the replicate table
        :param marker_id: id of the marker
        :param df: dataframe containing the information
        :param min_count: minimal number of replicates in which the variant must be present
        :return: None
        """
        replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
            replicate_model.marker_id == marker_id)
        replicate_obj = engine.execute(replicate_select)
        for replicate in replicate_obj:
            df2 = self.df.loc[self.df['sample'] == replicate.biosample_name]
            df3 = df2['sequence'].value_counts().to_frame()
            df3.columns = ['sequence_count']
            df2 = df2.merge(df3, left_on='sequence', right_index=True)
            index_to_drop = list(df2.ix[df2.sequence_count < min_count].index)
            self.df.drop(index_to_drop, inplace=True)

    def min_replp(self, engine, replicate_model, marker_id, replicate_count):
        """
        Filter watching if a variant is present the more more than 1/3 of the number of sample replicates of a sample
        :param engine: engine of the database
        :param replicate_model: model of the replicate table
        :param marker_id: id of the marker
        :param df: data frame containing the information
        :param replicate_count: number of sample_replicate in each
        :return: None
        """
        replicate_select = select([replicate_model.biosample_name, replicate_model.name]).where(
            replicate_model.marker_id == marker_id)
        replicate_obj = engine.execute(replicate_select)
        for replicate in replicate_obj:
            df2 = self.df.loc[self.df['sample'] == replicate.biosample_name]
            df3 = df2['sequence'].value_counts().to_frame()
            df3.columns = ['sequence_count']
            df2 = df2.merge(df3, left_on='sequence', right_index=True)
            index_to_drop = list(df2.ix[df2.sequence_count < ((1 / 3) * replicate_count)].index)
            self.df.drop(index_to_drop, inplace=True)

    @staticmethod
    def filter_fasta(variant_list, sample_fasta, chimera):
        """
        Function used to create fasta file from the dataframe
        :param variant_list: list of the variant present in the data frame
        :param sample_fasta: file name choose to name the fasta file
        :param chimera: choose if the function is used for the pcr_error or the chimera step
        :return: None
        """
        with open(sample_fasta, 'w') as fout:
            for variant in variant_list:
                if chimera is True:
                    fout.write(
                        ">{};size={};\n".format(variant, str(len(variant)))
                    )
                else:
                    fout.write(">{}\n".format(variant))
                fout.write(variant + '\n')

    def pcr_error(self, engine, replicate_model, marker_id, var_prop, pcr_error_by_sample):
        """
        Function used to eliminate variants from the data frame in which a pcr error is spotted
        :param engine: engine of the database
        :param replicate_model: model of the replicate table
        :param variant_model: model of the variants table
        :param df: data frame containing the information
        :param marker_id: id of the marker
        :param var_prop: threshold choose by the user
        :param pcr_error_by_sample: boolean used to choose if the analysis is made by sample or sample_replicate
        :return: None
        """
        select_replicate = select([replicate_model.biosample_name, replicate_model.name]) \
            .where(replicate_model.marker_id == marker_id)
        replicate_obj = engine.execute(select_replicate)

        for replicate in replicate_obj:

            if pcr_error_by_sample is True:
                data_variant = self.df.loc[self.df['sample'] == replicate.biosample_name]
                sample_fasta_name = 'data/output/Filter/{}_{}.fasta'.format(replicate.biosample_name, marker_id)

            else:
                sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
                data_variant = self.df.loc[self.df['sample_replicate'] == sample_replicate]
                sample_fasta_name = 'data/output/Filter/{}_{}.fasta'.format(sample_replicate, marker_id)

            variant_list_series = data_variant['sequence']
            if variant_list_series.empty is False:
                variant_list = sorted(set(list(variant_list_series)))
                self.filter_fasta(variant_list, sample_fasta_name, False)
                shortest_sequence = min(variant_list)
                L = len(shortest_sequence)
                id_raw = (L - 1) / L
                id_rounded = round(id_raw, 2)
                sample_tsv_name = sample_fasta_name.replace('.fasta', '.tsv')
                vsearch_usearch_global_args = {'db': sample_fasta_name,
                                               'usearch_global': sample_fasta_name,
                                               'id': str(id_rounded),
                                               'maxrejects': 0,
                                               'maxaccepts': 0,
                                               'userout': sample_tsv_name,
                                               'userfields': "query+target+alnlen+ids+mism+gaps",
                                               }
                vsearch_usearch_global = VSearch1(**vsearch_usearch_global_args)
                vsearch_usearch_global.run()
                column_names = ['query', 'target', 'alnlen', 'ids', 'mism', 'gaps']
                df2 = pandas.read_csv(sample_tsv_name, sep='\t', names=column_names)
                false_df2 = list(df2.ix[(df2.mism + df2.gaps) != 1].index)
                df2.drop(false_df2, inplace=True)
                df3 = df2[['query', 'target']]
                df4 = data_variant[['sequence', 'count']]
                df3 = df3.merge(df4, left_on=['query'], right_on=['sequence'])
                df3 = df3.merge(df4, left_on=['target'], right_on=['sequence'])
                df3['count_ratio'] = df3.count_x / df3.count_y
                df3 = df3[['query', 'count_ratio']]
                df3.index = df3['query']
                data_variant = data_variant.merge(df3, left_on='sequence', right_index=True)
                index_to_drop = list(data_variant.ix[data_variant.count_ratio < var_prop].index)
                # delete_variant = self.df(index_to_drop)
                self.df.drop(index_to_drop, inplace=True)

    def chimera(self, replicate_obj_list, marker_id, chimera_by):
        """
        Function using Vsearch to identify chimera or borderline variant
        :param replicate_obj_list: list of the sample and the corresponding sample replicate
        :param df: data frame containing the information
        :param marker_id: id of the marker
        :param chimera_by: parameter to choose between sample, sample replicate or all
        :return: None
        """

        for replicate in replicate_obj_list:

            if chimera_by == "sample":
                df2 = self.df.loc[self.df['sample'] == replicate.biosample_name]
                filename = '{}_{}_repl.fasta'.format(replicate.biosample_name, marker_id)
                repl_fasta_name = os.path.join(tempdir, filename)

            elif chimera_by == "sample_replicate":
                sample_replicate = '{}_{}'.format(replicate.biosample_name, replicate.name)
                df2 = self.df.loc[self.df['sample_replicate'] == sample_replicate]
                # repl_fasta_name = 'data/output/Filter/{}_{}_repl.fasta'.format(sample_replicate, marker_id)
                filename = '{}_{}_repl.fasta'.format(sample_replicate, marker_id)
                repl_fasta_name = os.path.join(tempdir, filename)

            sorted_repl_fasta_name = repl_fasta_name.replace('.fasta', '_sorted.fasta')
            variant_list_series = df2['sequence']
            if variant_list_series.empty is False:
                variant_list_series = df2['sequence']
                variant_list = sorted(set(list(variant_list_series)))
                self.filter_fasta(variant_list, repl_fasta_name, True)
                vsearch_sortbysize_args = {"sortbysize": repl_fasta_name, "output": sorted_repl_fasta_name}
                vsearch_sortbysize = Vsearch2(**vsearch_sortbysize_args)
                vsearch_sortbysize.run()
                del vsearch_sortbysize
                borderline_repl_filename = sorted_repl_fasta_name.replace('sorted', 'borderline')
                nonchimera_repl_filename = sorted_repl_fasta_name.replace('sorted', 'correct')
                chimera_repl_filename = sorted_repl_fasta_name.replace('sorted', 'chimera')
                vsearch_chimera_args = {
                    "uchime_denovo": sorted_repl_fasta_name,
                    "borderline": borderline_repl_filename,
                    "nonchimeras": nonchimera_repl_filename,
                    "chimeras": chimera_repl_filename
                }
                vsearch_chimera = Vsearch3(**vsearch_chimera_args)
                vsearch_chimera.run()
                del vsearch_chimera
                with open(chimera_repl_filename, 'r') as fin:
                    for line in fin:

                        if ">" in line:
                            line = line.replace('>', '')
                            variant_sequence = line.strip().split(';')[0]
                            df_variant = df2.loc[df2['sequence'] == variant_sequence]

                            if df_variant.empty is False:
                                indexs_to_drop = df_variant.index.tolist()
                                self.df.drop(indexs_to_drop, inplace=True)
                            else:
                                pass
                borderline_variants = []
                for record in SeqIO.parse(borderline_repl_filename, 'fasta'):
                    borderline_variants.append(record.description.strip().split(';')[0])
                self.df['is_borderline'] = (self.df['sequence'] in borderline_variants)

    def renkonen(self, number_of_replicate, renkonen_tail):
        """
        Function calculating distance between sample replicate of the sample to determine if the sample replicate
        must be deleted or not
        :param engine: engine of the database
        :param replicate_model: model of the replicate table
        :param df: data frame containing the information
        :param marker_id: id of the marker
        :param number_of_replicate: Number of replicate by sample
        :return: None
        """
        columns_name = ['repl_i', 'repl_j', 'distance']
        df_read_count_per_sample_replicate = self.df.groupby(by=['sample_replicate'])['count'].sum()
        df_read_count_per_sample_replicate = df_read_count_per_sample_replicate.to_frame()
        df_read_count_per_sample_replicate.columns = ['replicate_count']
        df_read_count_per_sample_replicate = self.df.merge(df_read_count_per_sample_replicate, left_on='sample_replicate', right_index=True)
        df_read_count_per_sample_replicate['proportion'] = df_read_count_per_sample_replicate['count'] / df_read_count_per_sample_replicate['replicate_count']
        # df_replicate = df_read_count_per_sample_replicate.groupby(by=['sample'])['sample_replicate'].to_frame()
        samples = df_read_count_per_sample_replicate['sample']
        samples = list(set(samples.tolist()))
        for sample in samples:
            df_permutation_distance = pandas.DataFrame(columns=columns_name)
            df_replicate = df_read_count_per_sample_replicate.loc[df_read_count_per_sample_replicate['sample'] == sample]
            replicates = list(set(df_replicate['sample_replicate'].tolist()))
            for combi in itertools.permutations(replicates, 2):
                combi = list(combi)
                df_repli = df_replicate.loc[df_replicate['sample_replicate'] == combi[0]]
                data_repli = df_repli[['sequence', 'sample_replicate', 'proportion']]
                df_replj = df_replicate.loc[df_replicate['sample_replicate'] == combi[1]]
                data_replj = df_replj[['sequence', 'sample_replicate', 'proportion']]
                df_replij = data_repli.append(data_replj)
                group_repl = df_replij.groupby(by=['sequence'])['proportion'].min()
                distance = 1 - group_repl.sum()
                query = [combi[0], combi[1], distance]
                df_permutation_distance.loc[len(df_permutation_distance)] = query
            # df_calc = df_permutation_distance.loc[df_permutation_distance['repl_i'] == combi[0]]
            indices_to_drop = list(
                df_permutation_distance.loc[df_permutation_distance.distance > renkonen_tail].index.tolist()
            )
            df_permutation_distance.drop(indices_to_drop, inplace=True)
            repl_list = list(set(df_permutation_distance['repl_i'].tolist()))
            for repl_i in repl_list:
                df_calc = df_permutation_distance.loc[df_permutation_distance['repl_i'] == repl_i]
                if len(df_calc) > ((number_of_replicate -1)/2):
                    indices_to_drop = self.df.loc[self.df['sample_replicate'] == repl_i].index.tolist()
                    self.df.drop(indices_to_drop, inplace=True)

        # TODO Do not use database if data is already in dataframe
        # Todo Create df_with_permutations: sample, replicate_i, replicate_j
        # TODO Create permutations [1,2], [1,3], ...
        # TODO For each sample and permutation combination [1,2] in df_with_permutations
        # TODO Create df1 <- sample, repl_i, sequence, prop_seq
        # TODO Create df2 <- sample, repl_j, sequence, prop_seq
        # TODO Merge df1.merge(df2, left_on=[sample, sequence], right_on=[sample, sequence])

    def indel(self, delete_var):
        """

        :param delete_var:
        :return:
        """
        df2 = self.df.copy()
        df2['len_sequence'] = len(df2.sequence)
        df2['modulo3'] = len(df2.sequence)%3
        modulo3 = df2['modulo3'].tolist()
        sequence_length_max = max(modulo3, key=modulo3.count)
        if delete_var:
            indices_to_drop = df2.loc[df2.modulo3 != sequence_length_max].index.tolist()
            self.df.drop(indices_to_drop, inplace=True)
        else:
            df2['is_pseudogene_indel'] = (df2.modulo3 != sequence_length_max)
            is_pseudogene = df2['is_pseudogene_indel']
            self.df = self.df.merge(is_pseudogene.to_frame(), left_index=True, right_index=True)

    def codon_stop(self, df_codon_stop_per_genetic_code, genetic_code, delete_var):
        """
        Function searching stop codon the different reading frame of a sequence and tag/delete variant if there is
        stop codon in the 3 reading frames
        :param df_codon_stop_per_genetic_code: data frame which contains all codon stop classified by genetic code
        :param genetic_code: genetic code to search stop codon in the df_codon_stop_per_genetic_code dataframe
        :param delete_var: option which define if the variants must be deleted or tagged
        :return: void
        """
        sequences = self.df['sequence'].tolist()
        # Get the list of codon stop according to the genetic code given by the user
        df2 = df_codon_stop_per_genetic_code.loc[df_codon_stop_per_genetic_code['genetic_code'] == genetic_code]
        codon_stop_list = df2['codon'].tolist()
        indices_to_drop = []
        # Take sequence 1 by 1
        for sequence in sequences:
            # Keep the entire sequence for select in the dataframe
            entire_sequence = sequence
            # define the current reading_frame
            reading_frame = 1
            # Number of times when a codon stop is fin in the sequence
            fail = 0
            for i in range(3):
                # if a codon stop is not find in the first reading frame no need to search the others
                if reading_frame > 1 and fail == 0:
                    break
                # Get a new sequence eliminating the first nucleotide
                if reading_frame == 2:
                    sequence = sequence[1:]
                # Get a new sequence eliminating the first and the second nucleotides
                elif reading_frame == 3:
                    sequence = sequence[2:]
                # Divide the sequence in length 3 codons
                codon_list = re.findall('...', sequence)
                codon_stop_count = 0
                # Count the number of codon stop in the sequence length
                for codon_stop in codon_stop_list:
                    if codon_stop in codon_list:
                        codon_stop_count += 1
                # If a codon stop or more are found then 1 fail is added
                if codon_stop_count > 0:
                    print(entire_sequence)
                    fail += 1
                reading_frame += 1
            # If a codon stop is find in every reading frames, all the indexes where the variant sequence appeared are
            # kept and will be tagged or deleted at the loop end
            if fail == 3:
                indices_to_drop.append(self.df.loc[self.df['sequence'] == entire_sequence].index().tolist())
        # Depending on the user choice the variant will be tagged in dataframe or removed from it
        if delete_var:
            self.df.drop(indices_to_drop, inplace=True)
        else:
            self.df['is_pseudogene_codon_stop'] = (self.df.index in indices_to_drop)

    def consensus(self):
        """
        Function used to display the read average of the remaining variant
        :return:
        """
        variants_sequences = self.df["sequence"]
        variants_sequences =list(set(variants_sequences))
        read_average_columns = ['variant', 'read_average']
        read_average_df = pandas.DataFrame(columns=read_average_columns)
        for variant in variants_sequences:
            variant_df = self.df.loc[self.df["sequence"] == variant]
            read_average = round(variant_df["count"].sum()/len(variant_df['count']), 0)
            read_average_df.loc[len(read_average_df)] = [variant, read_average]
        self.df = self.df.merge(read_average_df, left_on='sequence', right_on='variant')
        self.df = self.df.drop(columns=['variant'])

    def filtered_variants(self):
        return self.df


































