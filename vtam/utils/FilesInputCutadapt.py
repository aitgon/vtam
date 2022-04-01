import os
import pandas as pd
try:
    from Bio.Alphabet import generic_dna
except ImportError:
    generic_dna = None

from Bio.Seq import Seq


class FilesInputCutadapt(object):
    """ make fasta files to use as input for cutadapt demultiplexing """

    def __init__(self, file_path, mergedFasta, no_reverse, tag_to_end):

        self.file_path = file_path
        self.mergedFasta = mergedFasta
        
        self.tag_to_end = tag_to_end
        self.no_reverse = no_reverse

        self.df = pd.read_csv(file_path, sep='\t', header=0)
        self.df.columns = self.df.columns.str.lower().copy()

        self.selected_df = self.df.loc[self.df['mergedfasta'] == self.mergedFasta].copy()
        self.length = self.selected_df.shape[0]
        self.dict = self.selected_df.to_dict(orient='list')
        self.mergedfasta_list = self.dict["mergedfasta"]


    def tags_file(self):
        
        relPath, _ = os.path.split(self.file_path)
        self.tagsFile = os.path.join(relPath, 'tagsFile.fasta')
        sample_names = self.get_sample_names()

        with open(self.tagsFile, 'at') as tags:

            for sample_name in sample_names:

                name, _, _, fwd, rev = sample_name

                
                if not "_reversed" in name:
                    if generic_dna:  # Biopython <1.78
                        tagRevRC = str(Seq(rev, generic_dna).reverse_complement())
                    else:  # Biopython =>1.78
                        tagRevRC = str(Seq(rev).reverse_complement())
                    
                    if not self.tag_to_end:
                        tags.write(f">{name}\n^{fwd}...{tagRevRC}$\n")
                    else:
                        tags.write(f">{name}\n{fwd};min_overlap={str(len(fwd))}...{tagRevRC};min_overlap={str(len(tagRevRC))}\n")
                    
                else:
                    if generic_dna:  # Biopython <1.78
                        tagFwdRC = str(Seq(fwd, generic_dna).reverse_complement())
                    else:  # Biopython =>1.78
                        tagFwdRC = str(Seq(fwd).reverse_complement())
                
                    if not self.tag_to_end:
                        tags.write(f">{name}\n^{rev}...{tagFwdRC}$\n")
                    else:
                        tags.write(f">{name}\n{rev};min_overlap={str(len(rev))}...{tagFwdRC};min_overlap={str(len(tagFwdRC))}\n")

        return self.tagsFile


    def primers(self):
                
        self.primers_result = []

        for i in range(self.length):
            primerFwd =f'{self.dict["primerfwd"][i]}'
            primerRev =f'{self.dict["primerrev"][i]}'
            marker =f'{self.dict["marker"][i]}'            
            primer = (marker, primerFwd, primerRev, len(primerFwd), len(primerRev))
            if not primer in self.primers_result:
                self.primers_result.append(primer)

        return self.primers_result


    def get_sample_names(self):
        
            sample_names = []
            for i in range(self.length):
                sample = f'{self.dict["sample"][i]}'
                marker = f'{self.dict["marker"][i]}'
                fwd = f'{self.dict["tagfwd"][i]}'
                rev = f'{self.dict["tagrev"][i]}'

                name = (f"marker[{marker}]sample[{sample}]fwd[{fwd}]rev[{rev}]", marker, sample, fwd, rev)

                if name not in sample_names:
                    sample_names.append(name)

            if self.no_reverse:
                samples_names_revered = []
                for name in sample_names:
                    samples_names_revered.append(name)
                    name_reversed = (name[0] + "_reversed",  name[1], name[2],  name[3], name[4])
                    samples_names_revered.append(name_reversed)
                return samples_names_revered
                
            return sample_names


    def get_df_info(self):
        return self.dict

    def remove_tags_file(self):
        if os.path.exists(self.tagsFile):
            os.remove(self.tagsFile)

        