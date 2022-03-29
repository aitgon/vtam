import os
import pandas as pd
try:
    from Bio.Alphabet import generic_dna
except ImportError:
    generic_dna = None

from Bio.Seq import Seq


class FilesInputCutadapt(object):
    """ make fasta files to use as input for cutadapt demultiplexing """

    def __init__(self, file_path, mergedFasta, no_reverse, tag_to_end, primer_to_end):

        self.file_path = file_path
        self.mergedFasta = mergedFasta
        
        self.tag_to_end = tag_to_end
        self.primer_to_end = primer_to_end
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

            for i, name in enumerate(sample_names):

                if self.no_reverse:
                    y = i // 2
                else:
                    y = i
                
                if not "_reversed" in name:
                    if generic_dna:  # Biopython <1.78
                        tagRevRC = str(Seq(self.dict['tagrev'][y], generic_dna).reverse_complement())
                    else:  # Biopython =>1.78
                        tagRevRC = str(Seq(self.dict['tagrev'][y]).reverse_complement())
                    
                    if not self.tag_to_end:
                        tags.write(f">{name}\n^{self.dict['tagfwd'][y]}...{tagRevRC}$\n")
                    else:
                        tags.write(f">{name}\n{self.dict['tagfwd'][y]};min_overlap={str(len(self.dict['tagfwd'][y]))}...{tagRevRC};min_overlap={str(len(tagRevRC))}\n")
                    
                else:
                    if generic_dna:  # Biopython <1.78
                        tagFwdRC = str(Seq(self.dict['tagfwd'][y], generic_dna).reverse_complement())
                    else:  # Biopython =>1.78
                        tagFwdRC = str(Seq(self.dict['tagfwd'][y]).reverse_complement())
                
                    if not self.tag_to_end:
                        tags.write(f">{name}\n^{self.dict['tagrev'][y]}...{tagFwdRC}$\n")
                    else:
                        tags.write(f">{name}\n{self.dict['tagrev'][y]};min_overlap={str(len(self.dict['tagrev'][y]))}...{tagFwdRC};min_overlap={str(len(tagFwdRC))}\n")

        return self.tagsFile


    def primers(self):
                
        self.primers_result = []

        for i in range(self.length):
            primerFwd = self.dict["primerfwd"][i] 
            primerRev = self.dict["primerrev"][i] 
            marker = self.dict["marker"][i]
            primer = (marker, primerFwd, primerRev, str(len(primerFwd)), str(len(primerRev)))
            if not  primer in self.primers_result:
                self.primers_result.append(primer)

        return self.primers_result


    def get_sample_names(self):
        
            sample_names = []
            for i in range(self.length):
                name = f'{self.dict["sample"][i]}'
                if name not in sample_names:
                    sample_names.append(name)
                

            if self.no_reverse:
                samples_names_revered = []
                for name in sample_names:
                    samples_names_revered.append(name)
                    samples_names_revered.append(name + "_reversed")
                return samples_names_revered
                
            return sample_names


    def get_df_info(self):
        return self.dict

    def remove_tags_file(self):
        if os.path.exists(self.tagsFile):
            os.remove(self.tagsFile)

        