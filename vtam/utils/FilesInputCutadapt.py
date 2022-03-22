import os
import pandas as pd
try:
    from Bio.Alphabet import generic_dna
except ImportError:
    generic_dna = None

from Bio.Seq import Seq


class FilesInputCutadapt(object):
    """ make fasta files to use as input for cutadapt demultiplexing """

    def __init__(self, file_path, mergedFasta, no_reverse=False, tag_to_end=False, primer_to_end=False):

        self.file_path = file_path
        self.mergedFasta = mergedFasta
        
        self.tag_to_end = tag_to_end
        self.primer_to_end = primer_to_end
        self.no_reverse = no_reverse

        self.df = pd.read_csv(file_path, sep='\t', header=0)
        self.columns = self.df.columns.str.lower()
        self.length = self.df.shape[0]

        self.sample = [self.df.iloc[i]['sample'] for i in range(self.length) if self.df.iloc[i].mergedfasta == self.mergedFasta]

        self.tag_fwd = [self.df.iloc[i].tagfwd for i in range(self.length) if self.df.iloc[i].mergedfasta == self.mergedFasta]
        self.tag_rev = [self.df.iloc[i].tagrev for i in range(self.length) if self.df.iloc[i].mergedfasta == self.mergedFasta]

        self.primer_fwd = [self.df.iloc[i].primerfwd for i in range(self.length) if self.df.iloc[i].mergedfasta == self.mergedFasta]
        self.primer_rev = [self.df.iloc[i].primerrev for i in range(self.length) if self.df.iloc[i].mergedfasta == self.mergedFasta]

        self.mergedfasta_list = [self.df.iloc[i].mergedfasta for i in range(self.length)]

    def tags_file(self):
        
        relPath, _ = os.path.split(self.file_path)
        self.tagsFile = os.path.join(relPath + 'tagsFile.fasta')

        with open(self.tagsFile, 'at') as tags:
            for tagFwd, tagRev, sample  in zip(self.tag_fwd, self.tag_rev, self.sample):

                tagFwd = tagFwd.strip()
                tagRev = tagRev.strip()
                sample = sample.strip()

                print(f'len(self.tag_fwd) : {len(self.tag_fwd)}')

                if generic_dna:  # Biopython <1.78
                    tagRevRC = str(Seq(tagRev, generic_dna).reverse_complement())
                else:  # Biopython =>1.78
                    tagRevRC = str(Seq(tagRev).reverse_complement())
                
                if not self.tag_to_end:
                    print("self.tag_to_end 1", self.tag_to_end)
                    tags.write(f">{sample}\n^{tagFwd}...{tagRevRC}$\n")
                else:
                    print("self.tag_to_end 2", self.tag_to_end)
                    tags.write(f">{sample}\n{tagFwd}...{tagRevRC}\n")
                
                if self.no_reverse:
                    print("self.no_reverse", self.no_reverse)
                    if generic_dna:  # Biopython <1.78
                        tagFwdRC = str(Seq(tagFwd, generic_dna).reverse_complement())
                    else:  # Biopython =>1.78
                        tagFwdRC = str(Seq(tagFwd).reverse_complement())
                
                    if self.tag_to_end:
                        tags.write(f">{sample}_reversed\n ^{tagRev}...{tagFwdRC}$\n")
                    else:
                        tags.write(f">{sample}_reversed\n {tagRev}...{tagFwdRC}\n")

        return self.tagsFile


    def primers_file(self):
        
        relPath, _ = os.path.split(self.file_path)
        self.primersFile = os.path.join(relPath + 'primersFile.fasta')
    
        with open(self.primersFile, 'at') as primers:
            for primerFwd, primerRev, sample in zip(self.primer_fwd, self.primer_rev, self.sample):
                
                primerFwd = primerFwd.strip()
                primerRev = primerRev.strip()
                sample = sample.strip()

                print(f'len(self.primer_fwd) : {len(self.primer_fwd)}')

                if generic_dna:  # Biopython <1.78
                    primerRevRC = str(Seq(primerRev, generic_dna).reverse_complement())
                else:  # Biopython =>1.78
                    primerRevRC = str(Seq(primerRev).reverse_complement())

                if not self.primer_to_end:
                    primers.write(f">{sample}\n^{primerFwd}...{primerRevRC}$\n")
                else:
                    primers.write(f">{sample}\n{primerFwd}...{primerRevRC}\n")
                
                if self.no_reverse:
                    if generic_dna:  # Biopython <1.78
                        primersFwdRC = str(Seq(primerFwd, generic_dna).reverse_complement())
                    else:  # Biopython =>1.78
                        primersFwdRC = str(Seq(primerFwd).reverse_complement())

                    if self.primer_to_end:
                        primers.write(f">{sample}_reversed\n^{primerRev}...{primersFwdRC}$\n")
                    else:
                        primers.write(f">{sample}_reversed \n{primerRev}...{primersFwdRC}\n")
        
        return self.primersFile


    def get_sample_names(self):

        if self.no_reverse:
            samples = []
            for name in self.sample:
                samples.append(name)
                samples.append(name + "_reversed")
            return samples

        return self.sample

    def get_mergedfasta(self):

        return self.mergedfasta_list

    def remove_tags_file(self):
        if os.path.exists(self.tagsFile):
            os.remove(self.tagsFile)


    def remove_primers_file(self):
        if os.path.exists(self.primersFile):
            os.remove(self.primersFile)
        