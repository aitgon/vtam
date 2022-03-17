import os
import pandas as pd
try:
    from Bio.Alphabet import generic_dna
except ImportError:
    generic_dna = None

from Bio.Seq import Seq



class FilesInputCutadapt(object):
    """ make fasta files to use as input for cutadapt demultiplexing """

    def __init__(self, file_path, tag_to_end=False, primer_to_end=False):

        self.file_path = file_path
        
        self.tag_to_end = tag_to_end
        self.primer_to_end = primer_to_end

        self.df = pd.read_csv(file_path, sep='\t', header=0)
        self.columns = self.df.columns.str.lower()
        self.length = self.df.shape[0]

        self.sampleName = [self.df.iloc[i].sampleName for i in range(self.length)]

        self.tag_fwd = [self.df.iloc[i].tagfwd for i in range(self.length)]
        self.tag_rev = [self.df.iloc[i].tagrev for i in range(self.length)]

        self.primer_fwd = [self.df.iloc[i].primerfwd for i in range(self.length)]
        self.primer_rev = [self.df.iloc[i].primerrev for i in range(self.length)]

            
    def tags_file(self):
        
        relPath, _ = os.path.split(self.file_path)
        self.tagsFile = os.path.join(relPath + 'tagsFile.fasta')

        with open(self.tagsFile, 'at') as tags:
            for tagFwd, tagRev, sampleName  in zip(self.tag_fwd, self.tag_rev, self.sampleName):
                if generic_dna:  # Biopython <1.78
                    tagRevRC = str(Seq(tagRev, generic_dna).reverse_complement())
                else:  # Biopython =>1.78
                    tagRevRC = str(Seq(tagRev).reverse_complement())
                
                if self.tag_to_end:
                    tags.write(f">{sampleName}\n ^{tagFwd}...{tagRevRC}$\n")
                else:
                    tags.write(f">{sampleName}\n {tagFwd}...{tagRevRC}\n")
    
        self.tagsFile


    def primers_file(self):
        
        relPath, _ = os.path.split(self.file_path)
        self.primersFile = os.path.join(relPath + 'primersFile.fasta')
    
        with open('./primers.fasta', 'at') as primers:
            for primerFwd, primerRev, sampleName in zip(self.primer_fwd, self.primer_rev, self.sampleName):
                if generic_dna:  # Biopython <1.78
                    primerRevRC = str(Seq(primerRev, generic_dna).reverse_complement())
                else:  # Biopython =>1.78
                    primerRevRC = str(Seq(primerRev).reverse_complement())

                if self.primer_to_end:
                    primers.write(f">{sampleName}\n ^{primerFwd}...{primerRevRC}$\n")
                else:
                    primers.write(f">{sampleName}\n{primerFwd}...{primerRevRC}\n")
        
        self.primersFile


    def remove_tags_file(self):
        if os.path.exists(self.tagsFile):
            os.remove(self.tagsFile)


    def remove_primers_file(self):
        if os.path.exists(self.primersFile):
            os.remove(self.primersFile)
        