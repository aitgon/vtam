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
        self.columns = self.df.columns.str.lower()
        self.length = self.df.shape[0]

        self.run = [self.df.iloc[i]['run'] for i in range(self.length) if self.df.iloc[i].mergedfasta == self.mergedFasta]
        self.marker = [self.df.iloc[i]['marker'] for i in range(self.length) if self.df.iloc[i].mergedfasta == self.mergedFasta]
        self.sample = [self.df.iloc[i]['sample'] for i in range(self.length) if self.df.iloc[i].mergedfasta == self.mergedFasta]
        self.replicate = [self.df.iloc[i]['replicate'] for i in range(self.length) if self.df.iloc[i].mergedfasta == self.mergedFasta]

        self.sampleInfo = [(run.strip(), marker.strip(), sample.strip(), str(replicate)) for run, marker, sample, replicate 
            in zip(self.run, self.marker,  self.sample, self.replicate)]
        self.tag_fwd = [self.df.iloc[i].tagfwd for i in range(self.length) if self.df.iloc[i].mergedfasta == self.mergedFasta]
        self.tag_rev = [self.df.iloc[i].tagrev for i in range(self.length) if self.df.iloc[i].mergedfasta == self.mergedFasta]

        self.primer_fwd = [self.df.iloc[i].primerfwd for i in range(self.length) if self.df.iloc[i].mergedfasta == self.mergedFasta]
        self.primer_rev = [self.df.iloc[i].primerrev for i in range(self.length) if self.df.iloc[i].mergedfasta == self.mergedFasta]

        self.mergedfasta_list = [self.df.iloc[i].mergedfasta for i in range(self.length)]

    def tags_file(self):
        
        relPath, _ = os.path.split(self.file_path)
        self.tagsFile = os.path.join(relPath, 'tagsFile.fasta')

        with open(self.tagsFile, 'at') as tags:
            for sample_info, tagFwd, tagRev  in zip(self.sampleInfo, self.tag_fwd, self.tag_rev):
                
                if generic_dna:  # Biopython <1.78
                    tagRevRC = str(Seq(tagRev, generic_dna).reverse_complement())
                else:  # Biopython =>1.78
                    tagRevRC = str(Seq(tagRev).reverse_complement())
                
                if not self.tag_to_end:
                    tags.write(f">{sample_info}\n^{tagFwd}...{tagRevRC}$\n")
                else:
                    tags.write(f">{sample_info}\n{tagFwd};min_overlap={str(len(tagFwd))}...{tagRevRC};min_overlap={str(len(tagRevRC))}\n")
                
                if self.no_reverse:
                    if generic_dna:  # Biopython <1.78
                        tagFwdRC = str(Seq(tagFwd, generic_dna).reverse_complement())
                    else:  # Biopython =>1.78
                        tagFwdRC = str(Seq(tagFwd).reverse_complement())
                
                    if not self.tag_to_end:
                        tags.write(f">{sample_info}_reversed\n^{tagRev}...{tagFwdRC}$\n")
                    else:
                        tags.write(f">{sample_info}_reversed\n{tagRev};min_overlap={str(len(tagRev))}...{tagFwdRC};min_overlap={str(len(tagFwdRC))}\n")

        return self.tagsFile


    def primers(self):
        
        relPath, _ = os.path.split(self.file_path)
        self.primersFile = os.path.join(relPath + 'primersFile.fasta')
    
        # with open(self.primersFile, 'at') as primers:
        self.primers_result = []
        for primerFwd, primerRev, sample in zip(self.primer_fwd, self.primer_rev, self.sample):
            
            primerFwd = primerFwd.strip()
            primerRev = primerRev.strip()
            sample = sample.strip()

            if not (primerFwd, primerRev, str(len(primerFwd)), str(len(primerRev))) in self.primers_result:
                self.primers_result.append((primerFwd, primerRev, str(len(primerFwd)), str(len(primerRev))))
                
        return self.primers_result


    def get_sample_names(self):
        sample_names = [f"{info[0]}_{info[1]}_{info[2]}_{info[3]}" for info 
            in self.sampleInfo]

        if self.no_reverse:
            samples_names_revered = []
            for name in sample_names:
                samples_names_revered.append(name)
                samples_names_revered.append(name + "_reversed")
            return samples_names_revered
            
        return sample_names

    def get_sample_info(self):
        return self.sampleInfo

    def get_mergedfasta(self):
        return self.mergedfasta_list

    def get_df_info(self):
        return 

    def remove_tags_file(self):
        if os.path.exists(self.tagsFile):
            os.remove(self.tagsFile)

        