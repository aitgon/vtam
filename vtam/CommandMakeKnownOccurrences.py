import os
import numpy as np
import pandas as pd
import pathlib

class CommandMakeKnownOccurrences(object):
    """Class for the makeKnownOccurrences command"""

    @staticmethod
    def main(asvTable, sampleTypes, mockComposition, known_occurrences = './known_occurrences.tsv', missing_occurrences = './missing_occurrences.tsv', habitat_proportion = 0.5):

        asv = pd.read_csv(asvTable, sep='\t')
        mock = pd.read_csv(mockComposition, sep='\t')
        samples = pd.read_csv(sampleTypes, sep='\t')

        # create output directories if missing
        pathlib.Path(os.path.dirname(known_occurrences)).mkdir(parents=True, exist_ok=True)
        pathlib.Path(os.path.dirname(missing_occurrences)).mkdir(parents=True, exist_ok=True)
        
        # copy mock_composition file
        occurrences = mock[mock['action']!='tolerate'].replace(np.nan, '', regex=True).astype(str)

        # add mock, variant and action attributes if missing
        if 'mock' not in occurrences:
            occurrences.insert(3,'mock', '1')
        if 'variant' not in occurrences:
            occurrences.insert(4,'variant', '0')
        if 'action' not in occurrences:
            occurrences.insert(5,'action', 'keep')

        # finding missing occurrences and adding variant number for expected sequences of mock samples
        missing = pd.DataFrame(columns=['run','marker','sample','sequence','tax_name'])
        for i, moc in occurrences.iterrows():
            row = asv[(asv['sequence'] == moc['sequence']) & (asv['marker'] == moc['marker']) & (asv['run'] == moc['run'])][['marker','run','variant',moc['sample'],'sequence']].to_dict('records')
            if row:
                occurrences.at[i,'variant'] = row[0]['variant']
            else:
                occurrences.at[i,'variant'] = 0
            if not row or row[0][moc['sample']] == 0:
                dfToConcat = pd.DataFrame({'run':moc['run'],'marker':moc['marker'],'sample':moc['sample'],'sequence':moc['sequence'],'tax_name':moc['tax_name']}, index=[0])
                missing = pd.concat([missing, dfToConcat] ,ignore_index=True)
                #replaced append by concat; append method deprecated
                #missing = missing.append({'run':moc['run'],'marker':moc['marker'],'sample':moc['sample'],'sequence':moc['sequence'],'tax_name':moc['tax_name']},ignore_index=True)
        
        # inner function to avoid repetitions
        def addDelete(asvRow, sample, isMock):
            newOccurrence = {}
            nonlocal occurrences
            for key in ['run','marker','variant','sequence']:
                newOccurrence[key] = asvRow[key]
            newOccurrence['sample'] = sample
            newOccurrence['action'] = 'delete'
            if  'ltg_tax_name' in row:
                newOccurrence['tax_name'] = asvRow['ltg_tax_name']
            if isMock:
                newOccurrence['mock'] = '1'
            else :
                newOccurrence['mock'] = '0'
            occurrences = pd.concat([occurrences, pd.DataFrame(newOccurrence, index=[0])], ignore_index=True)
            
        # delete occurrences of unexpected sequences from mock samples and from negative samples
        for _,samp in samples[samples['sample_type'].isin(['mock','negative'])].iterrows():
            for _,row in asv[(asv[samp['sample']]!=0) & (asv['run'] == samp['run'])].iterrows():
                if samp['sample_type'] == 'negative':
                    addDelete(row,samp['sample'],False)
                elif row['sequence'] not in list(mock[(mock['sample'] == samp['sample']) & (mock['run']==samp['run'])]['sequence']):
                    addDelete(row,samp['sample'],True)
                
        # delete occurrences where N_ih/N_i < habitat_proportion
        listHabitats = samples['habitat'].dropna().unique().tolist()
        for  _, row in asv.iterrows():
            habitats = dict.fromkeys(listHabitats,0)
            total = 0
            for _, samp in samples.dropna().iterrows():
                habitats[samp['habitat']] += row[samp['sample']]
                total += row[samp['sample']]
            for h, count in habitats.items():
                if count > 0 and count/total < float(habitat_proportion):
                    for samp in list(samples[(samples['habitat'] == h) & (samples['sample_type'] == 'real')]['sample']):
                        if row[samp] != 0 :
                            addDelete(row,samp,False)
                            
        occurrences.to_csv(known_occurrences, sep="\t", header=True, index=False)

        if not missing.empty:
            missing.to_csv(missing_occurrences, sep="\t", header=True, index=False)
