import os
import pandas
import sqlalchemy
import subprocess
import sys

from vtam.utils.Logger import Logger
from vtam.utils.RunnerWopmars import RunnerWopmars
from vtam.utils.constants import FilterLFNreference_records
from vtam.utils.FileDecompression import FileDecompression
from vtam.utils.FileSampleInformation import FileSampleInformation

class CommandFilterOptimize(object):
    """Class for the Merge command"""

    @staticmethod
    def main(arg_parser_dic):

        ###################################################################
        #
        # Check if the given files are compressed and decompress them
        #
        ###################################################################
        
        path_to_folder = arg_parser_dic['sorteddir']
        sortedInfo_df = FileSampleInformation(arg_parser_dic['sortedinfo']).read_tsv_into_df()
        decompressed_files = []
        sortedInfoDecompressed_df = sortedInfo_df.copy()
        for sortedFasta in sortedInfo_df[['sortedfasta']].drop_duplicates().values:
            sortedFasta = sortedFasta[0]
            if sortedFasta.endswith(".gz"):
                path_to_gz_file = os.path.join(path_to_folder, sortedFasta)
                decompress = FileDecompression(path_to_gz_file)
                fileDecompressed = decompress.gzip_decompression()
                _, relPath = os.path.split(fileDecompressed)
                sortedInfoDecompressed_df.loc[sortedInfoDecompressed_df['sortedfasta'] == sortedFasta, 'sortedfasta'] = relPath
                ##keep the decompressed files to erase them later
                decompressed_files.append(relPath)
        if decompressed_files: 
            sortedInfoDecompressed_df.to_csv(arg_parser_dic['sortedinfo'], sep="\t", header=True, index=False)


        ###################################################################
        #
        # Create FilterLFNreference table and fill it
        #
        ###################################################################

        engine = sqlalchemy.create_engine('sqlite:///{}'.format(str(arg_parser_dic['db'])), echo=False)
        meta = sqlalchemy.MetaData()
        filter_lfn_reference = sqlalchemy.Table(
            'FilterLFNreference', meta,
            sqlalchemy.Column('filter_id', sqlalchemy.Integer, primary_key=True),
            sqlalchemy.Column('filter_name', sqlalchemy.String),
        )
        meta.create_all(engine)

        with engine.connect() as conn:
            for filter_rec in FilterLFNreference_records:
                filter_name = filter_rec['filter_name']
                select_row = conn.execute(
                    sqlalchemy.select([filter_lfn_reference.c.filter_id]).where(
                        filter_lfn_reference.c.filter_name == filter_name)).first()
                if select_row is None:  # variant_sequence IS NOT in the database, so INSERT it
                    conn.execute(
                        filter_lfn_reference.insert().values(
                            **filter_rec))

        wopmars_runner = RunnerWopmars(command=arg_parser_dic['command'], cli_args_dic=arg_parser_dic)
        wopmars_command = wopmars_runner.get_wopmars_command()

        ########################################################################################
        #
        # Run wopmars
        #
        ########################################################################################

        # Some arguments will be passed through environmental variables
        if 'threads' in arg_parser_dic:
            os.environ['VTAM_THREADS'] = str(arg_parser_dic['threads'])
        Logger.instance().info(wopmars_command)
        run_result = subprocess.run(wopmars_command, shell=True)

        ########################################################################################
        #
        # Delete decompressed files and restore original infofile
        #
        ########################################################################################

        if decompressed_files:
            #delete the decompressed files
            for file in decompressed_files:
                os.remove(os.path.join(arg_parser_dic['sorteddir'],file))
            #delete the csv info file with decompressed names and replace by the original one
            sortedInfo_df.to_csv(arg_parser_dic['sortedinfo'], sep="\t", header=True, index=False)


        sys.exit(run_result.returncode)
