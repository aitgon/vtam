import os
import sqlalchemy
import subprocess
import sys

from vtam.utils.Logger import Logger
from vtam.utils.RunnerWopmars import RunnerWopmars
from vtam.utils.constants import FilterLFNreference_records


class CommandFilterOptimize(object):
    """Class for the Merge command"""

    @staticmethod
    def main(arg_parser_dic):

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
        sys.exit(run_result.returncode)
