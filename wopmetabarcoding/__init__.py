#!/usr/bin/env python
import os

import jinja2

from wopmetabarcoding.utils.ArgParser import ArgParser
from wopmetabarcoding.utils.OptionManager import OptionManager

from wopmetabarcoding.utils.OptionManager import OptionManager
from wopmetabarcoding.utils.constants import parameters_numerical
from wopmetabarcoding.utils.utilities import tempdir


class VTAM(object):

    def __init__(self, sys_argv):
        #
        ##################################
        #
        # Load default numerical parameters
        #
        ##################################
        parameters = parameters_numerical.copy()
        #
        ##################################
        #
        # Read command-line interface arguments
        # Store arguments in OptionsManager
        #
        ##################################
        self.sys_argv = sys_argv
        parser = ArgParser.get_arg_parser()
        args = parser.parse_args(sys_argv)
        for k in vars(args):
            parameters[k] = vars(args)[k]
        #
        #############################################################
        #
        # Render Wopfile
        #
        #############################################################
        wopfile_file_name = 'Wopfile_merge.yml'
        wopfile_in_path = os.path.join(os.path.dirname(__file__), '../data', wopfile_file_name)
        wopfile_out_path = os.path.join(tempdir, wopfile_file_name)
        #
        with open(wopfile_in_path) as fin:
            template = jinja2.Template(fin.read())
        wopfile_rendered = template.render(parameters)
        with open(wopfile_out_path, "w") as fout:
            fout.write(wopfile_rendered)
        #
        #############################################################
        #
        # Run Wopfile
        #
        #############################################################
        # {'fastq_minovlen': 50, 'fastq_maxmergelen': 300, 'fastq_minmergelen': 100, 'fastq_minlen': 50, 'fastq_maxee': 1,
        #  'fastq_truncqual': 10, 'fastq_maxns': 0, 'threads': 8, 'fastq_ascii': 33, 'min_id': 0.8, 'minseqlength': 32,
        #  'overhang': 0, 'lfn_variant_threshold': 0.001, 'lfn_variant_replicate_threshold': 0.001,
        #  'lfn_biosample_replicate_threshold': 0.001, 'lfn_read_count_threshold': 10, 'min_replicate_number': 2,
        #  'pcr_error_var_prop': 0.1, 'renkonen_threshold': 0.5, 'genetic_table_number': 5, 'identity_threshold': 97,
        #  'include_prop': 90, 'min_number_of_taxa': 3,
        #  'db': '/home/gonzalez/Software/repositories/wopmetabarcodin/db.sqlite', 'dryrun': True, 'forceall': False,
        #  'log': None, 'params': None, 'targetrule': None, 'verbose': None,
        #  'fastqinfo': '/home/gonzalez/Software/repositories/wopmetabarcodin/test/test_parser.py',
        #  'fastainfo': '/home/gonzalez/Software/repositories/wopmetabarcodin/foo',
        #  'fastqdir': '/home/gonzalez/Software/repositories/wopmetabarcodin/test',
        #  'fastadir': '/home/gonzalez/Software/repositories/wopmetabarcodin/foo'}
        parameters['wopfile_out_path'] = wopfile_out_path
        cmd_wopmars = "wopmars -w {wopfile_out_path} --db {db} --fastqinfo {fastqinfo} --fastqdir {fastqdir} " \
                      "--fastainfo {fastainfo} --fastadir {fastqdir}".format(**parameters)
        import pdb; pdb.set_trace()
        # if not parameters['fastqinfo'] is None:
        #     cmd_wopmars += " --fastqinfo" + parameters['fastqinfo']
        # if not parameters['fastqdir'] is None:
        #     cmd_wopmars += " --fastqdir" + parameters['fastqdir']
        # if not parameters['fastainfo'] is None:
        #     cmd_wopmars += " --fastainfo" + parameters['fastainfo']
        # if not parameters['fastadir'] is None:
        #     cmd_wopmars += " --fastadir" + parameters['fastadir']
        #
        # wopmars options
        if parameters['dryrun']:
            cmd_wopmars += " -n"
        # cmd = "wopmars -w {wopfile_out_path} -D sqlite:///{db} -v -p".format(**args_dic)
        # logger.debug(
        #     "file: {}; line: {}; CMD: {}".format(__file__, inspect.currentframe().f_lineno, cmd))
        # if args_dic['dryrun']:
        #     cmd = cmd + " -n"
        # if args_dic['forceall']:
        #     cmd = cmd + " -F"
        # if 'targetrule' in args_dic:
        #     cmd = cmd + " -t {targetrule}".format(**args_dic)
        # if not args_dic['log'] is None:
        #     PathFinder.mkdir_p(os.path.dirname(args_dic['log']))
        #     Path(args_dic['log']).touch() # touch log
        #     cmd = cmd + " --log " + args_dic['log']
        # os.system(cmd)
        # sys.exit(0)


if __name__ == '__main__':
    VTAM()
