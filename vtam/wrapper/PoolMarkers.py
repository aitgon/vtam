import inspect
import pandas

from wopmars.models.ToolWrapper import ToolWrapper

from vtam.utils.Logger import Logger
from vtam.utils.OptionManager import OptionManager
from vtam.utils.PoolMarkerRunner import PoolMarkerRunner


class PoolMarkers(ToolWrapper):
    __mapper_args__ = {
        "polymorphic_identity": "vtam.wrapper.PoolMarkers"}

    # Input file
    __input_file_tax_assign = "ASVtable"

    # Output table
    __output_file_pooled_markers = "PooledMarkers"

    def specify_input_file(self):
        return [
            PoolMarkers.__input_file_tax_assign,
        ]

    def specify_input_table(self):
        return [
        ]

    def specify_output_file(self):
        return [
            PoolMarkers.__output_file_pooled_markers,
        ]

    def specify_params(self):
        return {
            # "log_verbosity": "int",
            # "log_file": "str"
        }

    def run(self):
        session = self.session()
        engine = session._WopMarsSession__session.bind
        # OptionManager.instance()['log_verbosity'] = int(self.option("log_verbosity"))
        # if not self.option("log_verbosity") is None:
        #     OptionManager.instance()['log_file'] = str(self.option("log_file"))

        #########################################################
        #
        # 1. Wrapper inputs, outputs and parameters
        #
        #########################################################
        Logger.instance().debug(
            "file: {}; line: {}; Wrapper inputs, outputs and parameters.".format(__file__,
                                                                                 inspect.currentframe().f_lineno, ))
        #
        # Input file
        input_file_asv_table = self.input_file(PoolMarkers.__input_file_tax_assign)

        # Output table models
        output_file_pooled_markers = self.output_file(PoolMarkers.__output_file_pooled_markers)

        input_file_asv_table_df = pandas.read_csv(input_file_asv_table, sep="\t", header=0)
        pool_marker_runner = PoolMarkerRunner(input_file_asv_table_df)
        pooled_marker_df = pool_marker_runner.get_pooled_marker_df()
        pooled_marker_df.to_csv(output_file_pooled_markers, sep="\t", index=False)

