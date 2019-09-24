from wopmetabarcoding.utils.Singleton import Singleton

class OptionManager(dict, Singleton):
    """
    The OptionManager contains command-line arguments and VTAM parameters_numerical

    OptionManager inherits from dict so it behaves exactly the same but takes only one possible argument.
    It is used like this: OptionManager.instance()['--log'] = "vtam.log"
    """

    def __init__(self, *args, **kwargs):
        """
        :param: dict_args: dictionnary containing the arguments passed to the command line.
        :return: void
        """
        super().__init__(*args, **kwargs)

