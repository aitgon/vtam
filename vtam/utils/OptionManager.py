import yaml

from vtam.utils.Singleton import Singleton

class OptionManager(dict, Singleton):
    """
    The OptionManager contains command-line arguments and VTAM parameters

    OptionManager inherits from dict so it behaves exactly the same but takes only one possible argument.
    It is used like this: OptionManager.instance()['--log'] = "vtam.log"
    """

    def __init__(self, *args, **kwargs):
        """
        :param: dict_args: dictionnary containing the arguments passed to the command line.
        :return: void
        """
        super().__init__(*args, **kwargs)

    def add_options(self, option_dic):
        """Takes a dictionnary with option:values and add it to the OptionManager"""
        for option_key in option_dic:
            if not option_key == 'params':
                option_value = option_dic[option_key]
                self[option_key] = option_value
            else:
                params_yml_path = option_dic[option_key]
                if not params_yml_path is None:
                    with open(params_yml_path, 'r') as fin_params_yml:
                        params_dic = yaml.load(fin_params_yml, Loader=yaml.SafeLoader)
                        for params_key in params_dic:
                            self[params_key] = params_dic[params_key]
