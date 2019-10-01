class VTAMexception(Exception):
    """
    VTAM expception
    """

    def __init__(self, message=""):
        """
        :param message: string that contains the precise error
        :return: void
        """
        self.__message = message

    def __str__(self):
        s = "{}".format(self.__message)
        return s
