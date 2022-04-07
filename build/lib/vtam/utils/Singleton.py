# Based on wopmars genetic_code
class Singleton(object):
    __instance = None

    @classmethod
    def instance(cls):
        if not cls.__instance:
            cls.__instance = cls()
        return cls.__instance

    @classmethod
    def _drop(cls):
        """Drop the instance (for testing purposes)."""
        if cls.__instance:
            if cls.__instance:
                cls.__instance = None
