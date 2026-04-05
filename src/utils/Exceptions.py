

class Error(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return f'{self.message}'


class WrongEnergyError(Error):
    pass


class WrongDistanceError(Error):
    pass
