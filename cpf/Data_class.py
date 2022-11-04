# Inherit these through detector classes
class CpfData:
    def __init__(self, data):
        """
        :param data:
        """
        self.x = data[0, :, :]
        self.y = data[0]
        self.I = data[0]
        return self.x

    def twoth(self):
        return convert()[:]

    def mask(self):
        # load mask file
        return masked_data
