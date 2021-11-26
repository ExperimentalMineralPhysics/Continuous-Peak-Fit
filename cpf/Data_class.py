



# Inherit these through detector classes
class CpfData:

    def __init__(self, data):
        """
        :param data:
        """
        self.x = data[0,:,:]
        self.y = data[0]
        self.I = data[0]

    def self.twoth:
        return self.convert()[:]

    def self.mask:
        # load mask file
        return masked_data