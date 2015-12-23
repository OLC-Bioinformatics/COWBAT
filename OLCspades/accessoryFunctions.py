#!/usr/bin/env python
__author__ = 'adamkoziol'


def make_path(inpath):
    """
    from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL
    :param inpath: string of the supplied path
    """
    import os
    import errno
    try:
        # os.makedirs makes parental folders as required
        os.makedirs(inpath)
    # Except os errors
    except OSError as exception:
        # If the os error is anything but directory exists, then raise
        if exception.errno != errno.EEXIST:
            raise


def make_dict():
    """Makes Perl-style dictionaries"""
    from collections import defaultdict
    return defaultdict(make_dict)


def printtime(string, start):
    """Prints a string in bold with the elapsed time
    :param string: a string to be printed in bold
    :param start: integer of the starting time
    """
    import time
    print('\n\033[1m' + "[Elapsed Time: {:.2f} seconds] {}".format(time.time() - start, string) + '\033[0m')


class GenObject(object):
    """Object to store static variables"""
    def __init__(self, x=None):
        start = (lambda y: y if y else {})(x)
        super(GenObject, self).__setattr__('datastore', start)

    def __getattr__(self, key):
        return self.datastore[key]

    def __setattr__(self, key, value):
        if value:
            self.datastore[key] = value
        else:
            self.datastore[key] = "NA"


class MetadataObject(object):
    """Object to store static variables"""
    def __init__(self):
        """Create datastore attr with empty dict"""
        super(MetadataObject, self).__setattr__('datastore', {})

    def __getattr__(self, key):
        """:key is retrieved from datastore if exists, for nested attr recursively :self.__setattr__"""
        if key not in self.datastore:
            self.__setattr__(key)
        return self.datastore[key]

    def __setattr__(self, key, value=GenObject(), **args):
        """Add :value to :key in datastore or create GenObject for nested attr"""
        if args:
            self.datastore[key].value = args
        else:
            self.datastore[key] = value

    def dump(self):
        """Print only the nested dictionary values; removes __methods__ and __members__ attributes"""
        metadata = {}
        for attr in self.datastore:
            metadata[attr] = {}
            if not attr.startswith('__'):
                if isinstance(self.datastore[attr], str):
                    metadata[attr] = self.datastore[attr]
                else:
                    metadata[attr] = self.datastore[attr].datastore
        return metadata
