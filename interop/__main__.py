import interop
import unittest
import argparse
from interop import py_interop_run,py_interop_metrics,py_interop_run_metrics,py_interop_comm,py_interop_table,py_interop_plot,py_interop_summary
from interop.CoreTests import CoreTests

def execute_from_commandline():
    """ Provide a command line interface to a general test script
    """

    parser = argparse.ArgumentParser(prog='interop_check', description='Test script for the InterOp Library Python Interface')
    parser.add_argument('--version', action='version', version=interop.__version__)
    parser.add_argument('--test', action='store_true', help='Run the unit tests')
    parser.set_defaults(test=False)
    param = parser.parse_args()
    if param.test:
        testsuite = unittest.makeSuite(CoreTests)
        unittest.TextTestRunner(verbosity=1).run(testsuite)

if __name__ == "__main__":
    execute_from_commandline();

