#!/usr/bin/env python

import glob


def test_imports():
    python_files = glob.glob('*/*.py')
    for python_file in python_files:
        if '__init__.py' not in python_file and 'test_' not in python_file:
            import_statement = 'import {}.{}'.format(python_file.split('/')[0], python_file.split('/')[1]
                                                     .replace('.py', ''))
            print(import_statement)
            exec(import_statement)
