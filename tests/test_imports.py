#!/usr/bin/env python
import os


def test_imports():
    """
    Test the imports
    """
    for root, dirs, files in os.walk('.'):
        for f in files:
            package = os.path.basename(root)
            module = os.path.splitext(f)[0]
            if f.endswith('.py') and '__' not in f and 'test_' not in f and root != '.' and package != 'get':
                if package != 'cowbat':
                    import_statement = 'import {package}.{module}'.format(package=package,
                                                                          module=module)
                    print(import_statement)
                    exec(import_statement)
                else:
                    with open(os.path.join(root, f), 'r') as python_file:
                        for line in python_file:
                            if line.startswith('from') or line.startswith('import'):
                                print(line.rstrip())
                                exec(line.rstrip())
