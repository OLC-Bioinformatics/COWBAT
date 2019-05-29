#!/usr/bin/env python
import os


def test_imports():
    """
    Test the imports
    """
    for root, dirs, files in os.walk('.'):
        for f in files:
            package = '.'.join(os.path.split(root.lstrip('./')))
            module = os.path.splitext(f)[0]
            if f.endswith('.py') and '__' not in f and '__' not in root and 'test_' not in f and root != '.' \
                    and f != 'setup.py':
                if 'cowbat' not in package:
                    import_statement = 'import {package}.{module}'.format(package=package,
                                                                          module=module)
                    print(import_statement.rstrip())
                    exec(import_statement.rstrip())
                # Try all the imports from within the file
                with open(os.path.join(root, f), 'r') as python_file:
                    for line in python_file:
                        if line.startswith('from') or line.startswith('import'):
                            statement = line.rstrip()
                            if len(statement.split(',')) == 1:
                                print(line.rstrip())
                                exec(line.rstrip())
                            else:
                                # from accessoryFunctions.accessoryFunctions import make_path, MetadataObject
                                if line.startswith('from'):
                                    if line.rstrip().endswith('\\'):
                                        line = line.rstrip().replace('\\', '') + next(python_file).rstrip().lstrip()
                                    import_base, import_line = line.split(' import ')
                                    imported_packages = import_line.split(',')
                                    for imported_package in imported_packages:
                                        # Don't check relative imports
                                        if '.' not in import_base:
                                            import_string = '{import_base} import {import_package}' \
                                                .format(import_base=import_base,
                                                        import_package=imported_package)
                                            print(import_string)
                                            exec(import_string)
                                # import sys, os
                                else:
                                    imported_packages = line.lstrip('import').split(', ')
                                    for imported_package in imported_packages:
                                        import_string = 'import {imp_package}' \
                                            .format(imp_package=imported_package.lstrip().rstrip())
                                        print(import_string)
                                        exec(import_string)
