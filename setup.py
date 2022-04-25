#!/usr/bin/env python
from setuptools import setup, find_packages
from distutils.util import convert_path
import os
__author__ = 'adamkoziol'

# Find the version
version = dict()
with open(convert_path(os.path.join('cowbat', 'version.py')), 'r') as version_file:
    exec(version_file.read(), version)

setup(
    name="COWBAT",
    version=version['__version__'],
    include_package_data=True,
    packages=find_packages(),
    scripts=[os.path.join('cowbat', 'assembly_pipeline.py'),
             os.path.join('cowbat', 'assembly_typing.py'),
             os.path.join('cowbat', 'validation', 'validate_cowbat.py')
             ],
    license='MIT',
    author='Adam Koziol',
    author_email='adam.koziol@inspection.gc.ca',
    description='CFIA OLC Workflow for Bacterial Assembly and Typing',
    url='https://github.com/OLC-Bioinformatics/COWBAT',
    long_description=open('README.md').read(),
)
