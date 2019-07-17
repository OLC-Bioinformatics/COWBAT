#!/usr/bin/env python
from setuptools import setup, find_packages
import os
__author__ = 'adamkoziol'
setup(
    name="COWBAT",
    version="0.5.0.7",
    include_package_data=True,
    packages=find_packages(),
    scripts=[os.path.join('cowbat', 'assembly_pipeline.py'),
             os.path.join('cowbat', 'assembly_typing.py'),
             os.path.join('cowbat', 'validation', 'validate_cowbat.py')
             ],
    license='MIT',
    author='Adam Koziol',
    author_email='adam.koziol@canada.ca',
    description='CFIA OLC Workflow for Bacterial Assembly and Typing',
    url='https://github.com/OLC-Bioinformatics/COWBAT',
    long_description=open('README.md').read(),
    install_requires=['interop']
)
