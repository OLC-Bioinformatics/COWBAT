#!/usr/bin/env python
from setuptools import setup, find_packages
__author__ = 'adamkoziol'
setup(
    name="COWBAT",
    version="0.2.0",
    include_package_data=True,
    license='MIT',
    author='Adam Koziol',
    author_email='adam.koziol@inspection.gc.ca',
    description='CFIA OLC Workflow for Bacterial Assembly and Typing',
    url='https://github.com/OLC-Bioinformatics/COWBAT',
    long_description=open('README.md').read()
)
