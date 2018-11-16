#!/usr/bin/env python
from setuptools import setup, find_packages
__author__ = 'adamkoziol'
setup(
    name="COWBAT",
    version="0.4.4.1",
    include_package_data=True,
    packages=find_packages(),
    scripts=['cowbat/assembly_pipeline.py', 'cowbat/assembly_typing.py', 'cowbat/validation/validate_data.py'],
    license='MIT',
    author='Adam Koziol',
    author_email='adam.koziol@canada.ca',
    description='CFIA OLC Workflow for Bacterial Assembly and Typing',
    url='https://github.com/OLC-Bioinformatics/COWBAT',
    long_description=open('README.md').read()
)
