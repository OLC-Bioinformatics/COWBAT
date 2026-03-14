#!/usr/bin/env python3

"""
Setup script for the COWBAT package.
"""

# Standard imports
import os
import re

# Third-party imports
from setuptools import setup, find_packages

__author__ = 'adamkoziol'

# Find the version
version = {}

# Set the path of the version file
version_script = os.path.join('cowbat', 'version.py')

with open(version_script, 'r', encoding='utf-8') as version_file:
    version_content = version_file.read()
    version_match = re.search(
        r"^__version__ = ['\"]([^'\"]*)['\"]",
        version_content,
        re.M
    )
    if version_match:
        version['__version__'] = version_match.group(1)
    else:
        raise RuntimeError(
            f"Unable to find version string in {version_script}.")

setup(
    name="COWBAT",
    version=version['__version__'],
    include_package_data=True,
    packages=find_packages(),
    scripts=[
        os.path.join('cowbat', 'assembly_pipeline.py'),
        os.path.join('cowbat', 'assembly_typing.py'),
        os.path.join('cowbat', 'validation', 'validate_cowbat.py'),
        os.path.join('cowbat', 'teacup.py'),
    ],
    license='MIT',
    author='Adam Koziol',
    author_email='adam.koziol@inspection.gc.ca',
    description='CFIA OLC Workflow for Bacterial Assembly and Typing',
    url='https://github.com/OLC-Bioinformatics/COWBAT',
    long_description=open('README.md', encoding='utf-8').read(),
)
