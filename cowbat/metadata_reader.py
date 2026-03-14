#!/usr/bin/env python3

"""
This script reads metadata from JSON files, processes it, and updates the
metadata objects using the CustomBox class.
"""

# Standard library imports
import json
import os
from typing import List

# Third-party imports
from olctools.accessoryFunctions.metadata import CustomBox

__author__ = 'adamkoziol'


def read_metadata(*, sample: CustomBox) -> List[CustomBox]:
    """
    Reads metadata from JSON files, processes it, and updates the metadata
    objects.

    :param sequence_path: Path containing sequencing run files
    :return: A list of CustomBox objects with updated metadata.
    """
    # Check if the metadata file exists and is not empty
    if os.path.isfile(sample.json_file) and os.stat(sample.json_file).st_size:
        try:
            with open(sample.json_file, encoding='utf-8') as metadata_report:
                json_data = json.load(metadata_report)

            # Create the metadata object
            metadata_obj = CustomBox()

            # Initialize the metadata categories as CustomBox objects
            for attr, value in json_data.items():
                if isinstance(value, dict):
                    metadata_obj[attr] = CustomBox(value)
                else:
                    metadata_obj[attr] = value

            return metadata_obj
        except ValueError:
            return sample
    else:
        return None
