#! /usr/env/python3

"""
Read in a sample sheet, determine if it is from a MiSeq or NextSeq, and
extract the necessary information
"""

# Standard imports
import copy
from glob import glob
import gzip
import logging
import os
from typing import Any, Dict, List, Optional

# Third-party imports
from genemethods.assemblypipeline.metadata_reader import read_metadata
from olctools.accessoryFunctions.accessoryFunctions import (
    filer,
    make_path,
    relative_symlink
)
from olctools.accessoryFunctions.metadata import CustomBox

# Constants
NUM_LINES_TO_READ = 1000


def determine_and_parse_sample_sheet(file_path: str) -> Dict[str, Any]:
    """
    Determines whether the sample sheet is for MiSeq or NextSeq based on its
    content and calls the appropriate parsing function.

    :param file_path: Path to the sample sheet file.
    :return: A dictionary containing parsed data from the sample sheet.
    """
    logging.info(
        "Determining the type of sample sheet for file: %s",
        file_path)

    # Open the sample sheet file and read the first few lines
    with open(file_path, "r", encoding='utf-8') as file:
        lines = file.readlines()

    # Check for specific keywords to determine the sample sheet type
    for line in lines:
        if "Local Run Manager Analysis Id" in line or "iemfileversion" in line:
            logging.info("MiSeq sample sheet detected for file: %s", file_path)
            return parse_miseq_sample_sheet(file_path=file_path)
        if "FileFormatVersion" in line or "InstrumentPlatform" in line:
            logging.info(
                "NextSeq sample sheet detected for file: %s",
                file_path)
            return parse_nextseq_sample_sheet(file_path=file_path)

    # If no specific keywords are found, raise an error
    logging.error(
        "Unable to determine the sample sheet type for file: %s",
        file_path)
    raise ValueError("Unable to determine the sample sheet type.")


def parse_miseq_sample_sheet(file_path: str) -> Dict[str, Any]:
    """
    Parses a MiSeq sample sheet and organizes the data into sections.

    :param file_path: Path to the MiSeq sample sheet file.
    :return: A dictionary containing parsed data from the sample sheet.
    """
    logging.info("Parsing MiSeq sample sheet: %s", file_path)

    # Open the sample sheet file and read all lines
    with open(file_path, "r", encoding='utf-8') as file:
        lines = file.readlines()

    # Initialize a dictionary to store parsed data
    data = {
        "Header": {},  # Stores key-value pairs from the [Header] section
        "ReadsList": [],   # Stores read lengths from the [Reads] section
        "Settings": {},  # Stores key-value pairs from the [Settings] section
        "Data": []     # Stores sample data from the [Data] section
    }

    current_section = None  # Tracks the current section being parsed
    headers = []  # Initialize headers

    # Iterate through each line in the file
    for line in lines:
        line = line.strip()  # Remove leading/trailing whitespace
        if not line:
            continue  # Skip empty lines

        # Check if the line indicates a new section
        if line.startswith("[") and line.endswith("]"):
            current_section = line[1:-1]  # Update the current section
            logging.debug("Entering section: %s", current_section)
            continue

        # Parse lines based on the current section
        if current_section == "Header" or current_section == "Settings":
            # Split the line into key-value pairs and store them
            key, value = line.split(',', 1)
            data[current_section][key.strip()] = value.strip()
            logging.debug("Parsed %s: %s", key.strip(), value.strip())
        elif current_section == "Reads":
            # Append read lengths to the Reads list
            data["ReadsList"].append(int(line))
            logging.debug("Parsed read length: %d", int(line))
        elif current_section == "Data":
            if not headers:
                # The first line in the Data section contains headers
                headers = line.split(',')
                logging.debug("Parsed headers: %s", headers)
            else:
                # Subsequent lines contain sample data
                values = line.split(',')

                sample_data = {
                    headers[i]: values[i] for i in range(len(headers))
                }
                data["Data"].append(sample_data)
                logging.debug("Parsed sample data: %s", sample_data)

    # Convert Reads list to a dictionary with forward and reverse read lengths
    if len(data["ReadsList"]) == 2:
        data["Reads"] = {
            "forward_read_length": data["ReadsList"][0],
            "reverse_read_length": data["ReadsList"][1]
        }
        logging.info(
            "Parsed read lengths: forward=%d, reverse=%d",
            data["Reads"]["forward_read_length"],
            data["Reads"]["reverse_read_length"]
        )

    logging.info("Finished parsing MiSeq sample sheet: %s", file_path)
    return data


def parse_nextseq_sample_sheet(file_path: str) -> Dict[str, Any]:
    """
    Parses a NextSeq sample sheet and organizes the data into sections.

    :param file_path: Path to the NextSeq sample sheet file.
    :return: A dictionary containing parsed data from the sample sheet.
    """
    logging.info("Parsing NextSeq sample sheet: %s", file_path)

    # Open the sample sheet file and read all lines
    with open(file_path, "r", encoding='utf-8') as file:
        lines = file.readlines()

    # Initialize a dictionary to store parsed data
    data = {
        "Header": {},  # Stores key-value pairs from the [Header] section
        "Reads": {},  # Stores read cycle counts from the [Reads] section
        "Sequencing_Settings": {},  # Stores key-value pairs from the
                                    # [Sequencing_Settings] section
        "BCLConvert_Settings": {},  # Stores key-value pairs from the
                                    # [BCLConvert_Settings] section
        "BCLConvert_Data": [],  # Stores sample data from the
                                # [BCLConvert_Data] section
        "Cloud_Settings": {},  # Stores key-value pairs from the
                               # [Cloud_Settings] section
        "Cloud_Data": []  # Stores sample data from the [Cloud_Data] section
    }

    current_section = None  # Tracks the current section being parsed

    # Initialize headers for data sections
    bcl_data_headers = []
    cloud_data_headers = []

    # Iterate through each line in the file
    for line in lines:

        line = line.strip()  # Remove leading/trailing whitespace
        if not line:
            continue  # Skip empty lines

        # Check if the line indicates a new section
        if line.startswith("[") and (
                line.endswith("]") or line.endswith("],")):

            # Update the current section. Remove square brackets
            current_section = line.replace(',', '').replace('[', '') \
                .replace(']', '')
            logging.debug("Entering section: %s", current_section)
            continue

        # Create a variable to store all the sections of the sample sheet
        sections = [
            "Header",
            "Sequencing_Settings",
            "BCLConvert_Settings",
            "Cloud_Settings"
        ]

        # Parse lines based on the current section
        if current_section in sections:

            # Split the line into key-value pairs and store them
            key, value = line.split(',', 1)
            data[current_section][key.strip()] = value.strip()
            logging.debug("Parsed %s: %s", key.strip(), value.strip())
        elif current_section == "Reads":

            # Split the line into key-value pairs and store them as integers
            key, value = line.split(',', 1)
            data["Reads"][key.strip()] = int(value.strip())

            logging.debug(
                "Parsed read cycle count: %s = %d",
                key.strip(), int(value.strip())
            )
        elif current_section == "BCLConvert_Data":
            if not bcl_data_headers:

                # The first line in the BCLConvert_Data section contains
                # headers
                bcl_data_headers = line.split(',')
                logging.debug("Parsed headers: %s", bcl_data_headers)
            else:

                # Subsequent lines contain sample data
                values = line.split(',')
                sample_data = {
                    bcl_data_headers[i]: values[i]
                    for i in range(len(bcl_data_headers))}
                data["BCLConvert_Data"].append(sample_data)
                logging.debug("Parsed BCLConvert_Data sample: %s", sample_data)
        elif current_section == "Cloud_Data":
            if not cloud_data_headers:

                # The first line in the Cloud_Data section contains headers
                cloud_data_headers = line.split(',')
                logging.debug("Parsed headers: %s", cloud_data_headers)
            else:

                # Subsequent lines contain sample data
                values = line.split(',')
                sample_data = {
                    cloud_data_headers[i]: values[i]
                    for i in range(len(cloud_data_headers))}

                # Add Sample_Plate attribute based on ProjectName
                sample_data["Sample_Plate"] = sample_data.get(
                    "ProjectName", "")
                data["Cloud_Data"].append(sample_data)
                logging.debug("Parsed Cloud_Data sample: %s", sample_data)

    # Adhere to the same naming scheme as the MiSeq Reads section
    data["Reads"]['forward_read_length'] = data["Reads"]["Read1Cycles"]
    data["Reads"]['reverse_read_length'] = data["Reads"]["Read2Cycles"]

    logging.info("Finished parsing NextSeq sample sheet: %s", file_path)
    return data


def process_sample(
    commit: str,
    data: Dict[str, str],
    path: str,
) -> List[CustomBox]:
    """
    Processes a single sample from the sample sheet and creates a metadata
    object.

    :param commit: String of the Teacup COWBAT commit
    :param data: Dictionary containing sample data.
    :param path: Path to the output directory.
    :return: List of CustomBox objects
    """
    # Initialise a list to store the sample metadata
    samples = []

    # Determine whether the sample is from a MiSeq or NextSeq run and set the
    # key to access the dictionary appropriately
    try:
        _ = data["Data"]
        key = "Data"
    except KeyError:
        key = "Cloud_Data"

    for i, sample in enumerate(data[key]):
        logging.info("Processing sample: %s", sample["Sample_ID"])

        # Try and replicate the Illumina rules to create file names from
        # "Sample_Name"
        sample_name = sample_namer(raw_sample_name=sample["Sample_ID"])
        logging.debug("Sanitized sample name: %s", sample_name)

        samples = populate_metadata(
            basic_assembly=False,
            commit=commit,
            fastq_name=sample_name,
            samples=samples,
            sequence_path=path,
            data=data,
            sample_number=i
        )

    # Return the list of samples
    return samples


def sample_namer(raw_sample_name: str) -> str:
    """
    Tries to replicate the Illumina rules to create names from 'Sample_Name'.

    This function replaces spaces and certain special characters in the sample
    name with hyphens or removes them entirely to create a valid file name.

    :param raw_sample_name: The raw sample name to be processed.
    :return: A sanitized sample name suitable for use as a file name.
    """
    logging.debug("Sanitizing sample name: %s", raw_sample_name)

    # Characters to replace with hyphens
    replace_with_hyphen = [" ", ".", "=", "/", "---", "--"]

    # Characters to remove
    remove_chars = ["+", "#"]

    # Strip trailing whitespace from the raw sample name
    sample_name = raw_sample_name.rstrip()

    # Replace specified characters with hyphens
    for char in replace_with_hyphen:
        sample_name = sample_name.replace(char, "-")
        logging.debug("Replaced '%s' with '-': %s", char, sample_name)

    # Remove specified characters
    for char in remove_chars:
        sample_name = sample_name.replace(char, "")
        logging.debug("Removed '%s': %s", char, sample_name)

    logging.debug("Final sanitized sample name: %s", sample_name)
    return sample_name


def populate_metadata(
    basic_assembly: bool,
    commit: str,
    fastq_name: str,
    samples: List[CustomBox],
    sequence_path: str,
    data: Dict = None,
    sample_number: int = None
) -> List[CustomBox]:
    """
    Populate metadata for a given FASTQ file and update the list of samples.

    This function creates a metadata object for a given FASTQ file, sets up the
    necessary directories and files, and updates the list of samples with the
    new metadata. If previous metadata exists, it appends it to the samples
        list.

    :param fastq_name: Name of the FASTQ file (without extension).
    :param samples: List of CustomBox objects representing samples.
    :param sequence_path: Path to the directory containing the sequence files.
    :return: Updated list of CustomBox objects representing samples.

    Example usage:
    >>> samples = []
    >>> updated_samples = populate_metadata('sample1', samples,
        '/path/to/sequences')
    >>> print(updated_samples)
    [<CustomBox object at 0x...>]
    """

    # Initialize a new CustomBox for metadata
    metadata = CustomBox()
    metadata.name = fastq_name

    # Create the .general attribute
    metadata.general = CustomBox()

    # Create the .commands attribute
    metadata.commands = CustomBox()
    metadata.commands.test = 'test'

    # Set the destination folder
    output_dir = os.path.join(sequence_path, fastq_name)

    # Make the destination folder
    make_path(output_dir)
    logging.debug("Created output directory: %s", output_dir)

    # Add the output directory to the metadata
    metadata.general.output_directory = output_dir

    # Add the path to the JSON output file
    metadata.json_file = os.path.join(
        metadata.general.output_directory,
        f'{metadata.name}_metadata.json'
    )

    # Grab metadata from previous runs
    previous_metadata = read_metadata(sample=metadata)

    logging.debug("Previous metadata: %s", previous_metadata)

    # Update samples (if required)
    if previous_metadata:
        samples.append(previous_metadata)
        logging.info("Updated samples with previous metadata")
    # Otherwise, create new metadata using the appropriate method based on if
    # the sample sheet was provided
    else:
        if basic_assembly:
            samples.append(
                basic_metadata_populate(
                    commit=commit,
                    fastq_name=fastq_name,
                    metadata=metadata,
                    output_dir=output_dir,
                    sequence_path=sequence_path
                )
            )
        else:
            samples.append(
                sample_sheet_metadata_populate(
                    commit=commit,
                    data=data,
                    sample_name=fastq_name,
                    sequence_path=sequence_path,
                    sample_number=sample_number
                )
            )
    return samples


def basic_metadata_populate(
    commit: str,
    fastq_name: str,
    metadata: CustomBox,
    output_dir: str,
    sequence_path: str
) -> CustomBox:
    """
    Populate basic metadata for a given FASTQ file.

    This function populates a CustomBox object with metadata for a given FASTQ
    file, including linking FASTQ files to the output directory, setting log
    file paths, and adding commit information.

    :param commit: Commit information to be added to the metadata.
    :param fastq_name: Name of the FASTQ file (without extension).
    :param metadata: CustomBox object to be populated with metadata.
    :param output_dir: Path to the output directory.
    :param sequence_path: Path to the directory containing the sequence files.
    :return: Populated CustomBox object with metadata.

    Example usage:
    >>> metadata = CustomBox()
    >>> populated_metadata = basic_metadata_populate(
    ...     commit='abc123',
    ...     fastq_name='sample1',
    ...     metadata=metadata,
    ...     output_dir='/path/to/output',
    ...     sequence_path='/path/to/sequences'
    ... )
    >>> print(populated_metadata)
    <CustomBox object at 0x...>
    """
    logging.info("Added metadata for sample: %s", fastq_name)

    # Get the fastq files specific to the fastq_name
    specific_fastq = glob(
        os.path.join(sequence_path, f'{fastq_name}*.fastq*')
    )
    logging.debug("Specific FASTQ files: %s", specific_fastq)

    # Initialize the run category
    metadata.run = CustomBox()

    # Create the .fastq_files category of metadata
    metadata.general.fastq_files = []

    # Link the files to the output folder
    for fastq in sorted(specific_fastq):
        relative_symlink(fastq, output_dir)
        metadata.general.fastq_files.append(
            os.path.join(output_dir, os.path.basename(fastq))
        )
        logging.debug("Linked FASTQ file: %s", fastq)

    # Add the log locations to the metadata object
    metadata.general.logout = os.path.join(
        sequence_path, metadata.name, f'{metadata.name}_log_out.txt'
    )
    metadata.general.log_err = os.path.join(
        sequence_path, metadata.name, f'{metadata.name}_log_err.txt'
    )

    # Add the path to the JSON output file
    metadata.json_file = os.path.join(
        metadata.general.output_directory,
        f'{metadata.name}_metadata.json'
    )

    # Add the commit information
    metadata.general.commit = commit

    # Run the read length method
    read_length(sample=metadata)

    return metadata


def sample_sheet_metadata_populate(
    commit: str,
    data: Dict[str, Any],
    sample_name: str,
    sequence_path: str,
    sample_number: int
) -> CustomBox:
    """
    Populate metadata for a sample from a sample sheet.

    This function creates a CustomBox object for a sample, sets up the
    necessary directories and files, and populates the metadata with
    information from the sample sheet.

    :param commit: Commit information to be added to the metadata.
    :param data: Dictionary containing sample sheet data.
    :param sample_name: Name of the sample.
    :param sequence_path: Path to the directory containing the sequence files.
    :param sample_number: Sample number to be added to the metadata.
    :return: Populated CustomBox object with metadata.

    Example usage:
    >>> data = {'Header': {'Project': 'ExampleProject'}, 'Sample_ID': 'S1'}
    >>> metadata = sample_sheet_metadata_populate(
    ...     commit='abc123',
    ...     data=data,
    ...     sample_name='sample1',
    ...     sequence_path='/path/to/sequences',
    ...     sample_number=1
    ... )
    >>> print(metadata)
    <CustomBox object at 0x...>
    """
    # Create a CustomBox object for storing nested static variables
    strain_metadata = CustomBox(default_box=True)

    # Set the sample name in the object
    strain_metadata.name = sample_name

    # Create the commands attribute
    strain_metadata.commands = CustomBox()

    # Create the 'General' category for strain_metadata
    strain_metadata.general = CustomBox(
        {'output_directory': os.path.join(sequence_path, sample_name),
         'pipeline_commit': commit, 'log_out': os.path.join(
             sequence_path, sample_name, f'{sample_name} _log_out.txt'),
         'log_err': os.path.join(
             sequence_path, sample_name, f'{sample_name} _log_err.txt')})

    # Add the path to the JSON output file
    strain_metadata.json_file = os.path.join(
        strain_metadata.general.output_directory,
        f'{strain_metadata.name}_metadata.json'
    )
    logging.debug("Set general attributes for sample: %s", sample_name)

    # Add the header object to strain_metadata
    strain_metadata.run = CustomBox(copy.copy(data['Header']))

    # Capture Sample_ID, Sample_Name, I7_Index_ID, index1, I5_Index_ID,
    # index2, Sample_Project
    for key, value in data.items():
        strain_metadata.run[key] = value if value else "NA"
        logging.debug("Set run.%s = %s", key, value if value else "NA")

    # Add the sample number
    strain_metadata.run.SampleNumber = sample_number
    logging.debug(
        "Set run.SampleNumber = %d",
        strain_metadata.run.SampleNumber)

    return strain_metadata


def basic(commit: str, sequence_path: str) -> List[CustomBox]:
    """
    Processes FASTQ files, extracts metadata, and sets up the output
    directories.

    :param commit: String of the commit
    :param sequence_path: Path to the directory containing FASTQ files.
    :return: A list of CustomBox objects containing sample metadata.
    """
    logging.info("Processing FASTQ files in directory: %s", sequence_path)

    # Initialize a list of samples
    samples: List[CustomBox] = []

    # Grab any .fastq files in the path
    fastq_files = glob(os.path.join(sequence_path, '*.fastq*'))
    logging.debug("Found FASTQ files: %s", fastq_files)

    # Extract the base name of the globbed name + path provided
    fastq_names = map(lambda x: os.path.split(x)[1], filer(fastq_files))

    # Iterate through the names of the fastq files
    for fastq_name in sorted(fastq_names):
        logging.info("Processing FASTQ file: %s", fastq_name)
        samples = populate_metadata(
            basic_assembly=True,
            commit=commit,
            fastq_name=fastq_name,
            samples=samples,
            sequence_path=sequence_path
        )

    return samples


def read_length(sample: CustomBox) -> None:
    """
    Calculates the read length of the FASTQ files. Short reads will not be
    able to be assembled properly with the default parameters used for
    SKESA.

    :param sample: CustomBox object containing sample metadata.
    """
    logging.info('Estimating read lengths of FASTQ files')

    # Populate empty attributes
    sample.run.Date = 'NA'
    sample.run.Sample_Plate = 'NA'

    # Only perform this step if the forward and reverse lengths have
    # not been loaded into the metadata
    if not hasattr(sample.run, 'forward_length') and \
            not hasattr(sample.run, 'reverse_length'):

        # Initialize the .header and .commands attributes for each sample
        sample.header = CustomBox()
        sample.commands = CustomBox()

        # Only process the samples if the file type is a list
        if isinstance(sample.general.fastq_files, list):
            # Set the forward fastq to be the first entry in the list
            forward_fastq = sorted(sample.general.fastq_files)[0]
            logging.debug(
                "Processing forward FASTQ file: %s",
                forward_fastq)

            # Read the first 1000 lines of the FASTQ file
            forward_reads = read_first_1000_lines(forward_fastq)

            # Decode the bytes object to a string
            forward_reads_str = decode_bytes(
                forward_reads, sample, 'forward_length')

            # Calculate the length of the forward reads
            sample.run.forward_length = calculate_read_length(
                forward_reads_str)
            logging.debug(
                "Forward read length: %d",
                sample.run.forward_length)

            # For paired end analyses, also calculate the length of the
            # reverse reads
            if len(sample.general.fastq_files) == 2:
                reverse_fastq = sorted(sample.general.fastq_files)[1]
                logging.debug(
                    "Processing reverse FASTQ file: %s",
                    reverse_fastq)
                reverse_reads = read_first_1000_lines(reverse_fastq)

                # Decode the bytes object to a string
                reverse_reads_str = decode_bytes(
                    reverse_reads, sample, 'reverse_length')

                # Calculate the length of the reverse reads
                sample.run.reverse_length = calculate_read_length(
                    reverse_reads_str)
                logging.debug(
                    "Reverse read length: %d",
                    sample.run.reverse_length)

            # Populate metadata of single end reads with 'NA'
            else:
                sample.run.reverse_length = 0
                logging.debug(
                    "Single end read detected, reverse length set to 0")


def read_first_1000_lines(file_path: str) -> bytes:
    """
    Reads the first 1000 lines of a FASTQ file.

    :param file_path: Path to the FASTQ file.
    :return: The first 1000 lines of the file as bytes.
    """
    logging.debug("Reading first 1000 lines of file: %s", file_path)
    lines = []

    # Determine the appropriate open function based on file extension
    open_func = gzip.open if file_path.endswith('.gz') else open

    try:
        # Open the file and read the first 1000 lines
        with open_func(file_path, 'rb') as f:
            for _ in range(NUM_LINES_TO_READ):
                line = f.readline()
                if not line:
                    break  # Stop if end of file is reached
                lines.append(line)
    except (OSError, IOError) as exc:
        logging.error("Error reading file %s: %s", file_path, exc)
        return b''

    # Join the list of lines into a single bytes object
    logging.debug("Successfully read first 1000 lines of file: %s", file_path)
    return b''.join(lines)


def decode_bytes(data: bytes, sample: CustomBox, attr: str) -> Optional[str]:
    """
    Decodes a bytes object to a string and handles UnicodeDecodeError.

    :param data: The bytes object to decode.
    :param sample: The sample object to update in case of error.
    :param attr: The attribute to set to 0 in case of error.
    :return: The decoded string or None if decoding fails.
    """
    logging.debug("Decoding bytes for attribute: %s", attr)
    try:
        return data.decode('utf-8')
    except UnicodeDecodeError:
        setattr(sample.run, attr, 0)
        logging.error(
            "UnicodeDecodeError: Setting %s to 0 for sample: %s",
            attr,
            sample.name)
        return None


def calculate_read_length(reads: Optional[str]) -> int:
    """
    Calculates the length of the reads from the FASTQ file content.

    :param reads: The FASTQ file content as a string.
    :return: The maximum read length or 0 if reads is None.
    """
    logging.debug("Calculating read length")
    if reads is None:
        return 0

    try:
        return max(
            len(sequence) for iterator, sequence in
            enumerate(reads.split('\n'))
            if iterator % 4 == 1
        )
    except (ValueError, TypeError):
        logging.error("Error calculating read length")
        return 0
