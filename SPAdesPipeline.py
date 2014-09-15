__author__ = 'akoziol'

import re
import os
from collections import defaultdict
import json
import runMetadataOptater

# The path is still hardcoded as, most of the time, this script is run from within Pycharm.
os.chdir("/media/nas1/akoziol/Pipeline_development/SPAdesPipelineSandbox")
path = os.getcwd()

runMetadata = runMetadataOptater.functionsGoNOW()

# print json.dumps(runMetadata, sort_keys=True, indent=4)

