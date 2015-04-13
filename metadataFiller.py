__author__ = 'blais'


def filler(metadata, geneSeekrList):
    """Properly populates the metadata dictionary - when I tried to populate the metadata dictionary within the
     multi-processed functions, it was returned as a list (likely due to how the files were returned. This essentially
     iterates through the list, and populates a dictionary appropriately"""
    # Make a copy of the metadata dictionary
    geneSeekrMetadata = metadata
    # Iterate through geneSeekrList
    for item in geneSeekrList:
        # each item in geneSeekrList is a dictionary entry
        # iterate through all the dictionaries
        for name in item:
            # The way the dictionaries were created, they should have the format:
            # metadataOutput[name]["1.General"]["geneSeekrProfile"] = geneSeekrOutput[name], so
            # e.g. "1.General"
            for generalCategory in item[name]:
                # e.g. "geneSeekrProfile"
                for specificCategory in item[name][generalCategory]:
                    # Populate the dictionary
                    if specificCategory not in geneSeekrMetadata[name][generalCategory]:
                        geneSeekrMetadata[name][generalCategory][specificCategory] = str(item[name][generalCategory][specificCategory])
    # Return the beautifully-populated dictionary
    return geneSeekrMetadata