# Formats

This folder contains the modules for reading and writing different spectroscopy 
file formats. Each module for a given format should define functions for:
1. Reading a file into a spectrum object.
2. Writing a file from a spectrum object.
3. Providing format specific metadata for reading/writing
4. Setting format specific metadata for reading/writing

The __init__.py file should define the supported formats

Currently and tobe supported formats are:
1. JCAMP DX format
2. csv format and associated tsv etc TODO (to find)
3. excel format in some form TODO (to find)
4. ool format in some form TODO (to find)
