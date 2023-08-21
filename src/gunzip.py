import gzip
import os

# Gunzip file from a path to another
def gunzip(src_path, dst_path):
    # Open input, output file stream
    with gzip.open(src_path, 'rb') as src_file, open(dst_path, 'wb') as dst_file:
        # Just write source file content to destination file
        dst_file.write(src_file.read())
