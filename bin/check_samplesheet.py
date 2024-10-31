#!/usr/bin/env python

##############################################################################################
# Copyright 2022 The Johns Hopkins University Applied Physics Laboratory LLC
# All rights reserved.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this
# software and associated documentation files (the "Software"), to deal in the Software
# without restriction, including without limitation the rights to use, copy, modify,
# merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
# OR OTHER DEALINGS IN THE SOFTWARE.
#

"""Provide a command line tool to validate and transform tabular samplesheets."""

import os
import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path
import gzip

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FORMATS = (
        ".fq.gz",
        ".fastq.gz",
        ".fastq.clean.gz",
        ".fq.clean.gz",
        ".fastq",
        ".fq"
    )

    def __init__(
        self,
        sample_col="sample",
        first_col="fastq_1",
        second_col="fastq_2",
        single_col="single_end",
        dir_col="directory",
        needscompressing="needscompressing",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "sample").
            first_col (str): The name of the column that contains the first (or only)
                FASTQ file path (default "fastq_1").
            second_col (str): The name of the column that contains the second (if any)
                FASTQ file path (default "fastq_2").
            single_col (str): The name of the new column that will be inserted and
                records whether the sample contains single- or paired-end sequencing
                reads (default "single_end").
            dir_col (str): The name of the new column that will be inserted and
                records whether the sample is a directory of fastq files or not.

        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._first_col = first_col
        self._second_col = second_col
        self._dir_col = dir_col
        self._single_col = single_col
        self._needscompressing = needscompressing
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_sample(row)
        # check if the path of row['fastq_1'] is a directory or file
        # if directory, then row['fastq_1'] is the directory name
        # if file, then row['fastq_1'] is the file name , check if the path is a directory too
        if row['fastq_1'] and row['fastq_1'] != '' and not any(row['fastq_1'].endswith(extension) for extension in self.VALID_FORMATS):
            row[self._single_col] = True
            row[self._dir_col] = True
            self._seen.add((row[self._sample_col], row[self._first_col]))
        else:
            row[self._single_col] = True
            row[self._dir_col] = False
            self._validate_first(row)
            self._validate_second(row)
            self._validate_pair(row)
            self._seen.add((row[self._sample_col], row[self._first_col]))
        # for each attribute in row, strip spaces on either side
        for key in row:
            if isinstance(row[key], str):
                row[key] = row[key].strip()
                # if row[key] == 'trim':
                # uppercase the trim value
                # row[key] = row[key].upper()
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        assert len(row[self._sample_col]) > 0, "Sample input is required."
        # Sanitize samples slightly.
        row[self._dir_col] = False
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")
        for char in ['/', '\\', ':', '*', '?', '"', '<', '>', '|', '(', ')', '[', ']', '{', '}', '#', '%', '&', '+', '!', '@', '$', '^', '`', '~', ';', ',']:
            row[self._sample_col] = row[self._sample_col].replace(char, "")
    def compress_fastq(self,  fastqfile):
        # use gzip to compress the fastq file
        gzipfile = fastqfile + ".gz"
        with open(fastqfile, 'rb') as f_in:
            with gzip.open(gzipfile, 'wb') as f_out:
                f_out.writelines(f_in)
            f_out.close()
        f_in.close()
        return

    def _validate_first(self, row):
        """Assert that the first FASTQ entry is non-empty and has the right format."""
        assert len(row[self._first_col]) > 0, "At least the first FASTQ file is required."
        row[self._needscompressing] = self._validate_fastq_format(row[self._first_col])
        return row

    def _validate_second(self, row):
        """Assert that the second FASTQ entry has the right format if it exists."""

        if self._second_col not in row:
            row[self._second_col] = None
        if row[self._second_col] and len(row[self._second_col]) > 0:
            row[self._needscompressing] = self._validate_fastq_format(row[self._second_col])
        return row

    def _validate_pair(self, row):
        """Assert that read pairs have the same file extension. Report pair status."""
        if row[self._first_col] and row[self._second_col]:
            row[self._single_col] = False
            assert (
                Path(
                    row[self._first_col]).suffixes[-2:] == Path(row[self._second_col]).suffixes[-2:]
            ), "FASTQ pairs must have the same file extensions."
        else:
            row[self._single_col] = True
        return row

    def _validate_fastq_format(self, filename):
        """Assert that a given filename has one of the expected FASTQ extensions."""
        assert any(filename.endswith(extension) for extension in self.VALID_FORMATS), (
            f"The FASTQ file has an unrecognized extension: {filename}\n"
            f"It should be one of: {', '.join(self.VALID_FORMATS)}"
        )

        if filename.endswith('.gz'):
            return None
        elif filename.endswith('.fastq') or filename.endswith('.fq'):
            return True
        else:
            raise ValueError("The fastq file must be compressed with gzip and in fq or fastq format.")

    def validate_unique_samples(self):
        """
        Assert that the combination of sample name and FASTQ filename is unique.

        In addition to the validation, also rename the sample if more than one sample,
        FASTQ file combination exists.

        """
        assert len(self._seen) == len(
            self.modified), "The pair of sample name and FASTQ must be unique."
        if len({pair[0] for pair in self._seen}) < len(self._seen):
            counts = Counter(pair[0] for pair in self._seen)
            seen = Counter()
            for row in self.modified:
                sample = row[self._sample_col]
                seen[sample] += 1
                if counts[sample] > 1:
                    row[self._sample_col] = f"{sample}_T{seen[sample]}"


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    header = handle.readline().rstrip()
    # header = f"{header},single_end,directory"
    # peek = handle.read(2048)
    sniffer = csv.Sniffer()
    if not sniffer.has_header(header):
        logger.critical(
            f"The given sample sheet does not appear to contain a header.")
        # sys.exit(1)
    dialect = sniffer.sniff(header)
    handle.seek(0)
    return dialect


def check_samplesheet(file_in, file_out):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

    Validate the general shape of the table, expected columns, and each row. Also add
    an additional column which records whether one or two FASTQ reads were found.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following structure,
        see also the `viral recon samplesheet`_::

            sample,fastq_1,fastq_2
            SAMPLE_PE,SAMPLE_PE_RUN1_1.fastq.gz,SAMPLE_PE_RUN1_2.fastq.gz
            SAMPLE_PE,SAMPLE_PE_RUN2_1.fastq.gz,SAMPLE_PE_RUN2_2.fastq.gz
            SAMPLE_SE,SAMPLE_SE_RUN1_1.fastq.gz,



    """
    required_columns = {"sample", "fastq_1"}
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="", encoding='utf-8-sig') as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            logger.critical(
                f"The sample sheet **must** contain the column headers: {', '.join(required_columns)}.")
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_samples()
    header = list(reader.fieldnames)
    if "fastq_2" not in header:
        header.insert(1, "fastq_2")
    header.insert(1, "single_end")
    header.insert(1, "directory")
    header.insert(1, "needscompressing")
    # remove any header that is empty string
    header = [x for x in header if x != '']
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            # remove any key that is empty string
            row = {k: v for k, v in row.items() if k != ''}
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level,
                        format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
