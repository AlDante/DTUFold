"""Data cleansing for filenames and amino acid sequences.

  Typical usage example:

  cleaned_seq = clean_query_sequence(query_sequence)
  cleaned_job = clean_query_sequence(jobname)

  where jobname and query_sequence are strings.

"""

import re

def clean_query_sequence(seq):
    """Returns cleaned-up amino acid sequence.

        Args:
            seq: An unprocessed string of amino acids

        Returns:
            A contiguous string of uppercase characters with no whitespace

        Raises:
            None.
    """
    # remove whitespace
    query_sequence = "".join(seq.split())
    # remove non-alpha characters and uppercase the remainder
    query_sequence = re.sub(r'[^a-zA-Z]','', query_sequence).upper()
    return query_sequence

def clean_jobname(name):
    """Returns cleaned up job name.

        Args:
            name: A filename as a string.

        Returns:
            A cleaned filename as a contiguous string of characters with no whitespace

        Raises:
            None.
    """
    # remove whitespace
    jobname = "".join(name.split())
    jobname = re.sub(r'\W+', '', jobname)
    return jobname

