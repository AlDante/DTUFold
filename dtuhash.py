"""Hash function for ensuring uniqueness of output files.

  Typical usage example:

  hashed_name = add_hash(jobname, query_sequence)

  where jobname is the current job and query_sequence is the peptide being analysed. hashed_name is the jobname
  with a 10-character SHA1 _Hash appended.

"""

import hashlib


def add_hash(x, y):
    """Returns hash of x and y.

        Args:
            x: An string, typically a filename.
            y: A string, typically a peptide sequence

        Returns:
            A 10 character sha1 _Hash value appended to x (the filename)

        Raises:
            None.
    """
    return x + "_" + hashlib.sha1(y.encode()).hexdigest()[:10]
