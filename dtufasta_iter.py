"""Extract sequences from FASTA. Code taken from Biostars.

    Typical usage example:

    fiter = fasta_iter(fasta_name)

    where fasta_name is the name of the file containing the sequences in FASTA format.

    https://www.biostars.org/p/710/
    for python3


"""
from itertools import groupby


def fasta_iter(fasta_name):
    """Yields an iterator over a FASTA file.

    |
    Typical usage example: db_type, accession, id_name, headerStr, seq = fasta_iter("Toxins.fasta")

    |
    Arguments;
    fasta_name: Name of FASTA file to process.

    |
    Yields:
    An iterator with the following elements:

        db_type: Database ID ('sp' for SWISS-PROT/Uniprot, 'prf' for PRF etc.
                   (see https://en.wikipedia.org/wiki/FASTA_format)
        |
        accession:  Accession number for the sequence, e.g. 'P0CJ40'.
                (see https://www.uniprot.org/help/accession_numbers/)
        |
        id_name:    Entry name for the sequence, e.g. OXLA_BOTMA.
                    (see https://www.uniprot.org/help/entry_name)
        |
        headerStr3: Remaining header information, e.g.
                    L-amino-acid oxidase (Fragment) OS=Bothrops marajoensis OX=157554 PE=1 SV=1
        |
        seq:        Peptide sequence, e.g. AHDGNPLEECFREDDEEFFLEIAKNGLTATSNPKRVVIV

        |
        Modified from Brent Pedersen
        Correct Way To Parse A Fasta File In Python given a fasta file. yield tuples of header, sequence
    """
    "first open the FASTA file"
    fh = open(fasta_name)

    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()
        db_type, accession, headerStr2 = headerStr.split("|")
        id_name, headerStr3 = headerStr2.split(' ', 1)

        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())

        yield (db_type, accession, id_name, headerStr3, seq)
