"""Define templates for Alphafold

    Typical usage example:

    template_features = mk_mock_template(query_sequence)
    template_features = mk_template(jobname)

    where fasta_name is the name of the file containing the sequences in FASTA format.

    https://www.biostars.org/p/710/
    for python3

"""
import numpy as np
from alphafold.data import templates
from alphafold.data import pipeline
from alphafold.data.tools import hhsearch



# @title Function: Make Template
def mk_template(jobname):
    """creates a template for Alphafold

    :param int query_sequence: Peptide sequence.
    :return: template features
    :rtype: TemplateSearchResult

    :example:

    >>> template_features = mk_template(query_sequence)
    """

    template_featurizer = templates.TemplateHitFeaturizer(
        mmcif_dir="templates/",
        max_template_date="2100-01-01",
        max_hits=20,
        kalign_binary_path="kalign",
        release_dates_path=None,
        obsolete_pdbs_path=None)

    hhsearch_pdb70_runner = hhsearch.HHSearch(binary_path="hhsearch", databases=[jobname])

    a3m_lines = "\n".join(open(f"{jobname}.a3m", "r").readlines())
    hhsearch_result = hhsearch_pdb70_runner.query(a3m_lines)
    hhsearch_hits = pipeline.parsers.parse_hhr(hhsearch_result)
    templates_result = template_featurizer.get_templates(query_sequence=query_sequence,
                                                         query_pdb_code=None,
                                                         query_release_date=None,
                                                         hits=hhsearch_hits)
    return templates_result.features


def mk_mock_template(query_sequence):
    """creates a mock template for Alphafold with zero input, confidence -1

    Since Alphafold's model requires a template input, we create a blank example with zero input, confidence -1

    :param int query_sequence: Peptide sequence.
    :return: template features
    :rtype: dict

    :example:

    >>> template_features = mk_mock_template(query_sequence)
    """

    ln = len(query_sequence)
    output_templates_sequence = "-" * ln
    output_confidence_scores = np.full(ln, -1)
    templates_all_atom_positions = np.zeros((ln, templates.residue_constants.atom_type_num, 3))
    templates_all_atom_masks = np.zeros((ln, templates.residue_constants.atom_type_num))
    templates_aatype = templates.residue_constants.sequence_to_onehot(output_templates_sequence,
                                                                      templates.residue_constants.HHBLITS_AA_TO_ID)
    template_features = {'template_all_atom_positions': templates_all_atom_positions[None],
                         'template_all_atom_masks': templates_all_atom_masks[None],
                         'template_sequence': [f'none'.encode()],
                         'template_aatype': np.array(templates_aatype)[None],
                         'template_confidence_scores': output_confidence_scores[None],
                         'template_domain_names': [f'none'.encode()],
                         'template_release_date': [f'none'.encode()]}
    return template_features

