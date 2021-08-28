"""Predict peptide structures with Alphafold

    Typical usage example:

    set_bfactor(unrelaxed_pdb_path, plddts[r], idx_res, chains)
    template_features = mk_template(jobname)

    AlphaFold produces a per-residue estimate of its confidence on a scale from 0 - 100 . This confidence measure is
    called pLDDT and corresponds to the model’s predicted score on the lDDT-Cα metric. It is stored in the B-factor
    fields of the mmCIF and PDB files available for download (although unlike a B-factor, higher pLDDT is better).
    pLDDT is also used to colour-code the residues of the model in the 3D structure viewer. The following rules of
    thumb provide guidance on the expected reliability of a given region:

    - Regions with pLDDT > 90 are expected to be modelled to high accuracy. These should be suitable for any
      application that benefits from high accuracy (e.g. characterising binding sites).
    - Regions with pLDDT between 70 and 90 are expected to be modelled well (a generally good backbone prediction).
    - Regions with pLDDT between 50 and 70 are low confidence and should be treated with caution.
    - The 3D coordinates of regions with pLDDT < 50 often have a ribbon-like appearance and should not be interpreted.
      We show in our paper that pLDDT < 50 is a reasonably strong predictor of disorder, i.e. it suggests such a region
      is either unstructured in physiological conditions or only structured as part of a complex.
    - Structured domains with many inter-residue contacts are likely to be more reliable than extended linkers or
      isolated long helices.
    - Unphysical bond lengths and clashes do not usually appear in confident regions. Any part of a structure with
      several of these should be disregarded.

    Note that the PDB and mmCIF files contain coordinates for all regions, regardless of their pLDDT score.
    It is up to the user to interpret the model judiciously, in accordance with the guidance above.

    https://alphafold.ebi.ac.uk/faq

"""

# @title Procedure: Set B-Factor
from string import ascii_uppercase

import alphafold
import haiku
import numpy as np
from alphafold.common import protein
from alphafold.model import config
from alphafold.model import data
from alphafold.model import features as affeatures
from alphafold.model import model
from alphafold.relax import relax


def set_bfactor(pdb_filename, bfac: list[float], idx_res: list[str], chains: list[str]):
    """Store B-factor in PDB file.

    AlphaFold produces a per-residue estimate of its confidence on a scale from 0 - 100 . This confidence measure is
    called pLDDT and corresponds to the model’s predicted score on the lDDT-Cα metric. It is stored in the B-factor
    fields of the mmCIF and PDB files available for download (although unlike a B-factor, higher pLDDT is better).

    :param str pdb_filename: Name of PDB file.
    :param list[float] bfac: B-Factor.
    :param list[str] idx_res: List of residues.
    :param list[str] chains: Peptide sequence per chain.

    :example:

    >>> set_bfactor(pdb_path, plddts[r], idx_res, chains)

    """

    I = open(pdb_filename, "r").readlines()
    O = open(pdb_filename, "w")
    for line in I:
        if line[0:6] == "ATOM  ":
            seq_id = int(line[22:26].strip()) - 1
            seq_id = np.where(idx_res == seq_id)[0][0]
            O.write(f"{line[:21]}{chains[seq_id]}{line[22:60]}{bfac[seq_id]:6.2f}{line[66:]}")
    O.close()


# @title Function collect_model_weights()
def collect_model_weights(num_models):
    use_model = {}
    # Model_runner_1 is used for models 1 and 2
    # Model_runner_3 is used for models 3, 4 and 5
    model_runner_1 = None
    model_runner_3 = None
    # if "model_params" not in dir(): model_params = {}
    model_params = {}
    for model_name in ["model_1", "model_2", "model_3", "model_4", "model_5"][:num_models]:
        print(f"configuring {model_name}")
        use_model[model_name] = True
        if model_name not in model_params:
            model_params[model_name] = data.get_model_haiku_params(model_name=model_name + "_ptm", data_dir=".")
            if model_name == "model_1":
                model_config = config.model_config(model_name + "_ptm")
                model_config.data.eval.num_ensemble = 1
                model_runner_1 = model.RunModel(model_config, model_params[model_name])
            if model_name == "model_3":
                model_config = config.model_config(model_name + "_ptm")
                model_config.data.eval.num_ensemble = 1
                model_runner_3 = model.RunModel(model_config, model_params[model_name])
        return use_model, model_config, model_params, model_runner_1, model_runner_3


# @title Function: Predict Structure
def predict_structure(prefix, model_runner_1: alphafold.model.model.RunModel,
                      model_runner_3: alphafold.model.model.RunModel,
                      feature_dict, Ls: list[int], model_params: haiku.Params, use_model, do_relax=False,
                      random_seed=0):
    """Predicts structure using AlphaFold for the given sequence.

    :param str prefix: Job name.
    :param alphafold.model.model.RunModel model_runner_1: Model_runner_1 is used for models 1 and 2.
    :param alphafold.model.model.RunModel model_runner_3: Model_runner_3 is used for models 3, 4 and 5.
    :param dict feature_dict: Features.
    :param list[int] Ls: Number of residues in each chain. I.e. for two homooligomers of length 10, a list [10,10s]
    :param haiku.Params: model_params: Model parameters.
    :param bool use_model: Whether to use the model.
    :param bool do_relax: Whether to use Amber relaxation.
    :param int random_seed: Random seed.

    :return: template features
    :rtype: TemplateSearchResult

    :example:

    >>> outs = predict_structure(jobname, feature_dict,
                    Ls=[len(query_sequence)]*homooligomer,
                    model_params=model_params, use_model=use_model,
                    do_relax=use_amber)
    """

    # Minkyung's code
    # add big enough number to residue index to indicate chain breaks
    idx_res = feature_dict['residue_index']
    L_prev = 0
    # Ls: number of residues in each chain
    for L_i in Ls[:-1]:
        idx_res[L_prev + L_i:] += 200
        L_prev += L_i
    chains = list("".join([ascii_uppercase[n] * L for n, L in enumerate(Ls)]))
    feature_dict['residue_index'] = idx_res

    # Run the models.
    plddts, paes = [], []
    unrelaxed_pdb_lines = []
    relaxed_pdb_lines = []

    print(f"Use_model {use_model}")

    for model_name, params in model_params.items():
        if model_name in use_model:
            print(f"running {model_name}")
            # swap params to avoid recompiling
            # note: models 1,2 have diff number of params compared to models 3,4,5
            if any(str(m) in model_name for m in [1, 2]): model_runner = model_runner_1
            if any(str(m) in model_name for m in [3, 4, 5]): model_runner = model_runner_3
            model_runner.params = params

            processed_feature_dict: affeatures.FeatureDict = model_runner.process_features(feature_dict,
                                                                                           random_seed=random_seed)
            # prediction_result is a dictionary of NumPy feature arrays
            prediction_result: dict = model_runner.predict(processed_feature_dict)
            unrelaxed_protein: protein.Protein = protein.from_prediction(processed_feature_dict, prediction_result)
            unrelaxed_pdb_lines.append(protein.to_pdb(unrelaxed_protein))
            plddts.append(prediction_result['plddt'])
            paes.append(prediction_result['predicted_aligned_error'])

            if do_relax:
                # Relax the prediction.
                amber_relaxer: relax.AmberRelaxation = relax.AmberRelaxation(max_iterations=0, tolerance=2.39,
                                                                             stiffness=10.0, exclude_residues=[],
                                                                             max_outer_iterations=20)
                relaxed_pdb_str, _, _ = amber_relaxer.process(prot=unrelaxed_protein)
                relaxed_pdb_lines.append(relaxed_pdb_str)

    # rerank models based on predicted lddt
    lddt_rank = np.mean(plddts, -1).argsort()[::-1]
    out = {}
    print("reranking models based on avg. predicted lDDT")
    for n, r in enumerate(lddt_rank):
        print(f"model_{n + 1} {np.mean(plddts[r])}")

        unrelaxed_pdb_path = f'{prefix}_unrelaxed_model_{n + 1}.pdb'
        with open(unrelaxed_pdb_path, 'w') as f:
            f.write(unrelaxed_pdb_lines[r])
        set_bfactor(unrelaxed_pdb_path, plddts[r], idx_res, chains)

        if do_relax:
            relaxed_pdb_path = f'{prefix}_relaxed_model_{n + 1}.pdb'
            with open(relaxed_pdb_path, 'w') as f: f.write(relaxed_pdb_lines[r])
            set_bfactor(relaxed_pdb_path, plddts[r], idx_res, chains)

        out[f"model_{n + 1}"] = {"plddt": plddts[r], "pae": paes[r]}
    return out
