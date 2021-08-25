# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.


def fasta_iter(fasta_name):
    """
    modified from Brent Pedersen
    Correct Way To Parse A Fasta File In Python
    given a fasta file. yield tuples of header, sequence
    """
    "first open the file outside "
    fh = open(fasta_name)

    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()

        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())

        yield (headerStr, seq)

    # --------------
    def add_hash(x,y):
      return x+"_"+hashlib.sha1(y.encode()).hexdigest()[:5]

# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    # Start by connecting gdrive into the google colab
    from google.colab import drive

    drive.mount('/content/gdrive')

    # "/content/gdrive/MyDrive/Uniprot-All-Snakes-Toxins/uniprot-taxonomy_8570+keyword_toxin+reviewed_yes.fasta"

    import os

    fasta_dir = '/content/gdrive/MyDrive/Uniprot-All-Snakes-Toxins'
    fasta_filename = 'uniprot-taxonomy_8570+keyword_toxin+reviewed_yes.fasta'
    fasta_path = os.path.join(fasta_dir, fasta_filename)
    with open(fasta_path, 'r') as fasta_h:
        fasta_h.read(1)
        fasta_h.close()
    #!head    {fasta_path}

    """
    https://www.biostars.org/p/710/"
     for python3
    """
    from itertools import groupby

    fiter = fasta_iter(fasta_path)

    # Need to get Protein Sequence number here
    for ff in fiter:
        headerStr, seq = ff
        # print(headerStr)
        print(seq)

    #@title Input protein sequence, then hit `Runtime` -> `Run all`
    from google.colab import files
    import os
    import os.path
    import re
    import hashlib

    query_sequence = 'PIAQIHILEGRSDEQKETLIREVSEAISRSLDAPLTSVRVIITEMAKGHFGIGGELASK' #@param {type:"string"}
    # remove whitespaces
    query_sequence = "".join(query_sequence.split())
    query_sequence = re.sub(r'[^a-zA-Z]','', query_sequence).upper()

    jobname = 'test' #@param {type:"string"}
    # remove whitespaces
    jobname = "".join(jobname.split())
    jobname = re.sub(r'\W+', '', jobname)
    jobname = add_hash(jobname, query_sequence)


    with open(f"{jobname}.fasta", "w") as text_file:
        text_file.write(">1\n%s" % query_sequence)

    # number of models to use
    #@markdown ---
    #@markdown ### Advanced settings
    msa_mode = "MMseqs2 (UniRef+Environmental)" #@param ["MMseqs2 (UniRef+Environmental)", "MMseqs2 (UniRef only)","single_sequence","custom"]
    num_models = 5 #@param [1,2,3,4,5] {type:"raw"}
    use_msa = True if msa_mode.startswith("MMseqs2") else False
    use_env = True if msa_mode == "MMseqs2 (UniRef+Environmental)" else False
    use_custom_msa = True if msa_mode == "custom" else False
    use_amber = False #@param {type:"boolean"}
    use_templates = False #@param {type:"boolean"}
    #@markdown ---
    #@markdown ### Experimental options
    homooligomer = 1 #@param [1,2,3,4,5,6,7,8] {type:"raw"}
    save_to_google_drive = False #@param {type:"boolean"}
    #@markdown ---
    #@markdown Don't forget to hit `Runtime` -> `Run all` after updating form

    if homooligomer > 1:
      if use_amber:
        print("amber disabled: amber is not currently supported for homooligomers")
        use_amber = False
      if use_templates:
        print("templates disabled: templates are not currently supported for homooligomers")
        use_templates = False

    with open(f"{jobname}.log", "w") as text_file:
        text_file.write("num_models=%s\n" % num_models)
        text_file.write("use_amber=%s\n" % use_amber)
        text_file.write("use_msa=%s\n" % use_msa)
        text_file.write("msa_mode=%s\n" % msa_mode)
        text_file.write("use_templates=%s\n" % use_templates)
        text_file.write("homooligomer=%s\n" % homooligomer)

    # decide which a3m to use
    if use_msa:
      a3m_file = f"{jobname}.a3m"
    elif use_custom_msa:
      a3m_file = f"{jobname}.custom.a3m"
      if not os.path.isfile(a3m_file):
        custom_msa_dict = files.upload()
        custom_msa = list(custom_msa_dict.keys())[0]
        header = 0
        import fileinput
        for line in fileinput.FileInput(custom_msa,inplace=1):
          if line.startswith(">"):
             header = header + 1
          if line.startswith("#"):
            continue
          if line.rstrip() == False:
            continue
          if line.startswith(">") == False and header == 1:
             query_sequence = line.rstrip()
          print(line, end='')

        os.rename(custom_msa, a3m_file)
        print(f"moving {custom_msa} to {a3m_file}")
    else:
      a3m_file = f"{jobname}.single_sequence.a3m"
      with open(a3m_file, "w") as text_file:
        text_file.write(">1\n%s" % query_sequence)

    if save_to_google_drive == True:
      from pydrive.drive import GoogleDrive
      from pydrive.auth import GoogleAuth
      from google.colab import auth
      from oauth2client.client import GoogleCredentials
      auth.authenticate_user()
      gauth = GoogleAuth()
      gauth.credentials = GoogleCredentials.get_application_default()
      drive = GoogleDrive(gauth)
      print("You are logged into Google Drive and are good to go!")

    #------------------------------------------------
      # @title Install dependencies
      %% bash - s $use_amber $use_msa $use_templates

      USE_AMBER =$1
      USE_MSA =$2
      USE_TEMPLATES =$3

      if [ ! -f AF2_READY]; then
      # install dependencies
      pip - q
      install
      biopython
      pip - q
      install
      dm - haiku
      pip - q
      install
      ml - collections
      pip - q
      install
      py3Dmol

      # download model
      if [ ! -d "alphafold/"]; then
      git
      clone
      https: // github.com / deepmind / alphafold.git - -quiet
      (cd alphafold; git checkout 0bab1bf84d9d887aba5cfb6d09af1e8c3ecbc408 --quiet)
      mv
      alphafold
      alphafold_
      mv
      alphafold_ / alphafold.
      # remove "END" from PDBs, otherwise biopython complains
      sed - i
      "s/pdb_lines.append('END')//" / content / alphafold / common / protein.py
      sed - i
      "s/pdb_lines.append('ENDMDL')//" / content / alphafold / common / protein.py
    fi

    # download model params (~1 min)
    if [ ! -d "params/"]; then
    wget - qnc
    https: // storage.googleapis.com / alphafold / alphafold_params_2021 - 07 - 14.
    tar
    mkdir
    params
    tar - xf
    alphafold_params_2021 - 07 - 14.
    tar - C
    params /
    rm
    alphafold_params_2021 - 07 - 14.
    tar
fi
touch
AF2_READY
fi
# download libraries for interfacing with MMseqs2 API
if [ ${USE_MSA} == "True"] | |[${USE_TEMPLATES} == "True"]; then
if [ ! -f MMSEQ2_READY]; then
apt - get - qq - y
update
2 > & 1
1 > / dev / null
apt - get - qq - y
install
jq
curl
zlib1g
gawk
2 > & 1
1 > / dev / null
touch
MMSEQ2_READY
fi
fi
# setup conda
if [ ${USE_AMBER} == "True"] | |[${USE_TEMPLATES} == "True"]; then
if [ ! -f CONDA_READY]; then
wget - qnc
https: // repo.anaconda.com / miniconda / Miniconda3 - latest - Linux - x86_64.sh
bash
Miniconda3 - latest - Linux - x86_64.sh - bfp / usr / local
2 > & 1
1 > / dev / null
rm
Miniconda3 - latest - Linux - x86_64.sh
touch
CONDA_READY
fi
fi
# setup template search
if [ ${USE_TEMPLATES} == "True"] & &[! -f HH_READY]; then
conda
install - y - q - c
conda - forge - c
bioconda
kalign3 = 3.2
.2
hhsuite = 3.3
.0
python = 3.7
2 > & 1
1 > / dev / null
touch
HH_READY
fi
# setup openmm for amber refinement
if [ ${USE_AMBER} == "True"] & &[! -f AMBER_READY]; then
conda
install - y - q - c
conda - forge
openmm = 7.5
.1
python = 3.7
pdbfixer
2 > & 1
1 > / dev / null
(cd / usr / local / lib / python3.7 / site-packages; patch -s -p0 < / content / alphafold_ / docker / openmm.patch)
wget - qnc
https: // git.scicore.unibas.ch / schwede / openstructure / - / raw / 7102
c63615b64735c4941278d92b554ec94415f8 / modules / mol / alg / src / stereo_chemical_props.txt
mv
stereo_chemical_props.txt
alphafold / common /
touch
AMBER_READY
fi

# ------------------------------------------------
# @title Import libraries
# setup the model
if "model" not in dir():
    # hiding warning messages
    import warnings
    from absl import logging
    import os
    import tensorflow as tf

    warnings.filterwarnings('ignore')
    logging.set_verbosity("error")
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    tf.get_logger().setLevel('ERROR')

    import sys
    import numpy as np
    import pickle
    from alphafold.common import protein
    from alphafold.data import pipeline
    from alphafold.data import templates
    from alphafold.model import data
    from alphafold.model import config
    from alphafold.model import model
    from alphafold.data.tools import hhsearch

    # plotting libraries
    import py3Dmol
    import matplotlib.pyplot as plt
    import ipywidgets
    from ipywidgets import interact, fixed, GridspecLayout, Output

if use_amber and "relax" not in dir():
    sys.path.insert(0, '/usr/local/lib/python3.7/site-packages/')
    from alphafold.relax import relax


def mk_mock_template(query_sequence):
    # since alphafold's model requires a template input
    # we create a blank example w/ zero input, confidence -1
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


def mk_template(jobname):
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


def set_bfactor(pdb_filename, bfac, idx_res, chains):
    I = open(pdb_filename, "r").readlines()
    O = open(pdb_filename, "w")
    for line in I:
        if line[0:6] == "ATOM  ":
            seq_id = int(line[22:26].strip()) - 1
            seq_id = np.where(idx_res == seq_id)[0][0]
            O.write(f"{line[:21]}{chains[seq_id]}{line[22:60]}{bfac[seq_id]:6.2f}{line[66:]}")
    O.close()


def predict_structure(prefix, feature_dict, Ls, model_params, use_model, do_relax=False, random_seed=0):
    """Predicts structure using AlphaFold for the given sequence."""

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

    for model_name, params in model_params.items():
        if model_name in use_model:
            print(f"running {model_name}")
            # swap params to avoid recompiling
            # note: models 1,2 have diff number of params compared to models 3,4,5
            if any(str(m) in model_name for m in [1, 2]): model_runner = model_runner_1
            if any(str(m) in model_name for m in [3, 4, 5]): model_runner = model_runner_3
            model_runner.params = params

            processed_feature_dict = model_runner.process_features(feature_dict, random_seed=random_seed)
            prediction_result = model_runner.predict(processed_feature_dict)
            unrelaxed_protein = protein.from_prediction(processed_feature_dict, prediction_result)
            unrelaxed_pdb_lines.append(protein.to_pdb(unrelaxed_protein))
            plddts.append(prediction_result['plddt'])
            paes.append(prediction_result['predicted_aligned_error'])

            if do_relax:
                # Relax the prediction.
                amber_relaxer = relax.AmberRelaxation(max_iterations=0, tolerance=2.39,
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