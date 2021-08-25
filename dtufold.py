"""Alphafold processing for snake toxin sequences.

  Typical usage example:

  ... = func(...)

  where ...

"""
import os
import sys
import shutil

import traceback  # for error dumps
import json  # for error dumps
from pathlib import Path

import tensorflow as tf

import dtuhash
from dtucitations import write_citations
from dtuclean import clean_jobname, clean_query_sequence
from dtufasta_iter import fasta_iter
from dtufileutils import save_results, zip_job_results
from dtummseqs import run_mmseqs2, process_homooligomers
from dtuplots import make_plots, show_pdb, plot_confidence, plot_plddt_legend
from dtupredict_structure import predict_structure
from dtutemplates import mk_template, mk_mock_template
from dtutimecheck import time_check

# @markdown ### MSA mode, Amber and Templates settings
msa_mode = "MMseqs2 (UniRef+Environmental)"
num_models = 5  # number of models to use (1-5)
use_msa = True if msa_mode.startswith("MMseqs2") else False
use_env = True if msa_mode == "MMseqs2 (UniRef+Environmental)" else False
use_custom_msa = True if msa_mode == "custom" else False
use_amber = True
use_templates = True

# Do not process homooligomers
homooligomer = 1

# @title Open FASTA file. This contains the 2000-odd sequences marked as toxins.
fasta_dir = '/content/drive/MyDrive/Shared/DTU'
fasta_filename = 'uniprot-taxonomy_8570+keyword_toxin+reviewed_yes.fasta'
fasta_path = os.path.join(fasta_dir, fasta_filename)
results_dir = '/content/drive/MyDrive/Shared/DTU/results_new/'


def show_this_run():
    """ Display some useful information about the run"""
    print("Current directory: ", os.getcwd())

    print("use_amber: ", use_amber)
    print("use_msa: ", use_msa)
    print("use_templates: ", use_templates)
    print("use_env: ", use_env)

    print("num_models: ", num_models)
    print("homooligomer: ", homooligomer)

    # Print disk usage if running on Colab
    if 'COLAB_GPU' in os.environ:
        print("I'm running on Colab")

        total, used, free = shutil.disk_usage("/")

        print("Total: %d GiB" % (total // (2 ** 30)))
        print("Used: %d GiB" % (used // (2 ** 30)))
        print("Free: %d GiB" % (free // (2 ** 30)))


# @markdown ### Procedure: Write log file headers
def write_log_headers(text_file):
    text_file.write("num_models=%s\n" % num_models)
    text_file.write("use_amber=%s\n" % use_amber)
    text_file.write("use_msa=%s\n" % use_msa)
    text_file.write("msa_mode=%s\n" % msa_mode)
    text_file.write("use_templates=%s\n" % use_templates)
    text_file.write("homooligomer=%s\n" % homooligomer)


def write_err(jobname, query_sequence, feature_dict, use_model):
    # Here if error. Save as much as we can about it
    with open(f"{jobname}.err", "w") as err_file:
        err_file.write(traceback.format_exc())
        err_file.write(query_sequence)
        json.dump(list(feature_dict), err_file)
        json.dump(use_model, err_file)


# @title Import libraries
# setup the model
if "model" not in dir():
    # hiding warning messages
    import warnings
    from absl import logging

    warnings.filterwarnings('ignore')
    logging.set_verbosity("error")
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
    tf.get_logger().setLevel('ERROR')

    from alphafold.data import pipeline
    from alphafold.model import data
    from alphafold.model import config
    from alphafold.model import model

    # plotting libraries

if "relax" not in dir():
    sys.path.insert(0, '/usr/local/lib/python3.7/site-packages/')

# @title Main loop


# Start processing
model_params = {}

fiter = fasta_iter(fasta_path)
loop_count = 0  # for debugging

for ff in fiter:

    time_check("Starting sequence: ")

    """
    # Uncomment to run only two sequences
    loop_count = loop_count + 1
    if loop_count > 2:
       fiter.close()
       raise ValueError
    """

    #################################################
    # get sequence
    db_type, accession, id_name, headerStr, seq = ff

    query_sequence = clean_query_sequence(seq)
    cleaned_jobname = clean_jobname(accession + id_name)
    jobname = dtuhash.add_hash(cleaned_jobname, query_sequence)
    print("jobname: ", jobname)
    print("query_sequence: {0} : {1}".format(len(query_sequence), query_sequence))

    #################################################
    # Restartability - skip over existing results files
    results_base = os.path.join(results_dir, jobname)
    results_file = Path(results_base + ".result.zip")

    if results_file.is_file():
        print("Skipping ", results_file)
        continue  # file exists

    a3m_file = f"{jobname}.a3m"

    # first line of job fasta file is the query sequence
    with open(f"{jobname}.fasta", "w") as text_file:
        text_file.write(">1\n%s" % query_sequence)

    with open(f"{jobname}.log", "w") as log_file:
        write_log_headers(log_file)

    #################################################
    # @title Gather input features, predict structure

    # collect model weights
    # use_model, model_config, model_params, model_runner_1, model_runner_3 = collect_model_weights(num_models)
    use_model = {}
    if "model_params" not in dir(): model_params = {}
    for model_name in ["model_1", "model_2", "model_3", "model_4", "model_5"][:num_models]:
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

    time_check("Model weights collected: ")

    #################################################
    # parse TEMPLATES
    if use_templates and os.path.isfile(f"{jobname}_hhm.ffindex"):
        template_features = mk_template(jobname)
    else:
        template_features = mk_mock_template(query_sequence * homooligomer)

    #################################################
    # parse MSA
    prefix = dtuhash.add_hash('tmp', query_sequence)

    a3m_dir = prefix + "_env" if use_env else prefix
    a3m_lines = run_mmseqs2(query_sequence, prefix, use_env=True, filter=True)
    msa, deletion_matrix = pipeline.parsers.parse_a3m(a3m_lines)
    msas, deletion_matrices = process_homooligomers(msa, deletion_matrix, homooligomer)

    time_check("Gathering features: ")

    # gather features
    feature_dict = {
        **pipeline.make_sequence_features(sequence=query_sequence * homooligomer,
                                          description="none",
                                          num_res=len(query_sequence) * homooligomer),
        **pipeline.make_msa_features(msas=msas, deletion_matrices=deletion_matrices),
        **template_features
    }

    # predict_structure can throw exceptions
    # ValueError: Amber minimization can only be performed on proteins with well-defined residues.
    # This protein contains at least one residue with no atoms.
    try:
        outs = predict_structure(jobname, feature_dict,
                                 Ls=[len(query_sequence)] * homooligomer,
                                 model_params=model_params, use_model=use_model,
                                 do_relax=use_amber)
    except:

        # Here if error. Save as much as we can about it
        write_err(jobname, query_sequence, feature_dict, use_model)
        write_citations(use_msa, use_env, use_templates, use_amber, jobname)
        zip_job_results(jobname, results_dir, a3m_dir)

        # Traceback uses sys.exc_info() to get current exception so no need to pass it in
        traceback.print_exc()
        time_check("Sequence processed with value error")
        continue

    time_check("Structure predicted: ")

    # Here if no problems with predict_structure - save output
    save_results(jobname, outs, msas, deletion_matrices, feature_dict, template_features)

    ##################################################
    # @title Display 3D structure {run: "auto"}

    color = "lDDT"  # @param ["chain", "lDDT", "rainbow"]
    show_sidechains = False
    show_mainchains = False

    # Display the models with vw.show(), plconfi.show() and plddt.show() respectively
    make_plots(msa, jobname, query_sequence, outs, homooligomer, num_models)
    for model_num in range(1, num_models + 1):
        vw = show_pdb(model_num, show_sidechains, show_mainchains, color)
        #def plot_confidence(Ln, plddt, pae, model_num=1, homooligomer=1):
        plconfi = plot_confidence(len(query_sequence), outs[model_name]["plddt"], outs[model_name]["pae"],
                                  homooligomer)

    if color == "lDDT":
        plddt = plot_plddt_legend()

    # Save results
    write_citations(use_msa, use_custom_msa, use_env, use_templates, use_amber, jobname)
    zip_job_results(jobname, results_dir, a3m_dir)
    time_check("Sequence processed")
