"""Predict peptide structures with Alphafold

    Typical usage example:

    serialize_result(f"{jobname}.outs", outs)
    outs = deserialize_result(f"{jobname}.outs")
    save_results(jobname, outs, msas, deletion_matrices, feature_dict, template_features)
    outs, msas, deletion_matrices, feature_dict, template_features = recover_results(jobname,
      outs, msas, deletion_matrices, feature_dict, template_features)


    ###  Code to recover the dump as follows:
    ###  with open(f"Q9PRP8SL117_ECHCA_7f5e67318d.out", "r") as out_file:
    ###      frozen = json.load(outs_file)
    ###      thawed = jsonpickle.decode(json.dumps(frozen))
    ###      print(thawed)

    https://alphafold.ebi.ac.uk/faq

"""

import fnmatch
import json
import os
from zipfile import ZipFile

import jsonpickle


# @title Procedure serialize_result()

def serialize_result(filename, data):
    """ Serialize a results file as a jsonpickle
    :param str filename: Name of file to save results in.
    :param str data: Data to save.

    :example:

    >>> serialize_result(f"{jobname}.outs", outs)
    """

    with open(filename, "w") as results_file:
        serialized = jsonpickle.encode(data)
        results_file.write(json.dumps(json.loads(serialized), indent=2))


def deserialize_result(filename):
    """ Deserialize an individual results file from a jsonpickle

    :param str filename: Name of results file.

    :return: template features
    :rtype: TemplateSearchResult

    :example:

    >>> outs = deserialize_result(filename)
    """

    with open(filename, "r") as results_file:
        frozen = json.load(results_file)
        thawed = jsonpickle.decode(json.dumps(frozen))
        return thawed


def save_results(jobname, outs, msas, deletion_matrices, feature_dict, template_features):
    """ Save all results files as jsonpickles """
    serialize_result(f"{jobname}.outs", outs)
    serialize_result(f"{jobname}.msas", msas)
    serialize_result(f"{jobname}.dels", deletion_matrices)
    serialize_result(f"{jobname}.feats", feature_dict)
    serialize_result(f"{jobname}.tfeats", template_features)


def recover_results(jobname):
    """ Read all results files from their jsonpickles """

    outs = deserialize_result(f"{jobname}.outs")
    msas = deserialize_result(f"{jobname}.msas")
    deletion_matrices = deserialize_result(f"{jobname}.dels")
    feature_dict = deserialize_result(f"{jobname}.feats")
    template_features = deserialize_result(f"{jobname}.tfeats")
    return outs, msas, deletion_matrices, feature_dict, template_features


# @title Procedure zipFilesInDir()
# https://thispointer.com/python-how-to-create-a-zip-archive-from-multiple-files-or-directory/

def zipFilesInDir(dirName, zipObj, pattern):
    # Iterate over all the files in directory
    for folderName, subfolders, filenames in os.walk(dirName):
        for filename in fnmatch.filter(filenames, pattern):
            # create complete filepath of file in directory
            filePath = os.path.join(folderName, filename)
            # Add file to zip
            zipObj.write(filePath, os.path.basename(filePath))


# @title Procedure zip_job_results
def zip_job_results(jobname, results_dir, a3m_dir):
    zip_path = os.path.join(results_dir, jobname + ".result.zip")
    log_path = jobname + ".log"
    a3m_path = a3m_dir
    # pdb_path = os.path.join(jobname, ".pdb")
    ddt_path = jobname + "_coverage_lDDT.png"
    bib_path = jobname + ".bibtex"
    pae_path = jobname + "_PAE.png"

    outs_path = jobname + ".outs"
    msas_path = jobname + ".msas"
    dels_path = jobname + ".dels"
    feats_path = jobname + ".feats"
    tfeats_path = jobname + ".tfeats"

    err_path = jobname + ".err"

    # Create a ZipFile object on my Google Drive and zip the results
    with ZipFile(zip_path, 'w') as zipObj:
        # Add multiple files to the zip
        zipObj.write(log_path)

        zipFilesInDir(a3m_path, zipObj, "*.*")

        zipFilesInDir(os.getcwd(), zipObj, jobname + "_*relaxed_model_*.pdb")

        if os.path.isfile(ddt_path): zipObj.write(ddt_path)
        if os.path.isfile(bib_path): zipObj.write(bib_path)
        if os.path.isfile(pae_path): zipObj.write(pae_path)
        if os.path.isfile(outs_path): zipObj.write(outs_path)
        if os.path.isfile(msas_path): zipObj.write(msas_path)
        if os.path.isfile(dels_path): zipObj.write(dels_path)
        if os.path.isfile(feats_path): zipObj.write(feats_path)
        if os.path.isfile(tfeats_path): zipObj.write(tfeats_path)
        if os.path.isfile(err_path): zipObj.write(err_path)
