# @title Function run_mmseqs2 - Call MMseqs2 to get MSA/templates

##########################################
# call mmseqs2
##########################################

import os
import random
import tarfile
import time

import requests
import tqdm.notebook

TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'


def process_homooligomers(msa, deletion_matrix, homooligomer: int = 1):
    """

    :param msa: query sequence to align with
    :param deletion_matrix:
    :param homooligomer: number of homooligomers
    :return: msas and deletion matrices per homooligomer
    """

    # make multiple copies of msa for each copy
    # AAA------
    # ---AAA---
    # ------AAA
    #
    # note: if you concat the sequences (as below), it does NOT work
    # AAAAAAAAA
    if homooligomer == 1:
        new_msas = [msa]
        new_mtxs = [deletion_matrix]
    else:
        new_msas = []
        new_mtxs = []
        for o in range(homooligomer):
            for msa, mtx in zip(new_msas, new_mtxs):
                num_res = len(msa[0])
                L = num_res * o
                R = num_res * (homooligomer - (o + 1))
                new_msas.append(["-" * L + s + "-" * R for s in msa])
                new_mtxs.append([[0] * L + m + [0] * R for m in mtx])
    return new_msas, new_mtxs


def run_mmseqs2(x, prefix, use_env=True, use_filter=True,
                use_templates=False, filter=None):
    def submit(seqs, mode, N=101):

        n, query = N, ""
        for seq in seqs:
            query += f">{n}\n{seq}\n"
            n += 1

        res = requests.post('https://a3m.mmseqs.com/ticket/msa', data={'q': query, 'mode': mode})
        try:
            out = res.json()
        except ValueError:
            out = {"status": "UNKNOWN"}
        return out

    def status(ID):
        res = requests.get(f'https://a3m.mmseqs.com/ticket/{ID}')
        try:
            out = res.json()
        except ValueError:
            out = {"status": "UNKNOWN"}
        return out

    def download(ID, path):
        res = requests.get(f'https://a3m.mmseqs.com/result/download/{ID}')
        with open(path, "wb") as out: out.write(res.content)

    # process input x
    seqs = [x] if isinstance(x, str) else x

    # compatibility to old option
    if filter is not None:
        use_filter = filter

    # setup mode
    if use_filter:
        mode = "env" if use_env else "all"
    else:
        mode = "env-nofilter" if use_env else "nofilter"

    # define path
    path = f"{prefix}_{mode}"
    if not os.path.isdir(path): os.mkdir(path)

    # call mmseqs2 api
    tar_gz_file = f'{path}/out.tar.gz'
    N, REDO = 101, True
    if not os.path.isfile(tar_gz_file):
        TIME_ESTIMATE = 150 * len(seqs)
        with tqdm.notebook.tqdm(total=TIME_ESTIMATE, bar_format=TQDM_BAR_FORMAT) as pbar:
            while REDO:
                pbar.set_description("SUBMIT")

                # Resubmit job until it goes through
                out = submit(seqs, mode, N)
                while out["status"] in ["UNKNOWN", "RATELIMIT"]:
                    # resubmit
                    time.sleep(5 + random.randint(0, 5))
                    out = submit(seqs, mode, N)

                # wait for job to finish
                ID, TIME = out["id"], 0
                pbar.set_description(out["status"])
                while out["status"] in ["UNKNOWN", "RUNNING", "PENDING"]:
                    t = 5 + random.randint(0, 5)
                    time.sleep(t)
                    out = status(ID)
                    pbar.set_description(out["status"])
                    if out["status"] == "RUNNING":
                        TIME += t
                        pbar.update(n=t)
                    # if TIME > 900 and out["status"] != "COMPLETE":
                    #  # something failed on the server side, need to resubmit
                    #  N += 1
                    #  break

                if out["status"] == "COMPLETE":
                    if TIME < TIME_ESTIMATE:
                        pbar.update(n=(TIME_ESTIMATE - TIME))
                    REDO = False

                if out["status"] == "ERROR":
                    REDO = False
                    raise Exception(
                        f'MMseqs2 API is giving errors. Please confirm your input is a valid protein sequence. If error persists, please try again an hour later.')

            # Download results
            download(ID, tar_gz_file)

    # prep list of a3m files
    a3m_files = [f"{path}/uniref.a3m"]
    if use_env: a3m_files.append(f"{path}/bfd.mgnify30.metaeuk30.smag30.a3m")

    # extract a3m files
    if not os.path.isfile(a3m_files[0]):
        with tarfile.open(tar_gz_file) as tar_gz:
            tar_gz.extractall(path)

            # templates
    if use_templates:
        templates = {}
        print("seq\tpdb\tcid\tevalue")
        for line in open(f"{path}/pdb70.m8", "r"):
            p = line.rstrip().split()
            M, pdb, qid, e_value = p[0], p[1], p[2], p[10]
            M = int(M)
            if M not in templates: templates[M] = []
            templates[M].append(pdb)
            if len(templates[M]) <= 20:
                print(f"{int(M) - N}\t{pdb}\t{qid}\t{e_value}")

        template_paths = {}
        for k, TMPL in templates.items():
            TMPL_PATH = f"{prefix}_{mode}/templates_{k}"
            if not os.path.isdir(TMPL_PATH):
                os.mkdir(TMPL_PATH)
                TMPL_LINE = ",".join(TMPL[:20])
                os.system(f"curl -s https://a3m-templates.mmseqs.com/template/{TMPL_LINE} | tar xzf - -C {TMPL_PATH}/")
                os.system(f"cp {TMPL_PATH}/pdb70_a3m.ffindex {TMPL_PATH}/pdb70_cs219.ffindex")
                os.system(f"touch {TMPL_PATH}/pdb70_cs219.ffdata")
            template_paths[k] = TMPL_PATH

    # gather a3m lines
    a3m_lines = {}
    for a3m_file in a3m_files:
        update_M, M = True, None
        for line in open(a3m_file, "r"):
            if len(line) > 0:
                if "\x00" in line:
                    line = line.replace("\x00", "")
                    update_M = True
                if line.startswith(">") and update_M:
                    M = int(line[1:].rstrip())
                    update_M = False
                    if M not in a3m_lines: a3m_lines[M] = []
                a3m_lines[M].append(line)

    # return results
    Ms = sorted(list(a3m_lines.keys()))
    a3m_lines = ["".join(a3m_lines[n]) for n in Ms]
    if use_templates:
        template_paths_ = []
        for n in Ms:
            if n not in template_paths:
                template_paths_.append(None)
                print(f"{n - N}\tno_templates_found")
            else:
                template_paths_.append(template_paths[n])
        template_paths = template_paths_

    if isinstance(x, str):
        return (a3m_lines[0], template_paths[0]) if use_templates else a3m_lines[0]
    else:
        return (a3m_lines, template_paths) if use_templates else a3m_lines
