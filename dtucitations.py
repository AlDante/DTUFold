
#@title Procedure write_citations()

def write_citations(jobname:str, use_msa:bool=False, use_custom_msa:bool=False, use_env:bool=False, use_templates:bool=False, use_amber:bool=False):

  citations = {
    "Ovchinnikov2021":  """@software{Ovchinnikov2021,
  author = {Ovchinnikov, Sergey and Steinegger, Martin and Mirdita, Milot},
  title = {{ColabFold - Making Protein folding accessible to all via Google Colab}},
  year = {2021},
  publisher = {Zenodo},
  version = {v1.0-alpha},
  doi = {10.5281/zenodo.5123297},
  url = {https://doi.org/10.5281/zenodo.5123297},
  comment = {The AlphaFold notebook}
  }""",
    "LevyKarin2020": """@article{LevyKarin2020,
  author = {{Levy Karin}, Eli and Mirdita, Milot and S{\"{o}}ding, Johannes},
  doi = {10.1186/s40168-020-00808-x},
  journal = {Microbiome},
  number = {1},
  title = {{MetaEuk—sensitive, high-throughput gene discovery, and annotation for large-scale eukaryotic metagenomics}},
  volume = {8},
  year = {2020},
  comment = {MetaEuk database}
  }""",
    "Delmont2020": """@article{Delmont2020,
  author = {Delmont, Tom O. and Gaia, Morgan and Hinsinger, Damien D. and Fremont, Paul and Guerra, Antonio Fernandez and Eren, A. Murat and Vanni, Chiara and Kourlaiev, Artem and D'Agata, Leo and Clayssen, Quentin and Villar, Emilie and Labadie, Karine and Cruaud, Corinne and Poulain, Julie and da Silva, Corinne and Wessner, Marc and Noel, Benjamin and Aury, Jean Marc and de Vargas, Colomban and Bowler, Chris and Karsenti, Eric and Pelletier, Eric and Wincker, Patrick and Jaillon, Olivier and Sunagawa, Shinichi and Acinas, Silvia G. and Bork, Peer and Karsenti, Eric and Bowler, Chris and Sardet, Christian and Stemmann, Lars and de Vargas, Colomban and Wincker, Patrick and Lescot, Magali and Babin, Marcel and Gorsky, Gabriel and Grimsley, Nigel and Guidi, Lionel and Hingamp, Pascal and Jaillon, Olivier and Kandels, Stefanie and Iudicone, Daniele and Ogata, Hiroyuki and Pesant, St{\'{e}}phane and Sullivan, Matthew B. and Not, Fabrice and Karp-Boss, Lee and Boss, Emmanuel and Cochrane, Guy and Follows, Michael and Poulton, Nicole and Raes, Jeroen and Sieracki, Mike and Speich, Sabrina},
  journal = {bioRxiv},
  title = {{Functional repertoire convergence of distantly related eukaryotic plankton lineages revealed by genome-resolved metagenomics}},
  year = {2020},
  comment = {SMAG database}
  }""",
    "Mitchell2019": """@article{Mitchell2019,
  author = {Mitchell, Alex L and Almeida, Alexandre and Beracochea, Martin and Boland, Miguel and Burgin, Josephine and Cochrane, Guy and Crusoe, Michael R and Kale, Varsha and Potter, Simon C and Richardson, Lorna J and Sakharova, Ekaterina and Scheremetjew, Maxim and Korobeynikov, Anton and Shlemov, Alex and Kunyavskaya, Olga and Lapidus, Alla and Finn, Robert D},
  doi = {10.1093/nar/gkz1035},
  journal = {Nucleic Acids Res.},
  title = {{MGnify: the microbiome analysis resource in 2020}},
  year = {2019},
  comment = {MGnify database}
  }""",
    "Eastman2017": """@article{Eastman2017,
  author = {Eastman, Peter and Swails, Jason and Chodera, John D. and McGibbon, Robert T. and Zhao, Yutong and Beauchamp, Kyle A. and Wang, Lee-Ping and Simmonett, Andrew C. and Harrigan, Matthew P. and Stern, Chaya D. and Wiewiora, Rafal P. and Brooks, Bernard R. and Pande, Vijay S.},
  doi = {10.1371/journal.pcbi.1005659},
  journal = {PLOS Comput. Biol.},
  number = {7},
  title = {{OpenMM 7: Rapid development of high performance algorithms for molecular dynamics}},
  volume = {13},
  year = {2017},
  comment = {Amber relaxation}
  }""",
    "Jumper2021": """@article{Jumper2021,
  author = {Jumper, John and Evans, Richard and Pritzel, Alexander and Green, Tim and Figurnov, Michael and Ronneberger, Olaf and Tunyasuvunakool, Kathryn and Bates, Russ and {\v{Z}}{\'{i}}dek, Augustin and Potapenko, Anna and Bridgland, Alex and Meyer, Clemens and Kohl, Simon A. A. and Ballard, Andrew J. and Cowie, Andrew and Romera-Paredes, Bernardino and Nikolov, Stanislav and Jain, Rishub and Adler, Jonas and Back, Trevor and Petersen, Stig and Reiman, David and Clancy, Ellen and Zielinski, Michal and Steinegger, Martin and Pacholska, Michalina and Berghammer, Tamas and Bodenstein, Sebastian and Silver, David and Vinyals, Oriol and Senior, Andrew W. and Kavukcuoglu, Koray and Kohli, Pushmeet and Hassabis, Demis},
  doi = {10.1038/s41586-021-03819-2},
  journal = {Nature},
  pmid = {34265844},
  title = {{Highly accurate protein structure prediction with AlphaFold.}},
  year = {2021},
  comment = {AlphaFold2 + BFD Database}
  }""",
    "Mirdita2019": """@article{Mirdita2019,
  author = {Mirdita, Milot and Steinegger, Martin and S{\"{o}}ding, Johannes},
  doi = {10.1093/bioinformatics/bty1057},
  journal = {Bioinformatics},
  number = {16},
  pages = {2856--2858},
  pmid = {30615063},
  title = {{MMseqs2 desktop and local web server app for fast, interactive sequence searches}},
  volume = {35},
  year = {2019},
  comment = {MMseqs2 search server}
  }""",
    "Steinegger2019": """@article{Steinegger2019,
  author = {Steinegger, Martin and Meier, Markus and Mirdita, Milot and V{\"{o}}hringer, Harald and Haunsberger, Stephan J. and S{\"{o}}ding, Johannes},
  doi = {10.1186/s12859-019-3019-7},
  journal = {BMC Bioinform.},
  number = {1},
  pages = {473},
  pmid = {31521110},
  title = {{HH-suite3 for fast remote homology detection and deep protein annotation}},
  volume = {20},
  year = {2019},
  comment = {PDB70 database}
  }""",
    "Mirdita2017": """@article{Mirdita2017,
  author = {Mirdita, Milot and von den Driesch, Lars and Galiez, Clovis and Martin, Maria J. and S{\"{o}}ding, Johannes and Steinegger, Martin},
  doi = {10.1093/nar/gkw1081},
  journal = {Nucleic Acids Res.},
  number = {D1},
  pages = {D170--D176},
  pmid = {27899574},
  title = {{Uniclust databases of clustered and deeply annotated protein sequences and alignments}},
  volume = {45},
  year = {2017},
  comment = {Uniclust30/UniRef30 database},
  }""",
    "Berman2003": """@misc{Berman2003,
  author = {Berman, Helen and Henrick, Kim and Nakamura, Haruki},
  booktitle = {Nat. Struct. Biol.},
  doi = {10.1038/nsb1203-980},
  number = {12},
  pages = {980},
  pmid = {14634627},
  title = {{Announcing the worldwide Protein Data Bank}},
  volume = {10},
  year = {2003},
  comment = {templates downloaded from wwPDB server}
  }""",
  }

  to_cite = [ "Jumper2021", "Ovchinnikov2021" ]
  if use_msa:       to_cite += ["Mirdita2019"]
  if use_msa:       to_cite += ["Mirdita2017"]
  if use_env:       to_cite += ["Mitchell2019"]
  if use_env:       to_cite += ["Delmont2020"]
  if use_env:       to_cite += ["LevyKarin2020"]
  if use_templates: to_cite += ["Steinegger2019"]
  if use_templates: to_cite += ["Berman2003"]
  if use_amber:     to_cite += ["Eastman2017"]

  with open(f"{jobname}.bibtex", 'w') as writer:
    for i in to_cite:
      writer.write(citations[i])
      writer.write("\n")

  #DEJ: print(f"Found {len(to_cite)} citation{'s' if len(to_cite) > 1 else ''} for tools or databases.")
  if use_custom_msa:
    print("Don't forget to cite your custom MSA generation method.")
