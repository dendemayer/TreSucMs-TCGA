::: {.related role="navigation" aria-label="related navigation"}
### Navigation

-   [index](genindex.html "General Index")
-   [modules](py-modindex.html "Python Module Index") \|
-   [DESeq2_pipeline 1.0 documentation](#) »
-   [Documentation: DESeq2_pipeline for TCGA]()
:::

::: {.document}
::: {.documentwrapper}
::: {.bodywrapper}
::: {.body role="main"}
::: {#documentation-deseq2-pipeline-for-tcga .section}
# Documentation: DESeq2_pipeline for TCGA[¶](#documentation-deseq2-pipeline-for-tcga "Permalink to this headline"){.headerlink}

"DESeq2_pipeline" a tool to choose, harvest and analyse expression data
of the TCGA-projects with help of the DESeq2 R package

build and activate the provided conda env:

::: {.highlight-bash .notranslate}
::: {.highlight}
    $ conda env create -f deseq_env.yaml
    $ conda activate deseq_pipeline
:::
:::

call the script without any options to enter the interctive mode and set
each option step by step:

::: {.highlight-bash .notranslate}
::: {.highlight}
    $ python main_deseq.py
:::
:::

print help page:

::: {.highlight-bash .notranslate}
::: {.highlight}
    $ python main_deseq.py --help
:::
:::

Usage: main_deseq.py \[OPTIONS\]

Options:

:   

    [[-D]{.option}, [\--download_data]{.option}]{.kbd}

    :   perform download and merging steps, without the DESeq2 analysis

    [[-A]{.option}, [\--analyse_data]{.option}]{.kbd}

    :   perform the DESeq2 analysis

    [[-o]{.option}, [\--out_path `TEXT`{.variable}]{.option}]{.kbd}

    :   The Path, where the results are saved to \[default: CWD\]

    [[-s]{.option}, [\--script_path `TEXT`{.variable}]{.option}]{.kbd}

    :   The Path, where the DESeq2_pipeline script is located, if not
        supplied, the current working directory is assumed

    [[-d]{.option}, [\--drugs `TEXT`{.variable}]{.option}]{.kbd}

    :   drug(s), like: -d drug1 -d drug2 or drugcombination(s), like: -d
        drug1,drug2

    [[-p]{.option}, [\--project `TEXT`{.variable}]{.option}]{.kbd}

    :   Project, that shall be applied, choose between TCGA- CESC,
        TCGA-HNSC, TCGA-LUSC, TCGA-ESCA, or combinations out of them
        like -p project1 -p project2

    [[-f]{.option}, [\--function `INTEGER`{.variable}]{.option}]{.kbd}

    :   running functions seperately lookup the manpage for detailed
        description of each function

    [[-t]{.option}, [\--treshold `INTEGER`{.variable}]{.option}]{.kbd}

    :   treshold for a higher division of the kaplan meier plots
        \[default: 0\]

    [[\--help]{.option}]{.kbd}

    :   Show this message and exit.

```{=html}
<!-- -->
```

`main_deseq.`{.sig-prename .descclassname}`call_with_options`{.sig-name .descname}[(]{.sig-paren}*[[\*]{.pre}]{.o}[[args]{.pre}]{.n}*, *[[\*\*]{.pre}]{.o}[[kwargs]{.pre}]{.n}*[)]{.sig-paren}[¶](#main_deseq.call_with_options "Permalink to this definition"){.headerlink}

:   "DESeq2_pipeline" a tool to choose, harvest and analyse expression
    data of the TCGA-projects with help of the DESeq2 R package

    within this function, every parameter is set needed for the analysis

    Param

    :   out_path: path for DESeq2 pipeline outputs

    Type

    :   out_path: str

    Param

    :   script_path: path to the DESeq2_pipeline repo

    Type

    :   script_path: str

    Param

    :   function: apply, if a specific function should be executet
        solely (also multiple functions possible)

    Type

    :   function: int

    Param

    :   drugs: applied drug(s)

    Type

    :   drugs: list of str

    Param

    :   project: list of projects choosen

    Type

    :   project: list of str

    Param

    :   download_data: bool flag, whether raw data needs to be
        downloaded

    Type

    :   download_data: bool

    Param

    :   analyse_data: bool flag, whether the deseq analyses shall be
        started

    Type

    :   download_data: bool

    Param

    :   treshold: parameter for the lifeline plots helping for the
        classification of expression data

    Type

    :   treshold: int

[]{#module-choose_therapy .target}

`choose_therapy.`{.sig-prename .descclassname}`Choose_drugs`{.sig-name .descname}[(]{.sig-paren}*[[SCRIP_PATH]{.pre}]{.n}*, *[[PROJECTS]{.pre}]{.n}*[)]{.sig-paren}[¶](#choose_therapy.Choose_drugs "Permalink to this definition"){.headerlink}

:   

    Param

    :   SCRIPT_PATH: path to the DESeq2_pipeline repo

    Type

    :   SCRIPT_PATH: str

    Param

    :   PROJECTS: list of projects choosen

    Type

    :   PROJECTS: list of str

    interactively requesting the drugs which shall be applied to the
    deseq approach

```{=html}
<!-- -->
```

`choose_therapy.`{.sig-prename .descclassname}`Choose_path_and_option`{.sig-name .descname}[(]{.sig-paren}*[[OUTPUT_PATH]{.pre}]{.n}*, *[[PROJECTS]{.pre}]{.n}*, *[[DRUGS]{.pre}]{.n}*, *[[function]{.pre}]{.n}*, *[[SCRIPT_PATH]{.pre}]{.n}*[)]{.sig-paren}[¶](#choose_therapy.Choose_path_and_option "Permalink to this definition"){.headerlink}

:   

    Param

    :   OUTPUT_PATH: path for DESeq2 pipeline outputs

    Type

    :   OUTPUT_PATH: str

    Param

    :   PROJECTS: list of projects choosen

    Type

    :   PROJECTS: list of str

    Param

    :   DRUGS: applied drug(s)

    Type

    :   DRUGS: list of str

    Param

    :   function: applied functions

    Type

    :   function: int

    Param

    :   SCRIPT_PATH: path to the DESeq2_pipeline repo

    Type

    :   SCRIPT_PATH: str

    interactively requesting whether download steps, analysis or both
    should be performed

```{=html}
<!-- -->
```

`choose_therapy.`{.sig-prename .descclassname}`Choose_project`{.sig-name .descname}[(]{.sig-paren}[)]{.sig-paren}[¶](#choose_therapy.Choose_project "Permalink to this definition"){.headerlink}

:   interactively requesting the Projects that shall be applied to the
    deseq approach

[]{#module-create_matrix_new .target}

`create_matrix_new.`{.sig-prename .descclassname}`correct_drugs`{.sig-name .descname}[(]{.sig-paren}*[[PROJECT]{.pre}]{.n}*, *[[OUTPUT_PATH]{.pre}]{.n}*, *[[logger]{.pre}]{.n}*[)]{.sig-paren}[¶](#create_matrix_new.correct_drugs "Permalink to this definition"){.headerlink}

:   

    Param

    :   PROJECT: list of projects choosen

    Type

    :   PROJECT: list of str

    Param

    :   OUTPUT_PATH: path for DESeq2 pipeline outputs

    Type

    :   OUTPUT_PATH: str

    Param

    :   logger: the adjustet logger with the right filehandler

    Type

    :   logger: logging instance

    creating table: 'OUTPUT_PATH/DF_3t_both_with_DRUG_combi.tsv',
    drugcombination in field 'drugnames' as ordered set, comma
    seperated. header includes: UUID case_id gender vital_status
    drugnames survivaltime years_to_last_follow_up year_of_birth
    year_of_death age_at_diagnosis PROJECT

    ::: {.highlight-bash .notranslate}
    ::: {.highlight}
        # for executing this step via terminal, issue -f 6 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 6
    :::
    :::

```{=html}
<!-- -->
```

`create_matrix_new.`{.sig-prename .descclassname}`create_snake_config`{.sig-name .descname}[(]{.sig-paren}*[[OUTPUT_PATH]{.pre}]{.n}*, *[[PROJECT_title]{.pre}]{.n}*, *[[DRUGS_title]{.pre}]{.n}*, *[[project_list]{.pre}]{.n}*, *[[DRUGS]{.pre}]{.n}*, *[[SCRIPT_PATH]{.pre}]{.n}*[)]{.sig-paren}[¶](#create_matrix_new.create_snake_config "Permalink to this definition"){.headerlink}

:   

    Param

    :   OUTPUT_PATH: path for DESeq2 pipeline outputs

    Type

    :   OUTPUT_PATH: str

    Param

    :   PROJECT_title: concatenated str of multiple projects

    Type

    :   PROJECT_title: str

    Param

    :   DRUGS_title: concatenated str of multiple drugs

    Type

    :   DRUGS_title: str

    Param

    :   project_list: list of applied projects

    Type

    :   project_list: list of str

    Param

    :   DRUGS: applied drug(s)

    Type

    :   DRUGS: list of str

    Param

    :   SCRIPT_PATH: path to the DESeq2_pipeline repo

    Type

    :   SCRIPT_PATH: str

    out of the log files in PROJECT_title/DRUGS_title/test_log.log parse
    out all outputfiles of the applied run and create
    PROJECT_title/DRUGS_title/snakemake_config.yaml

    ::: {.highlight-bash .notranslate}
    ::: {.highlight}
        # for executing this step via terminal, issue -f 15 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 15
    :::
    :::

```{=html}
<!-- -->
```

`create_matrix_new.`{.sig-prename .descclassname}`create_statistics_from_DESeq2_tables`{.sig-name .descname}[(]{.sig-paren}*[[OUTPUT_PATH]{.pre}]{.n}*, *[[DRUGS]{.pre}]{.n}*, *[[SCRIPT_PATH]{.pre}]{.n}*, *[[PROJECT]{.pre}]{.n}*, *[[logger]{.pre}]{.n}*[)]{.sig-paren}[¶](#create_matrix_new.create_statistics_from_DESeq2_tables "Permalink to this definition"){.headerlink}

:   

    Param

    :   OUTPUT_PATH: path for DESeq2 pipeline outputs

    Type

    :   OUTPUT_PATH: str

    Param

    :   DRUGS: applied drug(s)

    Type

    :   DRUGS: list of str

    Param

    :   SCRIPT_PATH: path to the DESeq2_pipeline repo

    Type

    :   SCRIPT_PATH: str

    Param

    :   PROJECT: list of projects choosen

    Type

    :   PROJECT: list of str

    Param

    :   logger: the adjustet logger with the right filehandler

    Type

    :   logger: logging instance

    adding some statistics to the result output tables from DESeq2 table
    is saved in
    OUTPUT_PATH/PROJECT_title/DRUGS_title/DESeq2_out\*\*\*/results_statistics.tsv

    ::: {.highlight-bash .notranslate}
    ::: {.highlight}
        # for executing this step via terminal, issue -f 8 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 8

        # when choosing multiple projects, call:
        $ python main_deseq.py -p TCGA-CESC  -p TCGA-HNSC -d cisplatin -f 10
    :::
    :::

```{=html}
<!-- -->
```

`create_matrix_new.`{.sig-prename .descclassname}`create_summary_table`{.sig-name .descname}[(]{.sig-paren}*[[OUTPUT_PATH]{.pre}]{.n}*, *[[PROJECT]{.pre}]{.n}*, *[[logger]{.pre}]{.n}*[)]{.sig-paren}[¶](#create_matrix_new.create_summary_table "Permalink to this definition"){.headerlink}

:   

    Param

    :   OUTPUT_PATH: path for DESeq2 pipeline outputs

    Type

    :   OUTPUT_PATH: str

    Param

    :   PROJECT: list of projects choosen

    Type

    :   PROJECT: list of str

    Param

    :   logger: the adjustet logger with the right filehandler

    Type

    :   logger: logging instance

    get the raw data out of the gzip compressed files and merge them
    together, save complete merged table in
    OUTPUT_PATH/PROJECT/summary_DF.tsv

    []{#summary-df-example}

    +------------------+---------------------------------------------------+
    | gene identifier  | case_id with counts                               |
    +==================+===================================================+
    | genes            | 6ff12a54-10da-4941-bfea-7b66e19b4be9 ...          |
    +------------------+---------------------------------------------------+
    | ENSG00000000003  | 3423 ...                                          |
    +------------------+---------------------------------------------------+
    | ENSG00000000005  | 0 ...                                             |
    +------------------+---------------------------------------------------+

    : [example structure for
    summary_DF.tsv:]{.caption-text}[¶](#id1 "Permalink to this table"){.headerlink}

    ::: {.highlight-bash .notranslate}
    ::: {.highlight}
        # for executing this step via terminal, issue -f 5 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 5
    :::
    :::

```{=html}
<!-- -->
```

`create_matrix_new.`{.sig-prename .descclassname}`download_clinical_tables`{.sig-name .descname}[(]{.sig-paren}*[[UUID]{.pre}]{.n}*, *[[PROJECT]{.pre}]{.n}*, *[[OUTPUT_PATH]{.pre}]{.n}*, *[[logger]{.pre}]{.n}*[)]{.sig-paren}[¶](#create_matrix_new.download_clinical_tables "Permalink to this definition"){.headerlink}

:   

    Param

    :   UUID: unique file identifier of the meta table

    Type

    :   UUID: str

    Param

    :   PROJECT: list of projects choosen

    Type

    :   PROJECT: list of str

    Param

    :   OUTPUT_PATH: path for DESeq2 pipeline outputs

    Type

    :   OUTPUT_PATH: str

    Param

    :   logger: the adjustet logger with the right filehandler

    Type

    :   logger: logging instance

    with the UUID the clinical tables will be downloaded in the
    OUTPUT_PATH/PROJECT:

    -   nationwidechildrens.org_clinical_patient\_\*\*\*\*.txt

    -   nationwidechildrens.org_clinical_drug\_\*\*\*\*.txt

    ::: {.highlight-bash .notranslate}
    ::: {.highlight}
        # for executing this step via terminal, issue -f 2 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 2
    :::
    :::

```{=html}
<!-- -->
```

`create_matrix_new.`{.sig-prename .descclassname}`download_data_files`{.sig-name .descname}[(]{.sig-paren}*[[PROJECT]{.pre}]{.n}*, *[[FILE_TYPE]{.pre}]{.n}*, *[[OUTPUT_PATH]{.pre}]{.n}*, *[[logger]{.pre}]{.n}*[)]{.sig-paren}[¶](#create_matrix_new.download_data_files "Permalink to this definition"){.headerlink}

:   

    Param

    :   PROJECT: list of projects choosen

    Type

    :   PROJECT: list of str

    Param

    :   FILE_TYPE: type of raw data to download from TCGA

    Type

    :   FILE_TYPE: str

    Param

    :   OUTPUT_PATH: path for DESeq2 pipeline outputs

    Type

    :   OUTPUT_PATH: str

    Param

    :   logger: the adjustet logger with the right filehandler

    Type

    :   logger: logging instance

    dowloading all data files and the belonging manifest with UUID to
    each file creating a subdir in the OUTPUT_PATH/PROJECT/ folder named
    "{PROJECT}\_data_files" the belonging MANIFEST.txt is saved in
    OUTPUT_PATH/PROJECT

    ::: {.highlight-bash .notranslate}
    ::: {.highlight}
        # for executing this step via terminal, issue -f 3 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 3
    :::
    :::

```{=html}
<!-- -->
```

`create_matrix_new.`{.sig-prename .descclassname}`merging_meta_infos`{.sig-name .descname}[(]{.sig-paren}*[[OUTPUT_PATH]{.pre}]{.n}*, *[[PROJECT]{.pre}]{.n}*, *[[logger]{.pre}]{.n}*[)]{.sig-paren}[¶](#create_matrix_new.merging_meta_infos "Permalink to this definition"){.headerlink}

:   

    Param

    :   OUTPUT_PATH: path for DESeq2 pipeline outputs

    Type

    :   OUTPUT_PATH: str

    Param

    :   PROJECT: list of projects choosen

    Type

    :   PROJECT: list of str

    Param

    :   logger: the adjustet logger with the right filehandler

    Type

    :   logger: logging instance

    merging the infos out of the 3 tables
    together(nationewidechildrens.org.., manifest, metainfo) creating 3
    new tables, one with innerjoin(both), two with left and right outer
    join to save those files where information is missing: \*
    DF_3t_left: no information about therapeutic_agents \* DF_3t_right:
    missing filename (case_id present)

    -   DF_3t_both holds infos like UUID, filename, md5, size, state,
        case_id, gender, vital_status, year_of_birth, year_of_death and
        some more..

    every table saved in OUTPUT_PATH/PROJECT/

    ::: {.highlight-bash .notranslate}
    ::: {.highlight}
        # for executing this step via terminal, issue -f 4 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 4
    :::
    :::

```{=html}
<!-- -->
```

`create_matrix_new.`{.sig-prename .descclassname}`meta_filter`{.sig-name .descname}[(]{.sig-paren}*[[PROJECT]{.pre}]{.n}*, *[[FILE_TYPE]{.pre}]{.n}*, *[[OUTPUT_PATH]{.pre}]{.n}*, *[[logger]{.pre}]{.n}*[)]{.sig-paren}[¶](#create_matrix_new.meta_filter "Permalink to this definition"){.headerlink}

:   

    Param

    :   PROJECT: list of projects choosen

    Type

    :   PROJECT: list of str

    Param

    :   FILE_TYPE: type of raw data to download from TCGA

    Type

    :   FILE_TYPE: str

    Param

    :   OUTPUT_PATH: path for DESeq2 pipeline outputs

    Type

    :   OUTPUT_PATH: str

    Param

    :   logger: the adjustet logger with the right filehandler

    Type

    :   logger: logging instance

    creating the meta_info.dat file in your OUTPUT_PATH/PROJECT
    directory

    ::: {.highlight-bash .notranslate}
    ::: {.highlight}
        # for executing this step via terminal, issue -f 1 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 1
    :::
    :::

```{=html}
<!-- -->
```

`create_matrix_new.`{.sig-prename .descclassname}`provide_DESeq2_table`{.sig-name .descname}[(]{.sig-paren}*[[PROJECT]{.pre}]{.n}*, *[[OUTPUT_PATH]{.pre}]{.n}*, *[[DRUGS]{.pre}]{.n}*, *[[SCRIPT_PATH]{.pre}]{.n}*, *[[logger]{.pre}]{.n}*[)]{.sig-paren}[¶](#create_matrix_new.provide_DESeq2_table "Permalink to this definition"){.headerlink}

:   

    Param

    :   PROJECT: list of projects choosen

    Type

    :   PROJECT: list of str

    Param

    :   OUTPUT_PATH: path for DESeq2 pipeline outputs

    Type

    :   OUTPUT_PATH: str

    Param

    :   DRUGS: applied drug(s)

    Type

    :   DRUGS: list of str

    Param

    :   SCRIPT_PATH: path to the DESeq2_pipeline repo

    Type

    :   SCRIPT_PATH: str

    Param

    :   logger: the adjustet logger with the right filehandler

    Type

    :   logger: logging instance

    filtering case_id according to DRUG query providing new table for
    DESeq2 analysis dependend on distinguishable factors of the tables
    provided, a singlefactorial (at least differences in vital state) or
    a mutlifactorial run in DESeq2 is performed (gender, therapy or
    project )

    ::: {.highlight-bash .notranslate}
    ::: {.highlight}
        # for executing this step via terminal, issue -f 7 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 7

        # when choosing multiple projects, call:
        $ python main_deseq.py -p TCGA-CESC  -p TCGA-HNSC -d cisplatin -f 9
    :::
    :::

```{=html}
<!-- -->
```

`create_matrix_new.`{.sig-prename .descclassname}`set_logger`{.sig-name .descname}[(]{.sig-paren}*[[OUTPUT_PATH]{.pre}]{.n}*, *[[PROJECT_title]{.pre}]{.n}*, *[[DRUGS_title]{.pre}]{.n}*[)]{.sig-paren}[¶](#create_matrix_new.set_logger "Permalink to this definition"){.headerlink}

:   

    Param

    :   OUTPUT_PATH: path for DESeq2 pipeline outputs

    Type

    :   OUTPUT_PATH: str

    Param

    :   PROJECT_title: merged project title out of multiple projects

    Type

    :   PROJECT_title: str

    Param

    :   DRUGS_title: merged drug title out of multiple drugs

    Type

    :   DRUGS_title: str

    every paths and options are set, configure here the logfiles, with
    which the snakemake config files are going to be created we create 2
    loggers, in case just one project is applied the logs are written in
    PROJECT/DRUGS_title/test_log.log in case multi project is applied,
    the logs are written in PROJECT_title/DRUGS_title/test_log.log with
    that, it is clear which config file shall be createt out of the
    logfiles present in one outputpath (the drugs path must therefore be
    created from the first fct, to write the log file also, the dir of
    the logfile must be logged, s.t. snakemake knows where the input
    file for the final snakemake configuratioin file is located

[]{#module-lifeline_summary_test_2 .target}

*[class]{.pre}* `lifeline_summary_test_2.`{.sig-prename .descclassname}`Lifeline_plot`{.sig-name .descname}[(]{.sig-paren}*[[OUTPUT_PATH]{.pre}]{.n}*, *[[PROJECT_title]{.pre}]{.n}*, *[[DRUGS_title]{.pre}]{.n}*, *[[path_to_result]{.pre}]{.n}*, *[[prefix]{.pre}]{.n}*, *[[treshold]{.pre}]{.n}*, *[[multi_project]{.pre}]{.n}*[)]{.sig-paren}[¶](#lifeline_summary_test_2.Lifeline_plot "Permalink to this definition"){.headerlink}

:   returning the ENSG on which the gene specific KaplanMeier plots are
    created

```{=html}
<!-- -->
```

`lifeline_summary_test_2.`{.sig-prename .descclassname}`lifelines_ENSG`{.sig-name .descname}[(]{.sig-paren}*[[OUTPUT_PATH]{.pre}]{.n}*, *[[DRUGS]{.pre}]{.n}*, *[[PROJECT_DRUG_UUID]{.pre}]{.n}*, *[[treshold]{.pre}]{.n}*[)]{.sig-paren}[¶](#lifeline_summary_test_2.lifelines_ENSG "Permalink to this definition"){.headerlink}

:   

    Param

    :   OUTPUT_PATH: path for DESeq2 pipeline outputs

    Type

    :   OUTPUT_PATH: str

    Param

    :   DRUGS: applied drug(s)

    Type

    :   DRUGS: list of strings

    Param

    :   PROJECT_DRUG_UUID: hash holding a project to the UUID of the
        belonging drugtable

    Type

    :   PROJECT_DRUG_UUID: dict

    Param

    :   treshold: parameter for the lifeline plots helping for the
        classification of expression data

    Type

    :   treshold: int

    ::: {.line}
    script creates the plots:
    :::

    ::: {.line}
    lifelines_cumulative_density.svg,
    :::

    ::: {.line}
    lifelines_multiple_groups.svg,
    :::

    ::: {.line}
    lifelines_parametric_models_2.svg,
    :::

    ::: {.line}
    lifelines_survival_fct.svg,
    :::

    ::: {.line}
    lifelines_table.tsv in OUTPUT_PATH
    :::

    -   with the different expression ENSGs (normalized counts given
        from deseq, median from them and separated in UP and DOWN)
        groups plottet with lifeline

    **needed**:

    -   ALL_PROJECTS_summary_dead_alive_reduced_INFO.tsv as

    **count_DF_MI** (just the info for building the multiindex for
    **count_DF**):

    []{#multiindex-table}

    +------------------+---------------------------------------------------+
    | variable         | value                                             |
    +==================+===================================================+
    | vital_status     | alive                                             |
    +------------------+---------------------------------------------------+
    | gender           | female                                            |
    +------------------+---------------------------------------------------+
    | case_id          | 6ff12a54-10da-4941-bfea-7b66e19b4be9              |
    +------------------+---------------------------------------------------+
    | PROJECT          | TCGA-CESC                                         |
    +------------------+---------------------------------------------------+

    : [example multiindex
    table]{.caption-text}[¶](#id2 "Permalink to this table"){.headerlink}

    -   DESeq2_MF_normalized_counts_reduced.tsv as

    []{#ensg-table}

    +-----------------------------------+-----------------------------------+
    | variable                          | value                             |
    +===================================+===================================+
    | ENSG00000000003                   | 4109.0073147311                   |
    +-----------------------------------+-----------------------------------+
    | ENSG00000000005                   | 0                                 |
    +-----------------------------------+-----------------------------------+

    : [ENSG table, **count_DF** (normalized counts with help of
    DESeq2)]{.caption-text}[¶](#id3 "Permalink to this table"){.headerlink}

    -   OUTPUT_PATH/DRUGS/DESeq2_MF_results_reduced.tsv as **DF_res**

    creates the **ENSG_list** on the basis of the DESeq_results in
    dependence of the resulttables, we get the 10 highest and 10 lowest
    logfoldchange

    --\> here we sort log2fold change wise and merge then the most
    different INCREASING and DECREASING cases with the
    TCGA-\*/DF_3t_both_with_DRUG_combi.tsv, in this table we have the
    fields:

    UUID case_id gender, vital_status, drugnames, survivaltime,
    years_to_last_follow_up, year_of_birth, year_of_death,
    age_at_diagnosis, PROJECT,

    with it we create the table for lifeline

    []{#lifeline-table}

    +------+-----------------+-------+------------------------------------+
    | i    | T               | E     | case_id                            |
    | ndex |                 |       |                                    |
    +======+=================+=======+====================================+
    | 0    | 5.              | True  | 9f                                 |
    |      | 197260273972604 |       | fa79fa-d2d8-48e1-8fd6-4b020ecf357c |
    +------+-----------------+-------+------------------------------------+
    | 1    | 0.8             | False | 0d                                 |
    |      | 602739726027397 |       | e19185-3517-4e30-925b-7eb1f5079ec2 |
    +------+-----------------+-------+------------------------------------+

    : [example table for lifeline
    input]{.caption-text}[¶](#id4 "Permalink to this table"){.headerlink}

    -   the up and down separation depends on the median of the
        normalized count matrix

    -   setting the threshold based on the median of the logfoldchange,
        delete out 10 % around it if -f 10 is applied

    ::: {.highlight-bash .notranslate}
    ::: {.highlight}
        # for executing this step via terminal, issue -f 11 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 11

        # you can try out different tresholds, new directorys are
        # created therefore:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 11 -t 10

        # if you want to include the new created outputs in your snakemake
        # config file, call function 15 (calls create_snake_config()):
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 15
    :::
    :::

[]{#module-walk_all_drug_frequency .target}

`walk_all_drug_frequency.`{.sig-prename .descclassname}`drug_frequency`{.sig-name .descname}[(]{.sig-paren}*[[OUTPUT_PATH]{.pre}]{.n}*, *[[DRUGS]{.pre}]{.n}*, *[[SCRIPT_PATH]{.pre}]{.n}*, *[[PROJECT_DRUG_UUID]{.pre}]{.n}*[)]{.sig-paren}[¶](#walk_all_drug_frequency.drug_frequency "Permalink to this definition"){.headerlink}

:   

    Param

    :   OUTPUT_PATH: path for DESeq2 pipeline outputs

    Type

    :   OUTPUT_PATH: str

    Param

    :   DRUGS: applied drug(s)

    Type

    :   DRUGS: list of str

    Param

    :   SCRIPT_PATH: path to the DESeq2_pipeline repo

    Type

    :   SCRIPT_PATH: str

    Param

    :   PROJECT_DRUG_UUID: hash holding a project to the UUID of the
        belonging drugtable

    Type

    :   PROJECT_DRUG_UUID: dict

    create an overview of all available drugs, of the applied projects
    in DESeq2 output dir

    ::: {.highlight-bash .notranslate}
    ::: {.highlight}
        # for executing this step via terminal, issue -f 12 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 12
    :::
    :::

```{=html}
<!-- -->
```

`walk_all_drug_frequency.`{.sig-prename .descclassname}`drug_frequency_all_single_projects`{.sig-name .descname}[(]{.sig-paren}*[[OUTPUT_PATH]{.pre}]{.n}*, *[[DRUGS]{.pre}]{.n}*, *[[SCRIPT_PATH]{.pre}]{.n}*, *[[PROJECT_DRUG_UUID]{.pre}]{.n}*[)]{.sig-paren}[¶](#walk_all_drug_frequency.drug_frequency_all_single_projects "Permalink to this definition"){.headerlink}

:   

    Param

    :   OUTPUT_PATH: path for DESeq2 pipeline outputs

    Type

    :   OUTPUT_PATH: str

    Param

    :   DRUGS: applied drug(s)

    Type

    :   DRUGS: list of str

    Param

    :   SCRIPT_PATH: path to the DESeq2_pipeline repo

    Type

    :   SCRIPT_PATH: str

    Param

    :   PROJECT_DRUG_UUID: hash holding a project to the UUID of the
        belonging drugtable

    Type

    :   PROJECT_DRUG_UUID: dict

    take every drug frequency out of the single project folders,
    therefore walk in the OUTPUT_PATH/TCGA-\[2..4\] folders...

    ::: {.highlight-bash .notranslate}
    ::: {.highlight}
        # for executing this step via terminal, issue -f 14 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 14
    :::
    :::

[]{#module-create_report .target}

`create_report.`{.sig-prename .descclassname}`create_report_pdf`{.sig-name .descname}[(]{.sig-paren}*[[OUTPUT_PATH]{.pre}]{.n}*, *[[DRUGS]{.pre}]{.n}*, *[[SCRIPT_PATH]{.pre}]{.n}*, *[[PROJECT_DRUG_UUID]{.pre}]{.n}*, *[[treshold]{.pre}]{.n}*[)]{.sig-paren}[¶](#create_report.create_report_pdf "Permalink to this definition"){.headerlink}

:   

    Param

    :   OUTPUT_PATH: path for DESeq2 pipeline outputs

    Type

    :   OUTPUT_PATH: str

    Param

    :   DRUGS: applied drug(s)

    Type

    :   DRUGS: list of str

    Param

    :   SCRIPT_PATH: path to the DESeq2_pipeline repo

    Type

    :   SCRIPT_PATH: str

    Param

    :   PROJECT_DRUG_UUID: hash holding a project to the UUID of the
        belonging drugtable

    Type

    :   PROJECT_DRUG_UUID: dict

    Param

    :   treshold: parameter for the lifeline plots helping for the
        classification of expression data

    Type

    :   treshold: int

    creating a report file with all visual outputs created in an
    analysis, saved at OUTPUT_PATH/PROJECT_title/DRUGS_title/REPORT.pdf

    ::: {.highlight-bash .notranslate}
    ::: {.highlight}
        # for executing this step via terminal, issue -f 13 to your call,
        # example:
        $ python main_deseq.py -p TCGA-CESC -d cisplatin -f 13
    :::
    :::

::: {.toctree-wrapper .compound}
:::

::: {#tables .section}
## Tables[¶](#tables "Permalink to this headline"){.headerlink}

-   [[example structure for summary_DF.tsv:]{.std
    .std-ref}](#summary-df-example){.reference .internal}

-   [[example multiindex table]{.std
    .std-ref}](#multiindex-table){.reference .internal}

-   [[ENSG table, count_DF (normalized counts with help of DESeq2)]{.std
    .std-ref}](#ensg-table){.reference .internal}

-   [[example table for lifeline input]{.std
    .std-ref}](#lifeline-table){.reference .internal}
:::
:::

::: {.clearer}
:::
:::
:::
:::

::: {.sphinxsidebar role="navigation" aria-label="main navigation"}
::: {.sphinxsidebarwrapper}
### [Table of Contents](#)

-   [Documentation: DESeq2_pipeline for TCGA](#){.reference .internal}
    -   [Tables](#tables){.reference .internal}

::: {role="note" aria-label="source link"}
### This Page

-   [Show Source](_sources/index.rst.txt)
:::

::: {#searchbox style="display: none" role="search"}
### Quick search {#searchlabel}

::: {.searchformwrapper}
:::
:::
:::
:::

::: {.clearer}
:::
:::

::: {.related role="navigation" aria-label="related navigation"}
### Navigation

-   [index](genindex.html "General Index")
-   [modules](py-modindex.html "Python Module Index") \|
-   [DESeq2_pipeline 1.0 documentation](#) »
-   [Documentation: DESeq2_pipeline for TCGA]()
:::

::: {.footer role="contentinfo"}
© Copyright 2021, Gabor Balogh. Created using
[Sphinx](https://www.sphinx-doc.org/) 3.5.2.
:::
