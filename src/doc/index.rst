.. toctree::
   :maxdepth: 2
   :caption: Contents:

Documentation: TreMSuc for TCGA
***************************************

"TreMSuc" a tool to choose, harvest and analyse expression and methylation data
of the TCGA-projects 

Build and activate the provided conda env:
------------------------------------------

.. code-block:: bash

    $ conda env create -f deseq_env.yaml
    $ conda activate deseq_pipeline

call the script without any options to enter the interactive mode and set
each option step by step:

.. code-block:: bash

    $ python main_deseq.py

print help page:

.. code-block:: bash

    $ python main_deseq.py --help


Documentation of modules, classes and functions:
------------------------------------------------


.. click:: shared.modules.main:call_with_options
   :prog: TreMSuc

.. automodule:: tcga_metilene.modules.main_metilene
    :members: 

Short tutorial:
---------------
Performing an example analysis:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The easiest way of applying a run is entering the interactive mode (it is
supposed that you cloned the deseq_pipeline git repository and cd into that
dir):

.. code-block:: bash

  $ python main_deseq.py

With it, every needed parameter is offered for further analyses. First, the
available projects are presented, based on that selection, available drugs
or drug combinations can be chosen.
  
In contrast to that, the parameter needed could be applied via command line.
An example terminal call for the projects TCGA-CESC and TCGA-HNSC together with
the drug cisplatin and the combination carboplatin,paclitaxel would be:

.. code-block:: bash

  $ python main_deseq.py -p TCGA-CESC -p TCGA-HNSC
  -d cisplatin -d carboplatin,paclitaxel -o /OUTPUT_path -D -A

First of all, the needed data for the selected projects is loaded via the
TCGA API and stored in:   

* /OUTPUT_path/TCGA-CESC/TCGA-CESC_data_files/ and
* /OUTPUT_path/TCGA-HNSC/TCGA-HNSC_data_files/


Intermediate merged tables and additional meta_data tables are stored in the
project directories:   

* /OUTPUT_path/TCGA-CESC/
* /OUTPUT_path/TCGA-HNSC/


First, single project analyses are performed. The actual analysis is
determined by the project, and by the drugs combination. The directory for
the drugs combination is created out of the applied drugs, so here, the
DRUGS_title is 'carboplatin,paclitaxel_cisplatin'.  


Everything below that drugs directory, is restricted to the chosen drugs
s.t. the results of both single project analyses are placed in:  

* /OUTPUT_path/TCGA-CESC/carboplatin,paclitaxel_cisplatin/
* /OUTPUT_path/TCGA-HNSC/carboplatin,paclitaxel_cisplatin/


After the single project analysis, the projects are combined. Those results
are stored in an additional directory, composed out of the applied projects,
so here, the PROJECT_title is: 'TCGA-CESC_TCGA-HNSC', those results are saved
in the directory:  

* /OUTPUT_path/TCGA-CESC_TCGA-HNSC/carboplatin,paclitaxel_cisplatin/


Since the analysis is determined by the project and drug combination, results
for 3 different approaches are created, two for the single projects and one
for the aggregation of the two projects. For all of them, a respective
REPORT.pdf is created, containing a summarized representation of the most
important results and plots, along with some explanations to them. They are
stored at:

* /OUTPUT_path/TCGA-CESC/carboplatin,paclitaxel_cisplatin/REPORT.pdf
* /OUTPUT_path/TCGA-HNSC/carboplatin,paclitaxel_cisplatin/REPORT.pdf
* /OUTPUT_path/TCGA-CESC_TCGA-HNSC/carboplatin,paclitaxel_cisplatin/REPORT.pdf

Recreate the performed analysis:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To rerun the analysis and reproduce all the outputs and results created with
it, a single Snakemake configuration file is created. It is stored in the
cloned repository location under the 'Snakes' subdir.
Since the analysis is determined by the composition of projects and drugs, the
unique filename of this configuration file is composed out of it. For the
example with CESC and HNSC, together with cisplatin and carboplatin,paclitaxel,
that would be:

* SCRIPT_path/Snakes/snakemake_config_TCGA-CESC_TCGA-HNSC_carboplatin,paclitaxel_cisplatin.yaml

The Snakefile needed is also hold available at:

* SCRIPT_path/Snakes/Snakefile 

This file must be edited and the path to the config yaml file, the OUTPUT_path
and the SCRIPT_path must be inserted.

With that, the Snakefile is configured to run the analyses again. Change the
directory into the SCRIPT_path/Snakes/ path and run for example:

.. code-block:: bash

   $ snakemake --cores 7

This would use 7 cores of your machine if available and make use of
parallelisation of steps where it is feasible.
