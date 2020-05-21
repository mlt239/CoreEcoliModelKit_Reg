This repository contains a modified Core E.coli model from https://github.com/varnerlab/CoreEcoliModelKit that performs Boolean transcriptional regulation of FBA results. Many of the files are altered to account for the regulatory network solver, which is based on the regulatory network solver for the original Cobra Toolbox implementation of the Core E.coli model:

Orth JD, Fleming RMT, and B.O. Palsson (2010) Reconstruction and Use of Microbial Metabolic Networks: the Core Escherichia coli Metabolic Model as an Educational Guide. EcoSal Plus 2013; doi:10.1128/ ecosalplus.10.2.1

To perform a regulated and non-regulated FBA, use the function included in ``run_flux.jl``, noting that the outputs are regulated_FBA, non-regulated FBA, and non-regulated flux bounds. To initially run these, I had to make sure the path in ``Include.jl`` was set to my current ``pwd()`` which is ``\src``, and I initially disabled Gurobi. I have changed these back to default for this upload.

A few new files are added to perform the regulation

``optimizeRegModel.jl`` Will iteratively calculate the regulatory network until a steady-state is reached by calling ``solveBooleanReg()``. It will then go on to perform an FBA and output the altered data dictionary, fluxes, deleted genes, constrained fluxes, and regulatory network states. The function ``optimizeRegModel()`` must be given a data_dictionary as input, while an initial regulatory state is optional (default is all 0's.)

``solveBooleanReg.jl`` contains the ``solveBooleanReg()`` function which calculates the new regulatory network with each iteration. This takes the data_dictionary and initial states that are used in ``optimizeRegModel()`` as inputs. 4 things happen here: Generating new Regulatory Inputs 1 states based on certain metabolite exchanges, performing an FBA after deleting genes that are deleted by the initial regulatory state (using ``deleteModelGenes()``), generating new Regulatory Inputs 2 states based on whether certain fluxes are non-zero, and generating new Regulatory Rules states. Note, functionas contained in ``Rules.jl`` are used to evaluate boolean rules, as ``eval()`` cannot be called locally.

``cobrafunctions.jl`` only contains one function right now, ``deleteModelGenes()``, which takes data_dictionary and a list of deleted genes (determined in ``solveBooleanReg()``) as input. ``deleteModelGenes()`` will constrain fluxes that are dependent on genes that were deleted/turned off as a result of the initial regulatory state. Other functions are planned to go here, such as Cobra Toolbox's ``buildRxnGeneMat`` which is necessary to generate the gene-reaction matrix that is used here. Without these, the model will need to be updated manually to account for new reactions/regulatory rules. Finally, ``deleteModelGenes()`` also used a function in ``Rules.jl`` to evaluate boolean rules.

Other accessory functions: ``RulesGenerator.jl`` generates the ``Rules.jl`` functions that are used to evaluate the Boolean rules. This translates RegulatoryRules and rules into functions that allow them to be evaluated.
