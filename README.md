==================================
CtxMod_BCVTB_EnergyPlus Dataset
==================================

This repository contains simulation files used for Context-based Thermodynamic Modeling of Buildings Spaces.

==================================
Dataset download
==================================

This dataset can also be downloaded from:

http://datasets.isr.tecnico.ulisboa.pt/BCVTB/

==================================
Citation
==================================

If you plan to use the dataset, please cite it as follows:

@PHDTHESIS {PedroFazendaphd,
    author  = "Pedro Fazenda",
    title   = "Reinforcement Learning to Optimize Occupant Comfort and Energy Usage in {HVAC} Systems and 
    Context-based Thermodynamic Modeling of Building Spaces",
    school  = "Instituto Superior T\'ecnico, Universidade T\'ecnica de Lisboa",
    year    = "2016",
    address = "Instituto Superior T\'ecnico, Av. Rovisco Pais 1, 1049-001 Lisboa, Portugal",
    month   = "",
    note    = "\"Instituto Superior T\'ecnico\" was formerly within  \"Universidade T\'ecnica de Lisboa\" and is now within \"Universidade de Lisboa\""
}

==================================
Brief Description
==================================

These simulation files concern the experimental setup and application example used for
Context-based Thermodynamic Modeling of Building Spaces. The repository contains simulation files for a single Thermal Zone
using the Buildings Controls Virtual Test Bed (BCVTB) and EnergyPlus. It also contains the Matlab simulation of
the same building, using a Context-based model. The repository includes the following files:


Simple.idf: The EnergyPlus input building used in the simulation.
system.xml: The BCVTB system description. 
Model.m: Simulation of the Context-based model.
