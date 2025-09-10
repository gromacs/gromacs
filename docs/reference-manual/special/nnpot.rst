.. _nnpot:

Neural Network Potentials
------------------------------------

The NNPot module allows the use of potentials based on neural networks or
machine learning in general. These models can be trained to reproduce forces
and energies at *ab initio* levels of accuracy, based on training data from
electronic structure calculations with e.g. DFT or CCSD(T). They take atomic
descriptors as their input, most commonly the atom positions and atomic numbers,
and output a value for the potential energy of the given conformation. The forces
can then be calculated by backpropagating the gradients through the machine
learning model (e.g. PyTorchs AutoGrad engine), w.r.t. the inputs.
Currently, only models trained in `PyTorch <https://pytorch.org/>`_,
exported using its TorchScript functionality, are supported.

Hybrid NNP/MM Simulations
^^^^^^^^^^^^^^^^^^^^^^^^^

The NNPot interface also supports hybrid NNP/MM simulations, where only a small
subsystem is modeled with the neural network potential, and the rest of the system
is modeled using regular force fields. The embedding scheme used is analoguous
to the additive mechanical embedding scheme in QMMM simulations, where the total
energy of the system is expressed as :math:`E_{tot} = E_{NNP} + E_{MM} + E_{NNP-MM}`,
where :math:`E_{NNP}` describes the energy predicted by the model for the NNP subsystem,
:math:`E_{MM}` describes all other interactions calculated using the classical MM
force field. The coupling term :math:`E_{NNP-MM}` term includes non-bonded electrostatic
and LJ interactions between the atoms in the NNP and MM regions, calculated as usual
in GROMACS. Bonded interactions are also described on the MM level, removing terms
consisting of bonds containing 2 NNP atoms, angles and settles containing 2 or 3 NNP atoms,
and dihedrals containing 3 or 4 NNP atoms. Broken chemical bonds between NNP and MM atoms
are currently *not* capped with a link atom, as is usual in QMMM simulations.
This might lead to unusual local chemical environments the model is not familiar with,
leading to unexpected behaviour. NNP/MM simulations with NNP subsystems that cut through
chemical bonds are therefore discouraged as of now. All necessary modifications
to the system topology are performed automatically during :ref:`gmx grompp` preprocessing.

Software Prerequisites
^^^^^^^^^^^^^^^^^^^^^^

To perform simulations with the NNPot module, GROMACS needs to be linked with
a `LibTorch installation <https://pytorch.org/get-started/locally/>`_ (version
2.0 or higher). For specific installations instructions please see
`this section <installing with Neural Network potential support>` of the install guide.

Usage
^^^^^

Simulations using the NNPot interface are controlled by setting :ref:`mdp` file options.
In particular, the relevant options, specifying the behaviour of a simulation
with ``nnpot-active`` set to ``true``, are:

-  ``nnpot-modelfile``: Specifies a path to a TorchScript-compiled model, either absolute
   or relative to the simulation directory. If not provided, the interface will look for
   a file named ``model.pt`` in the current working directory.
-  ``nnpot-input-group``: Specifies an [index group] defining the input atoms for
   the NNP subsystem. Defaults to ``System``, which performs a pure NNP simulation.
-  ``nnpot-model-input[1-4]``: These options can be used to specify the inputs
   for the model. Supported options are ``atom-positions``, a vector containing the input
   atom positions; ``atom-numbers``, a vector containing atomic numbers; ``box``, the unit
   vectors of the simulation box; ``pbc``, a boolean vector specifying PBC type. 

The inputs are passed to the model in order of their occurence in the mdp file. Note
that there are no default values for the model input, so not specifying the model
input will lead to errors. The model is expected to return a tensor containing the energy
of the system and, optionally, a tensor containing the forces on the input atoms.
This option can be useful for cases in which the forces can be computed by the model
by some technique that is more efficient than via backpropagation. If the model does not
provide its own forces, they are calculated by the interface as gradients
w.r.t. the *first* input tensor. \
You can specify the device on which you wish to run model inference using the
environment variable ``GMX_NN_DEVICE``. For now, only ``cpu`` and ``cuda`` are supported.
For ``cuda``, a CUDA-aware LibTorch installation should be installed, and the corresponding
CUDA installation should be available in your ``PATH``. 

Find below an example Python code snippet to export a pretrained model in
PyTorch using TorchScript. As an example, we use the popular ANI models
available from `TorchANI <https://github.com/aiqm/torchani/>`_. 

::

    import torch
    from torch import nn
    from torchani.models import ANI2x
    from typing import Optional

    class GmxNNPotModelWrapper(nn.Module):
        def __init__(self):
            super().__init__()

            # Load a pre-trained ANI-2x model
            self.model = ANI2x(periodic_table_index=True)

            # GROMACS and TorchANI use different unit conventions
            self.length_conversion = 10.0   # nm --> Å
            self.energy_conversion = 2625.5 # Hartree --> kJ/mol

        def forward(self, positions, atomic_numbers, 
                    box: Optional[torch.Tensor]=None, pbc: Optional[torch.Tensor]=None):
            
            # Prepare the inputs for the model
            atomic_numbers = atomic_numbers.unsqueeze(0)
            positions = positions.unsqueeze(0) * self.length_conversion
            if box is not None:
                box *= self.length_conversion

            # Forward pass
            result = self.model((atomic_numbers, positions), box, pbc)

            energy = result.energies[0] * self.energy_conversion

            return energy
        
    model = GmxNNPotModelWrapper()

    save_path = 'ani2x.pt'
    torch.jit.script(model).save(save_path)

The model can then be used in |Gromacs| by specifying the path to the saved model.
Take care that the LibTorch version linked to |Gromacs| matches the one that
was used to train/export the model.