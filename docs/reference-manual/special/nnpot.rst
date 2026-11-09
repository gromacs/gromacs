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
to the additive mechanical embedding scheme in QM/MM simulations, where the total
energy of the system is expressed as :math:`E_{tot} = E_{NNP} + E_{MM} + E_{NNP-MM}`,
where :math:`E_{NNP}` describes the energy predicted by the model for the NNP subsystem,
:math:`E_{MM}` describes all other interactions calculated using the classical MM
force field. The coupling term :math:`E_{NNP-MM}` term includes non-bonded electrostatic
and LJ interactions between the atoms in the NNP and MM regions, calculated as usual
in |Gromacs|. An exception to this is when the NNP model itself calculates the electrostatic
interactions, in which case the coupling term only includes the LJ interactions.
Bonded interactions are also described on the MM level, removing terms
consisting of bonds containing 2 NNP atoms, angles and settles containing 2 or 3 NNP atoms,
and dihedrals containing 3 or 4 NNP atoms. Broken chemical bonds between NNP and MM atoms
are capped with a link atom, a common practice in QM/MM simulations. The link atom is placed at
a fixed distance from the NNP atom along the vector pointing towards the MM atom, and
the force acting on it is redistributed to the NNP and MM atoms according to the chain rule,
as if it was a 2FD-type virtual site. The type of link atom and the distance from the NNP atom
can be specified via :ref:`mdp` options, see below. All necessary modifications to the system
topology are performed automatically during :ref:`gmx grompp` preprocessing.

Software Prerequisites
^^^^^^^^^^^^^^^^^^^^^^

To perform simulations with the NNPot module, |Gromacs| needs to be linked with
a `LibTorch installation <https://pytorch.org/get-started/locally/>`_ (version
2.0 or higher). For specific installations instructions please see
`this section <installing with Neural Network potential support>` of the install guide.

Usage
^^^^^

Simulations using the NNPot interface are controlled by setting :ref:`mdp` file options.
For a detailed description of all available options, please see the
:ref:`mdp options documentation<mdp-nnpot>`.
Since the 2027 version, the inputs are passed to the model as a dictionary of tensors.
That is, the model is expected to take a single dictionary as argument to its ``forward()``
method, e.g. ``forward(input_dict: Dict[str, torch.Tensor]) -> Dict[str, torch.Tensor]``.
Available inputs, and their corresponding shapes and dictionary keys, are as follows:

* ``positions``: A tensor of shape ``(N_ml, 3)`` containing the positions of the :math:`N_{\text{ml}}` ML atoms.
* ``atomic-numbers``: A tensor of shape ``(N_ml,)`` containing the atomic numbers of the ML atoms.
* ``pairs``: A tensor of shape ``(M, 2)`` containing the indices of ``M`` nearest-neighbor pairs.
    of ML atoms with in the ``pair-cutoff`` specified by the user. Available only if the ``pair-cutoff`` option is set
    to a non-zero value in the :ref:`mdp` file.
* ``shifts``: A tensor of shape ``(M, 3)`` containing periodic shift vectors to apply to the
    neighbor pairs in order to get the minimum image distance between them. Available only if the ``pair-cutoff`` option
    is set to a non-zero value in the :ref:`mdp` file.
* ``mm-positions``: A tensor of shape ``(N_mm, 3)`` containing the positions of the :math:`N_{\text{mm}}` MM atoms surrounding
    the NNP atoms. Available only if the :mdp-value:`nnpot-embedding=electrostatic-model` option is set in the :ref:`mdp` file.
* ``mm-charges``: A tensor of shape ``(N_mm,)`` containing the charges of the :math:`N_{\text{mm}}` MM atoms surrounding
    the NNP atoms. Available only if the :mdp-value:`nnpot-embedding=electrostatic-model` option is set in the :ref:`mdp` file.
* ``total-charge``: A tensor of shape ``(1,)`` containing the total charge of the NNP subsystem, as set by the user with the
    :mdp-value:`nnpot-nnp-charge` option.
* ``box``: A tensor of shape ``(3, 3)`` containing the vectors of the simulation box.
* ``pbc``: A boolean tensor of shape ``(3,)`` indicating whether PBCs are used along this dimension.

The precision of real-valued input tensors is depends on the precision |Gromacs| was compiled with,
and can be either 32-bit (single precision) or 64-bit (double precision) floating point. Atom indices
are passed as 32-bit integers. 
The model is expected to return a dictionary with a ``(1,)``-shape tensor containing the energy of the system under
the key ``energy`` and, optionally, a ``(N_ml, 3)`` tensor containing the forces on the input atoms (``forces``).
This option can be useful for cases in which the forces can be computed by the model
by some technique that is more efficient than via backpropagation. If the model does not
provide its own forces, they are calculated by the interface as gradients w.r.t. the position tensor.
In the case of :mdp-value:`nnpot-embedding=electrostatic-model`, forces on MM atoms
should be additionally returned by the model as a ``(N_mm, 3)`` tensor, under the key ``forces-mm``. \
By default, |Gromacs| will run inference on the GPU, if available. You can also explicitly
specify the device on which you wish to run model inference using the environment
variable ``GMX_NN_DEVICE``. The options ``cpu`` and ``gpu`` are supported. Option ``cuda``
is deprecated, and is treated as ``gpu``.
If multiple GPUs are available, the active device will be inherited from mdrun.
For ``gpu``, a CUDA/HIP-aware LibTorch installation should be installed, and the
corresponding CUDA/HIP installation should be available in your ``PATH``.
On Intel GPUs, a LibTorch with XPU support (PyTorch 2.5 or later) is needed. Unlike
for CUDA and HIP, the device is not inherited from mdrun: the model runs on the first
Intel GPU LibTorch enumerates, which can be restricted with ``ONEAPI_DEVICE_SELECTOR``
or ``ZE_AFFINITY_MASK``. When |Gromacs| itself is built with SYCL, two things apply at
run time. First, |Gromacs| and LibTorch must use one and the same SYCL runtime
(``libsycl``): build |Gromacs| with the oneAPI release whose runtime the LibTorch
package bundles, or make sure the runtime bundled with LibTorch is the one loaded;
with two different runtimes in the process either side can crash. Second, expose a
single SYCL backend, e.g. ``ONEAPI_DEVICE_SELECTOR=level_zero:gpu``, as |Gromacs|
itself recommends when it detects a GPU through more than one backend.

Find below an example Python code snippet to export a pretrained model in
PyTorch using TorchScript. As an example, we use the popular ANI models
available from `TorchANI <https://github.com/aiqm/torchani/>`_.

::

    import torch
    from torch import nn
    from torchani.models import ANI2x
    from typing import Dict

    class GmxNNPotModelWrapper(nn.Module):
        def __init__(self):
            super().__init__()

            # Load a pre-trained ANI-2x model
            self.model = ANI2x(periodic_table_index=True)

            # GROMACS and TorchANI use different unit conventions
            self.length_conversion = 10.0   # nm --> Å
            self.energy_conversion = 2625.5 # Hartree --> kJ/mol

        def forward(self, input_data: Dict[str, torch.Tensor]) -> Dict[str, torch.Tensor]:
            # Prepare the inputs for the model
            atomic_numbers = input_data["atomic-numbers"].unsqueeze(0)
            positions = input_data["positions"].unsqueeze(0) * self.length_conversion
            # positions = positions.detach().requires_grad_(True)
            if "box" in input_data and input_data["box"] is not None:
                cell = input_data["box"] * self.length_conversion
            else:
                cell = None
            if "pbc" in input_data and input_data["pbc"] is not None:
                pbc = input_data["pbc"]
            else:
                pbc = None

            # Forward pass
            result = self.model((atomic_numbers, positions), box, pbc)
            output = {}
            output["energy"] = result.energies[0] * self.energy_conversion

            return output

    model = GmxNNPotModelWrapper()

    save_path = 'ani2x.pt'
    torch.jit.script(model).save(save_path)

The model can then be used in |Gromacs| by specifying the path to the saved model.
Take care that the LibTorch version linked to |Gromacs| matches the one that
was used to train/export the model.
In the :ref:`mdp` file, the relevant options could look like:

::

    nnpot-active          = true
    nnpot-modelfile       = ani2x.pt
    nnpot-input-group     = System
