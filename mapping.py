"""
mapping.py - Map a Hamiltonian to a LinearCombinaisonPauliString

Copyright 2020-2021 Maxime Dion <maxime.dion@usherbrooke.ca>
This file has been modified by <Your,Name> during the
QSciTech-QuantumBC virtual workshop on gate-based quantum computing.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

from pauli_string import PauliString, LinearCombinaisonPauliString
import numpy as np


class Mapping:

    def fermionic_creation_annihilation_operators(self, n_qubits: int) -> tuple[list[LinearCombinaisonPauliString],
                                                                                list[LinearCombinaisonPauliString]]:
        raise NotImplementedError("abstract base class implementation")
        
    def fermionic_hamiltonian_to_qubit_hamiltonian(self, fermionic_hamiltonian)->LinearCombinaisonPauliString:
        """
        Do the mapping of a FermionicHamiltonian. First generates the LCPS representation of the creation/annihilation
        operators for the specific mapping. Uses the 'to_pauli_string_linear_combinaison' of the FermionicHamiltonian
        to generate the complete LCPS.

        Args:
            fermionic_hamiltonian (FermionicHamiltonian): A FermionicHamiltonian that provided a 
                'to_pauli_string_linear_combinaison' method.

        Returns:
            LinearCombinaisonPauliString: The LCPS representing the FermionicHamiltonian
        """

        creation_operators, annihilation_operators = self.fermionic_creation_annihilation_operators(fermionic_hamiltonian.number_of_orbitals())
        qubit_hamiltonian = fermionic_hamiltonian.to_linear_combinaison_pauli_string(creation_operators, annihilation_operators)
        return qubit_hamiltonian


class JordanWigner(Mapping):
    def __init__(self):
        """
        The Jordan-Wigner mapping
        """

        self.name = 'jordan-wigner'

    def fermionic_creation_annihilation_operators(self, n_qubits: int) -> tuple[list[LinearCombinaisonPauliString],
                                                                                list[LinearCombinaisonPauliString]]:
        """
        Build the LCPS reprensetations for the creation/annihilation operator for each qubit following 
        Jordan-Wigner mapping.

        Args:
            n_qubits (int): The number of orbitals to be mapped to the same number of qubits.

        Returns:
            list<LinearCombinaisonPauliString>, list<LinearCombinaisonPauliString>: Lists of the creation/annihilation
                operators for each orbital in the form of LinearCombinaisonPauliString.
        """

        creation_operators = list()
        annihilation_operators = list()

        eye = np.eye(n_qubits, dtype=bool)
        zeros = np.zeros(n_qubits, dtype=bool)

        for i in range(n_qubits):  # iterate over states
            c_i = 0.5 * PauliString.from_zx_bits(np.concatenate((zeros, eye[i]))) + (-0.5j) * PauliString.from_zx_bits(np.concatenate((eye[i], eye[i])))
            a_i = 0.5*PauliString.from_zx_bits(np.concatenate((zeros, eye[i]))) + 0.5j*PauliString.from_zx_bits(np.concatenate((eye[i], eye[i])))
            for j in range(i):
                a_i *= (1*PauliString.from_zx_bits(np.concatenate((eye[j], zeros))))
                c_i *= (1*PauliString.from_zx_bits(np.concatenate((eye[j], zeros))))
            annihilation_operators.append(a_i)
            creation_operators.append(c_i)

        return creation_operators, annihilation_operators


class Parity(Mapping):
    def __init__(self):
        """
        The Parity mapping
        """

        self.name = 'parity'

    def fermionic_creation_annihilation_operators(self, n_qubits: int) -> tuple[list[LinearCombinaisonPauliString],
                                                                                list[LinearCombinaisonPauliString]]:
        """
        Build the LCPS reprensetations for the creation/annihilation operator for each qubit following 
        Parity mapping.

        Args:
            n_qubits (int): The number of orbtials to be mapped to the same number of qubits.

        Returns:
            list<LinearCombinaisonPauliString>, list<LinearCombinaisonPauliString>: Lists of the creation/annihilation
                operators for each orbital in the form of LinearCombinaisonPauliString
        """

        creation_operators = list()
        annihilation_operators = list()
        
        ################################################################################################################
        # YOUR CODE HERE
        # OPTIONAL
        ################################################################################################################

        raise NotImplementedError()

        return creation_operators, annihilation_operators
