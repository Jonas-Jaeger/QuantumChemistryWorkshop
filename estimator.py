"""
estimator.py - To estimate expectation value of observables

Copyright 2020-2023 Maxime Dion <maxime.dion@usherbrooke.ca>
This file has been modified by **Jonas Jaeger** <jojaeger@cs.ubc.ca> during the
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

import time

import numpy as np

from qiskit import QuantumCircuit, execute
from qiskit.providers import Backend

from pauli_string import PauliString, LinearCombinaisonPauliString

from typing import Union, Optional
from numpy.typing import NDArray


class Estimator:
    def __init__(self, varform: QuantumCircuit, backend: Backend, execute_opts={}, record: Optional[object] = None):
        """
        An Estimator allows to transform an observable into a callable function. The observable is not set at the 
        initialization. The estimator will build the QuantumCircuit necessary to estimate the expected value of the 
        observable. Upon using the 'eval' method, it will execute these circuits and interpret the results to return an
        estimate of the expected value of the observable.

        Args:
            varform (QuantumCircuit): A paramatrized QuantumCircuit.
            backend (Backend): A qiskit backend. Could be a simulator are an actual quantum computer.
            execute_opts (dict, optional): Optional arguments to be passed to the qiskit.execute function.
                                           Defaults to {}.

            record (object, optional): And object that could be called on each evaluation to record the results. 
                Defaults to None.
        """

        self.varform = varform
        self.backend = backend
        self.execute_opts = execute_opts

        self.record = record

        # To be set attributes
        self.n_qubits = varform.num_qubits
        self.diagonalizing_circuits = list()
        self.diagonal_observables = list()

    def set_observable(self, observable: LinearCombinaisonPauliString) -> None:
        """
        Set the observable which the expectation value will be estimated.
        This sets the value of the attribute 'n_qubits'.
        The observable is converted into a list of diagonal observables along the circuit which performs this diagonalization.
        This is done using the 'diagonal_observables_and_circuits' method (defined at the subclass level).

        Args:
            observable (LinearCombinaisonPauliString): The observable to be evaluated.
        """

        self.diagonal_observables, self.diagonalizing_circuits = self.diagonal_observables_and_circuits(observable)

    def eval(self, params: Union[NDArray, list]) -> float:
        """
        Evaluate an estimation of the expectation value of the set observable.

        Args:
            params (list or NDArray): Parameter values at which the expectation value should be evaluated.
                Will be fed to the 'varform' paramatrized QuantumCircuit.

        Returns:
            float: The estimated expectation value of the observable.
        """

        t0 = time.time()
        state_circuit = self.prepare_state_circuit(params)
        circuits = self.assemble_circuits(state_circuit)

        expectation_value = 0.

        if isinstance(self.backend, Backend):
            job = execute(circuits, backend=self.backend, **self.execute_opts)
            results = job.result()

            for diag_observable, counts in zip(self.diagonal_observables, results.get_counts()):
                expectation_value += Estimator.estimate_diagonal_observable_expectation_value(diag_observable, counts)
        else:
            from qiskit.quantum_info import SparsePauliOp
            observables = [SparsePauliOp(str(diag_obs.pauli_strings[0])) for diag_obs in self.diagonal_observables]
            results = self.backend.run(circuits=circuits, observables=observables, run_options=self.execute_opts).result()
            for diag_observable, value in zip(self.diagonal_observables, results.values):
                assert np.isreal(diag_observable.coefs[0])
                expectation_value += np.real(diag_observable.coefs[0]) * value


        eval_time = time.time()-t0
        print("Eval time (s): ", eval_time)

        return expectation_value

    def prepare_state_circuit(self, params: Union[NDArray, list]) -> QuantumCircuit:
        """
        Assign parameter values to the variational circuit (varfom) to prepare the quantum state.

        Args:
            params (list or NDArray): Params to be assigned to the 'varform' QuantumCircuit.

        Returns:
            QuantumCircuit: The quantum state circuit
        """

        state_circuit = self.varform.bind_parameters(dict(zip(self.varform.parameters, params)))

        return state_circuit

    def assemble_circuits(self, state_circuit: QuantumCircuit) -> list[QuantumCircuit]:
        """
        For every diagonal observable, assemble the complete circuit with:
        - State preparation
        - Measurement circuit
        - Measurements

        Args:
            state_circuit (QuantumCircuit): The quantum state circuit

        Returns:
            list<QuantumCircuit>: The quantum circuits to be executed.
        """

        circuits = list()

        for diag_lcps, diag_circuit in zip(self.diagonal_observables, self.diagonalizing_circuits):
            assert len(diag_lcps.pauli_strings) == 1
            diag_pauli_string = diag_lcps.pauli_strings[0]

            assembled_qc = state_circuit + diag_circuit
            assembled_qc.measure_all()

            circuits.append(assembled_qc)


        return circuits

    @staticmethod
    def diagonal_pauli_string_eigenvalue(diagonal_pauli_string: PauliString, state: str) -> float:
        """
        Computes the eigenvalue (+1 or -1) of a diagonal pauli string for a basis state.

        Args:
            diagonal_pauli_string (PauliString): A diagonal pauli string
            state (str): a basis state (ex : '1011')

        Returns:
            float: The eigenvalue
        """

        state_array = -2*np.array(list(map(int, state))) + 1
        comp_basis_meas_indices = np.where(np.logical_not(diagonal_pauli_string.ids()))
        eigenvalue = np.prod(state_array[comp_basis_meas_indices]).item()

        return eigenvalue

    @staticmethod
    def estimate_diagonal_pauli_string_expectation_value(diagonal_pauli_string: PauliString, counts: dict) -> float:
        """
        Estimate the expectation value for a diagonal pauli string based on counts.

        Args:
            diagonal_pauli_string (PauliString): The diagonal pauli string (must be only I and Z)
            counts (dict): Contains the number of times each basis state was obtained from a measurement

        Returns:
            float: The expectation value of the Pauli string
        """

        expectation_value = 0.

        for basis_state, n in counts.items():
            eigenvalue = Estimator.diagonal_pauli_string_eigenvalue(diagonal_pauli_string, basis_state)
            expectation_value += n * eigenvalue

        n_total = sum(counts.values())
        expectation_value /= n_total # normalize

        return expectation_value

    @staticmethod
    def estimate_diagonal_observable_expectation_value(diagonal_observable: LinearCombinaisonPauliString,
                                                       counts: dict) -> float:
        """
        Estimate the expectation value for a diagonal observable (linear combinaison of pauli strings) based on counts.

        Args:
            diagonal_observable (LinearCombinaisonPauliString): The observable (must be only I and Z)
            counts (dict): Contains the number of times each basis state was obtained from a measurement

        Returns:
            float: The expectation value of the Observable
        """

        expectation_value = 0.

        for pauli_string, coef in zip(diagonal_observable.pauli_strings, diagonal_observable.coefs):
            pauli_str_expect_val = Estimator.estimate_diagonal_pauli_string_expectation_value(pauli_string, counts)
            assert np.isreal(coef)
            expectation_value += np.real(coef) * pauli_str_expect_val

        return expectation_value

    @staticmethod
    def diagonalizing_pauli_string_circuit(pauli_string: PauliString) -> tuple[QuantumCircuit, PauliString]:
        """
        Builds the circuit representing the transformation which diagonalizes the given PauliString.

        Args:
            pauli_string (PauliString): The pauli string to be diagonalized

        Returns:
            QuantumCircuit: A quantum circuit representing the transformation which diagonalizes the given PauliString
            PauliString: The diagonal PauliString                   
        """
        
        n_qubits = len(pauli_string)
        diagonalizing_circuit = QuantumCircuit(n_qubits)
        diagonal_pauli_string = list()

        for i, letter in enumerate(str(pauli_string)):
            if letter == "X":
                diagonalizing_circuit.h(n_qubits - i - 1)
                diagonal_pauli_string.append("Z")
            elif letter == "Y":
                diagonalizing_circuit.sdg(n_qubits - i - 1)
                diagonalizing_circuit.h(n_qubits - i - 1)
                diagonal_pauli_string.append("Z")
            else:
                diagonal_pauli_string.append(letter)
        diagonal_pauli_string = PauliString.from_str("".join(diagonal_pauli_string))


        return diagonalizing_circuit, diagonal_pauli_string


class BasicEstimator(Estimator):
    """
    The BasicEstimator should build 1 quantum circuit and 1 interpreter for each PauliString.
    The interpreter should be 1d array of size 2**number of qubits.
    It does not exploit the fact that commuting PauliStrings can be evaluated from a common circuit.
    """
    
    @staticmethod
    def diagonal_observables_and_circuits(observable: LinearCombinaisonPauliString) -> tuple[list[LinearCombinaisonPauliString],
                                                                                             list[QuantumCircuit]]:
        """
        This method converts each PauliString in the observable into :
        - A diagonal observable including the associated coefficient
        - The quantum circuit that performs this diagonalization.
        The diagonal observables and the quantum circuits are return as two list.

        Args:
            observable (LinearCombinaisonPauliString): An observable.

        Returns:
            list<LinearCombinaisonPauliString> : The diagonal observables (one PauliString long).
            list<list of QuantumCircuit> : The diagonalizing quantum circuits
        """
        
        diagonal_observables = list()
        diagonalizing_circuits = list()

        for coef, pauli_string in zip(observable.coefs, observable.pauli_strings):
            diag_circuit, diag_pauli_string = Estimator.diagonalizing_pauli_string_circuit(pauli_string)
            diagonal_observables.append(coef * diag_pauli_string)
            diagonalizing_circuits.append(diag_circuit)

        return diagonal_observables, diagonalizing_circuits


class BitwiseCommutingCliqueEstimator(Estimator):
    """
    The BitwiseCommutingCliqueEstimator should build 1 quantum circuit and 1 interpreter for each clique of PauliStrings.
    The interpreter should be 2d array of size (number of cliques ,2**number of qubits).
    It does exploit the fact that commuting PauliStrings can be evaluated from a common circuit.
    """

    @staticmethod
    def diagonal_observables_and_circuits(observable: LinearCombinaisonPauliString) -> tuple[list[LinearCombinaisonPauliString],
                                                                                             list[QuantumCircuit]]:
        """
        This method first divide the observable into bitwise commuting cliques. Each commuting clique is then converted
        into :
        - A diagonal observable including the associated coefficients
        - The quantum circuit that performs this diagonalization.
        The diagonal observables and the quantum circuits are return as two list.

        Args:
            observable (LinearCombinaisonPauliString): An observable.

        Returns:
            list<LinearCombinaisonPauliString> : The diagonal observables
            list<list of QuantumCircuit> : The diagonalizing quantum circuits
        """

        cliques = observable.divide_in_bitwise_commuting_cliques()

        diagonal_observables = list()
        diagonalizing_circuits = list()
        
        ################################################################################################################
        # YOUR CODE HERE
        # TO COMPLETE (after lecture on VQE)
        # Hint : the next method does the work for 1 PauliString + coef
        ################################################################################################################
        
        raise NotImplementedError()
            
        return diagonal_observables, diagonalizing_circuits


