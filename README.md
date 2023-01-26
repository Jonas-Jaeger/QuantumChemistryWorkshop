# molecule-incomplete-suggested-solution
Incomplete version of the suggested solution to find the ground state of a molecule using quantum computing.

Description of the files :
- hamiltonian.py : This files defines the FermionicHamiltonian class and subclasses. You should be able to partially complete it after the lecture on second quantization. The 'to_linear_combinaison_pauli_string' methods can be completed after the lecture on mapping.
- pauli_string.py : Defines PauliString and LinearCombinaisonPauliString class. You should be able to complete it after the lecture on mapping. The 'to_matrix' method is optional.
- mapping.py : Defines the JordanWigner mapping. You should be able to complete it after the lecture on mapping.
- estimator.py : Defines the abstract class Estimator and the BasicEstimator class. You should be able to complete it after the lecture on VQE.
- solver.py : Defines VQESolver and ExactSolver. You should be able to complete it after the lecture on VQE. The ExactSolver is optionnal.

Other files :
- Integrals_sto-3g_H2_d_0.7350_no_spin.npz : Contains the one body and two body integrals (no spin) for a H2 molecule with d=0.735. The two body is given in the physicist order.
- activity_mppaing_.ipynb and activity_vqe_.ipynb : Tutorial Jupyter notebooks to help you code the concepts seen in the respective activities.
