# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# [<img src="images/quantum_algorithms_tutorials.png" alt="drawing" width="100" align="left"/>][6]
#
# <h2 align="center">
# Estimation of the Bloch sphere angles
# </h2>
#
# **[[Homepage][6]]**
# **[[Open with the nbviewer][7]]**
#
# <h3 align="left">
# Introduction
# </h3>
#
# Sometimes in Quantum Computing one needs to check or test the output of the quantum algorithms. In this tutorial, we are going to show how one can estimate Bloch sphere angles that can be used for example for testing quantum state teleportation protocol or the simplest and smallest implementations of the HHL algorithm. For the same proposes one can use the Swap test, that is also introduced here. Actually, we have introduced the Swap test for testing the Bloch sphere angles estimation procedure. Here we talk mainly about quantum pure states, but at the end, we have added what one can do for a more general case (how to implement the quantum state tomography for a mixed state).
#
# In Quantum Computing one of the most important and basic ideas is the concept of the one-qubit quantum state. In general one-qubit state can be described with two complex numbers [[1]]:
#
# \begin{equation*}
# \left|\psi \right\rangle = \alpha \left|0 \right\rangle + \beta \left|1 \right\rangle
# \end{equation*}
#
# where $\alpha$ and $\beta$ are complex numbers. Although $2$ complex numbers can be described with $4$ real numbers, we will show that we need only two numbers in order to describe one-qubit state. $\left| \alpha \right|^2$ and $\left| \beta \right|^2$ are probabilities of measuring the qubit in respectively $\left| 0 \right\rangle$ and $\left| 1 \right\rangle$ states. The sum of the probabilities should be equal to unity (normalization):
#
# $$\left| \alpha \right|^2 + \left| \beta \right|^2 = 1$$
#
# Also, we have the global phase ambiguity. The idea of the global phase ambiguity is in the impossibility to distinguish between two states that differ from each other only by $e^{i \varphi_g}$ global phase [[1]]:
#
# \begin{equation*}
# \left|\psi_1 \right\rangle = \alpha \left| 0 \right\rangle + \beta \left| 1 \right\rangle
# \qquad
# \qquad
# \left|\psi_2 \right\rangle  = e^{i \varphi_g} \left( \alpha \left| 0 \right\rangle + \beta \left| 1 \right\rangle \right)
# \end{equation*}
#
# In other words, there is no measurement procedure that will let us find $e^{i \varphi_g}$. So, experimentally/practically $\left|\psi_1 \right\rangle$ and $\left|\psi_2 \right\rangle$ can be regarded as the same states [[1]]. That is why:
#
# \begin{align*}
# \left|\psi \right\rangle &= \alpha \left| 0 \right\rangle + \beta \left| 1 \right\rangle = e^{i\varphi_g}|\alpha| \left| 0 \right\rangle + \beta \left| 1 \right\rangle =
# \\
# &= e^{i\varphi_g} \left( |\alpha| \left| 0 \right\rangle + e^{-i\varphi_g} \beta  \left| 1 \right\rangle \right) \stackrel{\text{up to }\varphi_g}{=} \alpha' \left| 0 \right\rangle + \beta' \left| 1 \right\rangle
# \end{align*}
#
# where $\alpha' = |\alpha|$ is a real positive number and $\beta' = e^{-i\varphi_g} \beta$ is a complex number. So, we are left with one real number and one complex number that describe one-qubit state (overall $3$ real numbers). If we will add the normalization $|\alpha'|^2 + |\beta'|^2 = 1$ we can reduce the needed real numbers to describe one-qubit state:
#
# \begin{align*}
# \left|\psi \right\rangle &= \alpha' \left| 0 \right\rangle + \beta' \left| 1 \right\rangle = \alpha' \left| 0 \right\rangle + e^{i\varphi}|\beta'| \left| 1 \right\rangle =
# \\
# &=\cos\left(\frac{\theta}{2}\right) \left| 0 \right\rangle + e^{i\varphi} \sin\left(\frac{\theta}{2}\right) \left| 1 \right\rangle
# \end{align*}
#
# where $\theta \in [0, \pi ]$ is chosen such that $\cos\left(\frac{\theta}{2}\right) = \alpha'$ and $\sin\left(\frac{\theta}{2}\right) = |\beta'|$, because $\cos^2\left(\frac{\theta}{2}\right) + \sin^2\left(\frac{\theta}{2}\right) = 1$. $\phi \in [0, 2 \pi)$ is the phase of the $\beta'$ complex number. Thus we can describe one-qubit states with just two numbers: $\theta$ and $\phi$. The goal of this tutorial is to estimate $\theta$ and $\phi$ (aka Bloch sphere angles) with Qiskit. All possible values of the $\left(\theta,\phi \right)$ pair can be described geometrically with the Bloch sphere [[1]]:
#
# <img src="images/bloch_sphere.png" alt="drawing" width="300"/>
#
# Each point on the Bloch sphere represents a distinct quantum (we assume pure) state. If we want to obtain one-qubit quantum state with given $\theta$ and $\varphi$ we can apply Qiskit's $U(\theta, \phi, 0)$ gate to the qubit that is in the $\left| 0 \right\rangle$ state (the initial default quantum state for most QCs):
#
# \begin{align*}
# U(\theta, \phi, 0) \left| 0 \right\rangle &= 
# \begin{pmatrix}
# \cos\left(\frac{\theta}{2}\right) & -\sin\left(\frac{\theta}{2}\right) \\
# e^{i\varphi} \sin\left(\frac{\theta}{2}\right) & e^{i\varphi}\cos\left(\frac{\theta}{2}\right)
# \end{pmatrix}
# \begin{pmatrix}
# 1 \\
# 0
# \end{pmatrix}
# =
# \\
# &= \begin{pmatrix}
# \cos\left(\frac{\theta}{2}\right) \\
# e^{i\varphi} \sin\left(\frac{\theta}{2}\right)
# \end{pmatrix} = \cos\left(\frac{\theta}{2}\right) \left| 0 \right\rangle + e^{i\varphi} \sin\left(\frac{\theta}{2}\right) \left| 1 \right\rangle
# \end{align*}
#
# Now let's assume that after some quantum circuit we have a quantum state whose Bloch sphere angles ($\theta$ and $\varphi$) we want to estimate. In other words, we want to know what are the Bloch angles of the qubit. Here we are going to discuss one such procedure that finds the Bloch angles.
#
# In this procedure, we assume that we can recreate the $\left| \psi \right\rangle$ one-qubit state of interest as many times as we want. This means, that the quantum circuit that outputs $\left| \psi \right\rangle$ can be re-executed. Also, we assume that the qubit is not in the entangled state (it is **one** qubit state) and we don't have quantum errors (one-qubit **pure** state). The last two assumptions can be easily checked in our procedure, hence, we will check them with few code lines.
#
#
# For this procedure we need to calculate the expectation values of $X$, $Y$ and $Z$ operators. From these three real numbers (expectation values) we can determine the Bloch angles. Firstly let's mathematically calculate the $Z$ expectation value $\left\langle \psi \right| Z \left|\psi \right\rangle$:
#
# \begin{align*}
# \left\langle \psi \right| Z \left|\psi \right\rangle &= \left( \cos\left(\frac{\theta}{2}\right) \left\langle 0 \right| + e^{-i\varphi} \sin\left(\frac{\theta}{2}\right) \left\langle 1\right| \right) Z \left( \cos\left(\frac{\theta}{2}\right) \left| 0 \right\rangle + e^{i\varphi} \sin\left(\frac{\theta}{2}\right) \left| 1 \right\rangle \right) =
# \\
# &= \cos^2\left(\frac{\theta}{2}\right) - \sin^2\left(\frac{\theta}{2}\right)
# \end{align*}
#
# where we took into account that $Z\left| 0 \right\rangle = \left| 0 \right\rangle$, $Z \left| 1 \right\rangle = -\left| 1 \right\rangle$, $\left\langle i \right| \left| j \right\rangle = 0$ if $i \ne j$ and $\left\langle i \right| \left| j \right\rangle = 1$ if $i = j$. Note that $\cos^2\left(\frac{\theta}{2}\right) = P(0)$ is the probability of measuring the qubit in the $\left| 0 \right\rangle$ state and $\sin^2\left(\frac{\theta}{2}\right) = P(1)$ is the probability of measuring the qubit in the $\left| 1 \right\rangle$ state. So,
#
# $$\left\langle Z \right\rangle = \left\langle \psi \right| Z \left|\psi \right\rangle = P(0) - P(1) = 2 P(0) - 1$$
#
# and by taking into account that $\theta = 2 \arccos\left( \sqrt{ P(0)} \right)$:
#
# $$\theta =  2 \arccos\left( \sqrt{ \frac{\left\langle Z \right\rangle + 1}{2}} \right)$$
#
# As can be seen with $\left\langle Z \right\rangle = 2 P(0) - 1$ we can estimate the $\theta$. Therefore, half of our work can be completed by estimating the $\left\langle Z \right\rangle$ or equivalently by estimating the $P(0)$. But how to estimate $P(0)$ experimentally with the quantum computer? For that, we just need to re-execute the circuit $N$ times and count how many times we had measured $\left| 0 \right\rangle$ state ($N_0$). $\frac{N_0}{N}$ will be our approximation for $P(0)$:
#
# $$P(0) = \lim_{N \rightarrow \infty} \frac{N_0}{N}$$
#
# By increasing $N$ we can improve our estimation for $P(0)$ and thus improve our estimation for $\theta$.
#
# For estimating $\varphi$ we will need to calculate $\left\langle X \right\rangle$ and $\left\langle Y \right\rangle$ expectation values. Before proceeding, let's write the code that estimates the expectation value of the $Z$ operator. Note that we will use the same function for $\left\langle X \right\rangle$ and $\left\langle Y \right\rangle$, with some preprocessing. 
#
# First of all, here are the libraries that we are going to use:
#
# [1]: https://www.cambridge.org/am/academic/subjects/physics/quantum-physics-quantum-information-and-quantum-computation/quantum-computation-and-quantum-information-10th-anniversary-edition?format=HB
# [2]: https://en.wikipedia.org/wiki/Purity_(quantum_mechanics)
# [3]: https://en.wikipedia.org/wiki/Swap_test
# [4]: https://quantumcomputing.stackexchange.com/a/13055/9459
# [5]: https://quantumcomputing.stackexchange.com/a/14401/9459
# [6]: https://github.com/DavitKhach/quantum-algorithms-tutorials
# [7]: https://nbviewer.jupyter.org/github/DavitKhach/quantum-algorithms-tutorials/blob/master/bloch_angles_estimation.ipynb

# %%
from qiskit import *
import numpy as np
from random import random
import warnings


# %% [markdown]
# The function for $\left\langle Z \right\rangle$

# %%
def z_expectation_from_counts(counts, index, shots):
    """
    Calculate Z expectation value for one qubit. 
                    e.g. <Z>, <IZ>, <ZII>, <IZI>, <ZI>.
    
    :param shots:  The number of experiments/shots
    :param index:  index of the corresponding bit in the counts.
    :param counts: dict {'00': 358, '01': 311, '10': 109, '11': 246}
    :return:       Z expectation value <Z> = P(0) - P(1) = 2P(0) - 1, 
                   where P(i) is the probability for the |i> state
    """

    probability_of_0 = 0

    for key in counts.keys():
        if key[-index - 1] == '0':
            probability_of_0 += counts[key] / shots

    return 2 * probability_of_0 - 1


# %% [markdown]
# Now let's focus on $\left\langle X \right\rangle$ and $\left\langle Y \right\rangle$ and try to estimate $\varphi$ from them. The expectation value for $X$ operator:
#
# \begin{align*}
# \left\langle \psi \right| X \left| \psi \right\rangle &= \left( \cos\left(\frac{\theta}{2}\right) \left\langle 0 \right| + e^{-i\varphi} \sin\left(\frac{\theta}{2}\right) \left\langle 1 \right| \right) X \left( \cos\left(\frac{\theta}{2}\right) \left| 0 \right\rangle + e^{i\varphi} \sin\left(\frac{\theta}{2}\right) \left| 1 \right\rangle \right) =
# \\
# &= \frac{e^{i \varphi} + e^{-i \varphi}}{2} \sin{\theta} = \cos{\varphi} \sin{\theta}
# \end{align*}
#
# where we took into account that $X \left| 0 \right\rangle = \left| 1 \right\rangle$ and $X \left| 1 \right\rangle = \left| 0 \right\rangle$. As one can see, $\left\langle X \right\rangle$, in contrast to the $\left\langle Z \right\rangle$, has  dependence on $\varphi$ and we can obtain some information about $\varphi$ from it. For full estimation we still need the expectation value of the $Y$ operator:
#
# \begin{align*}
# \left\langle \psi \right| Y \left| \psi \right\rangle &= \left( \cos\left(\frac{\theta}{2}\right) \left\langle 0 \right| + e^{-i\varphi} \sin\left(\frac{\theta}{2}\right) \left\langle 1 \right| \right) Y \left( \cos\left(\frac{\theta}{2}\right) \left| 0 \right\rangle + e^{i\varphi} \sin\left(\frac{\theta}{2}\right) \left| 1 \right\rangle \right) =
# \\
# &= -i \frac{e^{i \varphi} - e^{-i \varphi}}{2} \sin{\theta} = \sin{\varphi} \sin{\theta}
# \end{align*}
#
# where we took into account that $Y \left| 0 \right\rangle = i\left| 1 \right\rangle$ and $Y \left| 1 \right\rangle = -i\left| 0 \right\rangle$. So, we have two equations:
#
# $$
# \begin{cases}
# \cos{\varphi} = \frac{\left\langle X \right\rangle}{\sin{\theta}} \\
# \sin{\varphi} = \frac{\left\langle Y \right\rangle}{\sin{\theta}}
# \end{cases}
# $$
#
# How to find $\varphi \in [0,2 \pi)$ from these two equations? Firstly we should take into account that the range of the usual principal value of the $\arccos$ function is in $[0,\pi]$. Also, if $\sin{\varphi} \geq 0$, then $0 \leq \varphi \leq \pi$ and if $\sin{\varphi} < 0$, then $\pi < \varphi < 2\pi$. Therefore:
#
# $$
# \begin{cases}
# \varphi = \arccos{\frac{\left\langle X \right\rangle}{\sin{\theta}}}, \qquad \quad \text{if} \quad  \sin{\varphi} \geq 0\\
# \varphi = 2\pi - \arccos{\frac{\left\langle X \right\rangle}{\sin{\theta}}}, \quad \text{if}\quad  \sin{\varphi} < 0
# \end{cases}
# $$
#
# So, for estimating $\theta$ we only need $\left\langle Z \right\rangle$, but for $\varphi$ we will need $\theta$, $\left\langle X \right\rangle$ and $\left\langle Y \right\rangle$. We already know how to estimate $\left\langle Z \right\rangle$ and now we are going to discuss the experimental estimation procedures for $\left\langle X \right\rangle$ and $\left\langle Y \right\rangle$. These procedures are very similar to the one that we used for $\left\langle Z \right\rangle$. Moreover, we will not write separate functions for $\left\langle X \right\rangle$ and $\left\langle Y \right\rangle$ operators and we will use the same `z_expectation_from_counts` function for them. 
#
# Note that $X = HZH$, where $H = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$ is the Hadamard gate. That is why
#
# $$\left\langle \psi \right| X \left| \psi \right\rangle = \left\langle \psi \right| HZH \left| \psi \right\rangle = \left\langle \psi' \right| Z \left| \psi' \right\rangle$$
#
# where $\left| \psi' \right\rangle = H \left| \psi \right\rangle$ and $\left\langle \psi' \right| = \left( H \left| \psi \right\rangle \right)^\dagger = \left\langle \psi \right| H$, because $H^{\dagger} = H$. In other words, $\left\langle \psi \right| X \left| \psi \right\rangle$ calculated for $\left| \psi \right\rangle$ is equal to $\left\langle \psi' \right| Z \left| \psi' \right\rangle$ calculated for $\left| \psi' \right\rangle$. Therefore we can calculate the $X$ expectation value by calculating the $Z$ expectation value after applying the $H$ gate in order to obtain $\left| \psi' \right\rangle$.
#
# A similar thing can be done for the $Y$ expectation value. Here, we should take into account that $Y = H_y Z H_y$, where  $H_y = \frac{1}{\sqrt{2}}\begin{pmatrix} 1 & -i \\ i & -1 \end{pmatrix}$:
#
# $$\left\langle \psi \right| Y \left| \psi \right\rangle = \left\langle \psi \right| H_y Z H_y \left| \psi \right\rangle = \left\langle \psi'' \right| Z \left| \psi'' \right\rangle$$
#
# where $\left| \psi'' \right\rangle = H_y \left| \psi \right\rangle$ and $\left\langle \psi'' \right| = \left( H_y \left| \psi \right\rangle \right)^\dagger = \left\langle \psi \right| H_y$, because $H_y^{\dagger} = H_y$. So, after applying $H_y$ to the qubit, we should calculate the $\left\langle \psi'' \right| Z \left| \psi'' \right\rangle$ that will give as the desired $\left\langle \psi \right| Y \left| \psi \right\rangle$.
#
# Now let's write a function that will calculate the expectation value for a given Pauli operator ($I$ is not included):
#
# [1]: https://www.cambridge.org/am/academic/subjects/physics/quantum-physics-quantum-information-and-quantum-computation/quantum-computation-and-quantum-information-10th-anniversary-edition?format=HB
# [2]: https://en.wikipedia.org/wiki/Purity_(quantum_mechanics)
# [3]: https://en.wikipedia.org/wiki/Swap_test
# [4]: https://quantumcomputing.stackexchange.com/a/13055/9459
# [5]: https://quantumcomputing.stackexchange.com/a/14401/9459
# [6]: https://github.com/DavitKhach/quantum-algorithms-tutorials
# [7]: https://nbviewer.jupyter.org/github/DavitKhach/quantum-algorithms-tutorials/blob/master/bloch_angles_estimation.ipynb

# %%
def one_qubit_pauli_expectation_value(qubit, circuit, backend, 
                                      shots, pauli_operator):
    """
    Expectation value of X or Y or Z operators
    
    :param pauli_operator: The given operators:
                    'x', 'y','z' or 'X', 'Y','Z'
    :param qubit:   the qubit whose Bloch angles we 
                    should estimate
    :param circuit: QuantumCircuit that prepares the 
                    qubit in a pure quantum state
    :param backend: Quantum Hardware or Quantum simulator
    :param shots:   the number of experiments/shots
    :return:        the expectation value
    """
    if len(circuit.cregs) != 0:
        raise NotImplementedError("The circuit should not have" 
                                  "classical registers, because" 
                                  "that case is not implemented.")

    classical_register = ClassicalRegister(1)
    circuit_expectation = QuantumCircuit(classical_register)
    
    # add the qregs of the circuit to the circuit_expectation
    for qreg in circuit.qregs:
        circuit_expectation.add_register(qreg)
    
    circuit_expectation += circuit

    if pauli_operator == 'x' or pauli_operator == 'X':
        circuit_expectation.h(qubit)
    elif pauli_operator == 'y' or pauli_operator == 'Y':
        circuit_expectation.u(np.pi / 2, np.pi / 2, np.pi / 2, qubit)  # H_y = UGate(np.pi / 2, np.pi / 2, np.pi / 2)
    elif pauli_operator != 'z' and pauli_operator != 'Z':
        # nothing should be done for the <z> case
        raise ValueError(f"Pauli operator should be equal" 
                         f"'x', 'y','z' or 'X', 'Y','Z'." 
                         f"It was given {pauli_operator}")

    circuit_expectation.measure(qubit, classical_register[0])
    counts = execute(circuit_expectation, 
                     backend, shots=shots).result().get_counts()

    return z_expectation_from_counts(counts, 0, shots)


# %% [markdown]
# We actually ready to write the main function for estimating the Bloch angles, but before that, we should note something. Here we are not using a generalization of the Bloch sphere for the mixed states. In order words, we assume that the given state is one qubit pure state. So, for example, if we have a Bell state $\frac{1}{\sqrt{2}}\left(\left|00\right\rangle + \left|11\right\rangle \right)$ and we want to run this procedure for estimating the Bloch angles for the first qubit the procedure should warn us that we have given not appropriate quantum state (the qubit is not in the one qubit pure state). Or if there exists some probability of quantum errors, the code should warn us about that as well. These cases can be checked with this condition:
#
# $$p = \frac{1 + \left\langle X \right\rangle^2 + \left\langle Y \right\rangle^2 + \left\langle Z \right\rangle^2}{2} \approx 1$$
#
# where $p = Tr{\rho^2}$ is the [purity][2], $\rho$ is the [density matrix][1] of the quantum state (for more look at the end). It can be proved that, if $p < 0$, then we don't have a pure state. We are not going to prove this. All we need from this is to write few lines of code that will check if $p \approx 1$. Otherwise, the quantum state is not a pure state. We are checking the approximate equality and not exact equality, because estimation of $P(0) = \lim_{N \rightarrow \infty} \frac{N_0}{N}$ is not exact in the first place and hence the accuracy of the whole procedure is dependent on the number of the experiments $N$.
#
# Here is the main function that is going to estimate $\theta$ and $\varphi$ Bloch angles:
#
# [1]: https://www.cambridge.org/am/academic/subjects/physics/quantum-physics-quantum-information-and-quantum-computation/quantum-computation-and-quantum-information-10th-anniversary-edition?format=HB
# [2]: https://en.wikipedia.org/wiki/Purity_(quantum_mechanics)
# [3]: https://en.wikipedia.org/wiki/Swap_test
# [4]: https://quantumcomputing.stackexchange.com/a/13055/9459
# [5]: https://quantumcomputing.stackexchange.com/a/14401/9459
# [6]: https://github.com/DavitKhach/quantum-algorithms-tutorials
# [7]: https://nbviewer.jupyter.org/github/DavitKhach/quantum-algorithms-tutorials/blob/master/bloch_angles_estimation.ipynb

# %%
def estimate_bloch_angles(qubit, circuit, backend, shots):
    """
    Estimates the Bloch angles theta and phi:
    
    |psi > = cos(theta / 2 ) |0> + e^{i phi} cos(theta / 2 ) |1>
    
    :param shots:   the number of experiments/shots
    :param backend: Quantum Hardware or Quantum simulator
    :param qubit:   the qubit whose Bloch angles we 
                    should estimate
    :param circuit: QuantumCircuit that prepares the qubit
                    in some state
    :return:        (theta, phi) tuple
    """

    x_expectation_value = one_qubit_pauli_expectation_value(qubit, 
                           circuit, backend, shots, "X")
    y_expectation_value = one_qubit_pauli_expectation_value(qubit, 
                           circuit, backend, shots, "Y")
    z_expectation_value = one_qubit_pauli_expectation_value(qubit, 
                           circuit, backend, shots, "Z")

    # check if it is a pure state
    purity = (1 + x_expectation_value**2 + \
              y_expectation_value**2 + z_expectation_value**2) / 2
    if not np.isclose(purity, 1, rtol=0, atol=1e-1):
        warnings.warn(f"For not pure one qubit states Bloch "
                      f"angles are not defined. The purity is "
                      f"equal to {purity}. The resulted "
                      f"estimations are not valid.")

    theta = 2 * np.arccos(np.sqrt((1 + z_expectation_value) / 2))
    # |0> or |1> state cases
    if np.isclose(theta, 0, rtol=0, atol=1e-2) or \
       np.isclose(theta, np.pi / 2, rtol=0, atol=1e-2):
        
        phi = 0
        return theta, phi
    
    # arccos_argument should be in [-1, 1]
    arccos_argument = x_expectation_value / np.sin(theta)
    if np.isclose(arccos_argument, 1, rtol=0, atol=1e-2):
        arccos_argument = 1
    elif np.isclose(arccos_argument, -1, rtol=0, atol=1e-2):
        arccos_argument = -1
    elif abs(arccos_argument) > 1:
        raise ValueError("Value error for arccos, try to "
                         "increase the shots in order to "
                         "improve estimation.")
    
    if (y_expectation_value / np.sin(theta)) > 0:  # = sin(phi) > 0
        phi = np.arccos(arccos_argument)
    elif (y_expectation_value / np.sin(theta)) < 0:  # = sin(phi) < 0
        phi = -np.arccos(x_expectation_value / np.sin(theta)) \
              + 2 * np.pi

    return theta, phi


# %% [markdown]
# Also, note that this procedure can be done slightly differently [[4]] (this is the QCSE answer of the author from which this tutorial is inspired) by taking into account that instead of $\langle Y \rangle$ we can only estimate the $sign \left(\left\langle Y \right\rangle \right)$ and in some cases, it will reduce the number of needed measurements. Nevertheless, let's apply our written function for a randomly generated quantum state.
#
# [1]: https://www.cambridge.org/am/academic/subjects/physics/quantum-physics-quantum-information-and-quantum-computation/quantum-computation-and-quantum-information-10th-anniversary-edition?format=HB
# [2]: https://en.wikipedia.org/wiki/Purity_(quantum_mechanics)
# [3]: https://en.wikipedia.org/wiki/Swap_test
# [4]: https://quantumcomputing.stackexchange.com/a/13055/9459
# [5]: https://quantumcomputing.stackexchange.com/a/14401/9459
# [6]: https://github.com/DavitKhach/quantum-algorithms-tutorials
# [7]: https://nbviewer.jupyter.org/github/DavitKhach/quantum-algorithms-tutorials/blob/master/bloch_angles_estimation.ipynb

# %%
backend = BasicAer.get_backend('qasm_simulator')
shots = 8192

bloch_theta = np.pi * random()
bloch_phi = 2 * np.pi * random()

qubit = QuantumRegister(1)
circuit_one_qubit = QuantumCircuit(qubit)

circuit_one_qubit.u(bloch_theta, bloch_phi, 0, qubit[0]) # creates cos(theta/2)|0> + e^{i phi}sin(theta/2)|1>

estimated_theta, estimated_phi = estimate_bloch_angles(qubit[0], 
                                  circuit_one_qubit, backend, shots)

print("theta = ", bloch_theta)
print("The estimated theta = ", estimated_theta)

print("phi = ", bloch_phi)
print("The estimated phi = ", estimated_phi)

# %% [markdown]
# Let's try to estimate the Bloch angles when the qubit is not in a pure state and see if the function will warn us about that. We will create one of the Bell states $\frac{1}{\sqrt{2}}\left(\left| 00 \right\rangle + \left| 11 \right\rangle \right)$ and try to run the code for the first qubit.

# %%
bloch_theta = np.pi * random()
bloch_phi = 2 * np.pi * random()

quantum_register = QuantumRegister(2)
circuit_bell = QuantumCircuit(quantum_register)

# create the Bell state
circuit_bell.h(quantum_register[0])
circuit_bell.cx(quantum_register[0], quantum_register[1])

estimated_theta, estimated_phi = estimate_bloch_angles(quantum_register[0], 
                                  circuit_bell, backend, shots=8192)

print("The estimated theta = ", estimated_theta)
print("The estimated phi = ", estimated_phi)


# %% [markdown]
# If you see the warning message then the code does what we wanted from it!
#
# <h3 align="left">
# The Swap test
# </h3>
#
# Now let's introduce one more interesting procedure that will help us to double-check the main procedure. In the first example, we chose Bloch angles randomly and estimated them with the `estimate_bloch_angles` function. At the end, we had the estimated Bloch angles and the true Bloch angles and it was easy to test the estimation procedure just by comparing them. What if we don't know the true Bloch angles and the comparison cannot be made. In this case, we can implement the [Swap test][3] method to determine the absolute value of the inner product between the input quantum state $\left|\psi_{in} \right\rangle$ and the created quantum state with the estimated Bloch angles $\left|\psi_{est} \right\rangle$. In other words, the output of the Swap test is equal to $\left|\left\langle \psi_{in} \right|\left| \psi_{est} \right\rangle \right|^2$, that is equal to $1$ if the quantum states are the same (or differ only by a global phase factor $\left| \psi_{est} \right\rangle = e^{i\phi}\left|\psi_{in} \right\rangle$). So, after running the Swap test and checking if $\left|\left\langle \psi_{in} \right|\left| \psi_{est} \right\rangle \right|^2 \approx 1$, we will know that the Bloch angles estimation procedure works.
#
# The Swap test for one qubit states requires $3$ qubits (for $n$ qubit states we will need $2 n + 1$ qubits). The first qubit is an auxiliary qubit on which we will apply measurement that eventually will estimate $\left|\left\langle \psi_{1} \right|\left| \psi_{2} \right\rangle \right|^2$. The other two qubits are in the given $\left| \psi_{1} \right\rangle$ and $\left| \psi_{2} \right\rangle$ (in our case $\left|\psi_{in} \right\rangle$ and $\left|\psi_{est} \right\rangle$) quantum states for which we want to calculate $\left|\left\langle \psi_{1} \right|\left| \psi_{2} \right\rangle \right|^2$. The circuit looks like this:
#
# <img src="images/SWAP_test.png" alt="drawing" width="400"/>
#
#
# where $\left| \Psi_{1} \right\rangle$ is the initial quantum state and $\left| \Psi_{2} \right\rangle$ is the final quantum state before measurement. Here we are going to prove that the expectation value of the $IIZ$ operator will give us the desired $\left|\left\langle \psi_{1} \right|\left| \psi_{2} \right\rangle \right|^2$. In other words, we should prove that:
#
# $$\left\langle \Psi_2 \right| IIZ \left| \Psi_2 \right\rangle = P(0) - P(1) = \left|\left\langle \psi_{1} \right|\left| \psi_{2} \right\rangle \right|^2$$
#
# where the $P(0)$ ($P(1)$) is the probability of measuring the auxiliary qubit in the $\left| 0 \right\rangle$ ($\left| 1 \right\rangle$) state.
#
# First of all, let's write down the initial quantum state:
#
# $$\left| \Psi_{1} \right\rangle = \left| \psi_{2} \right\rangle \left| \psi_{1} \right\rangle \left| 0 \right\rangle$$
#
# After the first Hadamard gate:
#
# $$ IIH \cdot \left| \Psi_{1} \right\rangle = \left| \psi_{2} \right\rangle \left| \psi_{1} \right\rangle \left| + \right\rangle = \frac{1}{\sqrt{2}} \left(\left| \psi_{2} \right\rangle \left| \psi_{1} \right\rangle \left| 0 \right\rangle + \left| \psi_{2} \right\rangle \left| \psi_{1} \right\rangle \left| 1 \right\rangle \right)$$
#
# Now let's apply the CSWAP gate:
#
# $$\text{CSWAP} \frac{1}{\sqrt{2}} \left(\left| \psi_{2} \right\rangle \left| \psi_{1} \right\rangle \left| 0 \right\rangle + \left| \psi_{2} \right\rangle \left| \psi_{1} \right\rangle \left| 1 \right\rangle \right) = \frac{1}{\sqrt{2}} \left(\left| \psi_{2} \right\rangle \left| \psi_{1} \right\rangle \left| 0 \right\rangle + \left| \psi_{1} \right\rangle \left| \psi_{2} \right\rangle \left| 1 \right\rangle \right)$$
#
# The CSWAP gate has changed (swapped) the quantum states when the control qubit was in the $\left| 1 \right\rangle$ state. Finally Let's apply the last Hadamard gate:
#
# $$IIH \cdot \frac{1}{\sqrt{2}} \left(\left| \psi_{2} \right\rangle \left| \psi_{1} \right\rangle \left| 0 \right\rangle + \left| \psi_{1} \right\rangle \left| \psi_{2} \right\rangle \left| 1 \right\rangle \right) = \frac{1}{2} \left( \left( \left| \psi_{2} \right\rangle \left| \psi_{1} \right\rangle + \left| \psi_{1} \right\rangle \left| \psi_{2} \right\rangle \right) \left| 0 \right\rangle + \left( \left| \psi_{2} \right\rangle \left| \psi_{1} \right\rangle - \left| \psi_{1} \right\rangle \left| \psi_{2} \right\rangle \right) \left| 1 \right\rangle \right)$$
#
# So, the final quantum state before the measurement:
#
# $$\left| \Psi_{2} \right\rangle  = \frac{1}{2} \left( \left( \left| \psi_{2} \right\rangle \left| \psi_{1} \right\rangle + \left| \psi_{1} \right\rangle \left| \psi_{2} \right\rangle \right) \left| 0 \right\rangle + \left( \left| \psi_{2} \right\rangle \left| \psi_{1} \right\rangle - \left| \psi_{1} \right\rangle \left| \psi_{2} \right\rangle \right) \left| 1 \right\rangle \right)$$
#
# Now we can calculate the expectation value of the $IIZ$ operator for $\left| \Psi_{2} \right\rangle$ state:
#
# $$\left\langle \Psi_2 \right| IIZ \left| \Psi_2 \right\rangle = \frac{1}{4} \left(\left| \left| \psi_{2} \right\rangle \left| \psi_{1} \right\rangle + \left| \psi_{1} \right\rangle \left| \psi_{2} \right\rangle \right|^2 -  \left| \left| \psi_{2} \right\rangle \left| \psi_{1} \right\rangle - \left| \psi_{1} \right\rangle \left| \psi_{2} \right\rangle \right|^2\right) = \left|\left\langle \psi_{1} \right|\left| \psi_{2} \right\rangle \right|^2$$
#
# Note that associated with the $IIZ$ operator we have two projector operators: $Pr_0 = II\left|0 \right\rangle \left\langle 0 \right|$ and $Pr_1 = II\left|1 \right\rangle \left\langle 1 \right|$. According to the definition of the projective measurement (page 87 of the [[1]]) the probabilities of measuring $\left|0 \right\rangle$ or $\left|1\right\rangle$ can be found the following way:
#
# $$
# P(0) = \left\langle \Psi_2 \right| Pr_0 \left| \Psi_2 \right\rangle  = \frac{1}{4} \left| \left| \psi_{2} \right\rangle \left| \psi_{1} \right\rangle + \left| \psi_{1} \right\rangle \left| \psi_{2} \right\rangle \right|^2\\
# P(1) = \left\langle \Psi_2 \right| Pr_1 \left| \Psi_2 \right\rangle = \frac{1}{4} \left| \left| \psi_{2} \right\rangle \left| \psi_{1} \right\rangle - \left| \psi_{1} \right\rangle \left| \psi_{2} \right\rangle \right|^2
# $$
#
# With this, we prove that:
#
# $$\frac{1}{4} \left(\left| \left| \psi_{2} \right\rangle \left| \psi_{1} \right\rangle + \left| \psi_{1} \right\rangle \left| \psi_{2} \right\rangle \right|^2 -  \left| \left| \psi_{2} \right\rangle \left| \psi_{1} \right\rangle - \left| \psi_{1} \right\rangle \left| \psi_{2} \right\rangle \right|^2\right) = P(0) - P(1)$$
#
# But this was equal also to the $\left|\left\langle \psi_{1} \right|\left| \psi_{2} \right\rangle \right|^2$. So:
#
# $$\left|\left\langle \psi_{1} \right|\left| \psi_{2} \right\rangle \right|^2 = P(0) - P(1) = 2 P(0) - 1$$
#
# Hence, we can use the same `z_expectation_from_counts` in order to calculate the expectation value for the $IIZ$ operator and consequently the desired square of the absolute value of the inner product between the quantum states $\left|\left\langle \psi_{1} \right|\left| \psi_{2} \right\rangle \right|^2$.
#
# Here is the function that implements the Swap test procedure:
#
#
# [1]: https://www.cambridge.org/am/academic/subjects/physics/quantum-physics-quantum-information-and-quantum-computation/quantum-computation-and-quantum-information-10th-anniversary-edition?format=HB
# [2]: https://en.wikipedia.org/wiki/Purity_(quantum_mechanics)
# [3]: https://en.wikipedia.org/wiki/Swap_test
# [4]: https://quantumcomputing.stackexchange.com/a/13055/9459
# [5]: https://quantumcomputing.stackexchange.com/a/14401/9459
# [6]: https://github.com/DavitKhach/quantum-algorithms-tutorials
# [7]: https://nbviewer.jupyter.org/github/DavitKhach/quantum-algorithms-tutorials/blob/master/bloch_angles_estimation.ipynb

# %%
def one_qubit_state_swap_test(qubit_1, qubit_2, auxiliary_qubit, circuit, backend, shots):
    """
    The Swap test on two one qubit states
    
    :param qubit_1: the first qubit in the |psi_1> (pure) state
    :param qubit_2: the second qubit in the |psi_2> (pure) state
    :param auxiliary_qubit: the auxiliary qubit for Swap test
    :param circuit: The circuit where all mentioned above 
                    qubits are defined
    :return:        returns |<psi_1|psi_2>|^2
    """
    if len(circuit.cregs) != 0:
        raise NotImplementedError("The circuit should not have "
                                  "classical registers, because "
                                  "that case is not implemented.")

    classical_register = ClassicalRegister(1)
    swap_circuit = QuantumCircuit(classical_register)

    for qreg in circuit.qregs:
        swap_circuit.add_register(qreg)

    swap_circuit += circuit
    
    # the main part
    swap_circuit.h(auxiliary_qubit)
    swap_circuit.cswap(auxiliary_qubit, qubit_1, qubit_2)
    swap_circuit.h(auxiliary_qubit)
    swap_circuit.measure(auxiliary_qubit, classical_register[0])

    counts = execute(swap_circuit, 
                     backend, shots=shots).result().get_counts()

    return z_expectation_from_counts(counts, 0, shots)


# %% [markdown]
# Now let's apply a sequence of random gates to a qubit and measure its Bloch angles. After finding the Bloch angles we will use the same state obtained with the same gate sequence and compare its state with another qubit that is going to be initialized by the estimated Bloch angles. Random gates will not allow us straightforwardly estimate the Bloch angles and hence the Bloch angles are not given to us in this problem as it was previously. We can either multiply the matrices that correspond to the gates to each other and try to estimate mathematically the Bloch angles or we can do it in a "quantum way" with the Swap test. If Bloch angles estimation is correct then the output of the Swap test should be $\approx 1$.

# %%
random_gate_number = 20

main_qubit = QuantumRegister(1)
test_qubit = QuantumRegister(1)
auxiliary_qubit = QuantumRegister(1)
circuit_random_gates = QuantumCircuit(main_qubit,
                        test_qubit, auxiliary_qubit)

for _ in range(random_gate_number):
    (param1, param2, param3) = (2 * np.pi * random(), 
                                2 * np.pi * random(), 
                                2 * np.pi * random())
    circuit_random_gates.u(param1, param2, param3, main_qubit)

estimated_theta, estimated_phi = estimate_bloch_angles(main_qubit[0], 
                                  circuit_random_gates, backend, shots)

# recreate the state
circuit_random_gates.u(estimated_theta, estimated_phi, 0, 
                         test_qubit[0])
swap_test_output = one_qubit_state_swap_test(main_qubit[0], 
                         test_qubit[0], auxiliary_qubit[0],
                         circuit_random_gates, backend, shots)

if np.isclose(swap_test_output, 1, rtol=0, atol=1e-1):
    print(f"The estimated theta and phi angles are correct. " 
          f"The output of the Swap test is equal to "
          f"|<psi_in|psi_est>|^2 = {swap_test_output}")
else:
    print(f"The estimated theta and phi angles are not correct. "
          f"The output of the Swap test is equal to "
          f"|<psi_in|psi_est>|^2 = {swap_test_output}")

# %% [markdown]
# <h3 align="left">
# The density matrix
# </h3>
#
# As was said this procedure works for one qubit pure state. For not pure states we should use density matrices. If the reader is not familiar with the density matrix then this part of the tutorial can be skipped. The goal of this tutorial is not to dive deep into explanations about density matrices, but here is a quick introduction [[5]] (for more info [[1]]):
#
# If we have a pure state:
#
# $$\left| \psi \right\rangle = \alpha \left| 0 \right\rangle + \beta \left| 1 \right\rangle$$
#
# The corresponding density matrix is:
#
# $$\rho = \left| \psi \right\rangle \langle \psi| = (\alpha \left| 0 \right\rangle + \beta \left| 1 \right\rangle)(\alpha^* \left\langle 0 \right| + \beta^* \left\langle 1 \right|) = \\= |\alpha|^2\left| 0 \right\rangle \langle 0| + \alpha \beta^*\left| 0 \right\rangle \left\langle 1 \right| + \alpha^* \beta|1 \rangle \langle 0|+ |\beta|^2 |1 \rangle \left\langle 1 \right| = 
# \begin{pmatrix} 
# |\alpha|^2 & \alpha \beta^* \\ 
# \alpha^* \beta & |\beta|^2 
# \end{pmatrix}$$
#
# $\rho$ has the same descriptive abilities for the quantum state as the $\left|\psi \right\rangle$. If the state is not pure and let's say with $p_1$ probability we have $\left| \psi_1 \right\rangle$ state and with $p_2$ probability we have $\left| \psi_2 \right\rangle$ state (and this can be extended to $n$ cases: $p_n$ probability for $\left| \psi_n \right\rangle$ state):
#
# \begin{equation*}
# \rho = p_1 \left| \psi_1 \right\rangle \left\langle \psi_1 \right| + p_2 \left| \psi_2 \right\rangle \left\langle \psi_2 \right| = 
# \\
# = p_1 
# \begin{pmatrix} 
# |\alpha_1|^2 & \alpha_1 \beta_1^* \\ 
# \alpha_1^* \beta_1 & |\beta_1|^2 
# \end{pmatrix} + p_2 
# \begin{pmatrix} 
# |\alpha_2|^2 & \alpha_2 \beta_2^* \\ 
# \alpha_2^* \beta_2 & |\beta_2|^2 \end{pmatrix} = 
# \\
# =\begin{pmatrix} 
# p_1|\alpha_1|^2 + p_2|\alpha_2|^2 & p_1 \alpha_1 \beta_1^* + p_2 \alpha_2 \beta_2^* \\ 
# p_1 \alpha_1^* \beta_1 + p_2 \alpha_2^* \beta_2 & p_1|\beta_1|^2 + p_2 |\beta_2|^2
# \end{pmatrix}
# \end{equation*}
#
# This is the density matrix that describes the quantum state of one qubit with the given probabilities mentioned above. For this case, we can't attribute to the quantum state one $\left|\psi \right\rangle$, but we can define a density matrix that will give us all needed measurement properties for the quantum state. In this case, when we don't have a pure state we should do the quantum state tomography for estimating the density matrix. Fortunately, it turns out that for the quantum state tomography all we need is the same set of expectation value measurements ($\langle X \rangle$, $\langle Y \rangle$ and $\langle Z \rangle$). The entries $\rho_{ij}$ of a one qubit density matrix can be expressed with the expectation values of the Pauli operators:
#
# $$
# \begin{align*}
# \rho = \begin{pmatrix}
# \rho_{11} & \rho_{12} \\
# \rho_{21} & \rho_{22}
# \end{pmatrix}=
# \frac{1}{2}
# \begin{pmatrix}
# 1 - \left\langle Z \right\rangle & \left\langle X \right\rangle - i \left\langle Y \right\rangle \\
# \left\langle X \right\rangle + i \left\langle Y \right\rangle & 1 + \left\langle Z \right\rangle
# \end{pmatrix}
# \end{align*}
# $$
#
# This can be proved by taking into account that for density matrices the expectation values are equal to $\langle X \rangle = Tr(\rho X)$, $\langle Y \rangle = Tr(\rho Y)$ and $\langle Z \rangle= Tr(\rho Z) $. Also from this one can derive the formula for the purity mentioned above (the purity is defined as $p = Tr(\rho^2)$).
#
# $$p = Tr(\rho^2) = \frac{1 + \left\langle X \right\rangle^2 + \left\langle Y \right\rangle^2 + \left\langle Z \right\rangle^2}{2}$$
#
# As one can see although we have restricted ourselves with pure states it is very easy to use the same ideas to implement the `estimate_one_qubit_density_matrix` function that will generalize our approach for dealing with mixed one qubit states.
#
# [1]: https://www.cambridge.org/am/academic/subjects/physics/quantum-physics-quantum-information-and-quantum-computation/quantum-computation-and-quantum-information-10th-anniversary-edition?format=HB
# [2]: https://en.wikipedia.org/wiki/Purity_(quantum_mechanics)
# [3]: https://en.wikipedia.org/wiki/Swap_test
# [4]: https://quantumcomputing.stackexchange.com/a/13055/9459
# [5]: https://quantumcomputing.stackexchange.com/a/14401/9459
# [6]: https://github.com/DavitKhach/quantum-algorithms-tutorials
# [7]: https://nbviewer.jupyter.org/github/DavitKhach/quantum-algorithms-tutorials/blob/master/bloch_angles_estimation.ipynb

# %% [markdown]
# **[[Homepage][6]]**
#
# <h3 align="left">
# References
# </h3>
#
# [[1]] [M.A. Nielsen, I.L. Chuang, Cambridge University Press New York, "Quantum Computation and Quantum Information: 10th Anniversary Edition
# 10th" (2011)][1]
#
# [[2]] [Purity (quantum mechanics) from Wikipedia, the free encyclopedia][2]
#
# [[3]] [Swap test from Wikipedia, the free encyclopedia][3]
#
# [[4]] [The auther's QCSE answer about estimating the the Bloch angles][4]
#
# [[5]] [The auther's QCSE answer about density matrices][5]
#
# [1]: https://www.cambridge.org/am/academic/subjects/physics/quantum-physics-quantum-information-and-quantum-computation/quantum-computation-and-quantum-information-10th-anniversary-edition?format=HB
# [2]: https://en.wikipedia.org/wiki/Purity_(quantum_mechanics)
# [3]: https://en.wikipedia.org/wiki/Swap_test
# [4]: https://quantumcomputing.stackexchange.com/a/13055/9459
# [5]: https://quantumcomputing.stackexchange.com/a/14401/9459
# [6]: https://github.com/DavitKhach/quantum-algorithms-tutorials
# [7]: https://nbviewer.jupyter.org/github/DavitKhach/quantum-algorithms-tutorials/blob/master/bloch_angles_estimation.ipynb
