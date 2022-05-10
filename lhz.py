import hspace
import networkx as nx
from qiskit import QuantumCircuit, Aer, execute


def edges_from_constraints(constraint_list):
    edges = []
    for constraint in constraint_list:
        for edge in zip(constraint, displace_list(constraint)):
            edges.append(tuple(sorted(edge)))
        edges.sort()
    return edges


def displace_list(a):
    displaced = a[:-1]
    displaced.insert(0, a[-1])
    return displaced


def LHZ_compiled(n):
    k = int(n * (n - 1) / 2)
    G_LHZ = nx.Graph()
    G_LHZ.add_nodes_from(range(k))
    T_constraints = hspace.generate_T_constraints(n)
    S_constraints = hspace.generate_S_constraints(n)
    T_list = [hspace.to_parity_basis(T, n) for T in T_constraints]
    S_list = [hspace.to_parity_basis(S, n) for S in S_constraints]
    T_edges = edges_from_constraints(T_list)
    S_edges = edges_from_constraints(S_list)
    edges = list(set(T_edges + S_edges))
    G_LHZ.add_edges_from(edges)
    return G_LHZ, T_list, S_list


def add_plaquette_term(qc, sites, omega=3):
    connections = []
    i = 0
    while i < len(sites) - 1:
        connections.append((sites[i], sites[i + 1]))
        i += 1
    for connection in connections:
        qc.cx(*connection)
    qc.rz(omega, sites[-1])
    for connection in connections[::-1]:
        qc.cx(*connection)


def add_quench_term(qc, theta, qubit):
    qc.rx(theta, qubit)
    pass


def add_local_terms(qc, gamma, k):
    for i in range(k):
        qc.rz(gamma, i)
    pass


def get_driver_step(G, theta):
    n_qubits = len(G.nodes)
    driver_unitary = QuantumCircuit(n_qubits, n_qubits)
    for qubit in range(n_qubits):
        add_quench_term(driver_unitary, theta, qubit)
    return driver_unitary


def get_phase_separation_LHZ(G, gamma, omega=3):
    n_logical = len(G.nodes)
    n_physical = len(G.edges)
    phase_unitary = QuantumCircuit(n_physical, n_physical)
    G_LHZ, T_terms, S_terms = LHZ_compiled(n_logical)
    constraints = T_terms + S_terms
    add_local_terms(phase_unitary, gamma, n_physical)
    for c in constraints:
        add_plaquette_term(phase_unitary, c, omega)
    return phase_unitary


def LHZ_QAOA_step(G, gamma, theta, omega=3):
    n_logical = len(G.nodes)
    phase_separation_unitary = get_phase_separation_LHZ(G, gamma, omega)
    driver_unitary = get_driver_step(LHZ_compiled(n_logical)[0], theta)
    return phase_separation_unitary.compose(driver_unitary)


def get_initial_state(G):
    n_qubits = len(G.nodes)
    initial_state = QuantumCircuit(n_qubits, n_qubits)
    for qubit in range(n_qubits):
        initial_state.h(qubit)
    return initial_state


def LHZ_QAOA(G, p, params):
    n_logical = len(G.nodes)
    n_physical = int(n_logical * (n_logical - 1) / 2)
    qc = QuantumCircuit(n_physical, n_physical)
    qc = qc.compose(get_initial_state(LHZ_compiled(n_logical)[0]))
    gamma = params[:p]
    theta = params[p : 2 * p]
    for param in zip(gamma, theta):
        qc = qc.compose(LHZ_QAOA_step(G, *param))
    qc.barrier(range(n_physical))
    qc.measure(range(n_physical), range(n_physical))
    return qc
