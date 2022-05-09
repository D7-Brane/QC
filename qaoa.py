def add_ising_term(circ, gamma, edge):
    i, j = edge
    circ.cx(i, j)
    circ.rz(gamma, j)
    circ.cx(i, j)
    pass


def add_quench_term(circ, theta, qubit):
    circ.rx(theta, qubit)
    pass


def get_phase_separation_step(G, gamma):
    n_qubits = len(G.nodes)
    phase_unitary = QuantumCircuit(n_qubits, n_qubits)
    ordered_edges = list(G.edges)
    ordered_edges.sort()
    for edge in ordered_edges:
        add_ising_term(phase_unitary, gamma, edge)
    return phase_unitary


def get_driver_step(G, theta):
    n_qubits = len(G.nodes)
    driver_unitary = QuantumCircuit(n_qubits, n_qubits)
    for qubit in range(n_qubits):
        add_quench_term(driver_unitary, theta, qubit)
    return driver_unitary


def get_QAOA_step(G, gamma, theta):
    phase_unitary = get_phase_separation_step(G, gamma)
    driver_unitary = get_driver_step(G, theta)
    return phase_unitary.compose(driver_unitary)


def get_initial_state(G):
    n_qubits = len(G.nodes)
    initial_state = QuantumCircuit(n_qubits, n_qubits)
    for qubit in range(n_qubits):
        initial_state.h(qubit)
    return initial_state


def invert_counts(counts):
    return {k[::-1]: v for k, v in counts.items()}


def objective(state, G):
    vev = 0
    for edge in G.edges:
        if state[edge[0]] != state[edge[1]]:
            vev += 1
    return -vev


def sample_vev(result, G):
    counts = invert_counts(result.get_counts())
    acc = 0
    for state in counts.keys():
        acc += objective(state, G) * counts[state]
    counts.values()
    total_counts = 0
    for count in counts.values():
        total_counts += count
    return acc / total_counts


def QAOA_circ(G, p, params):
    qc = QuantumCircuit(len(G.nodes), len(G.nodes))
    qc = qc.compose(get_initial_state(G))
    gamma = params[:p]
    theta = params[p:]
    for param in zip(gamma, theta):
        qc = qc.compose(get_QAOA_step(G, *param))
    qc.barrier(range(len(G.nodes)))
    qc.measure(range(len(G.nodes)), range(len(G.nodes)))
    return qc


def circuit_objective(G, p, params):
    backend = Aer.get_backend("qasm_simulator")
    qc = QAOA_circ(G, p, params)
    job = execute(qc, backend)
    result = job.result()
    return sample_vev(result, G)


def black_box_objective(G, p):
    backend = Aer.get_backend("qasm_simulator")

    def f(params):
        qc = QAOA_circ(G, p, params)
        job = execute(qc, backend)
        result = job.result()
        return sample_vev(result, G)

    return f
