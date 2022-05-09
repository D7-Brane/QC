def generate_basis(n):
    def standardize_state(state, length):
        off_length = length - len(state)
        if off_length > 0:
            standard_state = [0 for i in range(off_length)] + state
            return standard_state
        else:
            return state

    basis = []
    for element in range(2 ** n):
        basis.append([int(x) for x in list(bin(element)[2:])])
    generated_basis = [standardize_state(state, len(basis[-1])) for state in basis]
    return generated_basis


def log_to_phys(state):
    pos = 0
    phys_state = []
    while pos < len(state) - 1:
        phys_state.append((state[pos] + state[pos + 1] + 1) % 2)
        pos += 1
    phys_state.append((state[0] + state[-1] + 1) % 2)
    for i in range(len(state)):
        for j in range(len(state)):
            if i < j and (abs(i-j)>1 and abs(i-j)<len(state)-1):
                phys_state.append((state[i] + state[j] + 1) % 2)
    return phys_state


def remove_duplicates(induced_states):
    unique_states = []
    for state in induced_states:
        if state not in unique_states:
            unique_states.append(state)
    return unique_states


def sub_state(sites, state):
    reduced_state = []
    for site in sites:
        reduced_state.append(state[site])
    return reduced_state


def check_even(state):
    if state.count(0) % 2 == 0:
        return True
    else:
        return False


def generate_edges(n):
    edges = []
    for i in range(n):
        for j in range(n):
            if j > i:
                edges.append((i, j))
    return edges


def displace_edges(edges, displacement):
    displaced = []
    for edge in edges:
        displaced.append((edge[0] + displacement[0], edge[1] + displacement[1]))
    return displaced


def generate_T_constraints(n):
    T_num = n - 3
    init_edges = [(0, 1), (0, 2), (1, 2)]
    T_edges = [init_edges]
    i = 1
    for count in range(T_num):
        T_edges.append(displace_edges(init_edges, [i, i]))
        i += 1
    return T_edges


def generate_S_constraints(n):
    init_edges = [(0, 2), (1, 2), (0, 3), (1, 3)]
    S_edges = [init_edges]
    i = 1
    while i <= n - 4:
        j = 0
        while j <= i:
            S_edges.append(displace_edges(init_edges, [j, i]))
            j += 1
        i += 1
    return S_edges


def rules_edges(rules):
    edges = []
    for rule in rules:
        edges.append(rule[0])
    return edges


def generate_rules(n):
    sub_rules = []
    basis_counter = 0
    for i in range(n - 1):
        sub_rules.append([(i, i + 1), basis_counter])
        basis_counter += 1
    sub_rules.append([(n - 1, 0), basis_counter])
    basis_counter += 1
    for i in range(n):
        for j in range(n):
            if i != j:
                if ((i, j) not in rules_edges(sub_rules)) and (
                    (j, i) not in rules_edges(sub_rules)
                ):
                    sub_rules.append([(i, j), basis_counter])
                    basis_counter += 1
    return sub_rules


def get_rule_index(rules, edge):
    index = 0
    while index < len(rules):
        if set(rules[index][0]) == set(edge):
            return index
        index += 1


def to_parity_basis(list_of_edges, n):
    rules = generate_rules(n)
    indices = [get_rule_index(rules, edge) for edge in list_of_edges]
    indices.sort()
    return indices


def k_body_fix(positions, logical_basis):
    spurious_states = []
    filtered_states = []

    def reduced_state(positions, state):
        reduced_state = []
        for position in positions:
            reduced_state.append(state[position])
        return reduced_state

    for state in logical_basis:
        if reduced_state(positions, state).count(0) % 2 != 0:
            spurious_states.append(state)
        else:
            filtered_states.append(state)
    return spurious_states, filtered_states
