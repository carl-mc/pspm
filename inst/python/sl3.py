from typing import Tuple, Dict, List, cast, OrderedDict
from abc import ABC, abstractmethod
import collections as cl
import networkx as nx
import numpy as np
from scipy.special import logsumexp


class SpatialLattice:  
    # Encapsulates lattice and predictor info
    # Does not change during sampling or estimation

    N: int  # number of nodes
    g: nx.Graph  # the lattice of size N
    P: int  # number of node predictors
    X: np.ndarray  # a NxP matrix of node predictors
    L: int  # number of edge predictors
    Z: List[np.ndarray]  # list of NxN matrices with edge predictors

    def __init__(self, A: np.ndarray, Z: List[np.ndarray], X: np.ndarray):
        self.N = A.shape[0]
        self.g = nx.from_numpy_array(A)
        self.L = len(Z)
        self.Z = Z
        self.P = X.shape[1]
        self.X = X

        # Store edge predictors as edge attributes in graph
        edge_attr_dict = {}
        for i in range(self.N - 1):
            for j in range(i + 1, self.N, 1):
                key: Tuple[int, int] = (i, j)
                val: np.ndarray = np.zeros((1, self.L), dtype=float)
                for predictor_idx in range(self.L):
                    val[0, predictor_idx] = self.Z[predictor_idx][i, j]
                edge_attr_dict[key] = val
        nx.set_edge_attributes(self.g, values=edge_attr_dict, name='x')

    def get_adjacency_matrix(self) -> np.ndarray:
        # Returns the adjacency matrix for the lattice
        mat: np.ndarray = nx.adjacency_matrix(self.g)
        return mat

    def neighbors(self, node) -> any:
        return self.g.neighbors(node)


class NetworkStatistic(ABC):
    # Class for computing a network statistic on a SpatialLattice given a partitioning

    name: str  # the name of the statistic
    contributes_energy: bool  # whether this NS contributes to the energy of the partition model
    spatial_lattice: SpatialLattice  # the spatial lattice
    partitioning: List[int]  # a local copy of the partitioning
    statistic: np.ndarray  # the network statistic
    statistic_shape: Tuple[int, int]  # the shape of the network statistic (typically a Kx1 vector)

    def __init__(self, spatial_lattice: SpatialLattice, partitioning: List[int], contributes_energy: bool = True):
        self.spatial_lattice = spatial_lattice
        self.contributes_energy = contributes_energy
        self.partitioning = partitioning.copy()
        self.compute_statistic()
        self.statistic_shape = self.get_statistic_shape()
        self.name = self.get_name()

    def set_partitioning(self, partitioning: List[int]):
        # resets the partitioning and recomputes the statistic from scratch
        self.partitioning = partitioning.copy()
        self.compute_statistic()

    def update_partitioning(self, node: int, new_partition: int) -> None:
        # updates the partitioning and the statistic if node has new partition
        self.update_statistic(node, new_partition)
        self.partitioning[node] = new_partition

    @abstractmethod
    def get_name(self) -> str:
        # Returns the name of the network statistic
        pass

    @abstractmethod
    def get_statistic_shape(self) -> Tuple[int, int]:
        # Returns the shape of the network statistic (typically a Kx1 vector)
        pass

    @abstractmethod
    def compute_statistic(self) -> None:
        # computes the statistic
        pass

    @abstractmethod
    def update_statistic(self, node: int, new_partition: int) -> None:
        # updates the statistic given node has new partition
        pass

    @abstractmethod
    def get_deltas(self, node: int, candidate_partitions: List[int]) -> List[np.ndarray]:
        # gets changes in statistic if given node joins any of the candidate partitions
        pass


class LatticeNeighborsNS(NetworkStatistic):
    # NS that sums edge predictors for lattice edges within the same partition

    def __init__(self, spatial_lattice, partitioning, contributes_energy=True):
        super().__init__(spatial_lattice, partitioning, contributes_energy)

    def get_name(self) -> str:
        return 'lattice_neighbors'

    def get_statistic_shape(self) -> Tuple[int, int]:
        return self.spatial_lattice.L, 1

    def compute_statistic(self) -> None:
        self.statistic: np.ndarray = np.zeros((1, self.spatial_lattice.L))
        for u, v, d in self.spatial_lattice.g.edges(data=True):
            if self.partitioning[u] == self.partitioning[v]:
                self.statistic += d['x']

    def update_statistic(self, node: int, new_partition: int) -> None:
        delta: np.ndarray = np.zeros((1, self.spatial_lattice.L))
        current_partition: int = self.partitioning[node]
        for nb in self.spatial_lattice.g.neighbors(node):
            nb_partition: int = self.partitioning[nb]
            d = self.spatial_lattice.g.get_edge_data(node, nb)
            if nb_partition == current_partition:
                delta -= d['x']
            if nb_partition == new_partition:
                delta += d['x']
        self.statistic += delta

    def get_deltas(self, node: int, candidate_partitions: List[int]) -> List[np.ndarray]:
        n_candidates: int = len(candidate_partitions)
        deltas: List[np.ndarray] = [np.zeros((1, self.spatial_lattice.L)) for _ in range(n_candidates)]
        current_partition: int = self.partitioning[node]
        for nb in self.spatial_lattice.g.neighbors(node):
            nb_partition: int = self.partitioning[nb]
            d = self.spatial_lattice.g.get_edge_data(node, nb)
            for idx in range(n_candidates):
                if nb_partition == current_partition:
                    deltas[idx] -= d['x']
                if nb_partition == candidate_partitions[idx]:
                    deltas[idx] += d['x']
        return deltas


class NonContiguousPartitionsNS(NetworkStatistic):

    # NS that counts non-contiguous partitions
    subgraphs: Dict[int, nx.Graph]  # subgraph for each partition
    connection_status: Dict[int, bool]  # is_connected status for each partition
    articulation_points: Dict[int, List[int]]  # articulation points for each partition

    def __init__(self, spatial_lattice, partitioning, contributes_energy=True):
        super().__init__(spatial_lattice, partitioning, contributes_energy)

    def get_name(self) -> str:
        return 'non_contiguous_partitions'

    def get_statistic_shape(self) -> Tuple[int, int]:
        return 1, 1

    @staticmethod
    def _make_subgraph(g: nx.Graph, nodes: List[int]) -> nx.Graph:
        # Makes a subgraph of g based on nodes
        sg: nx.Graph = nx.Graph()
        sg.add_nodes_from((n, g.nodes[n]) for n in nodes)
        sg.add_edges_from((n, nbr)
                          for n, nbrs in g.adj.items() if n in nodes
                          for nbr, d in nbrs.items() if nbr in nodes)
        return sg

    def compute_statistic(self) -> None:
        self.statistic: np.ndarray = np.zeros((1, 1))

        # Make subgraphs
        partition_nodes: Dict[int, List[int]] = {}
        for node, partition in enumerate(self.partitioning):
            if partition not in partition_nodes.keys():
                partition_nodes[partition] = [node]
            else:
                partition_nodes[partition].append(node)
        self.subgraphs: Dict[int, nx.Graph] = {}
        for partition in partition_nodes.keys():
            self.subgraphs[partition] = self._make_subgraph(self.spatial_lattice.g, partition_nodes[partition])

        # Determine articulation points
        # (not needed for computing the statistic, but speeds up get_delta queries)
        self.articulation_points: Dict[int, List[int]] = {}
        for partition, sg in self.subgraphs.items():
            self.articulation_points[partition] = list(nx.articulation_points(sg))

        # Count non-contiguous subgraphs
        self.connection_status = {}
        non_contiguous_subgraphs: int = 0
        for partition, sg in self.subgraphs.items():
            is_connected: bool = nx.is_connected(sg)
            self.connection_status[partition] = is_connected
            if not is_connected:
                non_contiguous_subgraphs += 1
        self.statistic[0, 0] = float(non_contiguous_subgraphs)

    def update_statistic(self, node: int, new_partition: int) -> None:

        statistic_change: float = 0.0  # change in count of non-contiguous partitions

        # handle the 'old' partition (partition from which node is removed)
        old_partition: int = self.partitioning[node]
        old_subgraph: nx.Graph = self.subgraphs[old_partition]
        old_partition_was_disconnected: bool = not self.connection_status[old_partition]

        old_subgraph.remove_node(node)
        if old_subgraph.number_of_nodes() == 0:
            del self.subgraphs[old_partition]
            del self.articulation_points[old_partition]
            del self.connection_status[old_partition]
        else:
            old_partition_is_disconnected: bool = not nx.is_connected(old_subgraph)
            self.connection_status[old_partition] = not old_partition_is_disconnected
            if old_partition_was_disconnected and not old_partition_is_disconnected:
                statistic_change -= 1.0
            if not old_partition_was_disconnected and old_partition_is_disconnected:
                statistic_change += 1.0
            self.articulation_points[old_partition] = list(nx.articulation_points(old_subgraph))

        # handle the new partition (partition to which node is added)
        if new_partition in self.subgraphs.keys():  # new partition already present
            new_subgraph: nx.Graph = self.subgraphs[new_partition]
            new_partition_was_disconnected: bool = not self.connection_status[new_partition]

            # update the subgraph
            new_subgraph.add_node(node)
            for nb in self.spatial_lattice.neighbors(node):
                if new_subgraph.has_node(nb):
                    new_subgraph.add_edge(node, nb)

            # update the statistic
            new_partition_is_disconnected: bool = not nx.is_connected(new_subgraph)
            self.connection_status[new_partition] = not new_partition_is_disconnected
            if new_partition_was_disconnected and not new_partition_is_disconnected:
                statistic_change -= 1.0
            if not new_partition_was_disconnected and new_partition_is_disconnected:
                statistic_change += 1.0

            # update articulation points
            self.articulation_points[new_partition] = list(nx.articulation_points(new_subgraph))
        else:  # new partition not yet present
            new_subgraph: nx.Graph = nx.Graph()
            new_subgraph.add_node(node)
            self.subgraphs[new_partition] = new_subgraph
            self.articulation_points[new_partition] = []
            self.connection_status[new_partition] = True

        # update the statistic
        self.statistic[0, 0] = self.statistic[0, 0] + statistic_change

    def get_deltas(self, node: int, candidate_partitions: List[int]) -> List[np.ndarray]:
        old_partition: int = self.partitioning[node]
        old_partition_is_connected: bool = self.connection_status[old_partition]

        # compute delta in old partition (partition node is removed from)
        old_partition_delta: float = 0.0  # change in statistic due to change in old partition
        old_partition_articulation_points: List[int] = self.articulation_points[old_partition]
        if old_partition_is_connected:
            # removing a node from a connected partition only leads to delta if node is articulation point
            if node in old_partition_articulation_points:
                old_partition_delta += 1.0
        else:
            old_subgraph_copy: nx.Graph = self.subgraphs[old_partition].copy()
            old_subgraph_copy.remove_node(node)
            if nx.is_connected(old_subgraph_copy):
                old_partition_delta -= 1.0

        # compute delta in new partition (partition node is added to)
        candidate_deltas: List[np.ndarray] = []
        node_neighbor_partitions: List[int] = []
        for nb in self.spatial_lattice.neighbors(node):
            node_neighbor_partitions.append(self.partitioning[nb])
        for candidate_partition in candidate_partitions:
            candidate_delta: float = 0.0  # change in statistic due to change in candidate partition
            if candidate_partition == old_partition:  # no change in partitions -> move ahead in loop
                candidate_deltas.append(np.zeros((1, 1)))
                continue
            if candidate_partition not in self.subgraphs.keys():
                # new partition is contiguous
                pass
            else:
                candidate_partition_was_connected: bool = self.connection_status[candidate_partition]
                if candidate_partition_was_connected:
                    if candidate_partition not in node_neighbor_partitions:
                        # node has no neighbor with same partition -> is disconnected from remaining partition
                        candidate_delta += 1.0
                else:
                    candidate_subgraph_copy: nx.Graph = self.subgraphs[candidate_partition].copy()
                    candidate_subgraph_copy.add_node(node)
                    for nb in self.spatial_lattice.neighbors(node):
                        if candidate_subgraph_copy.has_node(nb):
                            candidate_subgraph_copy.add_edge(node, nb)
                    if nx.is_connected(candidate_subgraph_copy):
                        candidate_delta -= 1.0

            delta = np.array(old_partition_delta + candidate_delta).reshape(1, 1)
            candidate_deltas.append(delta)

        return candidate_deltas




class PartitionModel:
    # Class for sampling and estimating partitionings on a SpatialLattice object

    spatial_lattice: SpatialLattice
    partitioning: List[int]
    partition_ids: List[int]  # unique partition ids (unique elements in partitioning vector)
    network_statistics: OrderedDict[str, NetworkStatistic]
    parameters: Dict[str, np.ndarray]  # one parameter vector per network statistic
    force_contiguous: bool  # whether we force the partitioning to be contiguous
    likelihood_cache: Dict  # a cache of deltas and candidate partitions for speeding up the likelihood computation
    temperature: float  # temperature for sampling

    def __init__(self, spatial_lattice: SpatialLattice, partitioning: List[int],
                 force_contiguous: bool = False, temperature: float = 1.0,
                 net_stats: List = ["LatticeNeighborsNS"]):
        self.spatial_lattice = spatial_lattice
        self.partitioning = partitioning
        self.partition_ids = cast(List[int], np.unique(partitioning).tolist())
        self.network_statistics = cl.OrderedDict()  # ordered so we can assign/query parameters based on position
        self.force_contiguous = force_contiguous
        self.likelihood_cache = {}
        self.temperature = temperature

        # Set up network statistics()
        #  This can be expanded to add more supra-edge level stats
  
        if "LatticeNeighborsNS" in net_stats:
          lattice_neighbors_ns: LatticeNeighborsNS = LatticeNeighborsNS(spatial_lattice, partitioning)
          self.network_statistics[lattice_neighbors_ns.get_name()] = lattice_neighbors_ns
      
        # Set up force contiguous option
        if self.force_contiguous:
            # Add a network statistic that counts non-contiguous partitions, but does not contribute to the energy
            non_contiguous_ns: NonContiguousPartitionsNS = NonContiguousPartitionsNS(spatial_lattice,
                                                                                     partitioning,
                                                                                     contributes_energy=False)
            self.network_statistics[non_contiguous_ns.get_name()] = non_contiguous_ns
            # Throw error if initial partitioning is non_contiguous
            if non_contiguous_ns.statistic > 0:
                raise ValueError('Initial partitioning may not be non-contiguous when force_contiguous option is '
                                 'enabled.')

        # Initialize the parameter vectors
        self.parameters = {}
        for name, ns in self.network_statistics.items():
            if ns.contributes_energy:  # ignore network stats that don't contribute energy
                parameter_vector = np.zeros(ns.get_statistic_shape())
                self.parameters[name] = parameter_vector

    def set_partitioning(self, partitioning: List[int]) -> None:
        # Resets the partitioning and recomputes the statistics
        for ns in self.network_statistics.values():
            ns.set_partitioning(partitioning)
        self.partitioning = partitioning
        self.partition_ids = cast(List[int], np.unique(partitioning).tolist())

        # Wipe the likelihood cache
        self.likelihood_cache = {}

    def update_partitioning(self, node: int, new_partition: int) -> None:
        # Updates a single node's partition membership (and updates statistics)
        old_partition: int = self.partitioning[node]

        # update stats
        for ns in self.network_statistics.values():
            ns.update_partitioning(node, new_partition)

        # update partitioning
        self.partitioning[node] = new_partition

        # update partition ids
        if new_partition not in self.partition_ids:
            self.partition_ids.append(new_partition)
        if old_partition not in self.partitioning:
            self.partition_ids.remove(old_partition)

        # Wipe the likelihood cache
        self.likelihood_cache = {}

    def set_temperature(self, temperature: float):
        self.temperature = temperature

    def get_statistics(self) -> Dict[str, np.ndarray]:
        stats: Dict[str, np.ndarray] = {}
        for name, ns in self.network_statistics.items():
            stats[name] = ns.statistic
        return stats

    def get_energies(self) -> Dict[str, float]:
        # Computes the energy contribution of each statistic
        energies: Dict[str, float] = {}
        for name, ns in self.network_statistics.items():
            if ns.contributes_energy:
                parameter_vector: np.ndarray = self.parameters[name]
                statistic: np.ndarray = ns.statistic
                energy: float = np.matmul(statistic, parameter_vector)[0, 0]
                energies[name] = energy
        return energies

    def get_total_energy(self) -> float:
        # Computes the total energy of the partitioning
        energies: Dict[str, float] = self.get_energies()
        total_energy: float = 0.0
        for energy in energies.values():
            total_energy += energy
        return total_energy

    def get_parameter_lengths(self) -> Dict[str, int]:
        # Returns required parameter length for each network statistic
        lengths: Dict[str, int] = {}
        for name, ns in self.network_statistics.items():
            if ns.contributes_energy:
                K: int = ns.get_statistic_shape()[1]
                lengths[name] = K
        return lengths

    def update_parameters(self, new_parameters: Dict[str, np.ndarray]):
        # Updates the parameter vectors
        for name, params in new_parameters.items():
            new_parameter: np.ndarray = np.asarray(params)
            required_shape: Tuple[int, int] = self.parameters[name].shape
            new_parameter = new_parameter.reshape(required_shape)
            self.parameters[name] = new_parameter

    def update_parameters_theta(self, theta: List[float]):
        # Updates the parameter vectors using a single theta vector
        position: int = 0
        for name, ns in self.network_statistics.items():
            if ns.contributes_energy:
                K: int = self.parameters[name].shape[0]
                new_parameter: np.ndarray = np.asarray(theta[position:(position+K)]).reshape(K, 1)
                self.parameters[name] = new_parameter
                position = position + K

    def get_next_partition_id(self) -> int:
        # Returns a new (unused) partition id
        sorted_partition_ids: np.ndarray = np.sort(self.partition_ids)
        for idx, partition in enumerate(sorted_partition_ids):
            if idx < partition:
                return idx
        return np.max(sorted_partition_ids) + 1

    def is_node_isolated(self, node: int) -> bool:
        return self.partitioning.count(self.partitioning[node]) == 1

    def get_joinable_partitions(self, node: int) -> List[int]:
        # Returns partitions this node can join while only modifying incoming edges;
        # Used for conditional sampling
        joinable_partitions: List[int] = []

        if self.force_contiguous:
            # we only consider isolation partition...
            if self.is_node_isolated(node):
                joinable_partitions.append(self.partitioning[node])
            else:
                joinable_partitions.append(self.get_next_partition_id())
            # ...and partitions in neighborhood of node...
            for u in self.spatial_lattice.neighbors(node):
                neighbor_partition: int = self.partitioning[u]
                if neighbor_partition not in joinable_partitions:
                    joinable_partitions.append(self.partitioning[u])
            # ...then we remove those partitions that would lead to a non-contiguous partition
            non_contiguous_ns = self.network_statistics['non_contiguous_partitions']
            deltas: List[np.ndarray] = non_contiguous_ns.get_deltas(node, joinable_partitions)
            joinable_partitions: List[int] = [p for idx, p in enumerate(joinable_partitions)
                                              if deltas[idx][0, 0] < 1.0]
        else:
            joinable_partitions = self.partition_ids.copy()
            if not self.is_node_isolated(node):  # Add an isolation candidate
                joinable_partitions.append(self.get_next_partition_id())

        return joinable_partitions

    def get_conditional_distribution(self, node: int) -> Dict[str, np.ndarray]:
        # Computes conditional probability of node joining any partition

        # Get the candidate partitions
        candidate_partitions: List[int] = self.get_joinable_partitions(node)

        # For each candidate: Compute sum of delta energies
        n_candidates: int = len(candidate_partitions)
        delta_energies: np.ndarray = np.zeros((n_candidates,))
        for name, ns in self.network_statistics.items():
            if ns.contributes_energy:
                deltas: List[np.ndarray] = ns.get_deltas(node, candidate_partitions)
                parameter_vector: np.ndarray = self.parameters[name]
                for candidate_idx, delta in enumerate(deltas):
                    energy: float = np.matmul(delta, parameter_vector)[0, 0]
                    delta_energies[candidate_idx] -= energy

        # Make temperature adjustment
        if self.temperature != 1.0:
            delta_energies *= (1/self.temperature)

        # Compute log-probabilities
        log_probabilities: np.ndarray = delta_energies - logsumexp(delta_energies)
        
        # Return dict
        distribution: Dict[str, np.ndarray] = {'partition': np.asarray(candidate_partitions),
                                               'log_probability': log_probabilities}
        return distribution

    def sample(self, node: int) -> None:
        # Samples a new partition for node
        distribution: Dict[str, np.ndarray] = self.get_conditional_distribution(node)
        # probabilities: np.ndarray = np.exp(distribution['log_probability'])
        # partitions: np.ndarray = distribution['partition']
        new_partition: int = np.random.choice(distribution['partition'], 1, False, 
                                              np.exp(distribution['log_probability']))[0]
        if new_partition != self.partitioning[node]:
          self.update_partitioning(node, new_partition)

    def sample_all(self, M: int) -> None:
        # Samples all nodes M times
        for _ in range(M):
            for node in range(self.spatial_lattice.N):
                self.sample(node)

    def generate_likelihood_cache(self):
        # Makes a cache of deltas and joinable partitions to speed up repeated likelihood computation
        for node in range(self.spatial_lattice.N):
            # Get the candidate partitions
            candidate_partitions: List[int] = self.get_joinable_partitions(node)

            # Get the deltas for each (relevant) network statistic
            statistic_deltas: Dict[str, List[np.ndarray]] = {}
            for name, ns in self.network_statistics.items():
                if ns.contributes_energy:
                    deltas: List[np.ndarray] = ns.get_deltas(node, candidate_partitions)
                    statistic_deltas[name] = deltas

            # Add to cache
            self.likelihood_cache[node] = (np.asarray(candidate_partitions), statistic_deltas)

    def get_composite_log_likelihood(self) -> float:
        # Computes composite ll with current parameters
        ll: float = 0.0

        # Generate the cache (if it doesn't exist yet)
        if not self.likelihood_cache:
            self.generate_likelihood_cache()

        # Compute the log composite likelihood
        for node in range(self.spatial_lattice.N):

            # Compute the delta energies
            candidate_partitions, statistic_deltas = self.likelihood_cache[node]
            n_candidates: int = len(candidate_partitions)
            delta_energies: np.ndarray = np.zeros((n_candidates,))
            for name, deltas in statistic_deltas.items():
                parameter_vector: np.ndarray = self.parameters[name]
                for candidate_idx, delta in enumerate(deltas):
                    energy: float = np.matmul(delta, parameter_vector)[0, 0]
                    delta_energies[candidate_idx] -= energy

            # Compute log-probabilities
            log_probabilities: np.ndarray = delta_energies - logsumexp(delta_energies)

            # Get relevant log-probability
            observed_partition: int = self.partitioning[node]
            log_prob: float = log_probabilities[candidate_partitions == observed_partition][0]
            ll += log_prob
        return ll
