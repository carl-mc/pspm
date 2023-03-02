
## Herfindahl Index
def hhi(series):
    _, cnt = np.unique(series, return_counts=True)
    return np.square(cnt/cnt.sum()).sum()    
  
## Weighted Herfindahl Index across all units (g)
def hhi_all(g, x):
    x = np.array(x)
    g = np.array(g)
    gu, cnt = np.unique(g, return_counts=True)
    wgt = cnt
    result = 0
    cnt = 0
    for i in gu:
      result += hhi(x[g == i]) * wgt[cnt]
      cnt += 1
    return result

  
## Cracking
class VertexCrackNS(NetworkStatistic):
    # NS that evaluates how cracked groups in 
    groups: List
    weights: List
    stat: List
    
    def __init__(self, spatial_lattice, partitioning, contributes_energy=True):
        super().__init__(spatial_lattice, partitioning, contributes_energy)

    def get_name(self) -> str:
        return 'vertex_crack'
      
    def get_statistic_shape(self) -> Tuple[int, int]:
        return self.spatial_lattice.P, 1

    def compute_statistic(self) -> None:
        self.statistic: np.ndarray = np.zeros((1, self.spatial_lattice.P))
        self.groups: List(np.array) = []
        self.weights: List(np.array) = []
        self.stat: List(np.array) = []
        for c in range(self.spatial_lattice.P):
          g = np.array(self.spatial_lattice.X[:,c])
          gu, cnt = np.unique(g, return_counts=True)
          self.groups.append(gu.astype(int))
          self.weights.append(cnt)
          result = 0
          cnt = 0
          res = []
          for i in gu:
            res.append(1 - hhi(np.array(self.partitioning)[g == i]))
          self.stat.append(np.array(res))
          self.statistic[c] = sum(self.stat[c] * self.weights[c])


    def update_statistic(self, node: int, new_partition: int) -> None:
        # delta: np.ndarray = np.zeros((1, self.spatial_lattice.L))
        new_partitioning = self.partitioning.copy()
        new_partitioning[node] = new_partition
        for c in range(self.spatial_lattice.P):
          aff_group = self.spatial_lattice.X[node,c]
          g_part = np.array(new_partitioning)[np.array(self.spatial_lattice.X[:,c]) == aff_group]
          self.stat[c][self.groups[c] == aff_group] = (1 - hhi(g_part))
          self.statistic[c] = sum(self.stat[c] * self.weights[c])

    def get_deltas(self, node: int, candidate_partitions: List[int]) -> List[np.ndarray]:
        n_candidates: int = len(candidate_partitions)
        deltas: List[np.ndarray] = [np.zeros((1, self.spatial_lattice.P)) for _ in range(n_candidates)]
        new_partitioning = self.partitioning.copy()
        for idx in range(n_candidates):
          diff: np.ndarray = np.zeros((1, self.spatial_lattice.P))
          new_partitioning[node] = candidate_partitions[idx]
          for c in range(self.spatial_lattice.P):
            aff_group = self.spatial_lattice.X[node,c]
            g_part = np.array(new_partitioning)[self.spatial_lattice.X[:,c] == aff_group]
            hhi_diff = (1 - hhi(g_part)) - self.stat[c][self.groups[c] == aff_group]
            diff[c] =  hhi_diff * self.weights[c][self.groups[c] == aff_group]
          deltas[idx] = diff
          
        return deltas

# ## Cracking -- slow & safe
# class VertexCrackSlowNS(NetworkStatistic):
#     # NS that evaluates how cracked groups in 
#     
#     def __init__(self, spatial_lattice, partitioning, contributes_energy=True):
#         super().__init__(spatial_lattice, partitioning, contributes_energy)
# 
#     def get_name(self) -> str:
#         return 'vertex_crack_slow'
#       
#     def get_statistic_shape(self) -> Tuple[int, int]:
#         return self.spatial_lattice.P, 1
# 
#     def compute_statistic(self) -> None:
#         self.statistic: np.ndarray = np.zeros((1, self.spatial_lattice.P))
#         for c in range(self.spatial_lattice.P):
#           self.statistic[c] = 1 - hhi_all(self.spatial_lattice.X[:,c],self.partitioning)
# 
#     def update_statistic(self, node: int, new_partition: int) -> None:
#         new_partitioning = self.partitioning.copy()
#         new_partitioning[node] = new_partition
#         self.statistic: np.ndarray = np.zeros((1, self.spatial_lattice.P))
#         for c in range(self.spatial_lattice.P):
#           self.statistic[c] = 1 - hhi_all(self.spatial_lattice.X[:,c],new_partitioning)
# 
#     def get_deltas(self, node: int, candidate_partitions: List[int]) -> List[np.ndarray]:
#         n_candidates: int = len(candidate_partitions)
#         deltas: List[np.ndarray] = [np.zeros((1, self.spatial_lattice.P)) for _ in range(n_candidates)]
#         new_partitioning = self.partitioning.copy()
#         for idx in range(n_candidates):
#           new_stat: np.ndarray = np.zeros((1, self.spatial_lattice.P))
#           new_partitioning[node] = candidate_partitions[idx]
#           for c in range(self.spatial_lattice.P):
#             new_stat[c] = 1 - hhi_all(self.spatial_lattice.X[:,c],new_partitioning)
#           deltas[idx] = self.statistic - new_stat
#           
#         return deltas

## Packing - fast
class VertexPackNS(NetworkStatistic):
    # NS that evaluates how cracked groups in 
    groups: List
    weights: List
    stat: List
    
    def __init__(self, spatial_lattice, partitioning, contributes_energy=True):
        super().__init__(spatial_lattice, partitioning, contributes_energy)

    def get_name(self) -> str:
        return 'vertex_pack'
      
    def get_statistic_shape(self) -> Tuple[int, int]:
        return self.spatial_lattice.P, 1

    def compute_statistic(self) -> None:
        self.statistic: np.ndarray = np.zeros((1, self.spatial_lattice.P))
        gu, cnt = np.unique(self.partitioning, return_counts=True)
        self.groups = gu
        self.weights = cnt
        self.stat: List(np.array) = []
        for c in range(self.spatial_lattice.P):
          x = np.array(self.spatial_lattice.X[:,c])
          result = 0
          cnt = 0
          res = []
          for i in gu:
            res.append(hhi(x[self.partitioning == i]))
          self.stat.append(np.array(res))
          self.statistic[c] = sum(self.stat[c] * self.weights[c])


    def update_statistic(self, node: int, new_partition: int) -> None:
        # delta: np.ndarray = np.zeros((1, self.spatial_lattice.L))
        old_partition = self.partitioning[node]
        new_partitioning = self.partitioning.copy()
        new_partitioning[node] = new_partition

        # Adjust old weight
        old_idx = np.where(self.groups == old_partition)
        self.weights[old_idx] -= 1

        # Adjust new weight
        if all(self.groups != new_partition):
          self.groups= np.append(self.groups, new_partition)
          self.weights = np.append(self.weights, int(0))
          is_new = True
        else:
          is_new = False

        new_idx = np.where(self.groups == new_partition)
        self.weights[new_idx] += 1

        for c in range(self.spatial_lattice.P):
          if is_new:
            self.stat[c]= np.append(self.stat[c], int(0))
          
          for p in [old_idx,new_idx]:
            g_part = np.array(self.spatial_lattice.X[:,c])[np.array(new_partitioning) == self.groups[p]]
            self.stat[c][p] = hhi(g_part)
            
          self.statistic[c] = sum(self.stat[c] * self.weights)
            
          if self.weights[old_idx] == 0:
            self.stat[c] = np.delete(self.stat[c], old_idx)

        if self.weights[old_idx] == 0:
          self.weights = np.delete(self.weights, old_idx)
          self.groups = np.delete(self.groups, old_idx)
          
          
    def get_deltas(self, node: int, candidate_partitions: List[int]) -> List[np.ndarray]:
        n_candidates: int = len(candidate_partitions)
        deltas: List[np.ndarray] = [np.zeros((1, self.spatial_lattice.P)) for _ in range(n_candidates)]
        new_partitioning = self.partitioning.copy()
        old_partition = self.partitioning[node]
        
        # Adjust old weight
        old_idx = np.where(self.groups == old_partition)
        
        for idx in range(n_candidates):
          
          diff: np.ndarray = np.zeros((1, self.spatial_lattice.P))
          new_partitioning[node] = candidate_partitions[idx]
          
          temp_weights = self.weights.copy()
          temp_weights[old_idx] -= 1
          temp_groups = self.groups.copy()

          # Adjust new weight
          if all(temp_groups != candidate_partitions[idx]):
            is_new = True
            temp_weights = np.append(temp_weights, int(0))
            temp_groups = np.append(temp_groups, candidate_partitions[idx])
          else:
            is_new = False
          new_idx = np.where(temp_groups == candidate_partitions[idx])
          temp_weights[new_idx] += 1
          
          for c in range(self.spatial_lattice.P):
            temp_stat = self.stat[c].copy()
            if is_new:
              temp_stat = np.append(temp_stat, int(0))
            
            for p in [old_idx,new_idx]:
              g_part = np.array(self.spatial_lattice.X[:,c])[np.array(new_partitioning) == temp_groups[p]]
              temp_stat[p] = hhi(g_part)
              
            diff[c] = sum(temp_stat * temp_weights) - self.statistic
          deltas[idx] = diff
          
        return deltas

# 
# ## Packing / slow & safe
# class VertexPackSlowNS(NetworkStatistic):
#     # NS that evaluates how cracked groups in 
#     
#     def __init__(self, spatial_lattice, partitioning, contributes_energy=True):
#         super().__init__(spatial_lattice, partitioning, contributes_energy)
# 
#     def get_name(self) -> str:
#         return 'vertex_pack_slow'
#       
#     def get_statistic_shape(self) -> Tuple[int, int]:
#         return self.spatial_lattice.P, 1
# 
#     def compute_statistic(self) -> None:
#         self.statistic: np.ndarray = np.zeros((1, self.spatial_lattice.P))
#         for c in range(self.spatial_lattice.P):
#           self.statistic[c] = hhi_all(self.partitioning, self.spatial_lattice.X[:,c])
# 
#     def update_statistic(self, node: int, new_partition: int) -> None:
#         new_partitioning = self.partitioning.copy()
#         new_partitioning[node] = new_partition
#         self.statistic: np.ndarray = np.zeros((1, self.spatial_lattice.P))
#         for c in range(self.spatial_lattice.P):
#           self.statistic[c] = hhi_all(new_partitioning, self.spatial_lattice.X[:,c])
# 
#     def get_deltas(self, node: int, candidate_partitions: List[int]) -> List[np.ndarray]:
#         n_candidates: int = len(candidate_partitions)
#         deltas: List[np.ndarray] = [np.zeros((1, self.spatial_lattice.P)) for _ in range(n_candidates)]
#         new_partitioning = self.partitioning.copy()
#         for idx in range(n_candidates):
#           new_stat: np.ndarray = np.zeros((1, self.spatial_lattice.P))
#           new_partitioning[node] = candidate_partitions[idx]
#           for c in range(self.spatial_lattice.P):
#             new_stat[c] = hhi_all(new_partitioning, self.spatial_lattice.X[:,c])
#           deltas[idx] = self.statistic - new_stat
#           
#         return deltas
