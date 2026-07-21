#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace {

struct KdNode {
  int row;
  int axis;
  int left;
  int right;
};

class KdTree {
 public:
  KdTree(const NumericMatrix& x, const std::vector<int>& rows)
      : x_(x), dims_(x.ncol()), order_(rows) {
    nodes_.reserve(rows.size());
    root_ = build(0, static_cast<int>(order_.size()), 0);
  }

  void query(int row,
             int k,
             std::vector<int>& neighbors,
             std::vector<double>& distances) const {
    std::priority_queue<std::pair<double, int> > heap;
    search(root_, row, k, heap);
    const int size = static_cast<int>(heap.size());
    neighbors.resize(size);
    distances.resize(size);
    for (int i = size - 1; i >= 0; --i) {
      neighbors[i] = heap.top().second;
      distances[i] = heap.top().first;
      heap.pop();
    }
  }

 private:
  const NumericMatrix& x_;
  int dims_;
  int root_;
  std::vector<int> order_;
  std::vector<KdNode> nodes_;

  int build(int begin, int end, int depth) {
    if (begin >= end) return -1;
    const int axis = depth % dims_;
    const int middle = begin + (end - begin) / 2;
    std::nth_element(
      order_.begin() + begin,
      order_.begin() + middle,
      order_.begin() + end,
      [&](int a, int b) { return x_(a, axis) < x_(b, axis); }
    );
    const int node_index = static_cast<int>(nodes_.size());
    nodes_.push_back({order_[middle], axis, -1, -1});
    nodes_[node_index].left = build(begin, middle, depth + 1);
    nodes_[node_index].right = build(middle + 1, end, depth + 1);
    return node_index;
  }

  double squared_distance(int a, int b) const {
    double out = 0.0;
    for (int d = 0; d < dims_; ++d) {
      const double delta = x_(a, d) - x_(b, d);
      out += delta * delta;
    }
    return out;
  }

  void search(int node_index,
              int target,
              int k,
              std::priority_queue<std::pair<double, int> >& heap) const {
    if (node_index < 0) return;
    const KdNode& node = nodes_[node_index];
    const double delta = x_(target, node.axis) - x_(node.row, node.axis);
    const int near_child = delta <= 0.0 ? node.left : node.right;
    const int far_child = delta <= 0.0 ? node.right : node.left;

    search(near_child, target, k, heap);
    if (node.row != target) {
      const double distance = squared_distance(target, node.row);
      if (static_cast<int>(heap.size()) < k) {
        heap.push(std::make_pair(distance, node.row));
      } else if (distance < heap.top().first) {
        heap.pop();
        heap.push(std::make_pair(distance, node.row));
      }
    }

    const double limit = static_cast<int>(heap.size()) < k
      ? std::numeric_limits<double>::infinity()
      : heap.top().first;
    if (delta * delta <= limit) search(far_child, target, k, heap);
  }
};

struct FlowEdge {
  int to;
  int reverse;
  double capacity;
};

class Dinic {
 public:
  explicit Dinic(int vertices)
      : graph_(vertices), level_(vertices), next_(vertices) {}

  void add_edge(int from, int to, double capacity) {
    FlowEdge forward = {to, static_cast<int>(graph_[to].size()), capacity};
    FlowEdge backward = {from, static_cast<int>(graph_[from].size()), 0.0};
    graph_[from].push_back(forward);
    graph_[to].push_back(backward);
  }

  double max_flow(int source, int sink) {
    double total = 0.0;
    while (build_levels(source, sink)) {
      std::fill(next_.begin(), next_.end(), 0);
      while (true) {
        const double pushed = send_flow(source, sink,
                                        std::numeric_limits<double>::infinity());
        if (pushed <= kEpsilon) break;
        total += pushed;
      }
    }
    return total;
  }

  std::vector<bool> source_reachable(int source) const {
    std::vector<bool> reachable(graph_.size(), false);
    std::queue<int> pending;
    pending.push(source);
    reachable[source] = true;
    while (!pending.empty()) {
      const int vertex = pending.front();
      pending.pop();
      for (const FlowEdge& edge : graph_[vertex]) {
        if (edge.capacity > kEpsilon && !reachable[edge.to]) {
          reachable[edge.to] = true;
          pending.push(edge.to);
        }
      }
    }
    return reachable;
  }

 private:
  static constexpr double kEpsilon = 1e-10;
  std::vector<std::vector<FlowEdge> > graph_;
  std::vector<int> level_;
  std::vector<int> next_;

  bool build_levels(int source, int sink) {
    std::fill(level_.begin(), level_.end(), -1);
    std::queue<int> pending;
    pending.push(source);
    level_[source] = 0;
    while (!pending.empty()) {
      const int vertex = pending.front();
      pending.pop();
      for (const FlowEdge& edge : graph_[vertex]) {
        if (edge.capacity > kEpsilon && level_[edge.to] < 0) {
          level_[edge.to] = level_[vertex] + 1;
          pending.push(edge.to);
        }
      }
    }
    return level_[sink] >= 0;
  }

  double send_flow(int vertex, int sink, double available) {
    if (vertex == sink) return available;
    for (int& edge_index = next_[vertex];
         edge_index < static_cast<int>(graph_[vertex].size()); ++edge_index) {
      FlowEdge& edge = graph_[vertex][edge_index];
      if (edge.capacity <= kEpsilon || level_[edge.to] != level_[vertex] + 1) {
        continue;
      }
      const double pushed = send_flow(
        edge.to, sink, std::min(available, edge.capacity)
      );
      if (pushed > kEpsilon) {
        edge.capacity -= pushed;
        graph_[edge.to][edge.reverse].capacity += pushed;
        return pushed;
      }
    }
    return 0.0;
  }
};

struct PottsComponent {
  std::vector<int> rows;
  std::vector<std::pair<int, int> > edges;
  std::vector<int> classes;
  int neighbors;
};

double component_potts_energy(const std::vector<int>& current,
                              const std::vector<int>& observed,
                              const PottsComponent& component,
                              double unary) {
  double energy = 0.0;
  for (int local = 0; local < static_cast<int>(component.rows.size()); ++local) {
    const int row = component.rows[local];
    if (current[row] != observed[row]) energy += unary;
  }
  for (const std::pair<int, int>& edge : component.edges) {
    const int left = component.rows[edge.first];
    const int right = component.rows[edge.second];
    if (current[left] != current[right]) energy += 1.0;
  }
  return energy;
}

bool alpha_expand_component(std::vector<int>& current,
                            const std::vector<int>& observed,
                            const PottsComponent& component,
                            int alpha,
                            double unary) {
  const int count = static_cast<int>(component.rows.size());
  if (count < 2) return false;
  const int source = count;
  const int sink = count + 1;
  const double infinity = 1e12;
  std::vector<double> retain_cost(count, 0.0);
  std::vector<double> switch_cost(count, 0.0);

  for (int local = 0; local < count; ++local) {
    const int row = component.rows[local];
    const int label = current[row];
    retain_cost[local] = label == observed[row] ? 0.0 : unary;
    switch_cost[local] = alpha == observed[row] ? 0.0 : unary;
    if (label == alpha) {
      retain_cost[local] = infinity;
      switch_cost[local] = 0.0;
    }
  }

  Dinic graph(count + 2);
  for (const std::pair<int, int>& edge : component.edges) {
    const int left = edge.first;
    const int right = edge.second;
    const int left_label = current[component.rows[left]];
    const int right_label = current[component.rows[right]];
    const double e00 = left_label == right_label ? 0.0 : 1.0;
    const double e01 = left_label == alpha ? 0.0 : 1.0;
    const double e10 = alpha == right_label ? 0.0 : 1.0;
    const double pairwise = std::max(0.0, 0.5 * (e01 + e10 - e00));
    switch_cost[left] += e10 - e00 - pairwise;
    switch_cost[right] += e01 - e00 - pairwise;
    if (pairwise > 1e-12) {
      graph.add_edge(left, right, pairwise);
      graph.add_edge(right, left, pairwise);
    }
  }

  for (int local = 0; local < count; ++local) {
    const double shift = std::min(retain_cost[local], switch_cost[local]);
    retain_cost[local] -= shift;
    switch_cost[local] -= shift;
    graph.add_edge(source, local, std::max(0.0, switch_cost[local]));
    graph.add_edge(local, sink, std::max(0.0, retain_cost[local]));
  }

  graph.max_flow(source, sink);
  const std::vector<bool> source_side = graph.source_reachable(source);
  std::vector<int> changed_rows;
  std::vector<int> previous_labels;
  changed_rows.reserve(count);
  previous_labels.reserve(count);
  for (int local = 0; local < count; ++local) {
    if (!source_side[local]) {
      const int row = component.rows[local];
      if (current[row] != alpha) {
        previous_labels.push_back(current[row]);
        current[row] = alpha;
        changed_rows.push_back(row);
      }
    }
  }
  if (changed_rows.empty()) return false;

  const double after = component_potts_energy(current, observed, component, unary);
  for (std::size_t index = 0; index < changed_rows.size(); ++index) {
    current[changed_rows[index]] = previous_labels[index];
  }
  const double before = component_potts_energy(current, observed, component, unary);
  if (after > before + 1e-8) {
    return false;
  }
  for (const int row : changed_rows) current[row] = alpha;
  return true;
}

} // namespace

extern "C" SEXP _fibermargin_refine_spatial_graph_cpp(SEXP xy_s,
                                                           SEXP labels_s,
                                                           SEXP samples_s,
                                                           SEXP neighbors_s,
                                                           SEXP iterations_s,
                                                           SEXP consensus_s,
                                                           SEXP preserve_s,
                                                           SEXP margin_s,
                                                           SEXP weighted_s,
                                                           SEXP current_support_s) {
  BEGIN_RCPP
  NumericMatrix xy(xy_s);
  IntegerVector labels(labels_s);
  IntegerVector samples(samples_s);
  const int requested_k = as<int>(neighbors_s);
  const int iterations = std::max(1, as<int>(iterations_s));
  const double consensus = std::min(0.99, std::max(0.5, as<double>(consensus_s)));
  const double preserve = std::max(0.0, as<double>(preserve_s));
  const double required_margin = std::min(0.99, std::max(0.0, as<double>(margin_s)));
  const bool weighted = as<bool>(weighted_s);
  const double max_current_support = std::min(1.0, std::max(0.0, as<double>(current_support_s)));
  const int n = xy.nrow();

  std::vector<int> current(n);
  std::vector<int> next(n);
  std::vector<int> sample_codes(n);
  int n_classes = 0;
  for (int i = 0; i < n; ++i) {
    current[i] = labels[i];
    sample_codes[i] = samples[i];
    n_classes = std::max(n_classes, current[i]);
  }

  std::vector<int> sample_levels(sample_codes.begin(), sample_codes.end());
  std::sort(sample_levels.begin(), sample_levels.end());
  sample_levels.erase(std::unique(sample_levels.begin(), sample_levels.end()), sample_levels.end());

  std::vector<std::vector<int> > graph_rows(n);
  std::vector<std::vector<double> > graph_distances(n);
  IntegerVector sample_neighbors(sample_levels.size());
  for (std::size_t sample_index = 0; sample_index < sample_levels.size(); ++sample_index) {
    const int sample = sample_levels[sample_index];
    Rcpp::checkUserInterrupt();
    std::vector<int> rows;
    for (int i = 0; i < n; ++i) {
      if (sample_codes[i] == sample) rows.push_back(i);
    }
    if (rows.size() <= 1) continue;
    std::vector<bool> observed_classes(static_cast<std::size_t>(n_classes + 1), false);
    int sample_class_count = 0;
    for (int row : rows) {
      const int label = current[row];
      if (!observed_classes[static_cast<std::size_t>(label)]) {
        observed_classes[static_cast<std::size_t>(label)] = true;
        ++sample_class_count;
      }
    }
    const double class_scale = std::min(1.0, 5.0 / std::max(1, sample_class_count));
    const int automatic_k = std::min(
      31,
      std::max(11, static_cast<int>(std::round(
        1.6 * std::log2(static_cast<double>(rows.size())) * class_scale
      )))
    );
    const int chosen_k = requested_k > 0 ? requested_k : automatic_k;
    const int k = std::min(chosen_k, static_cast<int>(rows.size()) - 1);
    sample_neighbors[sample_index] = k;
    KdTree tree(xy, rows);
    for (int row : rows) tree.query(row, k, graph_rows[row], graph_distances[row]);
  }

  NumericVector confidence(n, 1.0);
  IntegerVector changed_per_iteration(iterations);
  std::vector<double> votes(static_cast<std::size_t>(n_classes + 1), 0.0);

  int completed = 0;
  for (int iteration = 0; iteration < iterations; ++iteration) {
    Rcpp::checkUserInterrupt();
    next = current;
    int changed = 0;
    for (int row = 0; row < n; ++row) {
      const std::vector<int>& nn = graph_rows[row];
      const std::vector<double>& dd = graph_distances[row];
      if (nn.empty()) continue;
      std::fill(votes.begin(), votes.end(), 0.0);
      const double bandwidth = std::max(dd.back(), 1e-12);
      double neighbor_weight = 0.0;
      for (std::size_t j = 0; j < nn.size(); ++j) {
        const double weight = weighted ? std::exp(-dd[j] / bandwidth) : 1.0;
        votes[current[nn[j]]] += weight;
        neighbor_weight += weight;
      }
      const double current_support = neighbor_weight > 0.0
        ? votes[current[row]] / neighbor_weight
        : 1.0;
      votes[current[row]] += preserve * neighbor_weight;
      const double total = neighbor_weight * (1.0 + preserve);
      int best_class = current[row];
      double best_vote = votes[best_class];
      double second_vote = 0.0;
      for (int cls = 1; cls <= n_classes; ++cls) {
        if (votes[cls] > best_vote) {
          second_vote = best_vote;
          best_vote = votes[cls];
          best_class = cls;
        } else if (cls != best_class && votes[cls] > second_vote) {
          second_vote = votes[cls];
        }
      }
      confidence[row] = total > 0.0 ? best_vote / total : 1.0;
      const double margin = total > 0.0 ? (best_vote - second_vote) / total : 0.0;
      if (best_class != current[row] && current_support <= max_current_support &&
          confidence[row] >= consensus && margin >= required_margin) {
        next[row] = best_class;
        ++changed;
      }
    }
    current.swap(next);
    changed_per_iteration[iteration] = changed;
    completed = iteration + 1;
    if (changed == 0) break;
  }

  IntegerVector output(n);
  for (int i = 0; i < n; ++i) output[i] = current[i];
  output.attr("confidence") = confidence;
  output.attr("changes") = changed_per_iteration[Range(0, completed - 1)];
  output.attr("neighbors") = sample_neighbors;
  return output;
  END_RCPP
}

extern "C" SEXP _fibermargin_direct_refiner_cpp(SEXP xy_s,
                                                        SEXP labels_s,
                                                        SEXP samples_s,
                                                        SEXP method_s,
                                                        SEXP neighbors_s) {
  BEGIN_RCPP
  NumericMatrix xy(xy_s);
  IntegerVector labels(labels_s);
  IntegerVector samples(samples_s);
  // 0 is the fixed-k unweighted neighbour-mode kernel. It is used both by the
  // GraphST correction protocol and by the explicitly labelled modal control.
  // 1 is the SpaGCN correction protocol.
  const int method = as<int>(method_s);
  const int requested_k = as<int>(neighbors_s);
  const int n = xy.nrow();
  int n_classes = 0;
  std::vector<int> sample_codes(n);
  for (int i = 0; i < n; ++i) {
    n_classes = std::max(n_classes, labels[i]);
    sample_codes[i] = samples[i];
  }

  std::vector<int> sample_levels(sample_codes.begin(), sample_codes.end());
  std::sort(sample_levels.begin(), sample_levels.end());
  sample_levels.erase(std::unique(sample_levels.begin(), sample_levels.end()), sample_levels.end());

  IntegerVector output = clone(labels);
  IntegerVector sample_neighbors(sample_levels.size());
  std::vector<int> counts(static_cast<std::size_t>(n_classes + 1), 0);
  for (std::size_t sample_index = 0; sample_index < sample_levels.size(); ++sample_index) {
    Rcpp::checkUserInterrupt();
    std::vector<int> rows;
    for (int i = 0; i < n; ++i) {
      if (sample_codes[i] == sample_levels[sample_index]) rows.push_back(i);
    }
    if (rows.size() <= 1) continue;
    const int k = std::min(requested_k, static_cast<int>(rows.size()) - 1);
    sample_neighbors[sample_index] = k;
    KdTree tree(xy, rows);
    std::vector<int> nearest;
    std::vector<double> distances;
    for (int row : rows) {
      tree.query(row, k, nearest, distances);
      std::fill(counts.begin(), counts.end(), 0);
      for (int neighbor : nearest) ++counts[labels[neighbor]];

      if (method == 0) {
        // Python's max(list, key=list.count) resolves ties by first occurrence.
        int best = labels[row];
        int best_count = -1;
        for (int neighbor : nearest) {
          const int candidate = labels[neighbor];
          if (counts[candidate] > best_count) {
            best = candidate;
            best_count = counts[candidate];
          }
        }
        output[row] = best;
      } else {
        ++counts[labels[row]];  // SpaGCN includes the focal spot in its vote.
        int best = labels[row];
        int best_count = counts[best];
        for (int cls = 1; cls <= n_classes; ++cls) {
          if (counts[cls] > best_count) {
            best = cls;
            best_count = counts[cls];
          }
        }
        if (counts[labels[row]] < k / 2.0 && best_count > k / 2.0) {
          output[row] = best;
        }
      }
    }
  }
  output.attr("neighbors") = sample_neighbors;
  return output;
  END_RCPP
}

extern "C" SEXP _fibermargin_alpha_expansion_potts_cpp(SEXP xy_s,
                                                        SEXP labels_s,
                                                        SEXP samples_s,
                                                        SEXP neighbors_s,
                                                        SEXP unary_s,
                                                        SEXP cycles_s) {
  BEGIN_RCPP
  NumericMatrix xy(xy_s);
  IntegerVector labels(labels_s);
  IntegerVector samples(samples_s);
  const int n = xy.nrow();
  const int requested_k = std::max(1, as<int>(neighbors_s));
  const double unary = std::max(1e-8, as<double>(unary_s));
  const int maximum_cycles = std::max(1, as<int>(cycles_s));
  if (labels.size() != n || samples.size() != n || xy.ncol() < 1) {
    stop("Incompatible coordinate, label, or specimen inputs.");
  }

  std::vector<int> observed(n);
  std::vector<int> current(n);
  std::vector<int> sample_codes(n);
  std::vector<int> sample_levels;
  sample_levels.reserve(n);
  for (int row = 0; row < n; ++row) {
    if (labels[row] == NA_INTEGER || labels[row] < 1 ||
        samples[row] == NA_INTEGER) {
      stop("Labels and specimen identifiers must be non-missing.");
    }
    observed[row] = labels[row];
    current[row] = labels[row];
    sample_codes[row] = samples[row];
    sample_levels.push_back(samples[row]);
  }
  std::sort(sample_levels.begin(), sample_levels.end());
  sample_levels.erase(
    std::unique(sample_levels.begin(), sample_levels.end()), sample_levels.end()
  );

  std::vector<PottsComponent> components;
  components.reserve(sample_levels.size());
  IntegerVector sample_neighbors(sample_levels.size());
  std::vector<int> local_index(n, -1);
  for (std::size_t sample_index = 0; sample_index < sample_levels.size(); ++sample_index) {
    Rcpp::checkUserInterrupt();
    PottsComponent component;
    const int sample = sample_levels[sample_index];
    for (int row = 0; row < n; ++row) {
      if (sample_codes[row] == sample) component.rows.push_back(row);
    }
    const int count = static_cast<int>(component.rows.size());
    if (count < 2) {
      sample_neighbors[sample_index] = 0;
      components.push_back(component);
      continue;
    }
    const int k = std::min(requested_k, count - 1);
    component.neighbors = k;
    sample_neighbors[sample_index] = k;
    for (int local = 0; local < count; ++local) {
      local_index[component.rows[local]] = local;
    }
    KdTree tree(xy, component.rows);
    std::vector<int> nearest;
    std::vector<double> distances;
    component.edges.reserve(static_cast<std::size_t>(count) * k);
    for (int local = 0; local < count; ++local) {
      const int row = component.rows[local];
      tree.query(row, k, nearest, distances);
      for (const int neighbor : nearest) {
        const int other = local_index[neighbor];
        if (other < 0) continue;
        component.edges.push_back(std::make_pair(
          std::min(local, other), std::max(local, other)
        ));
      }
    }
    std::sort(component.edges.begin(), component.edges.end());
    component.edges.erase(
      std::unique(component.edges.begin(), component.edges.end()), component.edges.end()
    );
    for (const int row : component.rows) component.classes.push_back(observed[row]);
    std::sort(component.classes.begin(), component.classes.end());
    component.classes.erase(
      std::unique(component.classes.begin(), component.classes.end()), component.classes.end()
    );
    for (const int row : component.rows) local_index[row] = -1;
    components.push_back(component);
  }

  IntegerVector changes_per_cycle(maximum_cycles);
  int completed_cycles = 0;
  for (int cycle = 0; cycle < maximum_cycles; ++cycle) {
    int changed = 0;
    for (const PottsComponent& component : components) {
      if (component.classes.size() < 2 || component.edges.empty()) continue;
      for (const int alpha : component.classes) {
        Rcpp::checkUserInterrupt();
        std::vector<int> before;
        before.reserve(component.rows.size());
        for (const int row : component.rows) before.push_back(current[row]);
        if (alpha_expand_component(current, observed, component, alpha, unary)) {
          for (std::size_t local = 0; local < component.rows.size(); ++local) {
            if (current[component.rows[local]] != before[local]) ++changed;
          }
        }
      }
    }
    changes_per_cycle[cycle] = changed;
    completed_cycles = cycle + 1;
    if (changed == 0) break;
  }

  double final_energy = 0.0;
  for (const PottsComponent& component : components) {
    final_energy += component_potts_energy(current, observed, component, unary);
  }
  IntegerVector output(n);
  for (int row = 0; row < n; ++row) output[row] = current[row];
  output.attr("neighbors") = sample_neighbors;
  output.attr("changes") = changes_per_cycle[Range(0, completed_cycles - 1)];
  output.attr("energy") = final_energy;
  output.attr("unary") = unary;
  return output;
  END_RCPP
}

extern "C" SEXP _fibermargin_fiber_margin_cpp(
  SEXP xy_s, SEXP labels_s, SEXP samples_s, SEXP control_s);

static const R_CallMethodDef CallEntries[] = {
  {"_fibermargin_refine_spatial_graph_cpp", (DL_FUNC) &_fibermargin_refine_spatial_graph_cpp, 10},
  {"_fibermargin_direct_refiner_cpp", (DL_FUNC) &_fibermargin_direct_refiner_cpp, 5},
  {"_fibermargin_alpha_expansion_potts_cpp", (DL_FUNC) &_fibermargin_alpha_expansion_potts_cpp, 6},
  {"_fibermargin_fiber_margin_cpp", (DL_FUNC) &_fibermargin_fiber_margin_cpp, 4},
  {NULL, NULL, 0}
};

extern "C" void R_init_fibermargin(DllInfo* dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
