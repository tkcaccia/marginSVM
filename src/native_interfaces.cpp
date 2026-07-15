#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <queue>
#include <random>
#include <unordered_map>
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

struct RffModel {
  int dims;
  int features;
  int classes;
  std::vector<int> class_codes;
  std::vector<double> center;
  std::vector<double> scale;
  std::vector<double> omega;
  std::vector<double> phase;
  std::vector<double> weights;
};

double dot_score(const RffModel& model, int cls, const double* phi) {
  const int stride = model.features + 1;
  const double* w = &model.weights[cls * stride];
  double score = w[model.features];
  for (int f = 0; f < model.features; ++f) {
    score += w[f] * phi[f];
  }
  return score;
}

void compute_phi(const RffModel& model,
                 const NumericMatrix& x,
                 int row,
                 std::vector<double>& phi) {
  const double scale = std::sqrt(2.0 / static_cast<double>(model.features));
  for (int f = 0; f < model.features; ++f) {
    double projection = model.phase[f];
    const int offset = f * model.dims;
    for (int d = 0; d < model.dims; ++d) {
      const double normalized = (x(row, d) - model.center[d]) / model.scale[d];
      projection += model.omega[offset + d] * normalized;
    }
    phi[f] = scale * std::cos(projection);
  }
}

std::vector<int> unique_sorted_codes(const std::vector<int>& labels,
                                     const std::vector<int>& rows) {
  std::vector<int> out;
  out.reserve(rows.size());
  for (int row : rows) {
    out.push_back(labels[row]);
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

RffModel train_model(const NumericMatrix& xy,
                     const std::vector<int>& y,
                     const std::vector<int>& rows,
                     const std::vector<int>& class_codes,
                     double gamma,
                     int n_features,
                     int epochs,
                     double lambda,
                     double learning_rate,
                     std::uint32_t seed) {
  RffModel model;
  model.dims = xy.ncol();
  model.features = std::max(4, n_features);
  model.classes = class_codes.size();
  model.class_codes = class_codes;
  model.center.assign(model.dims, 0.0);
  model.scale.assign(model.dims, 1.0);
  for (int d = 0; d < model.dims; ++d) {
    double low = std::numeric_limits<double>::infinity();
    double high = -std::numeric_limits<double>::infinity();
    for (int row : rows) {
      low = std::min(low, xy(row, d));
      high = std::max(high, xy(row, d));
    }
    model.center[d] = low;
    model.scale[d] = std::max(high - low, 1e-12);
  }
  model.omega.resize(static_cast<std::size_t>(model.features) * model.dims);
  model.phase.resize(model.features);
  model.weights.assign(static_cast<std::size_t>(model.classes) * (model.features + 1), 0.0);

  std::mt19937 rng(seed);
  std::normal_distribution<double> normal(0.0, std::sqrt(2.0 * gamma));
  std::uniform_real_distribution<double> uniform_phase(0.0, 2.0 * M_PI);

  for (double& value : model.omega) value = normal(rng);
  for (double& value : model.phase) value = uniform_phase(rng);

  std::unordered_map<int, int> class_index;
  for (int i = 0; i < model.classes; ++i) class_index[class_codes[i]] = i;

  std::vector<double> phi(model.features);
  std::vector<int> order(rows.begin(), rows.end());
  const int stride = model.features + 1;
  std::uint64_t iter = 0;

  for (int epoch = 0; epoch < epochs; ++epoch) {
    std::shuffle(order.begin(), order.end(), rng);
    for (int row : order) {
      compute_phi(model, xy, row, phi);
      const int positive = class_index[y[row]];
      const double epoch_progress = static_cast<double>(iter++) /
        std::max<std::size_t>(1, rows.size());
      const double eta = learning_rate / std::sqrt(1.0 + epoch_progress);

      double positive_score = dot_score(model, positive, phi.data());
      double negative_score = -std::numeric_limits<double>::infinity();
      int negative = -1;
      for (int cls = 0; cls < model.classes; ++cls) {
        if (cls == positive) continue;
        const double score = dot_score(model, cls, phi.data());
        if (score > negative_score) {
          negative_score = score;
          negative = cls;
        }
      }
      if (negative >= 0 && positive_score - negative_score < 1.0) {
        double* positive_weights = &model.weights[positive * stride];
        double* negative_weights = &model.weights[negative * stride];
        const double shrink = std::max(0.0, 1.0 - eta * lambda);
        for (int f = 0; f < model.features; ++f) {
          positive_weights[f] = shrink * positive_weights[f] + eta * phi[f];
          negative_weights[f] = shrink * negative_weights[f] - eta * phi[f];
        }
        positive_weights[model.features] += eta;
        negative_weights[model.features] -= eta;
      }
    }
  }

  return model;
}

int tile_id(const NumericMatrix& x,
            int row,
            const std::vector<int>& tiles,
            const std::vector<double>& mins,
            const std::vector<double>& maxs) {
  int id = 0;
  int multiplier = 1;
  const int dims = x.ncol();
  for (int d = 0; d < dims; ++d) {
    if (x(row, d) < mins[d] || x(row, d) > maxs[d]) return -1;
    const double span = maxs[d] - mins[d];
    int pos = 0;
    if (span > 0) {
      pos = static_cast<int>(std::floor((x(row, d) - mins[d]) / span * tiles[d]));
      if (pos < 0) pos = 0;
      if (pos >= tiles[d]) pos = tiles[d] - 1;
    }
    id += pos * multiplier;
    multiplier *= tiles[d];
  }
  return id;
}

void train_predict_block(const NumericMatrix& train_x,
                         const NumericMatrix& pred_x,
                         const std::vector<int>& y,
                         const std::vector<int>& train_rows,
                         const std::vector<int>& pred_rows,
                         IntegerVector& output,
                         NumericMatrix& scores,
                         double gamma,
                         int n_features,
                         int epochs,
                         double lambda,
                         double learning_rate,
                         std::uint32_t seed) {
  if (train_rows.empty() || pred_rows.empty()) return;
  std::vector<int> class_codes = unique_sorted_codes(y, train_rows);
  if (class_codes.empty()) return;
  if (class_codes.size() == 1) {
    for (int row : pred_rows) {
      output[row] = class_codes[0];
      scores(row, class_codes[0] - 1) = 1.0;
    }
    return;
  }

  RffModel model = train_model(
    train_x,
    y,
    train_rows,
    class_codes,
    gamma,
    n_features,
    epochs,
    lambda,
    learning_rate,
    seed
  );

  std::vector<double> phi(model.features);
  for (int row : pred_rows) {
    compute_phi(model, pred_x, row, phi);
    double best_score = -std::numeric_limits<double>::infinity();
    int best_class = model.class_codes[0];
    for (int cls = 0; cls < model.classes; ++cls) {
      const double score = dot_score(model, cls, phi.data());
      scores(row, model.class_codes[cls] - 1) = score;
      if (score > best_score) {
        best_score = score;
        best_class = model.class_codes[cls];
      }
    }
    output[row] = best_class;
  }
}

} // namespace

extern "C" SEXP _marginSVM_refine_spatial_graph_cpp(SEXP xy_s,
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

extern "C" SEXP _marginSVM_direct_refiner_cpp(SEXP xy_s,
                                                        SEXP labels_s,
                                                        SEXP samples_s,
                                                        SEXP method_s,
                                                        SEXP neighbors_s) {
  BEGIN_RCPP
  NumericMatrix xy(xy_s);
  IntegerVector labels(labels_s);
  IntegerVector samples(samples_s);
  const int method = as<int>(method_s);  // 0: GraphST, 1: SpaGCN
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

extern "C" SEXP _marginSVM_structured_spatial_svm_cpp(
  SEXP xy_s, SEXP labels_s, SEXP samples_s, SEXP control_s, SEXP verbose_s);

static const R_CallMethodDef CallEntries[] = {
  {"_marginSVM_refine_spatial_graph_cpp", (DL_FUNC) &_marginSVM_refine_spatial_graph_cpp, 10},
  {"_marginSVM_direct_refiner_cpp", (DL_FUNC) &_marginSVM_direct_refiner_cpp, 5},
  {"_marginSVM_structured_spatial_svm_cpp", (DL_FUNC) &_marginSVM_structured_spatial_svm_cpp, 5},
  {NULL, NULL, 0}
};

extern "C" void R_init_marginSVM(DllInfo* dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
