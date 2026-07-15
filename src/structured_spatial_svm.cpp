#include <Rcpp.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <queue>
#include <random>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace structured_svm {

struct Node {
  int row;
  int axis;
  int left;
  int right;
};

class Tree {
 public:
  Tree(const std::vector<double>& x, int n, int dims, const std::vector<int>& rows)
      : x_(x), n_(n), dims_(dims), order_(rows) {
    nodes_.reserve(rows.size());
    root_ = build(0, static_cast<int>(order_.size()), 0);
  }

  void query(int target, int k, std::vector<int>& out, std::vector<double>& distance) const {
    std::priority_queue<std::pair<double, int> > heap;
    search(root_, target, k, heap);
    out.resize(heap.size());
    distance.resize(heap.size());
    for (int i = static_cast<int>(heap.size()) - 1; i >= 0; --i) {
      out[i] = heap.top().second;
      distance[i] = heap.top().first;
      heap.pop();
    }
  }

 private:
  const std::vector<double>& x_;
  int n_;
  int dims_;
  int root_;
  std::vector<int> order_;
  std::vector<Node> nodes_;

  double value(int row, int dim) const { return x_[row + n_ * dim]; }

  int build(int begin, int end, int depth) {
    if (begin >= end) return -1;
    const int axis = depth % dims_;
    const int middle = begin + (end - begin) / 2;
    std::nth_element(order_.begin() + begin, order_.begin() + middle,
                     order_.begin() + end,
                     [&](int a, int b) { return value(a, axis) < value(b, axis); });
    const int index = nodes_.size();
    nodes_.push_back({order_[middle], axis, -1, -1});
    nodes_[index].left = build(begin, middle, depth + 1);
    nodes_[index].right = build(middle + 1, end, depth + 1);
    return index;
  }

  double squared_distance(int a, int b) const {
    double result = 0.0;
    for (int d = 0; d < dims_; ++d) {
      const double delta = value(a, d) - value(b, d);
      result += delta * delta;
    }
    return result;
  }

  void search(int node_index, int target, int k,
              std::priority_queue<std::pair<double, int> >& heap) const {
    if (node_index < 0) return;
    const Node& node = nodes_[node_index];
    const double delta = value(target, node.axis) - value(node.row, node.axis);
    const int near_child = delta <= 0 ? node.left : node.right;
    const int far_child = delta <= 0 ? node.right : node.left;
    search(near_child, target, k, heap);
    if (node.row != target) {
      const double distance = squared_distance(target, node.row);
      if (static_cast<int>(heap.size()) < k) {
        heap.push({distance, node.row});
      } else if (distance < heap.top().first) {
        heap.pop();
        heap.push({distance, node.row});
      }
    }
    const double limit = static_cast<int>(heap.size()) < k
      ? std::numeric_limits<double>::infinity() : heap.top().first;
    if (delta * delta <= limit) search(far_child, target, k, heap);
  }
};

bool invert_small_matrix(const std::vector<double>& matrix, int dims,
                         std::vector<double>& inverse) {
  std::vector<double> augmented(static_cast<std::size_t>(dims) * 2 * dims, 0.0);
  for (int i = 0; i < dims; ++i) {
    for (int j = 0; j < dims; ++j) augmented[i * 2 * dims + j] = matrix[i * dims + j];
    augmented[i * 2 * dims + dims + i] = 1.0;
  }
  for (int column = 0; column < dims; ++column) {
    int pivot = column;
    for (int row = column + 1; row < dims; ++row) {
      if (std::abs(augmented[row * 2 * dims + column]) >
          std::abs(augmented[pivot * 2 * dims + column])) pivot = row;
    }
    if (std::abs(augmented[pivot * 2 * dims + column]) < 1e-12) return false;
    if (pivot != column) {
      for (int j = 0; j < 2 * dims; ++j) {
        std::swap(augmented[column * 2 * dims + j], augmented[pivot * 2 * dims + j]);
      }
    }
    const double scale = augmented[column * 2 * dims + column];
    for (int j = 0; j < 2 * dims; ++j) augmented[column * 2 * dims + j] /= scale;
    for (int row = 0; row < dims; ++row) {
      if (row == column) continue;
      const double multiple = augmented[row * 2 * dims + column];
      for (int j = 0; j < 2 * dims; ++j) {
        augmented[row * 2 * dims + j] -= multiple * augmented[column * 2 * dims + j];
      }
    }
  }
  inverse.resize(dims * dims);
  for (int i = 0; i < dims; ++i) {
    for (int j = 0; j < dims; ++j) inverse[i * dims + j] = augmented[i * 2 * dims + dims + j];
  }
  return true;
}

void anisotropic_query(const std::vector<double>& x, int n, int dims, int target,
                       const std::vector<int>& candidates,
                       const std::vector<double>& euclidean,
                       int requested, double anisotropy,
                       std::vector<int>& output,
                       std::vector<double>& distance) {
  const int covariance_k = std::min(24, static_cast<int>(candidates.size()));
  std::vector<double> covariance(dims * dims, 0.0);
  for (int index = 0; index < covariance_k; ++index) {
    const int row = candidates[index];
    for (int a = 0; a < dims; ++a) {
      const double da = x[row + n * a] - x[target + n * a];
      for (int b = 0; b < dims; ++b) {
        covariance[a * dims + b] += da * (x[row + n * b] - x[target + n * b]);
      }
    }
  }
  double trace = 0.0;
  for (int d = 0; d < dims; ++d) trace += covariance[d * dims + d];
  const double ridge = std::max(trace / std::max(1, covariance_k * dims) * 0.12, 1e-12);
  for (double& value : covariance) value /= std::max(1, covariance_k);
  for (int d = 0; d < dims; ++d) covariance[d * dims + d] += ridge;
  std::vector<double> inverse;
  if (!invert_small_matrix(covariance, dims, inverse)) {
    output.assign(candidates.begin(), candidates.begin() + std::min(requested, (int)candidates.size()));
    distance.assign(euclidean.begin(), euclidean.begin() + output.size());
    return;
  }
  const double euclidean_scale = std::max(euclidean[std::min(covariance_k - 1,
    static_cast<int>(euclidean.size()) - 1)], 1e-12);
  std::vector<std::pair<double, int> > ranked;
  ranked.reserve(candidates.size());
  for (std::size_t index = 0; index < candidates.size(); ++index) {
    std::vector<double> delta(dims);
    for (int d = 0; d < dims; ++d) delta[d] = x[candidates[index] + n * d] - x[target + n * d];
    double mahalanobis = 0.0;
    for (int a = 0; a < dims; ++a) {
      for (int b = 0; b < dims; ++b) mahalanobis += delta[a] * inverse[a * dims + b] * delta[b];
    }
    const double mixed = (1.0 - anisotropy) * euclidean[index] / euclidean_scale +
      anisotropy * mahalanobis / std::max(1.0, static_cast<double>(dims));
    ranked.push_back({mixed, candidates[index]});
  }
  const int keep = std::min(requested, static_cast<int>(ranked.size()));
  std::partial_sort(ranked.begin(), ranked.begin() + keep, ranked.end());
  output.resize(keep);
  distance.resize(keep);
  for (int i = 0; i < keep; ++i) {
    output[i] = ranked[i].second;
    distance[i] = ranked[i].first;
  }
}

struct Edge {
  int a;
  int b;
  double spatial_weight;
};

struct Tile {
  std::vector<int> core;
  std::vector<int> halo;
  std::vector<double> low;
  std::vector<double> high;
  int depth;
  double impurity;
  int scale;
};

struct ModelResult {
  std::vector<int> rows;
  std::vector<double> probabilities;
  std::vector<double> taper;
  int classes;
  int landmarks;
};

double gini(const std::vector<int>& rows, const std::vector<int>& labels, int classes) {
  if (rows.empty()) return 0.0;
  std::vector<int> count(classes, 0);
  for (int row : rows) ++count[labels[row] - 1];
  double sum = 0.0;
  for (int value : count) {
    const double probability = static_cast<double>(value) / rows.size();
    sum += probability * probability;
  }
  return 1.0 - sum;
}

void adaptive_split(const std::vector<double>& x, int n, int dims,
                    const std::vector<int>& labels, int classes,
                    const std::vector<int>& rows,
                    const std::vector<double>& low,
                    const std::vector<double>& high,
                    int depth, int target, int minimum, int maximum_depth,
                    std::vector<Tile>& leaves) {
  const double parent_gini = gini(rows, labels, classes);
  int best_dimension = -1;
  double best_cut = 0.0;
  double best_gain = -1.0;
  std::vector<int> best_left;
  std::vector<int> best_right;

  for (int d = 0; d < dims; ++d) {
    std::vector<double> coordinate;
    coordinate.reserve(rows.size());
    for (int row : rows) coordinate.push_back(x[row + n * d]);
    const int middle = coordinate.size() / 2;
    std::nth_element(coordinate.begin(), coordinate.begin() + middle, coordinate.end());
    const double cut = coordinate[middle];
    std::vector<int> left;
    std::vector<int> right;
    for (int row : rows) {
      if (x[row + n * d] <= cut) left.push_back(row);
      else right.push_back(row);
    }
    if (static_cast<int>(left.size()) < minimum ||
        static_cast<int>(right.size()) < minimum) continue;
    const double child = (left.size() * gini(left, labels, classes) +
      right.size() * gini(right, labels, classes)) / rows.size();
    const double gain = parent_gini - child;
    if (gain > best_gain) {
      best_gain = gain;
      best_dimension = d;
      best_cut = cut;
      best_left.swap(left);
      best_right.swap(right);
    }
  }

  const bool mandatory = static_cast<int>(rows.size()) > target;
  const bool useful = best_gain >= 0.025 &&
    static_cast<int>(rows.size()) >= 2 * minimum;
  if (depth >= maximum_depth || best_dimension < 0 || (!mandatory && !useful)) {
    leaves.push_back({rows, std::vector<int>(), low, high, depth, parent_gini, 1});
    return;
  }

  std::vector<double> left_high(high);
  std::vector<double> right_low(low);
  left_high[best_dimension] = best_cut;
  right_low[best_dimension] = best_cut;
  adaptive_split(x, n, dims, labels, classes, best_left, low, left_high,
                 depth + 1, target, minimum, maximum_depth, leaves);
  adaptive_split(x, n, dims, labels, classes, best_right, right_low, high,
                 depth + 1, target, minimum, maximum_depth, leaves);
}

std::vector<int> balanced_counts(const std::vector<double>& low,
                                 const std::vector<double>& high,
                                 int desired) {
  const int dims = low.size();
  std::vector<int> count(dims, 1);
  while (std::accumulate(count.begin(), count.end(), 1, std::multiplies<int>()) < desired) {
    int best = 0;
    double best_width = -1.0;
    for (int d = 0; d < dims; ++d) {
      const double width = (high[d] - low[d]) / count[d];
      if (width > best_width) {
        best_width = width;
        best = d;
      }
    }
    ++count[best];
  }
  return count;
}

void add_regular_tiles(const std::vector<double>& x, int n, int dims,
                       const std::vector<int>& rows,
                       const std::vector<double>& low,
                       const std::vector<double>& high,
                       int target, std::vector<Tile>& tiles) {
  const int desired = std::max(1, static_cast<int>(std::ceil(
    static_cast<double>(rows.size()) / target)));
  const std::vector<int> count = balanced_counts(low, high, desired);
  int total = 1;
  for (int value : count) total *= value;
  for (int id = 0; id < total; ++id) {
    int remainder = id;
    std::vector<double> tile_low(dims);
    std::vector<double> tile_high(dims);
    for (int d = 0; d < dims; ++d) {
      const int position = remainder % count[d];
      remainder /= count[d];
      const double width = (high[d] - low[d]) / count[d];
      tile_low[d] = low[d] + position * width;
      tile_high[d] = position == count[d] - 1 ? high[d] : tile_low[d] + width;
    }
    std::vector<int> core;
    for (int row : rows) {
      bool inside = true;
      for (int d = 0; d < dims; ++d) {
        inside = inside && x[row + n * d] >= tile_low[d] && x[row + n * d] <= tile_high[d];
      }
      if (inside) core.push_back(row);
    }
    if (!core.empty()) tiles.push_back({core, std::vector<int>(), tile_low, tile_high, 0, 0.0, 0});
  }
}

void add_halos(std::vector<Tile>& tiles, const std::vector<double>& x,
               int n, int dims, const std::vector<int>& sample_rows, double overlap) {
  for (Tile& tile : tiles) {
    std::vector<double> width(dims);
    for (int d = 0; d < dims; ++d) width[d] = std::max(tile.high[d] - tile.low[d], 1e-12);
    for (int row : sample_rows) {
      bool inside = true;
      for (int d = 0; d < dims; ++d) {
        inside = inside && x[row + n * d] >= tile.low[d] - overlap * width[d] &&
          x[row + n * d] <= tile.high[d] + overlap * width[d];
      }
      if (inside) tile.halo.push_back(row);
    }
  }
}

double kernel(const std::vector<double>& x, int n, int dims, int a, int b,
              const std::vector<double>& low, const std::vector<double>& scale,
              double gamma) {
  double distance = 0.0;
  for (int d = 0; d < dims; ++d) {
    const double xa = (x[a + n * d] - low[d]) / scale[d];
    const double xb = (x[b + n * d] - low[d]) / scale[d];
    const double delta = xa - xb;
    distance += delta * delta;
  }
  return std::exp(-gamma * distance);
}

bool cholesky(std::vector<double>& matrix, int size) {
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sum = matrix[i * size + j];
      for (int k = 0; k < j; ++k) sum -= matrix[i * size + k] * matrix[j * size + k];
      if (i == j) {
        if (sum <= 1e-12) return false;
        matrix[i * size + j] = std::sqrt(sum);
      } else {
        matrix[i * size + j] = sum / matrix[j * size + j];
      }
    }
    for (int j = i + 1; j < size; ++j) matrix[i * size + j] = 0.0;
  }
  return true;
}

void forward_solve(const std::vector<double>& lower, int size,
                   const std::vector<double>& rhs, double* output) {
  for (int i = 0; i < size; ++i) {
    double value = rhs[i];
    for (int j = 0; j < i; ++j) value -= lower[i * size + j] * output[j];
    output[i] = value / lower[i * size + i];
  }
}

std::vector<int> select_landmarks(const Tile& tile, const std::vector<int>& labels,
                                  const std::vector<double>& support, int classes,
                                  int requested, std::mt19937& rng) {
  const int maximum = std::min(requested, static_cast<int>(tile.halo.size()));
  std::vector<std::vector<int> > by_class(classes);
  for (int row : tile.halo) by_class[labels[row] - 1].push_back(row);
  std::vector<int> selected;
  selected.reserve(maximum);
  std::unordered_set<int> used;
  used.reserve(static_cast<std::size_t>(maximum) * 2);
  const int per_class = std::max(1, maximum / std::max(1, classes * 3));
  for (int cls = 0; cls < classes; ++cls) {
    std::sort(by_class[cls].begin(), by_class[cls].end(),
              [&](int a, int b) { return support[a] < support[b]; });
    const int take = std::min(per_class, static_cast<int>(by_class[cls].size()));
    for (int i = 0; i < take; ++i) {
      selected.push_back(by_class[cls][i]);
      used.insert(by_class[cls][i]);
    }
  }
  std::vector<int> candidates(tile.halo);
  std::shuffle(candidates.begin(), candidates.end(), rng);
  std::stable_sort(candidates.begin(), candidates.end(),
                   [&](int a, int b) { return support[a] < support[b]; });
  for (int row : candidates) {
    if (static_cast<int>(selected.size()) >= maximum) break;
    if (used.insert(row).second) {
      selected.push_back(row);
    }
  }
  return selected;
}

ModelResult fit_tile(const Tile& tile, const std::vector<double>& x, int n, int dims,
                     const std::vector<int>& labels, const std::vector<double>& support,
                     int global_classes, int requested_landmarks, double gamma,
                     int epochs, double learning_rate, double lambda, double ramp,
                     bool cross_fitting,
                     bool pairwise_specialists, int pairwise_max,
                     std::uint32_t seed) {
  std::mt19937 rng(seed);
  const std::vector<int> landmark_rows = select_landmarks(
    tile, labels, support, global_classes, requested_landmarks, rng);
  const int m = landmark_rows.size();
  std::vector<int> class_codes;
  class_codes.reserve(global_classes);
  for (int row : tile.halo) class_codes.push_back(labels[row]);
  std::sort(class_codes.begin(), class_codes.end());
  class_codes.erase(std::unique(class_codes.begin(), class_codes.end()), class_codes.end());
  const int classes = class_codes.size();
  std::vector<int> class_index(global_classes + 1, -1);
  for (int c = 0; c < classes; ++c) class_index[class_codes[c]] = c;

  std::vector<double> low(dims);
  std::vector<double> scale(dims);
  for (int d = 0; d < dims; ++d) {
    low[d] = std::numeric_limits<double>::infinity();
    double high = -std::numeric_limits<double>::infinity();
    for (int row : tile.halo) {
      low[d] = std::min(low[d], x[row + n * d]);
      high = std::max(high, x[row + n * d]);
    }
    scale[d] = std::max(high - low[d], 1e-12);
  }

  std::vector<double> lower(m * m, 0.0);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= i; ++j) {
      lower[i * m + j] = kernel(x, n, dims, landmark_rows[i], landmark_rows[j],
                                low, scale, gamma);
    }
    lower[i * m + i] += 1e-4;
  }
  if (!cholesky(lower, m)) {
    for (int i = 0; i < m; ++i) lower[i * m + i] += 1e-2;
    cholesky(lower, m);
  }

  const int train_n = tile.halo.size();
  std::vector<double> phi(static_cast<std::size_t>(train_n) * m);
  std::vector<double> rhs(m);
  for (int i = 0; i < train_n; ++i) {
    for (int j = 0; j < m; ++j) {
      rhs[j] = kernel(x, n, dims, tile.halo[i], landmark_rows[j], low, scale, gamma);
    }
    forward_solve(lower, m, rhs, &phi[static_cast<std::size_t>(i) * m]);
  }

  const int stride = m + 1;
  // Two-fold cross-fitting prevents a point from validating its own noisy label.
  const int folds = cross_fitting && train_n >= 40 ? 2 : 1;
  std::vector<int> fold(train_n, 0);
  if (folds == 2) {
    for (int i = 0; i < train_n; ++i) {
      const std::uint32_t row = static_cast<std::uint32_t>(tile.halo[i]);
      fold[i] = static_cast<int>((row * 2654435761U + seed) & 1U);
    }
  }

  std::vector<double> weights(
    static_cast<std::size_t>(folds) * classes * stride, 0.0);
  std::vector<int> order(train_n);
  std::iota(order.begin(), order.end(), 0);
  std::vector<int> counts(classes, 0);
  for (int row : tile.halo) ++counts[class_index[labels[row]]];
  std::vector<double> class_weight(classes, 1.0);
  for (int c = 0; c < classes; ++c) {
    const double balanced = static_cast<double>(train_n) /
      (classes * std::max(1, counts[c]));
    class_weight[c] = std::min(2.0, std::max(0.75, std::sqrt(balanced)));
  }

  std::uint64_t iteration = 0;
  for (int epoch = 0; epoch < epochs; ++epoch) {
    std::shuffle(order.begin(), order.end(), rng);
    for (int position : order) {
      const int model = folds == 1 ? 0 : 1 - fold[position];
      const std::size_t model_offset =
        static_cast<std::size_t>(model) * classes * stride;
      const double* feature = &phi[static_cast<std::size_t>(position) * m];
      const int positive = class_index[labels[tile.halo[position]]];
      double positive_score = weights[model_offset + positive * stride + m];
      for (int f = 0; f < m; ++f) {
        positive_score += weights[model_offset + positive * stride + f] * feature[f];
      }
      int negative = -1;
      double negative_score = -std::numeric_limits<double>::infinity();
      for (int c = 0; c < classes; ++c) {
        if (c == positive) continue;
        double score = weights[model_offset + c * stride + m];
        for (int f = 0; f < m; ++f) {
          score += weights[model_offset + c * stride + f] * feature[f];
        }
        if (score > negative_score) {
          negative_score = score;
          negative = c;
        }
      }
      const double loss = 1.0 - positive_score + negative_score;
      if (negative < 0 || loss <= 0 || (epoch >= 2 && loss >= ramp)) continue;
      const double eta = learning_rate * class_weight[positive] /
        std::sqrt(1.0 + static_cast<double>(iteration++) / std::max(1, train_n));
      const double shrink = std::max(0.0, 1.0 - eta * lambda);
      for (int f = 0; f < m; ++f) {
        weights[model_offset + positive * stride + f] =
          shrink * weights[model_offset + positive * stride + f] + eta * feature[f];
        weights[model_offset + negative * stride + f] =
          shrink * weights[model_offset + negative * stride + f] - eta * feature[f];
      }
      weights[model_offset + positive * stride + m] += eta;
      weights[model_offset + negative * stride + m] -= eta;
    }
  }

  std::vector<double> probability(static_cast<std::size_t>(train_n) * global_classes, 0.0);
  std::vector<double> score(classes);
  for (int i = 0; i < train_n; ++i) {
    const std::size_t model_offset =
      static_cast<std::size_t>(fold[i]) * classes * stride;
    const double* feature = &phi[static_cast<std::size_t>(i) * m];
    double maximum = -std::numeric_limits<double>::infinity();
    for (int c = 0; c < classes; ++c) {
      score[c] = weights[model_offset + c * stride + m];
      for (int f = 0; f < m; ++f) {
        score[c] += weights[model_offset + c * stride + f] * feature[f];
      }
      maximum = std::max(maximum, score[c]);
    }
    double sum = 0.0;
    for (double& value : score) {
      value = std::exp(value - maximum);
      sum += value;
    }
    for (int c = 0; c < classes; ++c) {
      probability[static_cast<std::size_t>(i) * global_classes + class_codes[c] - 1] = score[c] / sum;
    }
  }

  // Boundary specialists revisit only ambiguous locally plausible class pairs.
  if (pairwise_specialists && classes > 1 && pairwise_max > 0) {
    std::unordered_map<std::uint64_t, int> pair_count;
    for (int i = 0; i < train_n; ++i) {
      int first = -1, second = -1;
      double first_p = -1.0, second_p = -1.0;
      for (int c = 0; c < classes; ++c) {
        const double value = probability[static_cast<std::size_t>(i) * global_classes + class_codes[c] - 1];
        if (value > first_p) {
          second_p = first_p; second = first; first_p = value; first = c;
        } else if (value > second_p) {
          second_p = value; second = c;
        }
      }
      if (second < 0 || first_p - second_p > 0.35) continue;
      const int a = std::min(first, second), b = std::max(first, second);
      ++pair_count[(static_cast<std::uint64_t>(a) << 32) | static_cast<std::uint32_t>(b)];
    }
    std::vector<std::pair<int, std::uint64_t> > ranked_pairs;
    for (const auto& item : pair_count) if (item.second >= 12) ranked_pairs.push_back({item.second, item.first});
    std::sort(ranked_pairs.begin(), ranked_pairs.end(),
      [](const std::pair<int, std::uint64_t>& a, const std::pair<int, std::uint64_t>& b) {
        return a.first > b.first;
      });
    if (static_cast<int>(ranked_pairs.size()) > pairwise_max) ranked_pairs.resize(pairwise_max);
    for (const auto& selected : ranked_pairs) {
      const int class_a = static_cast<int>(selected.second >> 32);
      const int class_b = static_cast<int>(selected.second & 0xffffffffU);
      std::vector<double> binary_weights(static_cast<std::size_t>(folds) * stride, 0.0);
      std::uint64_t binary_iteration = 0;
      for (int epoch = 0; epoch < std::max(4, epochs / 2); ++epoch) {
        std::shuffle(order.begin(), order.end(), rng);
        for (int position : order) {
          const int observed = class_index[labels[tile.halo[position]]];
          if (observed != class_a && observed != class_b) continue;
          const int model = folds == 1 ? 0 : 1 - fold[position];
          double* binary = &binary_weights[static_cast<std::size_t>(model) * stride];
          const double* feature = &phi[static_cast<std::size_t>(position) * m];
          const double sign = observed == class_a ? 1.0 : -1.0;
          double score = binary[m];
          for (int f = 0; f < m; ++f) score += binary[f] * feature[f];
          const double loss = 1.0 - sign * score;
          if (loss <= 0 || (epoch >= 2 && loss >= ramp)) continue;
          const double eta = learning_rate /
            std::sqrt(1.0 + static_cast<double>(binary_iteration++) / std::max(1, train_n));
          const double shrink = std::max(0.0, 1.0 - eta * lambda);
          for (int f = 0; f < m; ++f) binary[f] = shrink * binary[f] + eta * sign * feature[f];
          binary[m] += eta * sign;
        }
      }
      for (int i = 0; i < train_n; ++i) {
        double& pa = probability[static_cast<std::size_t>(i) * global_classes + class_codes[class_a] - 1];
        double& pb = probability[static_cast<std::size_t>(i) * global_classes + class_codes[class_b] - 1];
        const double pair_mass = pa + pb;
        if (pair_mass <= 1e-12 || std::abs(pa - pb) > 0.35) continue;
        const double* binary = &binary_weights[static_cast<std::size_t>(fold[i]) * stride];
        const double* feature = &phi[static_cast<std::size_t>(i) * m];
        double score = binary[m];
        for (int f = 0; f < m; ++f) score += binary[f] * feature[f];
        const double specialist = 1.0 / (1.0 + std::exp(-std::max(-20.0, std::min(20.0, score))));
        const double revised = 0.70 * pa / pair_mass + 0.30 * specialist;
        pa = pair_mass * revised;
        pb = pair_mass * (1.0 - revised);
      }
    }
  }

  std::vector<double> taper(train_n, 1.0);
  for (int i = 0; i < train_n; ++i) {
    double product = 1.0;
    for (int d = 0; d < dims; ++d) {
      const double center = 0.5 * (tile.low[d] + tile.high[d]);
      const double radius = std::max(0.5 * (tile.high[d] - tile.low[d]), 1e-12);
      const double normalized = std::abs(x[tile.halo[i] + n * d] - center) / (2.0 * radius);
      product *= std::max(1e-4, 1.0 - normalized);
    }
    taper[i] = std::pow(product, 1.0 / dims);
  }
  return {tile.halo, probability, taper, classes, m};
}

std::vector<double> two_component_trust(const std::vector<double>& evidence,
                                        const std::vector<int>& rows) {
  std::vector<double> posterior(rows.size(), 0.5);
  if (rows.size() < 8) return posterior;
  std::vector<double> sorted;
  sorted.reserve(rows.size());
  for (int row : rows) sorted.push_back(evidence[row]);
  std::sort(sorted.begin(), sorted.end());
  double mean0 = sorted[sorted.size() / 4];
  double mean1 = sorted[(3 * sorted.size()) / 4];
  double variance0 = 0.025;
  double variance1 = 0.025;
  double mixture = 0.5;
  for (int iteration = 0; iteration < 30; ++iteration) {
    double weight0 = 0.0;
    double weight1 = 0.0;
    double sum0 = 0.0;
    double sum1 = 0.0;
    for (std::size_t i = 0; i < rows.size(); ++i) {
      const double value = evidence[rows[i]];
      const double density0 = (1.0 - mixture) / std::sqrt(variance0) *
        std::exp(-0.5 * (value - mean0) * (value - mean0) / variance0);
      const double density1 = mixture / std::sqrt(variance1) *
        std::exp(-0.5 * (value - mean1) * (value - mean1) / variance1);
      posterior[i] = density1 / std::max(density0 + density1, 1e-15);
      weight0 += 1.0 - posterior[i];
      weight1 += posterior[i];
      sum0 += (1.0 - posterior[i]) * value;
      sum1 += posterior[i] * value;
    }
    mean0 = sum0 / std::max(weight0, 1e-8);
    mean1 = sum1 / std::max(weight1, 1e-8);
    variance0 = 0.0;
    variance1 = 0.0;
    for (std::size_t i = 0; i < rows.size(); ++i) {
      const double value = evidence[rows[i]];
      variance0 += (1.0 - posterior[i]) * (value - mean0) * (value - mean0);
      variance1 += posterior[i] * (value - mean1) * (value - mean1);
    }
    variance0 = std::max(variance0 / std::max(weight0, 1e-8), 1e-4);
    variance1 = std::max(variance1 / std::max(weight1, 1e-8), 1e-4);
    mixture = std::min(0.98, std::max(0.02, weight1 / rows.size()));
  }
  if (mean1 < mean0) for (double& value : posterior) value = 1.0 - value;
  return posterior;
}

struct DisjointSet {
  std::vector<int> parent;
  explicit DisjointSet(int n) : parent(n) { std::iota(parent.begin(), parent.end(), 0); }
  int find(int value) {
    while (parent[value] != value) {
      parent[value] = parent[parent[value]];
      value = parent[value];
    }
    return value;
  }
  void unite(int a, int b) {
    a = find(a);
    b = find(b);
    if (a != b) parent[b] = a;
  }
};

std::vector<char> protect_rare_components(const std::vector<int>& labels,
                                          const std::vector<int>& samples,
                                          const std::vector<std::vector<int> >& neighbors,
                                          const std::vector<double>& trust,
                                          const std::vector<double>& stability,
                                          const std::vector<double>& tile_disagreement,
                                          int classes, int graph_k) {
  const int n = labels.size();
  DisjointSet components(n);
  for (int i = 0; i < n; ++i) {
    const int limit = std::min(graph_k, static_cast<int>(neighbors[i].size()));
    for (int j = 0; j < limit; ++j) {
      const int other = neighbors[i][j];
      if (samples[i] == samples[other] && labels[i] == labels[other]) components.unite(i, other);
    }
  }
  std::vector<int> class_count(classes, 0);
  for (int label : labels) ++class_count[label - 1];
  const int maximum = *std::max_element(class_count.begin(), class_count.end());
  struct Summary { int count = 0; double trust = 0; double stability = 0; double agreement = 0; };
  std::unordered_map<int, Summary> summaries;
  for (int i = 0; i < n; ++i) {
    Summary& summary = summaries[components.find(i)];
    ++summary.count;
    summary.trust += trust[i];
    summary.stability += stability[i];
    summary.agreement += 1.0 - tile_disagreement[i];
  }
  std::vector<char> protect(n, 0);
  for (int i = 0; i < n; ++i) {
    const Summary& summary = summaries[components.find(i)];
    const double average_trust = summary.trust / summary.count;
    const double average_stability = summary.stability / summary.count;
    const double average_agreement = summary.agreement / summary.count;
    const bool rare = static_cast<double>(class_count[labels[i] - 1]) / std::max(1, maximum) < 0.12;
    protect[i] = rare && summary.count >= 4 && average_trust > 0.72 &&
      average_stability > 0.75 && average_agreement > 0.75;
  }
  return protect;
}

void project_simplex(double* values, int classes) {
  std::vector<double> sorted(values, values + classes);
  std::sort(sorted.begin(), sorted.end(), std::greater<double>());
  double cumulative = 0.0;
  int rho = 0;
  for (int j = 0; j < classes; ++j) {
    cumulative += sorted[j];
    const double threshold = (cumulative - 1.0) / (j + 1);
    if (sorted[j] > threshold) rho = j + 1;
  }
  cumulative = 0.0;
  for (int j = 0; j < rho; ++j) cumulative += sorted[j];
  const double threshold = (cumulative - 1.0) / rho;
  for (int c = 0; c < classes; ++c) values[c] = std::max(0.0, values[c] - threshold);
}

std::vector<int> decode_tv(const std::vector<double>& probability,
                           const std::vector<int>& labels,
                           const std::vector<double>& local_support,
                           const std::vector<std::vector<int> >& neighbors,
                           const std::vector<std::vector<double> >& distances,
                           int classes, double retention, double tv_strength,
                           double probability_floor, double coherence,
                           int iterations, int graph_k, bool experimental_v2,
                           const std::vector<double>& trust,
                           const std::vector<char>& protected_component,
                           std::vector<double>& final_probability) {
  const int n = labels.size();
  std::vector<int> class_count(classes, 0);
  for (int label : labels) ++class_count[label - 1];
  const int maximum_count = *std::max_element(class_count.begin(), class_count.end());
  std::vector<char> candidate(static_cast<std::size_t>(n) * classes, 0);
  for (int i = 0; i < n; ++i) {
    candidate[static_cast<std::size_t>(i) * classes + labels[i] - 1] = 1;
    const int limit = std::min(graph_k, static_cast<int>(neighbors[i].size()));
    for (int j = 0; j < limit; ++j) {
      candidate[static_cast<std::size_t>(i) * classes + labels[neighbors[i][j]] - 1] = 1;
    }
  }

  std::vector<double> unary(static_cast<std::size_t>(n) * classes);
  for (int i = 0; i < n; ++i) {
    const double rare = std::min(2.5, std::sqrt(static_cast<double>(maximum_count) /
      std::max(1, class_count[labels[i] - 1])));
    const double coherent = std::max(0.0, (local_support[i] - 0.65) / 0.35);
    const double base_keep = retention *
      (0.25 + local_support[i] * local_support[i] * rare) +
      coherence * coherent * coherent;
    const double keep = experimental_v2
      ? base_keep + (protected_component[i] ? retention : 0.0)
      : base_keep;
    for (int c = 0; c < classes; ++c) {
      const std::size_t index = static_cast<std::size_t>(i) * classes + c;
      unary[index] = candidate[index]
        ? -std::log(std::max(probability[index], probability_floor)) : 1e3;
      if (c != labels[i] - 1) unary[index] += keep;
    }
  }

  std::vector<Edge> edges;
  std::vector<int> degree(n, 0);
  for (int i = 0; i < n; ++i) {
    const int limit = std::min(graph_k, static_cast<int>(neighbors[i].size()));
    for (int j = 0; j < limit; ++j) {
      const int other = neighbors[i][j];
      if (i >= other) continue;
      double difference = 0.0;
      for (int c = 0; c < classes; ++c) {
        difference += std::abs(probability[static_cast<std::size_t>(i) * classes + c] -
                               probability[static_cast<std::size_t>(other) * classes + c]);
      }
      const double border = std::min(1.0, 0.5 * difference);
      const double spatial = std::exp(-distances[i][j] /
        std::max(distances[i][limit - 1], 1e-12));
      edges.push_back({i, other, spatial * (1.0 - border) * (1.0 - border)});
      ++degree[i];
      ++degree[other];
    }
  }
  const int maximum_degree = std::max(1, *std::max_element(degree.begin(), degree.end()));
  const double step = 0.9 / std::sqrt(2.0 * maximum_degree);

  std::vector<double> u(probability);
  std::vector<double> previous(u);
  std::vector<double> extrapolated(u);
  std::vector<double> dual(static_cast<std::size_t>(edges.size()) * classes, 0.0);
  std::vector<double> divergence(static_cast<std::size_t>(n) * classes);
  for (int iteration = 0; iteration < iterations; ++iteration) {
    for (std::size_t e = 0; e < edges.size(); ++e) {
      const double bound = tv_strength * edges[e].spatial_weight;
      for (int c = 0; c < classes; ++c) {
        const double gradient = extrapolated[static_cast<std::size_t>(edges[e].a) * classes + c] -
          extrapolated[static_cast<std::size_t>(edges[e].b) * classes + c];
        double value = dual[e * classes + c] + step * gradient;
        dual[e * classes + c] = std::max(-bound, std::min(bound, value));
      }
    }
    std::fill(divergence.begin(), divergence.end(), 0.0);
    for (std::size_t e = 0; e < edges.size(); ++e) {
      for (int c = 0; c < classes; ++c) {
        const double value = dual[e * classes + c];
        divergence[static_cast<std::size_t>(edges[e].a) * classes + c] += value;
        divergence[static_cast<std::size_t>(edges[e].b) * classes + c] -= value;
      }
    }
    previous = u;
    for (int i = 0; i < n; ++i) {
      double* row = &u[static_cast<std::size_t>(i) * classes];
      for (int c = 0; c < classes; ++c) {
        const std::size_t index = static_cast<std::size_t>(i) * classes + c;
        row[c] -= step * (divergence[index] + unary[index]);
      }
      project_simplex(row, classes);
    }
    for (std::size_t index = 0; index < u.size(); ++index) {
      extrapolated[index] = 2.0 * u[index] - previous[index];
    }
  }

  final_probability = u;
  std::vector<int> output(n);
  for (int i = 0; i < n; ++i) {
    output[i] = 1 + std::distance(
      u.begin() + static_cast<std::size_t>(i) * classes,
      std::max_element(u.begin() + static_cast<std::size_t>(i) * classes,
                       u.begin() + static_cast<std::size_t>(i + 1) * classes));
  }
  return output;
}

}  // namespace structured_svm

extern "C" SEXP _SpatialGraphRefine_structured_spatial_svm_cpp(SEXP xy_s,
                                                                SEXP labels_s,
                                                                SEXP samples_s,
                                                                SEXP control_s,
                                                                SEXP verbose_s) {
  BEGIN_RCPP
  NumericMatrix xy(xy_s);
  IntegerVector labels_r(labels_s);
  IntegerVector samples_r(samples_s);
  List control(control_s);
  const bool verbose = as<bool>(verbose_s);
  const int n = xy.nrow();
  const int dims = xy.ncol();
  const int classes = *std::max_element(labels_r.begin(), labels_r.end());
  const int neighbors_k = as<int>(control["neighbors"]);
  const int target = as<int>(control["target_tile_size"]);
  const double overlap = as<double>(control["overlap"]);
  const double gamma = as<double>(control["gamma"]);
  const int landmarks = as<int>(control["landmarks"]);
  const int epochs = as<int>(control["epochs"]);
  const double learning_rate = as<double>(control["learning_rate"]);
  const double lambda = as<double>(control["lambda"]);
  const double ramp = as<double>(control["ramp"]);
  const double retention = as<double>(control["retention"]);
  const double tv_strength = as<double>(control["tv_strength"]);
  const double graph_mix = as<double>(control["graph_mix"]);
  const double probability_floor = as<double>(control["probability_floor"]);
  const double coherence = as<double>(control["coherence"]);
  const double preserve_support = as<double>(control["preserve_support"]);
  const bool topology_abstention = as<double>(control["topology_abstention"]) > 0.5;
  const bool adaptive_tiles = as<double>(control["adaptive_tiles"]) > 0.5;
  const bool cross_fitting = as<double>(control["cross_fitting"]) > 0.5;
  const bool experimental_v2 = as<double>(control["experimental_v2"]) > 0.5;
  const int trust_neighbors = as<int>(control["trust_neighbors"]);
  const double anisotropy = as<double>(control["anisotropy"]);
  const bool pairwise_specialists = as<double>(control["pairwise_specialists"]) > 0.5;
  const int pairwise_max = as<int>(control["pairwise_max"]);
  const double change_threshold = as<double>(control["change_threshold"]);
  const double unresolved_threshold = as<double>(control["unresolved_threshold"]);
  const int tv_iterations = as<int>(control["tv_iterations"]);
  const int workers_requested = as<int>(control["workers"]);
  const std::uint32_t seed = static_cast<std::uint32_t>(as<int>(control["seed"]));

  std::vector<double> x(xy.begin(), xy.end());
  std::vector<int> labels(labels_r.begin(), labels_r.end());
  std::vector<int> samples(samples_r.begin(), samples_r.end());
  std::vector<double> support(n, 1.0);
  std::vector<double> support6;
  std::vector<double> support24;
  std::vector<double> support48;
  if (experimental_v2) {
    support6.assign(n, 1.0);
    support24.assign(n, 1.0);
    support48.assign(n, 1.0);
  }
  std::vector<std::vector<int> > neighbors(n);
  std::vector<std::vector<double> > distances(n);
  std::vector<double> blended(static_cast<std::size_t>(n) * classes, 0.0);
  std::vector<double> blended_square;
  if (experimental_v2) {
    blended_square.assign(static_cast<std::size_t>(n) * classes, 0.0);
  }
  std::vector<double> total_weight(n, 0.0);
  int tile_count = 0;

  std::vector<int> sample_levels(samples);
  std::sort(sample_levels.begin(), sample_levels.end());
  sample_levels.erase(std::unique(sample_levels.begin(), sample_levels.end()), sample_levels.end());
  int tile_index = 0;
  for (int sample : sample_levels) {
    if (verbose) Rcout << "sample " << sample << "\n";
    std::vector<int> rows;
    for (int i = 0; i < n; ++i) if (samples[i] == sample) rows.push_back(i);
    if (rows.size() < 2) continue;
    structured_svm::Tree tree(x, n, dims, rows);
    const int k = std::min(experimental_v2 ? trust_neighbors : neighbors_k,
                           static_cast<int>(rows.size()) - 1);
    for (int row : rows) {
      if (experimental_v2) {
        std::vector<int> candidates;
        std::vector<double> candidate_distance;
        const int candidate_k = std::min(std::max(96, 2 * trust_neighbors),
                                         static_cast<int>(rows.size()) - 1);
        tree.query(row, candidate_k, candidates, candidate_distance);
        structured_svm::anisotropic_query(x, n, dims, row, candidates, candidate_distance,
                                          k, anisotropy, neighbors[row], distances[row]);
      } else {
        tree.query(row, k, neighbors[row], distances[row]);
      }
      auto support_at = [&](int requested) {
        const int limit = std::min(requested, static_cast<int>(neighbors[row].size()));
        int same = 0;
        for (int j = 0; j < limit; ++j) same += labels[neighbors[row][j]] == labels[row];
        return limit > 0 ? static_cast<double>(same) / limit : 1.0;
      };
      support[row] = support_at(neighbors_k);
      if (experimental_v2) {
        support6[row] = support_at(6);
        support24[row] = support_at(24);
        support48[row] = support_at(48);
      }
    }

    std::vector<double> low(dims, std::numeric_limits<double>::infinity());
    std::vector<double> high(dims, -std::numeric_limits<double>::infinity());
    for (int row : rows) {
      for (int d = 0; d < dims; ++d) {
        low[d] = std::min(low[d], x[row + n * d]);
        high[d] = std::max(high[d], x[row + n * d]);
      }
    }
    std::vector<structured_svm::Tile> tiles;
    const int mandatory_depth = static_cast<int>(std::ceil(std::log2(
      std::max(1.0, static_cast<double>(rows.size()) / target))));
    if (adaptive_tiles) {
      structured_svm::adaptive_split(
        x, n, dims, labels, classes, rows, low, high, 0, target,
        std::max(250, target / 4), mandatory_depth + 3, tiles);
    }
    structured_svm::add_regular_tiles(x, n, dims, rows, low, high, target, tiles);
    structured_svm::add_halos(tiles, x, n, dims, rows, overlap);

    const int first_tile = tile_index;
    const int workers = std::max(1, std::min(workers_requested,
      static_cast<int>(tiles.size())));
    std::vector<structured_svm::ModelResult> results(tiles.size());
    std::atomic<std::size_t> next_tile(0);
    auto fit_next = [&]() {
      while (true) {
        const std::size_t index = next_tile.fetch_add(1);
        if (index >= tiles.size()) return;
        results[index] = structured_svm::fit_tile(
          tiles[index], x, n, dims, labels, support, classes, landmarks, gamma,
          epochs, learning_rate, lambda, ramp, cross_fitting,
          experimental_v2 && pairwise_specialists, pairwise_max,
          seed + 104729U * (first_tile + index + 1));
      }
    };
    if (workers == 1) {
      fit_next();
    } else {
      std::vector<std::thread> pool;
      pool.reserve(workers);
      for (int worker = 0; worker < workers; ++worker) pool.emplace_back(fit_next);
      for (std::thread& worker : pool) worker.join();
    }

    for (std::size_t index = 0; index < tiles.size(); ++index) {
      Rcpp::checkUserInterrupt();
      ++tile_index;
      if (verbose) Rcout << "  tile " << tile_index << " halo=" << tiles[index].halo.size() << "\n";
      structured_svm::ModelResult& result = results[index];
      for (std::size_t i = 0; i < result.rows.size(); ++i) {
        const int row = result.rows[i];
        const double weight = result.taper[i];
        total_weight[row] += weight;
        for (int c = 0; c < classes; ++c) {
          const double value = result.probabilities[i * classes + c];
          blended[static_cast<std::size_t>(row) * classes + c] +=
            weight * value;
          if (experimental_v2) {
            blended_square[static_cast<std::size_t>(row) * classes + c] +=
              weight * value * value;
          }
        }
      }
      ++tile_count;
    }
  }

  for (int i = 0; i < n; ++i) {
    if (total_weight[i] <= 0) {
      blended[static_cast<std::size_t>(i) * classes + labels[i] - 1] = 1.0;
      total_weight[i] = 1.0;
    }
    for (int c = 0; c < classes; ++c) {
      blended[static_cast<std::size_t>(i) * classes + c] /= total_weight[i];
      if (experimental_v2) {
        blended_square[static_cast<std::size_t>(i) * classes + c] /= total_weight[i];
      }
    }
  }

  std::vector<double> tile_disagreement(n, 0.0);
  std::vector<double> stability(n, 1.0);
  std::vector<double> trust(n, 0.5);
  std::vector<double> sample_mean_trust(n, 0.5);
  std::vector<char> protected_component(n, 0);
  if (experimental_v2) {
    std::vector<double> evidence(n, 0.0);
    for (int i = 0; i < n; ++i) {
      double variance = 0.0;
      for (int c = 0; c < classes; ++c) {
        const std::size_t index = static_cast<std::size_t>(i) * classes + c;
        variance += std::max(0.0, blended_square[index] - blended[index] * blended[index]);
      }
      tile_disagreement[i] = std::min(1.0, std::sqrt(std::max(0.0, variance)));
      const int limit = std::min(6, static_cast<int>(neighbors[i].size()));
      double perturbation = 0.0;
      for (int j = 0; j < limit; ++j) {
        const int other = neighbors[i][j];
        for (int c = 0; c < classes; ++c) {
          perturbation += 0.5 * std::abs(
            blended[static_cast<std::size_t>(i) * classes + c] -
            blended[static_cast<std::size_t>(other) * classes + c]);
        }
      }
      stability[i] = limit > 0 ? std::max(0.0, 1.0 - perturbation / limit) : 1.0;
      const double svm_label = blended[static_cast<std::size_t>(i) * classes + labels[i] - 1];
      evidence[i] = 0.28 * svm_label + 0.18 * support6[i] + 0.12 * support[i] +
        0.08 * support24[i] + 0.06 * support48[i] +
        0.14 * (1.0 - tile_disagreement[i]) + 0.14 * stability[i];
    }
    for (int sample : sample_levels) {
      std::vector<int> rows;
      for (int i = 0; i < n; ++i) if (samples[i] == sample) rows.push_back(i);
      const std::vector<double> posterior = structured_svm::two_component_trust(evidence, rows);
      for (std::size_t i = 0; i < rows.size(); ++i) {
        const double absolute = std::max(0.0, std::min(1.0,
          (evidence[rows[i]] - 0.35) / 0.50));
        trust[rows[i]] = std::max(absolute, posterior[i]);
      }
      double support_sum = 0.0;
      double discord_square = 0.0;
      double discord_cross = 0.0;
      for (int row : rows) {
        const double discord = 1.0 - support[row];
        const int limit = std::min(neighbors_k, static_cast<int>(neighbors[row].size()));
        double neighbor_discord = 0.0;
        for (int j = 0; j < limit; ++j) neighbor_discord += 1.0 - support[neighbors[row][j]];
        if (limit > 0) neighbor_discord /= limit;
        support_sum += support[row];
        discord_square += discord * discord;
        discord_cross += discord * neighbor_discord;
      }
      const double mean_support = support_sum / std::max<std::size_t>(1, rows.size());
      const double autocorrelation = discord_cross / std::max(discord_square, 1e-12);
      const auto logistic = [](double value) { return 1.0 / (1.0 + std::exp(-value)); };
      const double spatial_coherence = logistic(24.0 * (autocorrelation - 0.84)) *
        logistic(18.0 * (mean_support - 0.60));
      const double clean_coherence = logistic(35.0 * (mean_support - 0.90));
      const double sample_reliability = std::max(spatial_coherence, clean_coherence);
      for (int row : rows) {
        trust[row] = 1.0 - (1.0 - trust[row]) * (1.0 - 0.95 * sample_reliability);
      }
      double trust_sum = 0.0;
      for (int row : rows) trust_sum += trust[row];
      const double average_trust = trust_sum / std::max<std::size_t>(1, rows.size());
      for (int row : rows) sample_mean_trust[row] = average_trust;
    }
    protected_component = structured_svm::protect_rare_components(
      labels, samples, neighbors, trust, stability, tile_disagreement, classes, neighbors_k);
    for (int i = 0; i < n; ++i) if (protected_component[i]) trust[i] = std::max(trust[i], 0.70);
  }

  // Add graph evidence only where the cross-fitted SVM is uncertain. This
  // preserves learned nonlinear borders while stabilizing locally noisy labels.
  if (graph_mix > 0) {
    std::vector<double> graph_probability(classes);
    for (int i = 0; i < n; ++i) {
      std::fill(graph_probability.begin(), graph_probability.end(), 0.0);
      double first = 0.0;
      double second = 0.0;
      for (int c = 0; c < classes; ++c) {
        const double value = blended[static_cast<std::size_t>(i) * classes + c];
        if (value > first) {
          second = first;
          first = value;
        } else if (value > second) {
          second = value;
        }
      }
      double graph_total = 0.0;
      const int limit = std::min(neighbors_k, static_cast<int>(neighbors[i].size()));
      const double scale = limit == 0 ? 1.0 :
        std::max(distances[i][limit - 1], 1e-12);
      for (int j = 0; j < limit; ++j) {
        const double weight = std::exp(-distances[i][j] / scale);
        graph_probability[labels[neighbors[i][j]] - 1] += weight;
        graph_total += weight;
      }
      if (graph_total <= 0) continue;
      const double alpha = graph_mix * (1.0 - std::max(0.0, first - second));
      for (int c = 0; c < classes; ++c) {
        const std::size_t index = static_cast<std::size_t>(i) * classes + c;
        blended[index] = (1.0 - alpha) * blended[index] +
          alpha * graph_probability[c] / graph_total;
      }
    }
  }

  std::vector<double> decoded_probability;
  std::vector<int> output = structured_svm::decode_tv(
    blended, labels, support, neighbors, distances, classes, retention,
    tv_strength, probability_floor, coherence, tv_iterations, neighbors_k,
    experimental_v2, trust, protected_component,
    decoded_probability);
  std::vector<int> decision(n, 0);
  std::vector<double> selective_risk(n, 0.0);
  if (experimental_v2) {
    structured_svm::DisjointSet change_components(n);
    std::vector<char> proposed_change(n, 0);
    for (int i = 0; i < n; ++i) proposed_change[i] = output[i] != labels[i];
    for (int i = 0; i < n; ++i) {
      if (!proposed_change[i]) continue;
      const int limit = std::min(neighbors_k, static_cast<int>(neighbors[i].size()));
      for (int j = 0; j < limit; ++j) {
        const int other = neighbors[i][j];
        if (proposed_change[other] && samples[i] == samples[other]) {
          change_components.unite(i, other);
        }
      }
    }
    struct ChangeSummary { int count = 0; double support = 0; double stability = 0; double agreement = 0; };
    std::unordered_map<int, ChangeSummary> change_summary;
    std::unordered_map<std::uint64_t, int> change_label_count;
    std::unordered_map<int, int> sample_size;
    for (int sample : samples) ++sample_size[sample];
    for (int i = 0; i < n; ++i) {
      if (!proposed_change[i]) continue;
      ChangeSummary& summary = change_summary[change_components.find(i)];
      ++summary.count;
      summary.support += support[i];
      summary.stability += stability[i];
      summary.agreement += 1.0 - tile_disagreement[i];
      const std::uint64_t key = (static_cast<std::uint64_t>(change_components.find(i)) << 32) |
        static_cast<std::uint32_t>(labels[i] - 1);
      ++change_label_count[key];
    }
    std::vector<char> unresolved_component(n, 0);
    for (int i = 0; i < n; ++i) {
      if (!proposed_change[i]) continue;
      const ChangeSummary& summary = change_summary[change_components.find(i)];
      const int minimum = std::max(20, static_cast<int>(0.0002 * sample_size[samples[i]]));
      int dominant = 0;
      for (int c = 0; c < classes; ++c) {
        const std::uint64_t key = (static_cast<std::uint64_t>(change_components.find(i)) << 32) |
          static_cast<std::uint32_t>(c);
        const auto found = change_label_count.find(key);
        if (found != change_label_count.end()) dominant = std::max(dominant, found->second);
      }
      unresolved_component[i] = summary.count >= minimum &&
        (static_cast<double>(dominant) / summary.count > 0.45 ||
         sample_mean_trust[i] > 0.83) &&
        summary.stability / summary.count > 0.55 &&
        summary.agreement / summary.count > 0.60;
    }
    for (int i = 0; i < n; ++i) {
      double first = -1.0;
      double second = -1.0;
      for (int c = 0; c < classes; ++c) {
        const double value = decoded_probability[static_cast<std::size_t>(i) * classes + c];
        if (value > first) {
          second = first;
          first = value;
        } else if (value > second) second = value;
      }
      const double decoded_margin = first - std::max(0.0, second);
      const int proposal = output[i] - 1;
      const int original = labels[i] - 1;
      const double svm_advantage = std::max(0.0,
        blended[static_cast<std::size_t>(i) * classes + proposal] -
        blended[static_cast<std::size_t>(i) * classes + original]);
      const double correction_score = (1.0 - trust[i]) *
        std::sqrt(std::max(0.0, decoded_margin * svm_advantage)) *
        (0.55 + 0.45 * stability[i]) *
        (1.0 - support[i]) * (1.0 - support[i]);
      selective_risk[i] = 1.0 - correction_score;
      const double coherence_scale = 1.0 + 9.0 /
        (1.0 + std::exp(-60.0 * (sample_mean_trust[i] - 0.82)));
      const double effective_change_threshold = change_threshold * coherence_scale;
      if (output[i] == labels[i]) {
        decision[i] = 0;
      } else if (unresolved_component[i]) {
        decision[i] = 2;
        output[i] = labels[i];
      } else if (correction_score >= effective_change_threshold &&
                 (!protected_component[i] || correction_score >= 0.80)) {
        decision[i] = 1;
      } else if (correction_score >= unresolved_threshold) {
        decision[i] = 2;
        output[i] = labels[i];
      } else {
        decision[i] = 0;
        output[i] = labels[i];
      }
    }
  }
  std::vector<int> abstained_samples;
  if (!experimental_v2 && topology_abstention) {
    for (int sample : sample_levels) {
      double support_sum = 0.0;
      double discord_square = 0.0;
      double discord_cross = 0.0;
      int count = 0;
      for (int i = 0; i < n; ++i) {
        if (samples[i] != sample) continue;
        const double discord = 1.0 - support[i];
        double neighbor_discord = 0.0;
        for (int other : neighbors[i]) neighbor_discord += 1.0 - support[other];
        if (!neighbors[i].empty()) neighbor_discord /= neighbors[i].size();
        support_sum += support[i];
        discord_square += discord * discord;
        discord_cross += discord * neighbor_discord;
        ++count;
      }
      if (count == 0) continue;
      const double mean_support = support_sum / count;
      const double autocorrelation = discord_cross /
        std::max(discord_square, 1e-12);
      const bool coherent =
        (autocorrelation > 0.84 && mean_support > 0.60) || mean_support > 0.90;
      if (!coherent) continue;
      abstained_samples.push_back(sample);
      for (int i = 0; i < n; ++i) {
        if (samples[i] != sample) continue;
        output[i] = labels[i];
        for (int c = 0; c < classes; ++c) {
          decoded_probability[static_cast<std::size_t>(i) * classes + c] =
            c == labels[i] - 1 ? 1.0 : 0.0;
        }
      }
    }
  }
  for (int i = 0; !experimental_v2 && i < n; ++i) {
    if (support[i] < preserve_support) continue;
    output[i] = labels[i];
    for (int c = 0; c < classes; ++c) {
      decoded_probability[static_cast<std::size_t>(i) * classes + c] =
        c == labels[i] - 1 ? 1.0 : 0.0;
    }
  }

  NumericVector confidence(n);
  NumericVector margin(n);
  for (int i = 0; i < n; ++i) {
    double first = -1.0;
    double second = -1.0;
    for (int c = 0; c < classes; ++c) {
      const double value = decoded_probability[static_cast<std::size_t>(i) * classes + c];
      if (value > first) {
        second = first;
        first = value;
      } else if (value > second) {
        second = value;
      }
    }
    confidence[i] = first;
    margin[i] = first - std::max(0.0, second);
  }

  IntegerVector refined(output.begin(), output.end());
  NumericVector support_r(support.begin(), support.end());
  NumericVector trust_r(trust.begin(), trust.end());
  NumericVector disagreement_r(tile_disagreement.begin(), tile_disagreement.end());
  NumericVector stability_r(stability.begin(), stability.end());
  NumericVector risk_r(selective_risk.begin(), selective_risk.end());
  IntegerVector decision_r(decision.begin(), decision.end());
  LogicalVector protected_r(n);
  for (int i = 0; i < n; ++i) protected_r[i] = protected_component[i];
  IntegerVector abstained_r(abstained_samples.begin(), abstained_samples.end());
  return List::create(
    _["labels"] = refined,
    _["confidence"] = confidence,
    _["margin"] = margin,
    _["local_support"] = support_r,
    _["trust"] = trust_r,
    _["tile_disagreement"] = disagreement_r,
    _["perturbation_stability"] = stability_r,
    _["selective_risk"] = risk_r,
    _["decision"] = decision_r,
    _["protected_component"] = protected_r,
    _["abstained_samples"] = abstained_r,
    _["tiles"] = tile_count
  );
  END_RCPP
}
