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
                           int iterations, int graph_k,
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
    const double keep = retention *
      (0.25 + local_support[i] * local_support[i] * rare) +
      coherence * coherent * coherent;
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

extern "C" SEXP _marginSVM_structured_spatial_svm_cpp(SEXP xy_s,
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
  const int tv_iterations = as<int>(control["tv_iterations"]);
  const int workers_requested = as<int>(control["workers"]);
  const std::uint32_t seed = static_cast<std::uint32_t>(as<int>(control["seed"]));

  std::vector<double> x(xy.begin(), xy.end());
  std::vector<int> labels(labels_r.begin(), labels_r.end());
  std::vector<int> samples(samples_r.begin(), samples_r.end());
  std::vector<double> support(n, 1.0);
  std::vector<std::vector<int> > neighbors(n);
  std::vector<std::vector<double> > distances(n);
  std::vector<double> blended(static_cast<std::size_t>(n) * classes, 0.0);
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
    const int k = std::min(neighbors_k, static_cast<int>(rows.size()) - 1);
    for (int row : rows) {
      tree.query(row, k, neighbors[row], distances[row]);
      int same = 0;
      for (int other : neighbors[row]) same += labels[other] == labels[row];
      support[row] = neighbors[row].empty()
        ? 1.0 : static_cast<double>(same) / neighbors[row].size();
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
    }
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
    decoded_probability);
  std::vector<int> abstained_samples;
  if (topology_abstention) {
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
  for (int i = 0; i < n; ++i) {
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
  IntegerVector abstained_r(abstained_samples.begin(), abstained_samples.end());
  return List::create(
    _["labels"] = refined,
    _["confidence"] = confidence,
    _["margin"] = margin,
    _["local_support"] = support_r,
    _["abstained_samples"] = abstained_r,
    _["tiles"] = tile_count
  );
  END_RCPP
}
