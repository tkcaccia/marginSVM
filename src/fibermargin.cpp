#include <Rcpp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <numeric>
#include <thread>
#include <utility>
#include <vector>

using namespace Rcpp;

namespace fiber_margin {

constexpr double kGoldenRatioConjugate = 0.6180339887498949;
constexpr double kSqrtTwoConjugate = 0.4142135623730951;
constexpr double kSqrtThreeConjugate = 0.7320508075688772;
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr int kViews2D = 7;
constexpr int kViews3D = 9;
// One fixed central transport range keeps the enclosure rule single-scale and
// avoids a multiscale ladder in the public operator.
constexpr std::array<double, 1> kScales = {5.0};

struct Settings {
  int threads = 1;
};

int list_int(const List& control, const char* name, int fallback) {
  if (!control.containsElementNamed(name)) return fallback;
  return as<int>(control[name]);
}

Settings parse_settings(const List& control) {
  if (control.size() > 1 ||
      (control.size() == 1 && !control.containsElementNamed("threads"))) {
    stop("FiberMargin has no tuning controls; use `workers` to set CPU parallelism.");
  }
  Settings settings;
  settings.threads = list_int(control, "threads", settings.threads);
  if (settings.threads < 1 || settings.threads > 64) {
    stop("FiberMargin threads must be between 1 and 64.");
  }
  return settings;
}

double quantile_linear(std::vector<double> values, double probability) {
  if (values.empty()) return 0.0;
  std::sort(values.begin(), values.end());
  const double position = probability * static_cast<double>(values.size() - 1);
  const std::size_t lower = static_cast<std::size_t>(std::floor(position));
  const std::size_t upper = static_cast<std::size_t>(std::ceil(position));
  const double fraction = position - static_cast<double>(lower);
  return values[lower] * (1.0 - fraction) + values[upper] * fraction;
}

double median_inplace(std::vector<double>& values) {
  if (values.empty()) return 0.0;
  const std::size_t middle = values.size() / 2;
  std::nth_element(values.begin(), values.begin() + middle, values.end());
  const double upper = values[middle];
  if (values.size() % 2 == 1) return upper;
  const double lower = *std::max_element(values.begin(), values.begin() + middle);
  return 0.5 * (lower + upper);
}

double median(std::vector<double> values) {
  return median_inplace(values);
}

std::vector<double> robust_unit_coordinates(
    const NumericMatrix& xy, const std::vector<int>& rows) {
  const int dimensions = xy.ncol();
  const int size = static_cast<int>(rows.size());
  std::vector<double> output(static_cast<std::size_t>(size) * dimensions);
  for (int dimension = 0; dimension < dimensions; ++dimension) {
    std::vector<double> values(size);
    for (int i = 0; i < size; ++i) values[i] = xy(rows[i], dimension);
    const double low = quantile_linear(values, 0.01);
    const double high = quantile_linear(values, 0.99);
    const double span = std::max(high - low, 1e-12);
    for (int i = 0; i < size; ++i) {
      const double value = (xy(rows[i], dimension) - low) / span;
      output[static_cast<std::size_t>(i) * dimensions + dimension] =
        std::max(-0.1, std::min(1.1, value));
    }
  }
  return output;
}

std::vector<double> rotation_matrix(int dimensions, int view) {
  if (dimensions == 2) {
    const double angle = kPi * std::fmod(view * kGoldenRatioConjugate, 1.0);
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    return {cosine, -sine, sine, cosine};
  }

  const double u1 = std::fmod(0.5 + view * kSqrtTwoConjugate, 1.0);
  const double u2 = std::fmod(0.25 + view * kGoldenRatioConjugate, 1.0);
  const double u3 = std::fmod(0.75 + view * kSqrtThreeConjugate, 1.0);
  const double x = std::sqrt(1.0 - u1) * std::sin(2.0 * kPi * u2);
  const double y = std::sqrt(1.0 - u1) * std::cos(2.0 * kPi * u2);
  const double z = std::sqrt(u1) * std::sin(2.0 * kPi * u3);
  const double w = std::sqrt(u1) * std::cos(2.0 * kPi * u3);
  return {
    1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - z * w),
    2.0 * (x * z + y * w), 2.0 * (x * y + z * w),
    1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - x * w),
    2.0 * (x * z - y * w), 2.0 * (y * z + x * w),
    1.0 - 2.0 * (x * x + y * y)
  };
}

std::vector<double> path_unit(const std::vector<double>& coordinates,
                              int size, int dimensions, int view) {
  const std::vector<double> rotation = rotation_matrix(dimensions, view);
  std::vector<double> transformed(static_cast<std::size_t>(size) * dimensions);
  std::vector<double> low(dimensions, std::numeric_limits<double>::infinity());
  std::vector<double> high(dimensions, -std::numeric_limits<double>::infinity());
  for (int i = 0; i < size; ++i) {
    for (int destination = 0; destination < dimensions; ++destination) {
      double value = 0.0;
      for (int source = 0; source < dimensions; ++source) {
        value += coordinates[static_cast<std::size_t>(i) * dimensions + source] *
          rotation[static_cast<std::size_t>(source) * dimensions + destination];
      }
      transformed[static_cast<std::size_t>(i) * dimensions + destination] = value;
      low[destination] = std::min(low[destination], value);
      high[destination] = std::max(high[destination], value);
    }
  }
  for (int i = 0; i < size; ++i) {
    for (int dimension = 0; dimension < dimensions; ++dimension) {
      const std::size_t index = static_cast<std::size_t>(i) * dimensions + dimension;
      transformed[index] = std::max(0.0, std::min(
        1.0, (transformed[index] - low[dimension]) /
          std::max(high[dimension] - low[dimension], 1e-12)));
    }
  }
  return transformed;
}

std::uint64_t hilbert_code(const double* coordinate, int dimensions, int bits = 16) {
  std::array<std::uint64_t, 3> transformed = {0, 0, 0};
  const std::uint64_t maximum = (static_cast<std::uint64_t>(1) << bits) - 1;
  for (int dimension = 0; dimension < dimensions; ++dimension) {
    const double clipped = std::max(0.0, std::min(1.0 - 1e-12, coordinate[dimension]));
    transformed[dimension] =
      static_cast<std::uint64_t>(std::floor(clipped * maximum));
  }

  std::uint64_t q = static_cast<std::uint64_t>(1) << (bits - 1);
  while (q > 1) {
    const std::uint64_t p = q - 1;
    for (int axis = 0; axis < dimensions; ++axis) {
      if ((transformed[axis] & q) != 0) {
        transformed[0] ^= p;
      } else {
        const std::uint64_t exchange = (transformed[0] ^ transformed[axis]) & p;
        transformed[0] ^= exchange;
        transformed[axis] ^= exchange;
      }
    }
    q >>= 1;
  }
  for (int axis = 1; axis < dimensions; ++axis) {
    transformed[axis] ^= transformed[axis - 1];
  }
  std::uint64_t correction = 0;
  q = static_cast<std::uint64_t>(1) << (bits - 1);
  while (q > 1) {
    if ((transformed[dimensions - 1] & q) != 0) correction ^= q - 1;
    q >>= 1;
  }
  for (int axis = 0; axis < dimensions; ++axis) transformed[axis] ^= correction;

  std::uint64_t code = 0;
  for (int bit = 0; bit < bits; ++bit) {
    for (int axis = 0; axis < dimensions; ++axis) {
      code |= ((transformed[axis] >> bit) & 1ULL) <<
        (bit * dimensions + dimensions - 1 - axis);
    }
  }
  return code;
}

struct Geometry {
  std::vector<int> order;
  std::vector<double> isolation;
  std::vector<double> forward_action;
  std::vector<double> reverse_action;
};

Geometry path_geometry(const std::vector<double>& coordinates, int size,
                       int dimensions, int view) {
  Geometry geometry;
  geometry.order.resize(size);
  geometry.isolation.assign(size, 1.0);
  geometry.forward_action.assign(size, 0.0);
  geometry.reverse_action.assign(size, 0.0);
  std::iota(geometry.order.begin(), geometry.order.end(), 0);

  const std::vector<double> unit = path_unit(coordinates, size, dimensions, view);
  std::vector<std::uint64_t> code(size);
  for (int i = 0; i < size; ++i) {
    const double* row = &unit[static_cast<std::size_t>(i) * dimensions];
    code[i] = hilbert_code(row, dimensions);
  }
  std::sort(geometry.order.begin(), geometry.order.end(), [&](int left, int right) {
    return code[left] < code[right] || (code[left] == code[right] && left < right);
  });

  if (size < 2) return geometry;
  std::vector<double> step(size - 1);
  for (int position = 1; position < size; ++position) {
    const int left = geometry.order[position - 1];
    const int right = geometry.order[position];
    double squared = 0.0;
    for (int dimension = 0; dimension < dimensions; ++dimension) {
      const double delta =
        coordinates[static_cast<std::size_t>(left) * dimensions + dimension] -
        coordinates[static_cast<std::size_t>(right) * dimensions + dimension];
      squared += delta * delta;
    }
    step[position - 1] = std::sqrt(squared);
  }
  std::vector<double> positive;
  positive.reserve(step.size());
  for (double value : step) if (value > 1e-12) positive.push_back(value);
  const double reference = positive.empty() ? 1.0 : median(positive);

  // The incoming path edge is the local transport scale at its destination.
  std::vector<double> local_scale = step;
  for (double& value : local_scale) value = std::max(value, 1e-12);

  std::vector<double> point_scale(size, reference);
  point_scale[0] = local_scale.front();
  point_scale[size - 1] = local_scale.back();
  for (int position = 1; position < size - 1; ++position) {
    point_scale[position] = std::sqrt(local_scale[position - 1] * local_scale[position]);
  }

  // Transport always follows normalized geometric path length. A location is
  // isolated only when its adjacent path gaps exceed the median gap; that
  // factor is used solely to protect the final decision at the query.
  for (int position = 0; position < size; ++position) {
    const double gap = point_scale[position] / std::max(reference, 1e-12);
    geometry.isolation[geometry.order[position]] = std::max(1.0, gap);
  }
  for (int edge = 0; edge < size - 1; ++edge) {
    const double action = step[edge] / std::max(reference, 1e-12);
    geometry.forward_action[geometry.order[edge + 1]] = action;
    geometry.reverse_action[geometry.order[edge]] = action;
  }
  return geometry;
}

// A class is spatially enclosed at a path position when its transported mass
// arrives from both directions. The squared sum of directional amplitudes
// keeps unilateral support while giving an explicit bonus to balanced support:
// (sqrt(left) + sqrt(right))^2.
void two_sided_enclosure_field(const std::vector<double>& source_indicator,
                               int size, int classes,
                               const Geometry& geometry, double scale,
                               std::vector<double>& field) {
  const std::size_t width = static_cast<std::size_t>(size) * classes;
  std::vector<double> forward_field(width, 0.0);
  field.assign(width, 0.0);
  std::vector<double> state(classes, 0.0);
  for (int row : geometry.order) {
    const double decay = std::exp(-geometry.forward_action[row] / scale);
    const std::size_t offset = static_cast<std::size_t>(row) * classes;
    for (int cls = 0; cls < classes; ++cls) {
      state[cls] *= decay;
      forward_field[offset + cls] = state[cls];
    }
    for (int cls = 0; cls < classes; ++cls) {
      state[cls] += source_indicator[offset + cls];
    }
  }
  std::fill(state.begin(), state.end(), 0.0);
  for (auto iterator = geometry.order.rbegin(); iterator != geometry.order.rend(); ++iterator) {
    const int row = *iterator;
    const double decay = std::exp(-geometry.reverse_action[row] / scale);
    const std::size_t offset = static_cast<std::size_t>(row) * classes;
    for (int cls = 0; cls < classes; ++cls) {
      state[cls] *= decay;
      const double forward = forward_field[offset + cls];
      const double reverse = state[cls];
      field[offset + cls] = forward + reverse +
        2.0 * std::sqrt(forward * reverse);
    }
    for (int cls = 0; cls < classes; ++cls) {
      state[cls] += source_indicator[offset + cls];
    }
  }
}

bool has_spatial_extent(const std::vector<double>& coordinates, int size,
                        int dimensions) {
  if (size < 2) return false;
  for (int dimension = 0; dimension < dimensions; ++dimension) {
    const double first = coordinates[dimension];
    for (int row = 1; row < size; ++row) {
      if (coordinates[static_cast<std::size_t>(row) * dimensions + dimension] != first) {
        return true;
      }
    }
  }
  return false;
}

struct SampleResult {
  std::vector<int> labels;
  std::vector<int> candidate;
  std::vector<double> margin_score;
  std::vector<double> required;
  std::vector<double> dispersion;
  std::vector<double> isolation;
};

// Each spatial route supplies pointwise candidate evidence. The final label
// decision remains separate from its construction, so binary ballots and
// multiclass path transport share one selective margin rule.
struct MarginEvidence {
  std::vector<int> candidate;
  std::vector<double> margin_score;
  std::vector<double> required;
  std::vector<double> dispersion;
  std::vector<double> isolation;
};

struct NeighborNode {
  int row;
  int axis;
  int left;
  int right;
};

class NeighborTree {
 public:
  NeighborTree(const std::vector<double>& coordinates, int dimensions,
               const std::vector<int>& rows)
      : coordinates_(coordinates), dimensions_(dimensions), order_(rows) {
    nodes_.reserve(rows.size());
    root_ = build(0, static_cast<int>(order_.size()), 0);
  }

  void binary_counts(int target, int count, const std::vector<int>& labels,
                     std::vector<std::pair<double, int> >& heap,
                     int& count_zero, int& count_one,
                     int& nearest_label) const {
    heap.clear();
    search(root_, target, count, heap);
    count_zero = 0;
    count_one = 0;
    std::pair<double, int> nearest = heap.front();
    for (const std::pair<double, int>& entry : heap) {
      if (labels[entry.second] == 0) {
        ++count_zero;
      } else {
        ++count_one;
      }
      if (entry < nearest) nearest = entry;
    }
    nearest_label = labels[nearest.second];
  }

 private:
  const std::vector<double>& coordinates_;
  int dimensions_;
  int root_;
  std::vector<int> order_;
  std::vector<NeighborNode> nodes_;

  double coordinate(int row, int axis) const {
    return coordinates_[static_cast<std::size_t>(row) * dimensions_ + axis];
  }

  int build(int begin, int end, int depth) {
    if (begin >= end) return -1;
    const int axis = depth % dimensions_;
    const int middle = begin + (end - begin) / 2;
    std::nth_element(
      order_.begin() + begin, order_.begin() + middle, order_.begin() + end,
      [&](int left, int right) { return coordinate(left, axis) < coordinate(right, axis); }
    );
    const int node = static_cast<int>(nodes_.size());
    nodes_.push_back({order_[middle], axis, -1, -1});
    nodes_[node].left = build(begin, middle, depth + 1);
    nodes_[node].right = build(middle + 1, end, depth + 1);
    return node;
  }

  double squared_distance(int left, int right) const {
    double total = 0.0;
    for (int axis = 0; axis < dimensions_; ++axis) {
      const double delta = coordinate(left, axis) - coordinate(right, axis);
      total += delta * delta;
    }
    return total;
  }

  void search(int node_index, int target, int count,
              std::vector<std::pair<double, int> >& heap) const {
    if (node_index < 0) return;
    const NeighborNode& node = nodes_[node_index];
    const double difference = coordinate(target, node.axis) - coordinate(node.row, node.axis);
    const int near_child = difference <= 0.0 ? node.left : node.right;
    const int far_child = difference <= 0.0 ? node.right : node.left;
    search(near_child, target, count, heap);
    if (node.row != target) {
      const double distance = squared_distance(target, node.row);
      if (static_cast<int>(heap.size()) < count) {
        heap.push_back(std::make_pair(distance, node.row));
        std::push_heap(heap.begin(), heap.end());
      } else if (distance < heap.front().first) {
        std::pop_heap(heap.begin(), heap.end());
        heap.back() = std::make_pair(distance, node.row);
        std::push_heap(heap.begin(), heap.end());
      }
    }
    const double bound = static_cast<int>(heap.size()) < count ?
      std::numeric_limits<double>::infinity() : heap.front().first;
    if (difference * difference <= bound) search(far_child, target, count, heap);
  }
};

template <typename Function>
void parallel_ranges(int size, int threads, const Function& function) {
  const int workers = std::max(1, std::min(size, threads));
  if (workers == 1) {
    function(0, size);
    return;
  }
  std::vector<std::thread> pool;
  pool.reserve(workers);
  for (int worker = 0; worker < workers; ++worker) {
    const int begin = size * worker / workers;
    const int end = size * (worker + 1) / workers;
    pool.emplace_back([begin, end, &function]() {
      function(begin, end);
    });
  }
  for (std::thread& worker : pool) worker.join();
}

SampleResult apply_margin(const std::vector<int>& labels,
                          MarginEvidence evidence,
                          const Settings& settings) {
  const int size = static_cast<int>(labels.size());
  SampleResult result;
  result.labels = labels;
  result.candidate = std::move(evidence.candidate);
  result.margin_score = std::move(evidence.margin_score);
  result.required = std::move(evidence.required);
  result.dispersion = std::move(evidence.dispersion);
  result.isolation = std::move(evidence.isolation);
  parallel_ranges(size, settings.threads, [&](int begin, int end) {
    for (int row = begin; row < end; ++row) {
      if (result.candidate[row] != labels[row] &&
          result.margin_score[row] >= result.required[row]) {
        result.labels[row] = result.candidate[row];
      }
    }
  });
  return result;
}

MarginEvidence binary_ballot_evidence(
    const std::vector<double>& coordinates, const std::vector<int>& labels,
    int dimensions, const Settings& settings) {
  const int size = static_cast<int>(labels.size());
  MarginEvidence evidence;
  evidence.candidate = labels;
  evidence.margin_score.assign(size, 0.0);
  // Integer ballot scores require one net local-volume vote to displace the
  // supplied category. Exact ties therefore retain that category.
  evidence.required.assign(size, 1.0);
  evidence.dispersion.assign(size, 0.0);
  evidence.isolation.assign(size, 1.0);
  if (size < 3) return evidence;

  std::vector<int> rows(size);
  std::iota(rows.begin(), rows.end(), 0);
  const NeighborTree tree(coordinates, dimensions, rows);
  const int local_mass = std::min(size - 1, 1 << (dimensions + 4));
  parallel_ranges(size, settings.threads, [&](int begin, int end) {
    std::vector<std::pair<double, int> > heap;
    heap.reserve(local_mass);
    for (int row = begin; row < end; ++row) {
      int count_zero = 0;
      int count_one = 0;
      int nearest_label = labels[row];
      tree.binary_counts(
        row, local_mass, labels, heap, count_zero, count_one, nearest_label
      );
      // A tied local volume has no preferred rival. Retaining the observed
      // class makes the binary evidence obey the same candidate-first rule as
      // the multiclass branch and prevents diagnostic-only tie choices.
      const int rival = count_zero > count_one ? 0 :
        (count_one > count_zero ? 1 : labels[row]);
      const int rival_count = rival == 0 ? count_zero : count_one;
      const int observed_count = labels[row] == 0 ? count_zero : count_one;
      evidence.candidate[row] = rival;
      evidence.margin_score[row] = rival == labels[row] ? 0.0 :
        static_cast<double>(rival_count - observed_count);
    }
  });
  return evidence;
}

SampleResult refine_sample(const std::vector<double>& coordinates,
                           const std::vector<int>& labels, int dimensions,
                           int classes, const Settings& settings) {
  const int size = static_cast<int>(labels.size());
  MarginEvidence evidence;
  evidence.candidate = labels;
  evidence.margin_score.assign(size, 0.0);
  evidence.required.assign(size, 0.0);
  evidence.dispersion.assign(size, 0.0);
  evidence.isolation.assign(size, 1.0);
  // With no coordinate extent, a path order is arbitrary. Retaining the
  // supplied field avoids creating order-dependent spatial evidence.
  if (!has_spatial_extent(coordinates, size, dimensions)) {
    return apply_margin(labels, std::move(evidence), settings);
  }
  if (size < 20) return apply_margin(labels, std::move(evidence), settings);
  std::vector<double> count(classes, 0.0);
  int observed_classes = 0;
  for (int label : labels) {
    if (count[label] == 0.0) ++observed_classes;
    count[label] += 1.0;
  }
  if (observed_classes < 2) {
    return apply_margin(labels, std::move(evidence), settings);
  }
  if (classes == 2) {
    return apply_margin(
      labels, binary_ballot_evidence(coordinates, labels, dimensions, settings), settings);
  }

  const int views = dimensions == 2 ? kViews2D : kViews3D;
  const int charts = views;
  const std::size_t row_class = static_cast<std::size_t>(size) * classes;
  std::vector<double> observed_source(row_class, 0.0);
  for (int row = 0; row < size; ++row) {
    observed_source[static_cast<std::size_t>(row) * classes + labels[row]] = 1.0;
  }
  std::vector<Geometry> atlas(charts);
  checkUserInterrupt();
  parallel_ranges(charts, std::min(settings.threads, charts), [&](int begin, int end) {
    for (int view = begin; view < end; ++view) {
      atlas[view] = path_geometry(coordinates, size, dimensions, view);
    }
  });
  std::vector<float> chart_action;
  std::vector<float> chart_isolation(static_cast<std::size_t>(charts) * size);
  const auto compute_atlas = [&]() {
    chart_action.assign(static_cast<std::size_t>(charts) * row_class, 0.0F);
    const auto compute_view = [&](int view) {
      const Geometry& geometry = atlas[view];
      for (int row = 0; row < size; ++row) {
        chart_isolation[static_cast<std::size_t>(view) * size + row] =
          static_cast<float>(geometry.isolation[row]);
      }

      std::vector<double> mean_margin(row_class, 0.0);
      const auto add_enclosure_margin = [&]() {
        std::vector<double> field;
        for (double scale : kScales) {
          two_sided_enclosure_field(
            observed_source, size, classes, geometry, scale, field);
          for (int row = 0; row < size; ++row) {
            const std::size_t offset = static_cast<std::size_t>(row) * classes;
            double total = 0.0;
            for (int cls = 0; cls < classes; ++cls) total += field[offset + cls];
            total = std::max(total, 1e-12);
            const double observed = field[offset + labels[row]] / total;
            for (int cls = 0; cls < classes; ++cls) {
              const std::size_t index = offset + cls;
              mean_margin[index] += field[index] / total - observed;
            }
          }
        }
        const double scale_count = static_cast<double>(kScales.size());
        for (double& value : mean_margin) value /= scale_count;
      };

      add_enclosure_margin();
      for (std::size_t index = 0; index < row_class; ++index) {
        chart_action[static_cast<std::size_t>(view) * row_class + index] =
          static_cast<float>(mean_margin[index]);
      }
    };

    checkUserInterrupt();
    const int thread_count = std::min(settings.threads, views);
    if (thread_count == 1) {
      for (int view = 0; view < views; ++view) compute_view(view);
    } else {
      std::vector<std::thread> pool;
      std::vector<std::exception_ptr> errors(thread_count);
      pool.reserve(thread_count);
      for (int thread = 0; thread < thread_count; ++thread) {
        pool.emplace_back([&, thread]() {
          try {
            for (int view = thread; view < views; view += thread_count) {
              compute_view(view);
            }
          } catch (...) {
            errors[thread] = std::current_exception();
          }
        });
      }
      for (std::thread& worker : pool) worker.join();
      for (const std::exception_ptr& error : errors) {
        if (error) std::rethrow_exception(error);
      }
    }
  };

  const auto aggregate_atlas = [&](std::vector<double>& full_margin,
                                   std::vector<double>& full_dispersion) {
    std::vector<double> values(charts);
    full_margin.assign(row_class, 0.0);
    full_dispersion.assign(row_class, 0.0);
    for (std::size_t index = 0; index < row_class; ++index) {
      for (int atlas = 0; atlas < charts; ++atlas) {
        values[atlas] = chart_action[static_cast<std::size_t>(atlas) * row_class + index];
      }
      double mean = 0.0;
      for (int atlas = 0; atlas < charts; ++atlas) mean += values[atlas];
      mean /= static_cast<double>(charts);
      double variance = 0.0;
      for (int atlas = 0; atlas < charts; ++atlas) {
        const double delta = values[atlas] - mean;
        variance += delta * delta;
      }
      variance /= static_cast<double>(charts);
      full_margin[index] = mean;
      // This fixed two-sided chart-dispersion radius is twice the deterministic
      // atlas deviation sqrt(sum((m_a - mean)^2)) / A, not a sampling SE.
      full_dispersion[index] = 2.0 * std::sqrt(variance) /
        std::sqrt(static_cast<double>(charts));
    }
  };

  compute_atlas();
  std::vector<double> full_margin;
  std::vector<double> full_dispersion;
  aggregate_atlas(full_margin, full_dispersion);

  std::vector<double> positive;
  for (double value : count) if (value > 0.0) positive.push_back(value);
  const double reference = std::max(median(positive), 1.0);
  std::vector<double> deficit(classes, 1.0);
  for (int cls = 0; cls < classes; ++cls) {
    deficit[cls] = std::sqrt(reference / std::max(count[cls], 1.0));
  }

  const auto required_barrier = [&](int row, int rival, double dispersion) {
    return dispersion / deficit[rival] * evidence.isolation[row];
  };
  const auto best_rival = [&](const std::vector<double>& margin, int row) {
    const std::size_t offset = static_cast<std::size_t>(row) * classes;
    // Exact evidence ties retain the supplied label. This prevents a
    // zero-evidence change that would otherwise depend on class encoding.
    int rival = labels[row];
    for (int cls = 0; cls < classes; ++cls) {
      if (margin[offset + cls] > margin[offset + rival]) rival = cls;
    }
    return rival;
  };

  parallel_ranges(size, settings.threads, [&](int begin, int end) {
    std::vector<double> values(charts);
    for (int row = begin; row < end; ++row) {
      for (int atlas = 0; atlas < charts; ++atlas) {
        values[atlas] = chart_isolation[static_cast<std::size_t>(atlas) * size + row];
      }
      evidence.isolation[row] = median_inplace(values);

      const std::size_t offset = static_cast<std::size_t>(row) * classes;
      const int rival = best_rival(full_margin, row);
      const double required = required_barrier(
        row, rival, full_dispersion[offset + rival]);
      evidence.candidate[row] = rival;
      evidence.margin_score[row] = full_margin[offset + rival];
      evidence.required[row] = required;
      evidence.dispersion[row] = full_dispersion[offset + rival];
    }
  });
  return apply_margin(labels, std::move(evidence), settings);
}

}  // namespace fiber_margin

extern "C" SEXP _fibermargin_fiber_margin_cpp(SEXP xy_s, SEXP labels_s,
                                              SEXP samples_s, SEXP control_s) {
  BEGIN_RCPP
  NumericMatrix xy(xy_s);
  IntegerVector labels_r(labels_s);
  IntegerVector samples_r(samples_s);
  List control(control_s);
  const int n = xy.nrow();
  const int dimensions = xy.ncol();
  if (dimensions < 2 || dimensions > 3) stop("FiberMargin requires two or three coordinates.");
  if (labels_r.size() != n || samples_r.size() != n) stop("FiberMargin input lengths differ.");

  int classes = 0;
  for (int value : labels_r) {
    if (value == NA_INTEGER || value < 1) stop("FiberMargin labels must be observed factors.");
    classes = std::max(classes, value);
  }
  const fiber_margin::Settings settings = fiber_margin::parse_settings(control);
  IntegerVector output = clone(labels_r);
  IntegerVector candidate = clone(labels_r);
  NumericVector margin_score(n, 0.0);
  NumericVector required(n, 0.0);
  NumericVector dispersion(n, 0.0);
  NumericVector isolation(n, 1.0);

  std::vector<int> sample_levels;
  if (n > 0) {
    const int first_sample = samples_r[0];
    bool single_sample = true;
    for (int row = 1; row < n; ++row) {
      if (samples_r[row] != first_sample) {
        single_sample = false;
        break;
      }
    }
    if (single_sample) {
      sample_levels.push_back(first_sample);
    } else {
      sample_levels.assign(samples_r.begin(), samples_r.end());
      std::sort(sample_levels.begin(), sample_levels.end());
      sample_levels.erase(
        std::unique(sample_levels.begin(), sample_levels.end()), sample_levels.end());
    }
  }
  for (int sample : sample_levels) {
    checkUserInterrupt();
    std::vector<int> rows;
    for (int i = 0; i < n; ++i) if (samples_r[i] == sample) rows.push_back(i);
    std::vector<int> local_labels(rows.size());
    for (std::size_t i = 0; i < rows.size(); ++i) local_labels[i] = labels_r[rows[i]] - 1;
    const std::vector<double> coordinates =
      fiber_margin::robust_unit_coordinates(xy, rows);
    const fiber_margin::SampleResult result = fiber_margin::refine_sample(
      coordinates, local_labels, dimensions, classes, settings);
    for (std::size_t i = 0; i < rows.size(); ++i) {
      const int row = rows[i];
      output[row] = result.labels[i] + 1;
      candidate[row] = result.candidate[i] + 1;
      margin_score[row] = result.margin_score[i];
      required[row] = result.required[i];
      dispersion[row] = result.dispersion[i];
      isolation[row] = result.isolation[i];
    }
  }
  return List::create(
    _["labels"] = output,
    _["candidate"] = candidate,
    _["margin_score"] = margin_score,
    _["required"] = required,
    _["atlas_dispersion"] = dispersion,
    _["isolation"] = isolation,
    _["changed"] = output != labels_r
  );
  END_RCPP
}
