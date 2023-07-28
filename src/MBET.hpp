#include <limits>
#include <list>
#include <map>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include "common/bigraph.hpp"
#include "common/timer.hpp"

// #define MBETM
// #define ONLY_COUNTING
// #define OUTPUT_RESULT

static constexpr size_t NumBits = std::numeric_limits<size_t>::digits;
static constexpr size_t HighestBit = 1UL << (NumBits - 1);
using CompressedNeighbors = std::vector<std::pair<size_t, size_t>>;

struct ResultTreeNode {
  explicit ResultTreeNode() : is_end(false), has_prev(false) {}

  ~ResultTreeNode() {
    for (auto &[k, v] : children) v.reset();
  }

  std::map<vid_t, std::unique_ptr<ResultTreeNode>> children;
  bool is_end;
  bool has_prev;
};

struct NeighborTreeNode {
  explicit NeighborTreeNode(vid_t id, NeighborTreeNode *pa)
      : vid(id), is_end(false), is_res(false), is_use(false), sub_cnt(0), prev(nullptr), fork(nullptr), parent(pa) {}

  ~NeighborTreeNode() {
    for (auto &[k, v] : children) v.reset();
  }

  vid_t vid;
  std::map<vid_t, std::unique_ptr<NeighborTreeNode>> children;
  bool is_end;
  bool is_res;
  bool is_use;
  size_t sub_cnt;
  NeighborTreeNode *prev;
  std::pair<size_t, size_t> dfs_idx;
  NeighborTreeNode *fork;
  NeighborTreeNode *parent;

  friend std::ostream &operator<<(std::ostream &os, const NeighborTreeNode &n) {
    os << "[";
    os << "&=" << &n << ", ";
    os << "vid=" << n.vid << ", ";
    os << "children=(";
    for (const auto &[id, child] : n.children) {
      os << id << ",";
    }
    os << "), ";
    os << "is_end=" << n.is_end << ", ";
    os << "sub_cnt=" << n.sub_cnt << ", ";
    os << "prev=" << n.prev << ", ";
    os << "dfs_idx=" << n.dfs_idx.first << "|" << n.dfs_idx.second << ", ";
    return os;
  }
};

class MBET {
 public:
  explicit MBET(BiGraph bi_graph) : bi_graph_(std::move(bi_graph)), swap_vertices_(false), bic_cnt_(0) {}

  std::tuple<size_t, double, std::map<std::string, double>> Run() {
    main_watch_.Start();
    std::map<std::string, double> times;
    auto run_watch = RunTimer::StopWatch();

    {
      // renumber V
      std::vector<std::pair<vid_t, size_t>> Vs;
      for (auto &v : bi_graph_.VerticesV()) {
        Vs.emplace_back(v->vid, v->neighbors.size());
      }
      std::sort(Vs.begin(), Vs.end(), [](const auto &a, const auto &b) {
        return a.second > b.second || (a.second == b.second && a.first < b.first);
      });
      for (size_t i = 0; i < Vs.size(); i++) bi_graph_.VertexV(Vs[i].first)->vid = i;
      std::sort(bi_graph_.VerticesV().begin(), bi_graph_.VerticesV().end(),
                [](const std::shared_ptr<Vertex> &a, const std::shared_ptr<Vertex> &b) { return a->vid < b->vid; });
      for (auto &u : bi_graph_.VerticesU()) {
        std::sort(u->neighbors.begin(), u->neighbors.end(),
                  [](const std::shared_ptr<Vertex> &a, const std::shared_ptr<Vertex> &b) { return a->vid < b->vid; });
      }
    }
    {
      // renumber U
      std::vector<std::pair<vid_t, size_t>> Us;
      for (auto &u : bi_graph_.VerticesU()) {
        Us.emplace_back(u->vid, u->neighbors.size());
      }
      std::sort(Us.begin(), Us.end(), [](const auto &a, const auto &b) {
        return a.second > b.second || (a.second == b.second && a.first < b.first);
      });
      for (size_t i = 0; i < Us.size(); i++) bi_graph_.VertexU(Us[i].first)->vid = i;
      std::sort(bi_graph_.VerticesU().begin(), bi_graph_.VerticesU().end(),
                [](const std::shared_ptr<Vertex> &a, const std::shared_ptr<Vertex> &b) { return a->vid < b->vid; });
      for (auto &v : bi_graph_.VerticesV()) {
        std::sort(v->neighbors.begin(), v->neighbors.end(),
                  [](const std::shared_ptr<Vertex> &a, const std::shared_ptr<Vertex> &b) { return a->vid < b->vid; });
      }
    }
    if (std::getenv("ALLOW_SWAP") &&
        bi_graph_.VerticesV().front()->neighbors.size() < bi_graph_.VerticesU().front()->neighbors.size()) {
      std::cout << "[DBG] swap vertices U <=> V" << std::endl;
      std::swap(bi_graph_.VerticesU(), bi_graph_.VerticesV());
      swap_vertices_ = true;
    }

    times.emplace("sort", run_watch.Stop());
    run_watch.Report("MBET::Pre-sort");
    run_watch.Start();

    bic_cnt_ = 0;
    duplicated_cnt_ = 0;
    cand_tree_cnt_ = 0;
    result_tree_size_ = 0;
    result_tree_cnt_ = 0;
    cand_tree_size_cur_ = 0;
    cand_tree_size_max_ = 0;
    cand_tree_size_total_ = 0;
    result_tree_size_cur_ = 0;
    result_tree_size_max_ = 0;
    result_tree_size_total_ = 0;

    result_tree_roots_.resize(bi_graph_.VerticesV().size());
    traverse_order_u_.resize(bi_graph_.VerticesU().size());
    for (size_t i = 0; i < traverse_order_u_.size(); i++) traverse_order_u_[i] = i;
    traverse_order_u_tmp_.resize(bi_graph_.VerticesU().size());
    traverse_order_nbr_cnt_.resize(bi_graph_.VerticesV().size(), 0);
    traverse_order_v_range_.resize(bi_graph_.VerticesV().size());
    traverse_order_v_range_begin_.resize(bi_graph_.VerticesV().size(), 0);
    cand_tree_prev_node_.assign(bi_graph_.VerticesV().size(), {0, nullptr});
    result_tree_prev_dfs_idx_.assign(bi_graph_.VerticesV().size(), {0, 0});
    cache_inter_neighbors_.resize(bi_graph_.VerticesV().size());
    init_inter_neighbors_.resize(bi_graph_.VerticesV().size(), {(vid_t)-1, {}});

    std::vector<std::vector<size_t>> traverse_order_pre_pos(bi_graph_.VerticesV().size());
    for (auto &v : bi_graph_.VerticesV()) traverse_order_pre_pos[v->vid].reserve(v->neighbors.size());
    std::vector<std::pair<vid_t, NeighborTreeNode *>> traverse_order_last_node(bi_graph_.VerticesU().size(),
                                                                               {-1, nullptr});
    std::vector<size_t> traverse_order_depth_cnt(bi_graph_.VerticesU().size(), 0);

    build_traverse_order(0, traverse_order_u_.size(), 0);

    times.emplace("order", run_watch.Stop());
    run_watch.Report("MBET::Run-order");
    {
      std::cout << "[DBG] cnt=" << bic_cnt_ << " time=" << main_watch_.Stop() << std::endl;

      std::ifstream status("/proc/self/status");
      for (std::string line; std::getline(status, line);) {
        if (line.rfind("VmPeak:", 0) == 0) {
          std::cout << " [MEMORY] VmPeak: " << std::stod(line.substr(7)) / 1024.0 << "MB" << std::endl;
        }
        if (line.rfind("VmSize:", 0) == 0) {
          std::cout << " [MEMORY] VmSize: " << std::stod(line.substr(7)) / 1024.0 << "MB" << std::endl;
        }
        if (line.rfind("VmHWM:", 0) == 0) {
          std::cout << " [MEMORY] VmHWM: " << std::stod(line.substr(6)) / 1024.0 << "MB" << std::endl;
        }
        if (line.rfind("VmRSS:", 0) == 0) {
          std::cout << " [MEMORY] VmRSS: " << std::stod(line.substr(6)) / 1024.0 << "MB" << std::endl;
        }
      }
    }
    run_watch.Start();

    result_tree_root_vid_ = 0;
#ifndef ONLY_COUNTING
    VidList R;
    R.reserve(bi_graph_.VerticesU().front()->neighbors.size());
    std::vector<std::pair<vid_t, CompressedNeighbors>> CLs;
    auto update_compressed_root_neighbors = [&](vid_t root_id) {
      auto watch = DebugTimer<104>::AutoWatch();

      if (!CLs.empty() && CLs.front().first == root_id) {
        return;
      }
      CLs.resize(1);
      CLs.front().first = root_id;
      auto &CL = CLs.front().second;
      CL.clear();
      size_t num = bi_graph_.VertexV(root_id)->neighbors.size();
      size_t cnt = 0;
      while (num >= NumBits) {
        CL.emplace_back(cnt, -1);
        cnt++;
        num -= NumBits;
      }
      if (num > 0) {
        CL.emplace_back(cnt, 0);
        while (num > 0) {
          CL.back().second = (CL.back().second << 1) | 1;
          num--;
        }
      }
    };
#ifndef MBETM
    auto update_compressed_neighbors = [&](size_t pos, vid_t vid) {
      auto watch = DebugTimer<104>::AutoWatch();

      if (CLs.size() > pos && CLs[pos].first == vid) {
        return;
      }
      CLs.resize(pos + 1);
      CLs[pos].first = vid;
      auto &CPL = CLs[pos - 1].second;
      auto &CL = CLs[pos].second;
      CL.clear();
      if (init_inter_neighbors_[vid].first != result_tree_root_vid_) {
        init_inter_neighbors_[vid].first = result_tree_root_vid_;
        init_inter_neighbors_[vid].second.clear();
      }
      if (init_inter_neighbors_[vid].second.empty()) {
        auto auto_watch = DebugTimer<31>::AutoWatch();

        auto &root_nbrs = bi_graph_.VertexV(result_tree_root_vid_)->neighbors;
        auto &nbrs = bi_graph_.VertexV(vid)->neighbors;
        auto nit = nbrs.begin();
        for (auto [idx, cpn] : CPL) {
          size_t tmp = 0;
          while (cpn > 0) {
            auto ctz = __builtin_ctzl(cpn);
            auto nbr_vid = root_nbrs[idx * NumBits + ctz]->vid;
            while (nit != nbrs.end() && (*nit)->vid < nbr_vid) nit++;
            if (nit == nbrs.end()) break;
            auto low_bit = cpn & (-cpn);
            if ((*nit)->vid == nbr_vid) {
              tmp |= low_bit;
            }
            cpn ^= low_bit;
          }
          if (tmp != 0) {
            CL.emplace_back(idx, tmp);
          }
          if (nit == nbrs.end()) break;
        }
      } else {
        auto auto_watch = DebugTimer<32>::AutoWatch();

        auto &cnbrs = init_inter_neighbors_[vid].second;
        auto cnit = cnbrs.begin();
        for (const auto &[idx, cpn] : CPL) {
          while (cnit != cnbrs.end() && cnit->first < idx) cnit++;
          if (cnit == cnbrs.end()) break;
          if (cnit->first == idx) {
            size_t tmp = cnit->second & cpn;
            if (tmp > 0) {
              CL.emplace_back(idx, tmp);
            }
          }
        }
      }
    };
    update_compressed_root_neighbors(result_tree_root_vid_);
#endif
#endif
    for (traverse_order_u_pos_ = 0; traverse_order_u_pos_ < traverse_order_u_.size(); traverse_order_u_pos_++) {
      if (traverse_order_u_tmp_[traverse_order_u_pos_] ==
          static_cast<decltype(traverse_order_u_tmp_)::value_type>(-1)) {
        continue;
      }
      auto &u = bi_graph_.VertexU(traverse_order_u_[traverse_order_u_pos_]);
#ifndef MBETM
      while (result_tree_root_vid_ < u->neighbors.front()->vid) {
        if (result_tree_roots_[result_tree_root_vid_]) {
          result_tree_size_ = 0;
#ifndef ONLY_COUNTING
          {
            auto auto_watch = DebugTimer<102>::AutoWatch();

            result_tree_cnt_++;
            dfs_cnt_ = 0;
            dfs_result_tree(result_tree_roots_[result_tree_root_vid_].get(), result_tree_root_vid_, 0);
          }
          {
            auto auto_watch = DebugTimer<103>::AutoWatch();

            R.clear();
            Enumerate(result_tree_roots_[result_tree_root_vid_].get(), result_tree_root_vid_, CLs.front().second, R);
          }
#endif
          result_tree_roots_[result_tree_root_vid_].reset();
          result_tree_size_max_ = std::max(result_tree_size_max_, result_tree_size_cur_);
          result_tree_size_cur_ -= result_tree_size_;
          result_tree_size_total_ += result_tree_size_;
        }
        result_tree_root_vid_++;
#ifndef ONLY_COUNTING
        if (result_tree_root_vid_ < bi_graph_.VerticesV().size()) {
          update_compressed_root_neighbors(result_tree_root_vid_);
        }
#endif
      }
      if (result_tree_roots_[result_tree_root_vid_]) {
        auto rnode = result_tree_roots_[result_tree_root_vid_].get();
#ifndef ONLY_COUNTING
        R.clear();
        R.emplace_back(result_tree_root_vid_);
#endif
        for (size_t ni = 1; ni < u->neighbors.size(); ni++) {
          auto it = rnode->children.begin();
          while (it != rnode->children.end()) {
            if (it->first >= u->neighbors[ni]->vid) {
              break;
            }
            result_tree_size_ = 0;
#ifndef ONLY_COUNTING
            {
              auto auto_watch = DebugTimer<102>::AutoWatch();

              result_tree_cnt_++;
              dfs_cnt_ = 0;
              dfs_result_tree(it->second.get(), it->first, 0);
            }
            {
              auto auto_watch = DebugTimer<103>::AutoWatch();

              Enumerate(it->second.get(), it->first, CLs[ni - 1].second, R);
            }
#endif
            it = rnode->children.erase(it);
            result_tree_size_max_ = std::max(result_tree_size_max_, result_tree_size_cur_);
            result_tree_size_cur_ -= result_tree_size_;
            result_tree_size_total_ += result_tree_size_;
          }
          if (it == rnode->children.end() || it->first > u->neighbors[ni]->vid) {
            break;
          }
          rnode = it->second.get();
#ifndef ONLY_COUNTING
          update_compressed_neighbors(ni, it->first);
          R.emplace_back(it->first);
#endif
        }
      }
      {
        auto rnode = result_tree_roots_[u->neighbors.front()->vid].get();
        auto is_new = (rnode == nullptr);
        if (!is_new) {
          for (size_t ni = 1; ni < u->neighbors.size(); ni++) {
            if (rnode->children.empty() || rnode->children.begin()->first > u->neighbors[ni]->vid) {
              is_new = true;
              break;
            }
            rnode = rnode->children.begin()->second.get();
          }
          is_new |= !rnode->is_end;
        }
        if (!is_new) {
          for (auto &v : u->neighbors) {
            traverse_order_pre_pos[v->vid].emplace_back(traverse_order_u_pos_);
          }
          continue;
        }
      }
#endif
      auto cand_tree_root = std::make_unique<NeighborTreeNode>(-1, nullptr);
      size_t nbr_tree_size = 1;
      {
        auto auto_watch = DebugTimer<100>::AutoWatch();

        std::vector<vid_t> X;
        for (auto &v : u->neighbors) {
          for (auto &uu_pos : traverse_order_pre_pos[v->vid]) {
            auto uu = traverse_order_u_[uu_pos];
            NeighborTreeNode *parent_node = nullptr;
            if (traverse_order_last_node[uu].first != traverse_order_u_pos_) {
              traverse_order_last_node[uu] = {traverse_order_u_pos_, cand_tree_root.get()};
              traverse_order_depth_cnt[uu] = 0;
            } else {
              parent_node = traverse_order_last_node[uu].second;
              parent_node->sub_cnt--;
            }
            auto [it, res] = traverse_order_last_node[uu].second->children.try_emplace(v->vid, nullptr);
            if (res) {
              it->second = std::make_unique<NeighborTreeNode>(v->vid, parent_node);
              nbr_tree_size++;
            }
            it->second->sub_cnt++;
            traverse_order_last_node[uu].second = it->second.get();
            if (++traverse_order_depth_cnt[uu] == u->neighbors.size()) {
              X.emplace_back(uu);
              if (uu_pos > 0) {
                auto pos = uu_pos - 1;
                while (pos > 0) {
                  if (traverse_order_u_tmp_[pos] != static_cast<decltype(traverse_order_u_tmp_)::value_type>(-1)) {
                    break;
                  }
                  X.emplace_back(traverse_order_u_[pos]);
                  pos--;
                }
                if (traverse_order_u_tmp_[pos] == static_cast<decltype(traverse_order_u_tmp_)::value_type>(-1)) {
                  X.emplace_back(traverse_order_u_[pos]);
                }
              }
            }
          }
          traverse_order_pre_pos[v->vid].emplace_back(traverse_order_u_pos_);
        }
        if (!X.empty()) {
          auto node = traverse_order_last_node[X.front()].second;
          node->sub_cnt = 0;
          while (node->sub_cnt == 0) {
            auto np = node->parent;
            if (!np) {
              cand_tree_root->children.erase(node->vid);
              nbr_tree_size--;
              break;
            }
            auto next_it = np->children.erase(np->children.begin());
            nbr_tree_size--;
            if (next_it != np->children.end()) {
              break;
            }
            node = np;
          }
        }
#ifdef ONLY_COUNTING
        bic_cnt_++;
        if (bic_cnt_ % 1000000 == 0) {
          std::cout << "[DBG] cnt=" << bic_cnt_ << " time=" << main_watch_.Stop() << std::endl;

          std::ifstream status("/proc/self/status");
          for (std::string line; std::getline(status, line);) {
            if (line.rfind("VmPeak:", 0) == 0) {
              std::cout << " [MEMORY] VmPeak: " << std::stod(line.substr(7)) / 1024.0 << "MB" << std::endl;
            }
            if (line.rfind("VmSize:", 0) == 0) {
              std::cout << " [MEMORY] VmSize: " << std::stod(line.substr(7)) / 1024.0 << "MB" << std::endl;
            }
            if (line.rfind("VmHWM:", 0) == 0) {
              std::cout << " [MEMORY] VmHWM: " << std::stod(line.substr(6)) / 1024.0 << "MB" << std::endl;
            }
            if (line.rfind("VmRSS:", 0) == 0) {
              std::cout << " [MEMORY] VmRSS: " << std::stod(line.substr(6)) / 1024.0 << "MB" << std::endl;
            }
          }

          std::cout << " [STAT] duplicated_cnt=" << duplicated_cnt_ << std::endl;
          std::cout << " [STAT] cand_tree_size_cur=" << cand_tree_size_cur_ << std::endl;
          std::cout << " [STAT] cand_tree_size_max=" << cand_tree_size_max_ << std::endl;
          std::cout << " [STAT] cand_tree_size_total=" << cand_tree_size_total_ << std::endl;
          std::cout << " [STAT] result_tree_size_cur=" << result_tree_size_cur_ << std::endl;
          std::cout << " [STAT] result_tree_size_max=" << result_tree_size_max_ << std::endl;
          std::cout << " [STAT] result_tree_size_total=" << result_tree_size_total_ << std::endl;
        }
#else
        X.emplace_back(u->vid);
        VidList Y(u->neighbors.size());
        for (size_t i = 0; i < u->neighbors.size(); i++) Y[i] = u->neighbors[i]->vid;
        if (traverse_order_u_pos_ > 0) {
          auto pos = traverse_order_u_pos_ - 1;
          while (pos > 0) {
            if (traverse_order_u_tmp_[pos] != static_cast<decltype(traverse_order_u_tmp_)::value_type>(-1)) {
              break;
            }
            X.emplace_back(traverse_order_u_[pos]);
            pos--;
          }
          if (traverse_order_u_tmp_[pos] == static_cast<decltype(traverse_order_u_tmp_)::value_type>(-1)) {
            X.emplace_back(traverse_order_u_[pos]);
          }
        }
        Report(X, Y);
#endif
        cand_tree_size_total_ += nbr_tree_size;
        cand_tree_size_cur_ += nbr_tree_size;
        cand_tree_size_max_ = std::max(cand_tree_size_max_, cand_tree_size_cur_);
        for (auto &[id, child] : cand_tree_root->children) {
          build_end(child.get());
        }
      }
      {
        auto auto_watch = DebugTimer<101>::AutoWatch();

        Find(std::move(cand_tree_root));
        cand_tree_size_cur_ -= nbr_tree_size;
      }
#ifdef MBETM
      {
        for (auto &v : u->neighbors) {
          if (!result_tree_roots_[v->vid]) {
            continue;
          }
          result_tree_size_ = 0;
#ifndef ONLY_COUNTING
          result_tree_root_vid_ = v->vid;
          {
            auto auto_watch = DebugTimer<102>::AutoWatch();

            result_tree_cnt_++;
            dfs_cnt_ = 0;
            dfs_result_tree(result_tree_roots_[result_tree_root_vid_].get(), result_tree_root_vid_, 0);
          }
          update_compressed_root_neighbors(result_tree_root_vid_);
          {
            auto auto_watch = DebugTimer<103>::AutoWatch();

            R.clear();
            Enumerate(result_tree_roots_[result_tree_root_vid_].get(), result_tree_root_vid_, CLs.front().second, R);
          }
#endif
          result_tree_roots_[v->vid].reset();
          result_tree_size_max_ = std::max(result_tree_size_max_, result_tree_size_cur_);
          result_tree_size_cur_ -= result_tree_size_;
          result_tree_size_total_ += result_tree_size_;
        }
      }
#endif
    }
#ifndef MBETM
    while (result_tree_root_vid_ < bi_graph_.VerticesV().size()) {
      if (result_tree_roots_[result_tree_root_vid_]) {
        result_tree_size_ = 0;
#ifndef ONLY_COUNTING
        {
          auto auto_watch = DebugTimer<102>::AutoWatch();

          result_tree_cnt_++;
          dfs_cnt_ = 0;
          dfs_result_tree(result_tree_roots_[result_tree_root_vid_].get(), result_tree_root_vid_, 0);
        }
        {
          auto auto_watch = DebugTimer<103>::AutoWatch();

          R.clear();
          Enumerate(result_tree_roots_[result_tree_root_vid_].get(), result_tree_root_vid_, CLs.front().second, R);
        }
#endif
        result_tree_roots_[result_tree_root_vid_].reset();
        result_tree_size_max_ = std::max(result_tree_size_max_, result_tree_size_cur_);
        result_tree_size_cur_ -= result_tree_size_;
        result_tree_size_total_ += result_tree_size_;
      }
      result_tree_root_vid_++;
#ifndef ONLY_COUNTING
      if (result_tree_root_vid_ < bi_graph_.VerticesV().size()) {
        update_compressed_root_neighbors(result_tree_root_vid_);
      }
#endif
    }
#endif

    times.emplace("run", run_watch.Stop());
    run_watch.Report("MBET::Run");

    DebugTimer<100>::Report(" Run-find-build");
    DebugTimer<101>::Report(" Run-find-find");
    DebugTimer<11>::Report("  build prev/dfs_idx");
    DebugTimer<12>::Report("  result iter");
    DebugTimer<13>::Report("  add sub_cnt");
    DebugTimer<14>::Report("  tree fork");
    DebugTimer<15>::Report("  tree prune");
    DebugTimer<16>::Report("  build is_end");
    DebugTimer<17>::Report("  build is_use");
    DebugTimer<21>::Report("  build init_inter_neighbors");
    DebugTimer<102>::Report(" Run-enum-init");
    DebugTimer<103>::Report(" Run-enum-enum");
    DebugTimer<31>::Report("  intersect first-time");
    DebugTimer<32>::Report("  intersect from-cache");
    DebugTimer<33>::Report("  cache push");
    DebugTimer<34>::Report("  cache pop");
    DebugTimer<35>::Report("  extract report");
    DebugTimer<104>::Report(" Run-comp-nbr");
    std::cout << "[STAT] duplicated_cnt=" << duplicated_cnt_ << std::endl;
    std::cout << "[STAT] cand_tree_size_max=" << cand_tree_size_max_ << std::endl;
    std::cout << "[STAT] cand_tree_size_total=" << cand_tree_size_total_ << std::endl;
    std::cout << "[STAT] result_tree_size_max=" << result_tree_size_max_ << std::endl;
    std::cout << "[STAT] result_tree_size_total=" << result_tree_size_total_ << std::endl;

    return {bic_cnt_, main_watch_.Stop(), times};
  }

 private:
  void Report(const VidList &L, const VidList &R) {
    bic_cnt_++;
#ifdef OUTPUT_RESULT
    OutputResult(bi_graph_, L, R, swap_vertices_);
#else
    if (bic_cnt_ % 1000000 == 0) {
      std::cout << "[DBG] cnt=" << bic_cnt_ << " time=" << main_watch_.Stop() << std::endl;

      std::ifstream status("/proc/self/status");
      for (std::string line; std::getline(status, line);) {
        if (line.rfind("VmPeak:", 0) == 0) {
          std::cout << " [MEMORY] VmPeak: " << std::stod(line.substr(7)) / 1024.0 << "MB" << std::endl;
        }
        if (line.rfind("VmSize:", 0) == 0) {
          std::cout << " [MEMORY] VmSize: " << std::stod(line.substr(7)) / 1024.0 << "MB" << std::endl;
        }
        if (line.rfind("VmHWM:", 0) == 0) {
          std::cout << " [MEMORY] VmHWM: " << std::stod(line.substr(6)) / 1024.0 << "MB" << std::endl;
        }
        if (line.rfind("VmRSS:", 0) == 0) {
          std::cout << " [MEMORY] VmRSS: " << std::stod(line.substr(6)) / 1024.0 << "MB" << std::endl;
        }
      }

      std::cout << " [STAT] duplicated_cnt=" << duplicated_cnt_ << std::endl;
      std::cout << " [STAT] cand_tree_size_cur=" << cand_tree_size_cur_ << std::endl;
      std::cout << " [STAT] cand_tree_size_max=" << cand_tree_size_max_ << std::endl;
      std::cout << " [STAT] cand_tree_size_total=" << cand_tree_size_total_ << std::endl;
      std::cout << " [STAT] result_tree_size_cur=" << result_tree_size_cur_ << std::endl;
      std::cout << " [STAT] result_tree_size_max=" << result_tree_size_max_ << std::endl;
      std::cout << " [STAT] result_tree_size_total=" << result_tree_size_total_ << std::endl;
    }
#endif
  }

  void Report(const CompressedNeighbors &CNL, const VidList &R) {
    bic_cnt_++;

    VidList L;
    {
      auto watch = DebugTimer<35>::AutoWatch();

      auto &nbrs = bi_graph_.VertexV(result_tree_root_vid_)->neighbors;
      for (auto [idx, cpn] : CNL) {
        while (cpn > 0) {
          L.emplace_back(nbrs[idx * NumBits + __builtin_ctzl(cpn)]->vid);
          cpn &= cpn - 1;
        }
      }
    }
#ifdef OUTPUT_RESULT
    OutputResult(bi_graph_, L, R, swap_vertices_);
#else
    if (bic_cnt_ % 1000000 == 0) {
      std::cout << "[DBG] cnt=" << bic_cnt_ << " time=" << main_watch_.Stop() << std::endl;

      std::ifstream status("/proc/self/status");
      for (std::string line; std::getline(status, line);) {
        if (line.rfind("VmPeak:", 0) == 0) {
          std::cout << " [MEMORY] VmPeak: " << std::stod(line.substr(7)) / 1024.0 << "MB" << std::endl;
        }
        if (line.rfind("VmSize:", 0) == 0) {
          std::cout << " [MEMORY] VmSize: " << std::stod(line.substr(7)) / 1024.0 << "MB" << std::endl;
        }
        if (line.rfind("VmHWM:", 0) == 0) {
          std::cout << " [MEMORY] VmHWM: " << std::stod(line.substr(6)) / 1024.0 << "MB" << std::endl;
        }
        if (line.rfind("VmRSS:", 0) == 0) {
          std::cout << " [MEMORY] VmRSS: " << std::stod(line.substr(6)) / 1024.0 << "MB" << std::endl;
        }
      }

      std::cout << " [STAT] duplicated_cnt=" << duplicated_cnt_ << std::endl;
      std::cout << " [STAT] cand_tree_size_cur=" << cand_tree_size_cur_ << std::endl;
      std::cout << " [STAT] cand_tree_size_max=" << cand_tree_size_max_ << std::endl;
      std::cout << " [STAT] cand_tree_size_total=" << cand_tree_size_total_ << std::endl;
      std::cout << " [STAT] result_tree_size_cur=" << result_tree_size_cur_ << std::endl;
      std::cout << " [STAT] result_tree_size_max=" << result_tree_size_max_ << std::endl;
      std::cout << " [STAT] result_tree_size_total=" << result_tree_size_total_ << std::endl;
    }
#endif
  }

  void build_traverse_order(size_t range_begin, size_t range_end, size_t nbr_index) {
    std::vector<vid_t> range_nbr;
    size_t out_range_cnt = 0;
    for (size_t i = range_begin; i < range_end; i++) {
      traverse_order_u_tmp_[i] = traverse_order_u_[i];
      auto &nbrs = bi_graph_.VertexU(traverse_order_u_[i])->neighbors;
      if (nbrs.size() <= nbr_index) {
        out_range_cnt++;
        continue;
      }
      auto nbr = nbrs[nbr_index]->vid;
      if (traverse_order_nbr_cnt_[nbr]++ == 0) {
        range_nbr.emplace_back(nbr);
      }
    }
    if (out_range_cnt == 0 && range_nbr.size() == 1) {
#ifdef MBETM
      traverse_order_v_range_[range_nbr.front()].emplace_back(range_begin, range_end - 1);
#endif
      traverse_order_nbr_cnt_[range_nbr.front()] = 0;
      build_traverse_order(range_begin, range_end, nbr_index + 1);
      return;
    }
    if (range_nbr.size() > 1) {
      std::sort(range_nbr.begin(), range_nbr.end());
    }
    size_t sub_range_begin = range_begin;
    for (auto &nbr : range_nbr) {
      auto cnt = traverse_order_nbr_cnt_[nbr];
      traverse_order_nbr_cnt_[nbr] = sub_range_begin;
      sub_range_begin += cnt;
    }
    for (size_t i = range_begin; i < range_end; i++) {
      auto &nbrs = bi_graph_.VertexU(traverse_order_u_tmp_[i])->neighbors;
      if (nbrs.size() <= nbr_index) {
        traverse_order_u_[sub_range_begin++] = traverse_order_u_tmp_[i];
      } else {
        traverse_order_u_[traverse_order_nbr_cnt_[nbrs[nbr_index]->vid]++] = traverse_order_u_tmp_[i];
      }
    }
    for (auto &pivot_nbr : range_nbr) {
      auto nbr = pivot_nbr;
      pivot_nbr = traverse_order_nbr_cnt_[nbr];
      traverse_order_nbr_cnt_[nbr] = 0;
    }
    sub_range_begin = range_begin;
    for (auto sub_range_end : range_nbr) {
      if (sub_range_end - sub_range_begin > 1) {
#ifdef MBETM
        auto &nbrs = bi_graph_.VertexU(traverse_order_u_[sub_range_begin])->neighbors;
        traverse_order_v_range_[nbrs[nbr_index]->vid].emplace_back(sub_range_begin, sub_range_end - 1);
#endif
        build_traverse_order(sub_range_begin, sub_range_end, nbr_index + 1);
      } else {
#ifdef MBETM
        auto &nbrs = bi_graph_.VertexU(traverse_order_u_[sub_range_begin])->neighbors;
        for (size_t i = nbr_index; i < nbrs.size(); i++) {
          traverse_order_v_range_[nbrs[i]->vid].emplace_back(sub_range_begin, sub_range_begin);
        }
#endif
      }
      sub_range_begin = sub_range_end;
    }
    for (; sub_range_begin < range_end - 1; sub_range_begin++) {
      traverse_order_u_tmp_[sub_range_begin] = static_cast<decltype(traverse_order_u_tmp_)::value_type>(-1);
    }
  }

  void build_end(NeighborTreeNode *node) {
    if (node->sub_cnt > 0) {
      node->is_end = true;
      node->sub_cnt = 0;
    }
    for (auto &[id, child] : node->children) {
      build_end(child.get());
    }
  }

  void dfs_pre(NeighborTreeNode *node) {
    dfs_cnt_++;
    node->dfs_idx.first = dfs_cnt_;
    auto vid = node->vid;
    if (cand_tree_prev_node_[vid].first != cand_tree_cnt_) {
      cand_tree_prev_node_[vid].first = cand_tree_cnt_;
      cand_tree_prev_node_[vid].second = nullptr;
    }
    node->prev = cand_tree_prev_node_[vid].second;
    cand_tree_prev_node_[vid].second = node;
    for (auto &[id, child] : node->children) {
      dfs_pre(child.get());
    }
    node->dfs_idx.second = dfs_cnt_;
  }

  bool check_traverse_order_v_range(const std::pair<size_t, size_t> &pre_range, size_t idx) {
    if (idx >= traverse_order_v_range_idx_.size()) {
      return true;
    }
    const auto &ranges = traverse_order_v_range_[traverse_order_v_range_idx_[idx].first];
    for (auto &ri = traverse_order_v_range_idx_[idx].second; ri < ranges.size(); ri++) {
      if (ranges[ri].first >= pre_range.first && ranges[ri].second >= max_traverse_order_v_range_begin_) {
        if (ranges[ri].first > pre_range.second) {
          return false;
        }
        max_traverse_order_v_range_begin_ = std::max(max_traverse_order_v_range_begin_, ranges[ri].first);
        if (check_traverse_order_v_range(ranges[ri], idx + 1)) {
          return true;
        }
      }
    }
    return false;
  }

  void build_use(NeighborTreeNode *node, size_t depth, bool always_use = false) {
    if (always_use) {
      node->is_use = true;
    }
    if (node->is_end) {
      if (always_use) {
        node->is_res = true;
      } else {
        traverse_order_v_range_idx_.resize(depth);
        size_t d = depth;
        for (auto *s_node = node; s_node != nullptr; s_node = s_node->parent) {
          const auto &ranges = traverse_order_v_range_[s_node->vid];
          auto &ri = traverse_order_v_range_begin_[s_node->vid];
          while (ri < ranges.size() && ranges[ri].second <= traverse_order_u_pos_) {
            ri++;
          }
          traverse_order_v_range_idx_[--d] = {s_node->vid, ri};
        }
        max_traverse_order_v_range_begin_ = 0;
        if (!check_traverse_order_v_range({0, traverse_order_u_.size()}, 0)) {
          node->is_res = true;
          for (auto *s_node = node; s_node != nullptr; s_node = s_node->parent) {
            s_node->is_use = true;
          }
        }
      }
    }
    for (auto &[id, child] : node->children) {
      build_use(child.get(), depth + 1, always_use || node->is_res);
    }
  }

  void dfs_reduce(NeighborTreeNode *node) {
    size_t sum_child_sub_cnt = 0;
    for (const auto &[id, child] : node->children) sum_child_sub_cnt += child->sub_cnt;
    if (sum_child_sub_cnt < node->sub_cnt) {
      node->is_end = true;
    } else {
      node->is_end = false;
    }
    node->sub_cnt = 0;
    for (const auto &[id, child] : node->children) dfs_reduce(child.get());
  }

  void dfs_find(ResultTreeNode *r_node, NeighborTreeNode *node, size_t depth) {
#ifdef MBETM
    if (node->is_use)
#endif
    {
      auto watch = DebugTimer<12>::AutoWatch();

      if (r_node == nullptr) {
        if (!result_tree_roots_[node->vid]) {
          result_tree_size_cur_++;
          result_tree_roots_[node->vid] = std::make_unique<ResultTreeNode>();
        }
        r_node = result_tree_roots_[node->vid].get();
      } else {
        auto [it, res] = r_node->children.try_emplace(node->vid, nullptr);
        if (res) {
          result_tree_size_cur_++;
          it->second = std::make_unique<ResultTreeNode>();
        }
        r_node = it->second.get();
      }
    }

    if (node->is_end) {
      {
        auto watch = DebugTimer<13>::AutoWatch();

        for (auto *s_node = node; s_node != nullptr; s_node = s_node->parent) s_node->sub_cnt++;
      }
      if (
#ifdef MBETM
          node->is_res &&
#endif
          !r_node->is_end) {
        r_node->is_end = true;

#ifdef ONLY_COUNTING
        bic_cnt_++;
        if (bic_cnt_ % 1000000 == 0) {
          std::cout << "[DBG] cnt=" << bic_cnt_ << " time=" << main_watch_.Stop() << std::endl;

          std::ifstream status("/proc/self/status");
          for (std::string line; std::getline(status, line);) {
            if (line.rfind("VmPeak:", 0) == 0) {
              std::cout << " [MEMORY] VmPeak: " << std::stod(line.substr(7)) / 1024.0 << "MB" << std::endl;
            }
            if (line.rfind("VmSize:", 0) == 0) {
              std::cout << " [MEMORY] VmSize: " << std::stod(line.substr(7)) / 1024.0 << "MB" << std::endl;
            }
            if (line.rfind("VmHWM:", 0) == 0) {
              std::cout << " [MEMORY] VmHWM: " << std::stod(line.substr(6)) / 1024.0 << "MB" << std::endl;
            }
            if (line.rfind("VmRSS:", 0) == 0) {
              std::cout << " [MEMORY] VmRSS: " << std::stod(line.substr(6)) / 1024.0 << "MB" << std::endl;
            }
          }

          std::cout << " [STAT] duplicated_cnt=" << duplicated_cnt_ << std::endl;
          std::cout << " [STAT] cand_tree_size_cur=" << cand_tree_size_cur_ << std::endl;
          std::cout << " [STAT] cand_tree_size_max=" << cand_tree_size_max_ << std::endl;
          std::cout << " [STAT] cand_tree_size_total=" << cand_tree_size_total_ << std::endl;
          std::cout << " [STAT] result_tree_size_cur=" << result_tree_size_cur_ << std::endl;
          std::cout << " [STAT] result_tree_size_max=" << result_tree_size_max_ << std::endl;
          std::cout << " [STAT] result_tree_size_total=" << result_tree_size_total_ << std::endl;
        }
#endif

        size_t nbr_tree_size = 1;
        auto nbr_tree_root = std::make_unique<NeighborTreeNode>(-1, nullptr);

        {
          std::vector<NeighborTreeNode *> stack(depth, nullptr);
          {
            size_t d = depth;
            for (auto *s_node = node; s_node != nullptr; s_node = s_node->parent) {
              stack[--d] = s_node;
            }
          }
          {
            auto watch = DebugTimer<14>::AutoWatch();

            std::vector<NeighborTreeNode *> s_node_vec(depth, nullptr);
            for (size_t si = 0; si < depth; si++) {
              for (size_t sj = 0; sj < si; sj++) {
                s_node_vec[sj] = stack[sj];
              }
              auto *s_node = stack[si];
              auto vid = s_node->vid;
              {
                nbr_tree_size++;
                auto child = std::make_unique<NeighborTreeNode>(vid, nullptr);
                s_node->fork = child.get();
                child->sub_cnt = s_node->sub_cnt;
                if (si == 0) {
                  nbr_tree_root->children.emplace(vid, std::move(child));
                } else {
                  child->parent = stack[si - 1]->fork;
                  stack[si - 1]->fork->children.emplace(vid, std::move(child));
                }
              }
              for (s_node = s_node->prev; s_node != nullptr; s_node = s_node->prev) {
                bool flag = false;
                for (size_t tsj = si; tsj > 0; tsj--) {
                  size_t sj = tsj - 1;
                  while (s_node_vec[sj] != nullptr && s_node_vec[sj]->dfs_idx.first > s_node->dfs_idx.first)
                    s_node_vec[sj] = s_node_vec[sj]->prev;
                  if (s_node_vec[sj] == nullptr) continue;
                  if (s_node_vec[sj]->dfs_idx.second >= s_node->dfs_idx.first) {
                    flag = true;
                    auto *f_node = s_node_vec[sj]->fork;
                    auto [it, res] = f_node->children.try_emplace(vid, nullptr);
                    if (res) {
                      nbr_tree_size++;
                      it->second = std::make_unique<NeighborTreeNode>(vid, f_node);
                    }
                    it->second->sub_cnt += s_node->sub_cnt;
                    s_node->fork = it->second.get();
                    break;
                  }
                }
                if (!flag) {
                  auto [it, res] = nbr_tree_root->children.try_emplace(vid, nullptr);
                  if (res) {
                    nbr_tree_size++;
                    it->second = std::make_unique<NeighborTreeNode>(vid, nullptr);
                  }
                  it->second->sub_cnt += s_node->sub_cnt;
                  s_node->fork = it->second.get();
                }
              }
            }
            s_node_vec.clear();
            s_node_vec.shrink_to_fit();
          }

          {
            auto watch = DebugTimer<15>::AutoWatch();

            size_t same_cnt = stack.back()->fork->sub_cnt;
            for (const auto &s_node : stack) {
              s_node->fork->sub_cnt -= same_cnt;
            }
            for (size_t tsi = depth; tsi > 0; tsi--) {
              size_t si = tsi - 1;
              if (tsi == depth) {
                auto *n_node = stack[si]->fork;
                auto p_node = (si == 0) ? nbr_tree_root.get() : stack[si - 1]->fork;
                p_node->children.erase(n_node->vid);
                continue;
              }
              auto *n_node = stack[si]->fork;
              if (n_node->children.empty() && n_node->sub_cnt == 0) {
                auto p_node = (si == 0) ? nbr_tree_root.get() : stack[si - 1]->fork;
                p_node->children.erase(n_node->vid);
              }
            }
          }
          stack.clear();
          stack.shrink_to_fit();
        }

        cand_tree_size_total_ += nbr_tree_size;
        cand_tree_size_cur_ += nbr_tree_size;
        cand_tree_size_max_ = std::max(cand_tree_size_max_, cand_tree_size_cur_);

        if (!nbr_tree_root->children.empty()) {
          {
            auto watch = DebugTimer<16>::AutoWatch();

            dfs_reduce(nbr_tree_root.get());
          }
          Find(std::move(nbr_tree_root));
        }

        nbr_tree_root.reset();
        cand_tree_size_cur_ -= nbr_tree_size;
      } else {
        duplicated_cnt_++;
      }
    }

    for (auto &[id, child] : node->children) {
      dfs_find(r_node, child.get(), depth + 1);
    }
  }

  void Find(std::unique_ptr<NeighborTreeNode> &&nbr_tree_root) {
    {
      auto watch = DebugTimer<11>::AutoWatch();

      cand_tree_cnt_++;
      dfs_cnt_ = 0;
      for (auto &[id, child] : nbr_tree_root->children) {
        dfs_pre(child.get());
      }
    }
#ifdef MBETM
    {
      auto watch = DebugTimer<17>::AutoWatch();

      for (auto &[id, child] : nbr_tree_root->children) {
        build_use(child.get(), 1, false);
      }
    }
#endif
    for (auto &[id, child] : nbr_tree_root->children) {
      dfs_find(nullptr, child.get(), 1);
    }
  }

  void dfs_result_tree(ResultTreeNode *node, vid_t vid, size_t parent_dfs_idx) {
    result_tree_size_++;
    if (init_inter_neighbors_[vid].first != result_tree_root_vid_) {
      auto watch = DebugTimer<21>::AutoWatch();

      init_inter_neighbors_[vid].first = result_tree_root_vid_;
      init_inter_neighbors_[vid].second.clear();
      auto &root_nbrs = bi_graph_.VertexV(result_tree_root_vid_)->neighbors;
      auto &nbrs = bi_graph_.VertexV(vid)->neighbors;
      if (root_nbrs.size() / NumBits + 1 < nbrs.size()) {
        auto nit = nbrs.begin();
        size_t idx = 0;
        size_t tmp = 0;
        for (size_t i = 0, j = 1; i < root_nbrs.size(); i++) {
          while (nit != nbrs.end() && (*nit)->vid < root_nbrs[i]->vid) nit++;
          if (nit == nbrs.end()) break;
          if ((*nit)->vid == root_nbrs[i]->vid) {
            tmp |= j;
          }
          if (j == HighestBit) {
            if (tmp != 0) {
              init_inter_neighbors_[vid].second.emplace_back(idx, tmp);
              tmp = 0;
            }
            idx++;
            j = 1;
          } else {
            j <<= 1;
          }
        }
        if (tmp != 0) {
          init_inter_neighbors_[vid].second.emplace_back(idx, tmp);
        }
      }
    }
    auto dfs_idx = ++dfs_cnt_;
    if (result_tree_prev_dfs_idx_[vid].first != result_tree_cnt_) {
      result_tree_prev_dfs_idx_[vid] = {result_tree_cnt_, 0};
    } else if (result_tree_prev_dfs_idx_[vid].second > parent_dfs_idx) {
      node->has_prev = true;
    }
    result_tree_prev_dfs_idx_[vid].second = dfs_idx;
    for (auto &[id, child] : node->children) {
      dfs_result_tree(child.get(), id, dfs_idx);
    }
  }

  bool Enumerate(ResultTreeNode *node, vid_t vid, const CompressedNeighbors &PL, VidList &R) {
    R.emplace_back(vid);
    CompressedNeighbors intersected_neighbors;
    bool parent_need_pop_cache = false;
    if (vid == result_tree_root_vid_) {
      if (node->is_end) {
        Report(PL, R);
      }
    } else {
      if (cache_inter_neighbors_[vid].empty() && init_inter_neighbors_[vid].second.empty()) {
        auto watch = DebugTimer<31>::AutoWatch();

        auto &root_nbrs = bi_graph_.VertexV(result_tree_root_vid_)->neighbors;
        auto &nbrs = bi_graph_.VertexV(vid)->neighbors;
        auto nit = nbrs.begin();
        for (auto [idx, cpn] : PL) {
          size_t tmp = 0;
          while (cpn > 0) {
            auto ctz = __builtin_ctzl(cpn);
            auto nbr_vid = root_nbrs[idx * NumBits + ctz]->vid;
            while (nit != nbrs.end() && (*nit)->vid < nbr_vid) nit++;
            if (nit == nbrs.end()) break;
            auto low_bit = cpn & (-cpn);
            if ((*nit)->vid == nbr_vid) {
              tmp |= low_bit;
            }
            cpn ^= low_bit;
          }
          if (tmp != 0) {
            intersected_neighbors.emplace_back(idx, tmp);
          }
          if (nit == nbrs.end()) break;
        }
      } else {
        auto watch = DebugTimer<32>::AutoWatch();

        auto &cnbrs = cache_inter_neighbors_[vid].empty() ? init_inter_neighbors_[vid].second
                                                          : cache_inter_neighbors_[vid].back();
        auto cnit = cnbrs.begin();
        for (const auto &[idx, cpn] : PL) {
          while (cnit != cnbrs.end() && cnit->first < idx) cnit++;
          if (cnit == cnbrs.end()) break;
          if (cnit->first == idx) {
            size_t tmp = cnit->second & cpn;
            if (tmp > 0) {
              intersected_neighbors.emplace_back(idx, tmp);
            }
          }
        }
      }

      if (node->is_end) {
        Report(intersected_neighbors, R);
      }

      if (node->has_prev) {
        auto watch = DebugTimer<33>::AutoWatch();

        parent_need_pop_cache = true;
        cache_inter_neighbors_[vid].emplace_back(std::move(intersected_neighbors));
      }
    }
    bool need_pop_cache = false;
    for (auto rit = node->children.rbegin(); rit != node->children.rend(); rit++) {
      if (Enumerate(rit->second.get(), rit->first,
                    (vid == result_tree_root_vid_)
                        ? PL
                        : ((node->has_prev) ? cache_inter_neighbors_[vid].back() : intersected_neighbors),
                    R)) {
        need_pop_cache = true;
      }
    }
    if (need_pop_cache) {
      auto watch = DebugTimer<34>::AutoWatch();

      for (auto &[id, child] : node->children) {
        if (child->has_prev) {
          cache_inter_neighbors_[id].pop_back();
        }
      }
    }
    R.pop_back();
    return parent_need_pop_cache;
  }

  BiGraph bi_graph_;
  bool swap_vertices_;
  SteadyClockWatch main_watch_;

  size_t bic_cnt_;
  std::vector<std::unique_ptr<ResultTreeNode>> result_tree_roots_;

  std::vector<vid_t> traverse_order_u_;
  std::vector<vid_t> traverse_order_u_tmp_;
  std::vector<vid_t> traverse_order_nbr_cnt_;
  std::vector<std::vector<std::pair<size_t, size_t>>> traverse_order_v_range_;
  std::vector<size_t> traverse_order_v_range_begin_;
  std::vector<std::pair<vid_t, size_t>> traverse_order_v_range_idx_;
  size_t max_traverse_order_v_range_begin_;
  size_t traverse_order_u_pos_;

  size_t cand_tree_size_cur_;
  size_t cand_tree_size_max_;
  size_t cand_tree_size_total_;
  size_t result_tree_size_;
  size_t result_tree_size_cur_;
  size_t result_tree_size_max_;
  size_t result_tree_size_total_;

  size_t cand_tree_cnt_;
  std::vector<std::pair<size_t, NeighborTreeNode *>> cand_tree_prev_node_;
  std::vector<std::pair<size_t, size_t>> result_tree_prev_dfs_idx_;
  size_t result_tree_root_vid_;
  size_t result_tree_cnt_;
  std::vector<std::list<CompressedNeighbors>> cache_inter_neighbors_;
  std::vector<std::pair<vid_t, CompressedNeighbors>> init_inter_neighbors_;
  size_t dfs_cnt_;
  size_t duplicated_cnt_;
};