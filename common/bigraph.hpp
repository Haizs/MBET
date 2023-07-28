#pragma once

#include <algorithm>
#include <cassert>
#include <fstream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

using vid_t = size_t;
using label_t = vid_t;

class Vertex {
 public:
  explicit Vertex(vid_t vid) : vid(vid), label(vid) {}
  Vertex(vid_t vid, label_t label) : vid(vid), label(label) {}

  bool operator<(const Vertex& other) const { return vid < other.vid; }

  bool is_neighbor(const vid_t& v_id) const {
    return std::find_if(neighbors.begin(), neighbors.end(),
                        [&](const std::shared_ptr<Vertex>& v) { return v->vid == v_id; }) != neighbors.end();
  }

  vid_t vid;
  label_t label;
  std::vector<std::shared_ptr<Vertex>> neighbors;
};

class BiGraph {
 public:
  BiGraph() = default;

  BiGraph(const BiGraph& other) {
    U_.reserve(other.U_.size());
    for (auto& u : other.U_) {
      U_.emplace_back(std::make_shared<Vertex>(u->vid, u->label));
    }
    V_.reserve(other.V_.size());
    for (auto& v : other.V_) {
      V_.emplace_back(std::make_shared<Vertex>(v->vid, v->label));
    }
    for (auto& u : other.U_) {
      auto& uu = U_[u->vid];
      uu->neighbors.reserve(u->neighbors.size());
      for (auto& nbr : u->neighbors) {
        uu->neighbors.emplace_back(V_[nbr->vid]);
      }
    }
    for (auto& v : other.V_) {
      auto& vv = V_[v->vid];
      vv->neighbors.reserve(v->neighbors.size());
      for (auto& nbr : v->neighbors) {
        vv->neighbors.emplace_back(U_[nbr->vid]);
      }
    }
  }

  BiGraph(BiGraph&& other) {
    std::swap(U_, other.U_);
    std::swap(V_, other.V_);
    other.U_.clear();
    other.V_.clear();
  }

  ~BiGraph() {
    for (auto& vtx : U_) vtx->neighbors.clear();
    for (auto& vtx : V_) vtx->neighbors.clear();
  }

  auto& VerticesU() { return U_; }
  auto& VerticesU() const { return U_; }

  auto& VerticesV() { return V_; }
  auto& VerticesV() const { return V_; }

  auto& VertexU(vid_t id) { return U_[id]; }
  auto& VertexU(vid_t id) const { return U_.at(id); }

  auto& VertexV(vid_t id) { return V_[id]; }
  auto& VertexV(vid_t id) const { return V_.at(id); }

  static BiGraph LoadFromFile(const std::string& path) {
    std::ifstream input(path);
    assert(input.is_open());

    BiGraph bi_graph;
    std::unordered_map<vid_t, std::shared_ptr<Vertex>> map_u;
    std::unordered_map<vid_t, std::shared_ptr<Vertex>> map_v;

    std::string line;
    while (!input.eof()) {
      std::getline(input, line);
      vid_t x = 0;
      vid_t y = 0;
      auto it = line.begin();
      if (*it == '%' || *it == '#' || it == line.end()) {
        continue;
      }
      while (it != line.end() && (*it < '0' || *it > '9')) it++;
      while (it != line.end() && (*it >= '0' && *it <= '9')) x = x * 10 + (*it++) - '0';
      while (it != line.end() && (*it < '0' || *it > '9')) it++;
      while (it != line.end() && (*it >= '0' && *it <= '9')) y = y * 10 + (*it++) - '0';

      auto& vx = map_u.try_emplace(x, std::make_shared<Vertex>(x)).first->second;
      vx->vid = x;
      auto& vy = map_v.try_emplace(y, std::make_shared<Vertex>(y)).first->second;
      vy->vid = y;
      vx->neighbors.emplace_back(vy);
      vy->neighbors.emplace_back(vx);
    }
    bi_graph.U_.reserve(map_u.size());
    for (auto& [k, v] : map_u) {
      std::sort(v->neighbors.begin(), v->neighbors.end(), [](const auto& a, const auto& b) { return *a < *b; });
      bi_graph.U_.emplace_back(std::move(v));
    }
    std::sort(bi_graph.U_.begin(), bi_graph.U_.end(), [](const auto& a, const auto& b) { return *a < *b; });
    bi_graph.V_.reserve(map_v.size());
    for (auto& [k, v] : map_v) {
      std::sort(v->neighbors.begin(), v->neighbors.end(), [](const auto& a, const auto& b) { return *a < *b; });
      bi_graph.V_.emplace_back(std::move(v));
    }
    std::sort(bi_graph.V_.begin(), bi_graph.V_.end(), [](const auto& a, const auto& b) { return *a < *b; });
    if (bi_graph.U_.size() > bi_graph.V_.size()) {
      std::swap(bi_graph.U_, bi_graph.V_);
    }

    // make U_ and V_ are contiguous from 0
    for (vid_t i = 0; i < bi_graph.U_.size(); i++) bi_graph.U_[i]->vid = i;
    for (vid_t i = 0; i < bi_graph.V_.size(); i++) bi_graph.V_[i]->vid = i;

    // report statistics
    size_t sum_deg_U = 0;
    size_t sum_deg_V = 0;
    for (auto& v : bi_graph.U_) {
      sum_deg_U += v->neighbors.size();
    }
    for (auto& v : bi_graph.V_) {
      sum_deg_V += v->neighbors.size();
    }
    assert(sum_deg_U == sum_deg_V);
    std::cout << "Loaded BiGraph from: " << path << std::endl;
    std::cout << "|U|: " << bi_graph.U_.size() << " , avg_deg_U: " << (double)sum_deg_U / bi_graph.U_.size()
              << std::endl
              << "|V|: " << bi_graph.V_.size() << " , avg_deg_V: " << (double)sum_deg_V / bi_graph.V_.size()
              << std::endl
              << "|E|: " << sum_deg_U << std::endl;

    return bi_graph;
  }

 private:
  std::vector<std::shared_ptr<Vertex>> U_;
  std::vector<std::shared_ptr<Vertex>> V_;
};

using VidList = std::vector<vid_t>;
std::ostream& operator<<(std::ostream& output, const VidList& vec) {
  output << "[";
  auto it = vec.begin();
  while (it != vec.end()) {
    output << *it;
    it++;
    if (it != vec.end()) std::cout << ", ";
  }
  output << "]";
  return output;
}

void OutputResult(const BiGraph& bi_graph, const VidList& L, const VidList& R, bool swap = false) {
  std::vector<label_t> label_L;
  label_L.reserve(L.size());
  for (auto& v : L) {
    label_L.emplace_back(bi_graph.VertexU(v)->label);
  }
  std::sort(label_L.begin(), label_L.end());
  std::vector<label_t> label_R;
  label_R.reserve(R.size());
  for (auto& v : R) {
    label_R.emplace_back(bi_graph.VertexV(v)->label);
  }
  std::sort(label_R.begin(), label_R.end());
  if (!swap) {
    std::cout << label_L << label_R << std::endl;
  } else {
    std::cout << label_R << label_L << std::endl;
  }
}