// Exact event-driven Markovian SIR on a static network.
//
// The dust tau-leaping model (stoch_mod_adj.R) evaluates every node at every
// dt step -- n * max_degree work per step -- however few events occur. For a
// sparse outbreak that is ~10^5 operations per realised infection. This
// simulates the SAME continuous-time process exactly, as first-passage
// percolation over the contact graph:
//
//   * a node confirmed infected at time t draws an infectious period
//     Exp(gamma), and for each not-yet-infected neighbour a transmission
//     delay Exp(rate); the transmission happens only if it lands inside the
//     infectious period;
//   * candidates go in a min-heap keyed on time, and the FIRST pop for a node
//     is its true infection time (later pops are stale and discarded).
//
// Rates match stoch_mod_adj.R exactly:
//     foi[i]  = beta / N_total * sum_{k in nbr(i)} I[k] * transmisibility[k]
//     hazard  = susceptibility[i] * foi[i]
//   so the per-edge rate k -> i is  beta / n * trans[k] * sus[i].
// `beta` here is the value passed to run_stoch_adj (run_stoch_network already
// scales the user's beta by N / k_mean).
//
// Cost is O(E log V) with no dt discretisation bias.
// Validity: time-homogeneous rates (scalar beta, constant gamma) and waning=0.
//
// Graph is CSR (nbr/ptr, 0-based), so work scales with EDGES, not
// n * max_degree.

#include <cpp11.hpp>
#include <vector>
#include <queue>
#include <random>
#include <cmath>
#include <limits>
#include <cstdint>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace cpp11;

namespace {

struct Ev {
  double t;
  int node;
  bool operator>(const Ev& o) const { return t > o.t; }
};

// One realisation. inf[i] = infection time, or +Inf if never infected.
void one_sim(const int n,
             const std::vector<int>& nbr, const std::vector<int>& ptr,
             const std::vector<double>& sus, const std::vector<double>& trans,
             const double beta, const double gamma, const double t_max,
             const std::vector<int>& seeds,
             std::mt19937_64& rng,
             std::vector<double>& inf, std::vector<char>& done) {
  const double INF = std::numeric_limits<double>::infinity();
  std::fill(inf.begin(), inf.end(), INF);
  std::fill(done.begin(), done.end(), 0);
  std::exponential_distribution<double> exp1(1.0);
  std::priority_queue<Ev, std::vector<Ev>, std::greater<Ev>> pq;

  const double base = beta / static_cast<double>(n);   // matches foi's 1/N_total

  for (int s : seeds) if (s >= 0 && s < n) pq.push(Ev{0.0, s});

  while (!pq.empty()) {
    const Ev e = pq.top(); pq.pop();
    if (e.t > t_max) break;            // heap is time-ordered: nothing left in range
    const int i = e.node;
    if (done[i]) continue;             // stale candidate (infected earlier)
    done[i] = 1;
    inf[i]  = e.t;                     // first pop = true infection time

    // Infectious period, drawn once the infection is confirmed.
    const double rec_i = e.t + ((gamma > 0) ? exp1(rng) / gamma : INF);

    for (int k = ptr[i]; k < ptr[i + 1]; ++k) {
      const int j = nbr[k];
      if (done[j]) continue;
      const double rate = base * trans[i] * sus[j];
      if (rate <= 0) continue;
      const double t_ev = e.t + exp1(rng) / rate;      // first transmission i->j
      if (t_ev < rec_i && t_ev <= t_max) pq.push(Ev{t_ev, j});
    }
  }
}

}  // namespace

// n_sim x n matrix of infection times (Inf = never infected).
// nbr/ptr 0-based CSR; seeds 0-based node indices.
[[cpp11::register]]
doubles_matrix<> net_sir_event_times(int n,
                                     integers nbr, integers ptr,
                                     doubles susceptibility,
                                     doubles transmissibility,
                                     double beta, double gamma, double t_max,
                                     integers seeds, int n_sim, int seed,
                                     int n_threads) {
  const std::vector<int>    v_nbr(nbr.begin(), nbr.end());
  const std::vector<int>    v_ptr(ptr.begin(), ptr.end());
  const std::vector<double> v_sus(susceptibility.begin(), susceptibility.end());
  const std::vector<double> v_tr(transmissibility.begin(), transmissibility.end());
  const std::vector<int>    v_seed(seeds.begin(), seeds.end());

  // Simulate into a plain buffer so the loop can be threaded: each realisation
  // is independent, with its own RNG stream keyed on its index (so results do
  // NOT depend on the thread count or scheduling). R API calls are not
  // thread-safe, so the writable matrix is filled serially afterwards.
  std::vector<double> buf(static_cast<size_t>(n_sim) * n);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(n_threads)
#endif
  for (int s = 0; s < n_sim; ++s) {
    std::vector<double> inf(n);
    std::vector<char>   done(n);
    std::mt19937_64 rng(static_cast<uint64_t>(seed) * 1000003ULL +
                        static_cast<uint64_t>(s) * 7919ULL + 1ULL);
    one_sim(n, v_nbr, v_ptr, v_sus, v_tr, beta, gamma, t_max, v_seed, rng,
            inf, done);
    std::copy(inf.begin(), inf.end(), buf.begin() + static_cast<size_t>(s) * n);
  }

  writable::doubles_matrix<> out(n_sim, n);
  for (int s = 0; s < n_sim; ++s)
    for (int i = 0; i < n; ++i) out(s, i) = buf[static_cast<size_t>(s) * n + i];
  return out;
}
