#include "controller/MPC_diffDrive_fblin.h"

#include <cmath>
#include <fstream>
#include <unsupported/Eigen/MatrixFunctions>

static inline double sqr(double v) { return v*v; }

static inline double dist2(double x1, double y1, double x2, double y2) {
    return sqr(x1-x2) + sqr(y1-y2);
}

static inline double dist_point_segment2(double px, double py, const Segment& s) {
    double vx = s.x2 - s.x1;
    double vy = s.y2 - s.y1;
    double wx = px - s.x1;
    double wy = py - s.y1;
    double dv2 = vx*vx + vy*vy;
    if (dv2 <= std::numeric_limits<double>::epsilon()) {
        return dist2(px, py, s.x1, s.y1);
    }
    double t = (wx*vx + wy*vy) / dv2;
    if (t < 0.0) t = 0.0;
    else if (t > 1.0) t = 1.0;
    double projx = s.x1 + t * vx;
    double projy = s.y1 + t * vy;
    return dist2(px, py, projx, projy);
}

static inline double dist_point_segment(double px, double py, const Segment& s) {
    return std::sqrt(dist_point_segment2(px, py, s));
}

static inline bool segments_intersect(const Segment& a, const Segment& b) {
    auto orient = [](double ax, double ay, double bx, double by, double cx, double cy) {
        return (bx-ax)*(cy-ay) - (by-ay)*(cx-ax);
    };
    auto on_seg = [](double ax, double ay, double bx, double by, double px, double py) {
        return std::min(ax,bx) <= px + 1e-12 && px <= std::max(ax,bx) + 1e-12 &&
               std::min(ay,by) <= py + 1e-12 && py <= std::max(ay,by) + 1e-12;
    };

    double o1 = orient(a.x1,a.y1, a.x2,a.y2, b.x1,b.y1);
    double o2 = orient(a.x1,a.y1, a.x2,a.y2, b.x2,b.y2);
    double o3 = orient(b.x1,b.y1, b.x2,b.y2, a.x1,a.y1);
    double o4 = orient(b.x1,b.y1, b.x2,b.y2, a.x2,a.y2);

    if ( (o1>0 && o2<0 || o1<0 && o2>0) && (o3>0 && o4<0 || o3<0 && o4>0) ) return true;

    if (std::abs(o1) <= 1e-12 && on_seg(a.x1,a.y1,a.x2,a.y2, b.x1,b.y1)) return true;
    if (std::abs(o2) <= 1e-12 && on_seg(a.x1,a.y1,a.x2,a.y2, b.x2,b.y2)) return true;
    if (std::abs(o3) <= 1e-12 && on_seg(b.x1,b.y1,b.x2,b.y2, a.x1,a.y1)) return true;
    if (std::abs(o4) <= 1e-12 && on_seg(b.x1,b.y1,b.x2,b.y2, a.x2,a.y2)) return true;

    return false;
}

static inline double dist_segment_segment(const Segment& a, const Segment& b) {
    if (segments_intersect(a,b)) return 0.0;
    double d1 = dist_point_segment2(a.x1,a.y1, b);
    double d2 = dist_point_segment2(a.x2,a.y2, b);
    double d3 = dist_point_segment2(b.x1,b.y1, a);
    double d4 = dist_point_segment2(b.x2,b.y2, a);
    return std::sqrt(std::min(std::min(d1,d2), std::min(d3,d4)));
}

std::vector<std::vector<int>> cluster_circles_and_segments_long_separate(
    const std::vector<Circle>& circles,
    const std::vector<Segment>& segments,
    double tau,
    int max_cluster_size,
    double long_segment_threshold_m
) {
    const int n = static_cast<int>(circles.size());
    const int m = static_cast<int>(segments.size());

    // 1) separa segmenti "lunghi" e "corti" mantenendo gli indici originali
    std::vector<int> short_seg_idx;
    std::vector<int> long_seg_idx;
    short_seg_idx.reserve(m);
    long_seg_idx.reserve(m);

    for (int sj = 0; sj < m; ++sj) {
        const auto& s = segments[sj];
        double len = std::hypot(s.x2 - s.x1, s.y2 - s.y1);
        if (len > long_segment_threshold_m) long_seg_idx.push_back(sj);
        else short_seg_idx.push_back(sj);
    }

    // 2) DSU solo su (cerchi + segmenti corti)
    //    mapping: 0..n-1 cerchi
    //             n..n+ms-1 segmenti corti (ms = short_seg_idx.size())
    const int ms = static_cast<int>(short_seg_idx.size());
    const int total = n + ms;

    DSUCap dsu(total);

    auto global_idx_of_short_seg = [&](int local_short_k) -> int {
        // indice globale nel vettore "combined DSU"
        return n + local_short_k;
    };

    // 2.1) cerchio-cerchio (identico)
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double d = std::hypot(circles[i].x - circles[j].x, circles[i].y - circles[j].y);
            if (d <= circles[i].r + circles[j].r + tau) {
                dsu.unite(i, j, max_cluster_size);
            }
        }
    }

    // 2.2) cerchio-segmento (solo segmenti corti, identico criterio)
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < ms; ++k) {
            int sj = short_seg_idx[k];                 // indice reale in segments[]
            int j  = global_idx_of_short_seg(k);       // indice DSU
            double d_pt_seg = dist_point_segment(circles[i].x, circles[i].y, segments[sj]);
            if (d_pt_seg <= circles[i].r + tau) {
                dsu.unite(i, j, max_cluster_size);
            }
        }
    }

    // 2.3) segmento-segmento (solo segmenti corti, identico criterio)
    for (int a = 0; a < ms; ++a) {
        for (int b = a + 1; b < ms; ++b) {
            int sa = short_seg_idx[a];
            int sb = short_seg_idx[b];
            int i = global_idx_of_short_seg(a);
            int j = global_idx_of_short_seg(b);

            double dseg = dist_segment_segment(segments[sa], segments[sb]);
            if (dseg <= tau) {
                dsu.unite(i, j, max_cluster_size);
            }
        }
    }

    // 3) raggruppa per radice e converti indici DSU -> indici originali (cerchi/segmenti)
    std::unordered_map<int, std::vector<int>> by_root;
    by_root.reserve(total);

    for (int idx = 0; idx < total; ++idx) {
        by_root[dsu.find(idx)].push_back(idx);
    }

    // 4) prepara output: prima i cluster DSU (convertiti in indici "globali originali")
    //    Convenzione output (come prima):
    //      0..n-1 = circles
    //      n..n+m-1 = segments (indice globale = n + sj)
    std::vector<std::vector<int>> clusters;
    clusters.reserve(by_root.size() + long_seg_idx.size());

    // ordine deterministico: ordina le radici
    std::vector<int> roots;
    roots.reserve(by_root.size());
    for (auto& kv : by_root) roots.push_back(kv.first);
    std::sort(roots.begin(), roots.end());

    for (int root : roots) {
        std::vector<int> out;

        auto& members = by_root[root];
        out.reserve(members.size());

        for (int member : members) {
            if (member < n) {
                // cerchio: stesso indice
                out.push_back(member);
            } else {
                // segmento corto: converti DSU-index -> indice reale segmento sj -> output index = n + sj
                int k = member - n;
                int sj = short_seg_idx[k];
                out.push_back(n + sj);
            }
        }

        std::sort(out.begin(), out.end());
        clusters.push_back(std::move(out));
    }

    // 5) aggiungi i segmenti lunghi come cluster singoli (separati)
    //    Li metto in ordine crescente di indice segmento per determinismo.
    std::sort(long_seg_idx.begin(), long_seg_idx.end());
    for (int sj : long_seg_idx) {
        clusters.push_back(std::vector<int>{ n + sj });
    }

    // 6) opzionale: ordinamento globale dei cluster per avere output totalmente deterministico
    //    (prima ordina gli elementi interni, già fatto; qui ordino i cluster per il loro primo elemento)
    std::sort(clusters.begin(), clusters.end(),
              [](const std::vector<int>& a, const std::vector<int>& b) {
                  if (a.empty() || b.empty()) return a.size() < b.size();
                  return a[0] < b[0];
              });

    return clusters;
}


static double cross2(const point& O, const point& A, const point& B)
{
    return (A.x - O.x)*(B.y - O.y) - (A.y - O.y)*(B.x - O.x);
}

static double cross2(const P& O, const P& A, const P& B)
{
    return (A.x - O.x)*(B.y - O.y) - (A.y - O.y)*(B.x - O.x);
}

std::vector<point> convexHull(std::vector<point> P)
{
    if (P.size() < 3) return P;

    std::sort(P.begin(), P.end(), [](const point& a, const point& b){
      return a.x < b.x || (a.x == b.x && a.y < b.y);
    });

    // remove duplicates
    P.erase(std::unique(P.begin(), P.end(), [](const point& a, const point& b){
      return a.x == b.x && a.y == b.y;
    }), P.end());

    if (P.size() < 3) return P;

    std::vector<point> H;
    H.reserve(P.size()*2);

    // lower hull
    for (const auto& p : P)
    {
        while (H.size() >= 2 && cross2(H[H.size()-2], H[H.size()-1], p) < 0)
            H.pop_back();
        H.push_back(p);
    }

    // upper hull
    const size_t lower_size = H.size();
    for (int i = (int)P.size()-2; i >= 0; --i)
    {
        const auto& p = P[(size_t)i];
        while (H.size() > lower_size && cross2(H[H.size()-2], H[H.size()-1], p) < 0)
            H.pop_back();
        H.push_back(p);
    }

    if (!H.empty()) H.pop_back(); // remove duplicate start/end

    return H; // CCW
}

static void dedupKeepMaxParam(std::vector<P>& pts)
{
  std::sort(pts.begin(), pts.end(), [](const P& a, const P& b){
    if (a.x != b.x) return a.x < b.x;
    if (a.y != b.y) return a.y < b.y;
    return a.radius < b.radius; // così il max va in fondo
  });

  std::vector<P> out;
  out.reserve(pts.size());

  for (size_t i = 0; i < pts.size(); )
  {
    size_t j = i + 1;
    P best = pts[i];
    while (j < pts.size() && pts[j].x == pts[i].x && pts[j].y == pts[i].y)
    {
      if (pts[j].radius > best.radius) best = pts[j];
      ++j;
    }
    out.push_back(best);
    i = j;
  }

  pts.swap(out);
}

static std::vector<P> convexHull(std::vector<P> pts)
{
  if (pts.size() <= 2) return pts;

  std::sort(pts.begin(), pts.end(), [](const P& a, const P& b){
    return a.x < b.x || (a.x == b.x && a.y < b.y);
  });

  std::vector<P> H;
  H.reserve(pts.size()*2);

  // lower
  for (const auto& p : pts) {
    while (H.size() >= 2 && cross2(H[H.size()-2], H[H.size()-1], p) < 0)
      H.pop_back();
    H.push_back(p);
  }

  // upper
  size_t lower_size = H.size();
  for (int i = (int)pts.size()-2; i >= 0; --i) {
    const auto& p = pts[(size_t)i];
    while (H.size() > lower_size && cross2(H[H.size()-2], H[H.size()-1], p) < 0)
      H.pop_back();
    H.push_back(p);
  }

  if (!H.empty()) H.pop_back(); // ultimo = primo
  return H; // CCW, convessa
}

/*static LongestEdgeResult bestGapOnHull(const std::vector<P>& H)
{
  LongestEdgeResult out{};
  out.length = 0.0;
  out.clearance = 0.0;
  out.ia = out.ib = -1;

  const int n = (int)H.size();
  if (n < 2) return out;

  double best_clear = -1e300; // molto piccolo
  double best_d2 = -1.0;

  for (int i = 0; i < n; ++i) {
    int j = (i + 1) % n;

    const double dx = H[j].x - H[i].x;
    const double dy = H[j].y - H[i].y;
    const double d2 = dx*dx + dy*dy;
    const double d  = std::sqrt(d2);

    const double clear = d - H[i].radius - H[j].radius;

    // criterio principale: massima clearance
    // tie-break: distanza centri maggiore
    if (clear > best_clear + 1e-12 || (std::abs(clear - best_clear) <= 1e-12 && d2 > best_d2)) {
      best_clear = clear;
      best_d2 = d2;
      out.a = H[i];
      out.b = H[j];
      out.ia = i;
      out.ib = j;
      out.length = d;
      out.clearance = clear;
    }
  }

  // se tutti i clear sono negativi, comunque restituiamo il "meno peggio"
  return out;
}*/

static std::vector<LongestEdgeResult> gapsOnHullAboveDelta(const std::vector<P>& H, double delta)
{
  std::vector<LongestEdgeResult> outs;

  const int n = (int)H.size();
  if (n < 2) return outs;

  outs.reserve(n);

  for (int i = 0; i < n; ++i)
  {
    const int j = (i + 1) % n;

    const double dx = H[j].x - H[i].x;
    const double dy = H[j].y - H[i].y;
    const double d2 = dx*dx + dy*dy;
    const double d  = std::sqrt(d2);

    const double clear = d - H[i].radius - H[j].radius;

    if (clear + 1e-12 >= delta)   // clearance >= delta (tolleranza)
    {
      LongestEdgeResult out{};
      out.a = H[i];
      out.b = H[j];
      out.ia = i;
      out.ib = j;
      out.length = d;
      out.clearance = clear;
      outs.push_back(out);
    }
  }

  // Ordina: prima clearance, poi distanza tra centri
  std::sort(outs.begin(), outs.end(),
    [](const LongestEdgeResult& u, const LongestEdgeResult& v)
    {
      if (std::abs(u.clearance - v.clearance) > 1e-12)
        return u.clearance > v.clearance;
      return (u.length*u.length) > (v.length*v.length);
    });

  return outs;  // può essere vuoto se nessun gap >= delta
}

Eigen::Vector2d polygonAreaCentroid(const std::vector<point>& poly)
{
    const int n = (int)poly.size();
    if (n < 3) {
        // fallback: average of points
        Eigen::Vector2d c(0,0);
        for (auto& p : poly) c += Eigen::Vector2d(p.x, p.y);
        return (n>0) ? c / double(n) : Eigen::Vector2d(0,0);
    }

    double A2 = 0.0;  // 2*Area signed
    double Cx6A = 0.0;
    double Cy6A = 0.0;

    for (int i = 0; i < n; ++i)
    {
        const auto& p1 = poly[i];
        const auto& p2 = poly[(i+1) % n];
        const double cross = p1.x * p2.y - p2.x * p1.y;
        A2   += cross;
        Cx6A += (p1.x + p2.x) * cross;
        Cy6A += (p1.y + p2.y) * cross;
    }

    if (std::abs(A2) < 1e-12) {
        // degenerate polygon (almost line)
        Eigen::Vector2d c(0,0);
        for (auto& p : poly) c += Eigen::Vector2d(p.x, p.y);
        return c / double(n);
    }

    const double A = 0.5 * A2;
    return Eigen::Vector2d(Cx6A / (6.0 * A), Cy6A / (6.0 * A));
}

static Eigen::Vector2d closestPointOnSegment(
    const Eigen::Vector2d& P,
    const Eigen::Vector2d& A,
    const Eigen::Vector2d& B,
    double* t_out = nullptr)
{
  Eigen::Vector2d d = B - A;
  double denom = d.squaredNorm();

  double t = 0.0;
  if (denom > 1e-12) {
    t = (P - A).dot(d) / denom;
    t = std::clamp(t, 0.0, 1.0);
  } // se A==B, t resta 0 e Q=A

  if (t_out) *t_out = t;
  return A + t * d;
}

ClosestOnPolygonResult closestPointOnPolygonBoundary(
    const Eigen::Vector2d& P,
    const std::vector<Eigen::Vector2d>& poly)
{
  ClosestOnPolygonResult best;
  best.edge_index = -1;
  best.t = 0.0;
  best.dist2 = std::numeric_limits<double>::infinity();

  const int n = (int)poly.size();
  if (n < 2) return best;

  for (int i = 0; i < n; ++i) {
    const Eigen::Vector2d& A = poly[i];
    const Eigen::Vector2d& B = poly[(i + 1) % n];

    double t = 0.0;
    Eigen::Vector2d Q = closestPointOnSegment(P, A, B, &t);
    double d2 = (P - Q).squaredNorm();

    if (d2 < best.dist2) {
      best.dist2 = d2;
      best.closest_point = Q;
      best.edge_index = i;
      best.t = t;
    }
  }
  return best;
}

static inline double norm2(double x, double y) { return x*x + y*y; }

static double recomputeRadiusFromCentroid(const obstacle& o)
{
  if (o.polygon.empty()) return std::max(0.0, o.radius);

  double r2 = 0.0;
  for (const auto& p : o.polygon)
  {
    const double dx = p.x - o.x;
    const double dy = p.y - o.y;
    r2 = std::max(r2, dx*dx + dy*dy);
  }
  return std::sqrt(r2);
}

// ------------- Predizione -------------
static void predictTracks(std::vector<Track>& tracks, const AssocParams& p)
{
  for (auto& tr : tracks)
  {
    // predizione posizione
    tr.state.x += tr.state.vx * p.dt;
    tr.state.y += tr.state.vy * p.dt;

    // fai decadere la velocità verso zero (utile per “statici”)
    tr.state.vx *= p.vel_decay;
    tr.state.vy *= p.vel_decay;

    // (opzionale) se vuoi proprio forzare staticità:
    // if (tr.state.is_static) { tr.state.vx *= p.vel_decay; tr.state.vy *= p.vel_decay; }
  }
}

// ------------- Matrice costi N x M (flattened i + j*N) -------------
static distMatrix_t buildCostMatrix(const std::vector<Track>& tracks,
                                    const std::vector<obstacle>& dets)
{
  const size_t N = tracks.size();
  const size_t M = dets.size();
  distMatrix_t Cost(N * M, 0.0);

  for (size_t i = 0; i < N; ++i)
  {
    for (size_t j = 0; j < M; ++j)
    {
      const double dx = tracks[i].state.x - dets[j].x;
      const double dy = tracks[i].state.y - dets[j].y;
      Cost[i + j * N] = std::sqrt(dx*dx + dy*dy);
    }
  }
  return Cost;
}

// ------------- Hungarian solve -------------
static assignments_t solveHungarian(const distMatrix_t& Cost, size_t N, size_t M)
{
  AssignmentProblemSolver APS;
  assignments_t assignment;
  APS.Solve(Cost, N, M, assignment, AssignmentProblemSolver::optimal);
  return assignment; // size N, assignment[i] = j oppure -1
}

// ------------- Gating -------------
static void applyGating(assignments_t& assignment,
                        const distMatrix_t& Cost,
                        size_t N,
                        double gate_dist)
{
  for (size_t i = 0; i < assignment.size(); ++i)
  {
    const int j = assignment[i];
    if (j < 0) continue;

    const double d = Cost[i + size_t(j) * N];
    if (d > gate_dist)
      assignment[i] = -1;
  }
}

// ------------- Update track assegnati -------------
static void updateAssignedTracks(std::vector<Track>& tracks,
                                 const std::vector<obstacle>& dets,
                                 const assignments_t& assignment,
                                 const AssocParams& p)
{
  for (size_t i = 0; i < tracks.size(); ++i)
  {
    const int j = assignment[i];
    if (j < 0) continue;

    Track& tr = tracks[i];
    const obstacle& meas = dets[size_t(j)];

    // conserva ID del track (stabile) — NON usare quello della detection (spesso -1)
    const int track_id = tr.state.id;

    // salva vecchia posizione per stimare v
    const double x_prev = tr.state.x;
    const double y_prev = tr.state.y;

    // Update posizione (filtro esponenziale)
    tr.state.x = (1.0 - p.alpha_pos) * tr.state.x + p.alpha_pos * meas.x;
    tr.state.y = (1.0 - p.alpha_pos) * tr.state.y + p.alpha_pos * meas.y;

    // Stima velocità dal delta posizione
    const double dt = std::max(p.dt, 1e-6);
    tr.state.vx = (tr.state.x - x_prev) / dt;
    tr.state.vy = (tr.state.y - y_prev) / dt;

    // Forza velocità a tendere a zero (come vuoi tu)
    tr.state.vx *= p.vel_decay;
    tr.state.vy *= p.vel_decay;

    // COPIA attributi “semantici” dalla detection
    tr.state.is_person = meas.is_person;
    tr.state.is_segment = meas.is_segment;

    // poligono: aggiorna con quello della misura più recente
    tr.state.polygon = meas.polygon;

    // ricalcola radius usando il centroide del TRACK (tr.state.x,y) e il poligono aggiornato
    tr.state.radius = recomputeRadiusFromCentroid(tr.state);

    // ripristina ID
    tr.state.id = track_id;

    tr.missed = 0;
  }
}

// ------------- Gestione non assegnati e delete -------------
static void markAndDeleteLost(std::vector<Track>& tracks,
                              const assignments_t& assignment,
                              int max_missed)
{
  for (size_t i = 0; i < tracks.size(); ++i)
    if (assignment[i] < 0)
      tracks[i].missed++;

  tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                              [&](const Track& tr){ return tr.missed > max_missed; }),
               tracks.end());
}

// ------------- Detection non usate: spawn nuovi track -------------
static std::vector<bool> computeUsedDetections(const assignments_t& assignment, size_t M)
{
  std::vector<bool> used(M, false);
  for (int j : assignment)
    if (j >= 0) used[size_t(j)] = true;
  return used;
}

static void spawnNewTracks(std::vector<Track>& tracks,
                           const std::vector<obstacle>& dets,
                           const std::vector<bool>& used,
                           int& next_id)
{
  for (size_t j = 0; j < dets.size(); ++j)
  {
    if (used[j]) continue;

    Track tr;
    tr.state = dets[j];

    // ID stabile nuovo (ignora id della detection)
    tr.state.id = next_id++;

    // velocità iniziale
    tr.state.vx = 0.0;
    tr.state.vy = 0.0;

    // radius coerente col centroide del track (che coincide con dets[j].x,y)
    tr.state.radius = recomputeRadiusFromCentroid(tr.state);

    tr.missed = 0;
    tracks.push_back(std::move(tr));
  }
}

// ------------- Pipeline completa -------------
void associateAndTrack(std::vector<Track>& tracks,
                              int& next_id,
                              const std::vector<obstacle>& detections,
                              const AssocParams& p)
{
  // 1) predizione
  predictTracks(tracks, p);

  const size_t N = tracks.size();
  const size_t M = detections.size();

  // casi limite
  if (N == 0)
  {
    std::vector<bool> used(M, false);
    spawnNewTracks(tracks, detections, used, next_id);
    return;
  }
  if (M == 0)
  {
    for (auto& tr : tracks) tr.missed++;
    tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                                [&](const Track& tr){ return tr.missed > p.max_missed; }),
                 tracks.end());
    return;
  }

  // 2) cost matrix
  distMatrix_t Cost = buildCostMatrix(tracks, detections);

  // 3) Hungarian
  assignments_t assignment = solveHungarian(Cost, N, M);

  // 4) gating
  applyGating(assignment, Cost, N, p.gate_dist);

  // 5) update assegnati
  updateAssignedTracks(tracks, detections, assignment, p);

  // 6) missed + delete
  markAndDeleteLost(tracks, assignment, p.max_missed);

  // 7) spawn per detection non usate (dopo gating)
  std::vector<bool> used = computeUsedDetections(assignment, M);
  spawnNewTracks(tracks, detections, used, next_id);
}

static inline double clamp_m1_p1(double x) {
    if (x < -1.0) return -1.0;
    if (x >  1.0) return  1.0;
    return x;
}

// Spezza UN cluster (lista di indici globali) in rami.
// - Taglia sempre junction deg>=3
// - Se NON ci sono junction, taglia anche i gomiti (deg==2) con angolo < theta_cut_deg
//   (per V/W spesso è necessario)
// - keep_cut_nodes_as_own_cluster: mette i cut-nodes in un cluster a parte
/*std::vector<std::vector<int>> split_cluster_branches_mst(
    const std::vector<Circle>& circles,
    const std::vector<int>& cluster_idx,   // indici globali dei cerchi nel cluster
    double tau,
    bool keep_cut_nodes_as_own_cluster = true,
    int min_branch_size = 1,
    double theta_cut_deg = 150.0
) {
    const int m = (int)cluster_idx.size();
    if (m <= 2) return {cluster_idx};

    // 1) Build edges del grafo interno al cluster (stessa regola del DSU)
    std::vector<Edge> edges;
    edges.reserve((size_t)m * (m - 1) / 2);

    for (int ii = 0; ii < m; ++ii) {
        int gi = cluster_idx[ii];
        for (int jj = ii + 1; jj < m; ++jj) {
            int gj = cluster_idx[jj];
            const auto& A = circles[gi];
            const auto& B = circles[gj];
            double d = dist_centri(A, B);
            if (d <= A.r + B.r + tau) {
                edges.push_back({ii, jj, d}); // ii/jj sono indici locali 0..m-1
            }
        }
    }
    if (edges.empty()) return {cluster_idx};

    // 2) Kruskal MST su nodi locali 0..m-1
    std::sort(edges.begin(), edges.end(),
              [](const Edge& a, const Edge& b){ return a.w < b.w; });

    DSU dsu_mst(m);
    std::vector<std::vector<int>> adj(m);
    int added = 0;

    for (const auto& e : edges) {
        int ru = dsu_mst.find(e.u);
        int rv = dsu_mst.find(e.v);
        if (ru == rv) continue;

        dsu_mst.unite(e.u, e.v);
        adj[e.u].push_back(e.v);
        adj[e.v].push_back(e.u);
        added++;
        if (added == m - 1) break;
    }
    if (added != m - 1) {
        // grafo interno non connesso (caso raro se DSU era coerente)
        return {cluster_idx};
    }

    // 3) Decide i cut-nodes
    std::vector<char> is_cut(m, 0);

    int junction_count = 0;
    for (int u = 0; u < m; ++u) {
        if ((int)adj[u].size() >= 3) {
            is_cut[u] = 1;
            junction_count++;
        }
    }

    // Se non ci sono junction, allora siamo in forma "catena" (V/W/curve):
    // taglia tutti i gomiti con angolo < theta_cut
    if (junction_count == 0) {
        const double theta_cut = theta_cut_deg * M_PI / 180.0;

        for (int u = 0; u < m; ++u) {
            if ((int)adj[u].size() != 2) continue;

            int a = adj[u][0];
            int b = adj[u][1];

            const Circle& Cu = circles[ cluster_idx[u] ];
            const Circle& Ca = circles[ cluster_idx[a] ];
            const Circle& Cb = circles[ cluster_idx[b] ];

            double v1x = Ca.x - Cu.x, v1y = Ca.y - Cu.y;
            double v2x = Cb.x - Cu.x, v2y = Cb.y - Cu.y;

            double n1 = std::hypot(v1x, v1y);
            double n2 = std::hypot(v2x, v2y);
            if (n1 < 1e-9 || n2 < 1e-9) continue;

            double cosang = clamp_m1_p1((v1x*v2x + v1y*v2y) / (n1*n2));
            double theta = std::acos(cosang); // 0..pi

            if (theta < theta_cut) {
                is_cut[u] = 1; // gomito -> taglia (V/W)
            }
        }
    }

    bool any_cut = false;
    for (char c : is_cut) if (c) { any_cut = true; break; }
    if (!any_cut) return {cluster_idx};

    // 4) Componenti connesse ignorando i cut-nodes => rami
    std::vector<char> vis(m, 0);
    std::vector<std::vector<int>> branches;

    for (int i = 0; i < m; ++i) {
        if (is_cut[i] || vis[i]) continue;

        std::vector<int> comp_local;
        std::queue<int> q;
        q.push(i);
        vis[i] = 1;

        while (!q.empty()) {
            int u = q.front(); q.pop();
            comp_local.push_back(u);
            for (int v : adj[u]) {
                if (is_cut[v] || vis[v]) continue;
                vis[v] = 1;
                q.push(v);
            }
        }

        if ((int)comp_local.size() >= min_branch_size) {
            std::vector<int> out;
            out.reserve(comp_local.size());
            for (int li : comp_local) out.push_back(cluster_idx[li]); // locale->globale
            std::sort(out.begin(), out.end());
            branches.push_back(std::move(out));
        }
    }

    // 5) opzionale: cluster con i cut-nodes (giunzioni + gomiti)
    if (keep_cut_nodes_as_own_cluster) {
        std::vector<int> cuts;
        for (int i = 0; i < m; ++i) if (is_cut[i]) cuts.push_back(cluster_idx[i]);
        if (!cuts.empty()) {
            std::sort(cuts.begin(), cuts.end());
            branches.push_back(std::move(cuts));
        }
    }

    if (branches.empty()) return {cluster_idx};
    return branches;
}

// Wrapper: spezza tutti i cluster DSU in rami
std::vector<std::vector<int>> recluster_split_branches(
    const std::vector<Circle>& circles,
    const std::vector<std::vector<int>>& dsu_clusters,
    double tau,
    bool keep_cut_nodes_as_own_cluster = true,
    int min_branch_size = 1,
    double theta_cut_deg = 150.0
) {
    std::vector<std::vector<int>> out;
    for (const auto& cl : dsu_clusters) {
        auto pieces = split_cluster_branches_mst(
            circles, cl, tau,
            keep_cut_nodes_as_own_cluster,
            min_branch_size,
            theta_cut_deg
        );
        out.insert(out.end(), pieces.begin(), pieces.end());
    }

    // determinismo (opzionale)
    std::sort(out.begin(), out.end(), [](const auto& a, const auto& b){
        if (a.empty() || b.empty()) return a.size() < b.size();
        return a.front() < b.front();
    });

    return out;
}*/

MPC_diffDrive_fblin::MPC_diffDrive_fblin() {
    // Initialize class variables
    _robotParamsInitialized = false;
    _MPCparamsInitialized = false;
    _FBLINparamsInitialized = false;
    _controllerInitialized = false;
    _is_started = false;
    _is_started_complete = false;
    _init_P = false;
    _n_fail = 0;
    _is_blocked = false;
    _is_blocked_inside_hull = false;
    _is_closed = false;

    _linearVelocity = _angularVelocity = 0.0;

    // Initialize pointers
    _solver = NULL;
    _fblinController = _fblinSimController = NULL;
    _debug = _info = _error = NULL;
}

MPC_diffDrive_fblin::~MPC_diffDrive_fblin() {
    // Destroy solver object
    if (_solver)
    {
        delete _solver;
        _solver = NULL;
    }

    // Destroy linearization object
    if (_fblinController)
    {
        delete _fblinController;
        _fblinController = NULL;
    }
    if (_fblinSimController)
    {
        delete _fblinSimController;
        _fblinSimController = NULL;
    }
}

//void MPC_diffDrive_fblin::set_MPCparams(double samplingTime, int predictionHorizon, double q, double r, double a_max, double sp, int n_obs, double sa, double r_c, double r_red, double T_gp, double sensorRange, double e, int N_samples) {


    // Set MPC params
//    set_MPCparams(samplingTime, predictionHorizon, q, r, a_max, sp, n_obs, sa, r_c, r_red, T_gp, sensorRange, e, N_samples, lb, ub);
//}

void MPC_diffDrive_fblin::set_MPCparams(double samplingTime, int predictionHorizon, double q, double r, double a_max, double sp, double sa, double r_red, double T_gp, double sensorRange, double e, int N_samples) {
    // Set MPC parameters
    _q = q;
    _r = r;
    _a_max = a_max;
    _sp = sp;

    _sa = sa;
    //_r_c = r_c;
    _r_red = r_red;
    _T_gp = T_gp;
    _e = e;
    _N_samples = N_samples;
    _sensorRange = sensorRange;
    _MPC_Ts = samplingTime;
    _N = predictionHorizon;
    _dsl_a = 0.04;

    _optimVect_no_sl = Eigen::VectorXd::Zero(2*_N);

    // Set the initialization flag
    _MPCparamsInitialized = true;
}

void MPC_diffDrive_fblin::set_n_obs(const int n_obs) {
    // Set feedback linearization parameters
    _n_obs = n_obs;
    _n_groups = 0;
}

void MPC_diffDrive_fblin::set_r_c(const double& r_c) {
    // Set feedback linearization parameters
    _r_c = r_c;
}

void MPC_diffDrive_fblin::set_FBLINparams(double samplingTime, double pointPdistance) {
    // Set feedback linearization parameters
    _fblin_Ts = samplingTime;
    _Pdist = pointPdistance;

    // Set the initialization flag
    _FBLINparamsInitialized = true;
}

void MPC_diffDrive_fblin::set_robotParams(double wheelVelMax, double wheelVelMin, double wheelRadius, double track) {
    // Set robot parameters
    _wheelVelMax = wheelVelMax;
    _wheelVelMin = wheelVelMin;
    _wheelRadius = wheelRadius;
    _track = track;

    // Set the initialization flag
    _robotParamsInitialized = true;
}

bool MPC_diffDrive_fblin::initialize(const Eigen::Vector2d& final_ref_pos) {
    /** Preliminary checks */
    if (!_robotParamsInitialized)
    {
        errorMsg("[MPC_diffDrive_fblin.initialize] Call set_robotParams() before calling initialize()");
        return false;
    }
    if (!_FBLINparamsInitialized)
    {
        errorMsg("[MPC_diffDrive_fblin.initialize] Call set_FBLINparams() before calling initialize()");
        return false;
    }
    if (!_MPCparamsInitialized)
    {
        errorMsg("[MPC_diffDrive_fblin.initialize] Call set_MPCparams() before calling initialize()");
        return false;
    }
    if (std::remainder(_MPC_Ts, _fblin_Ts)>=1.0e-15)
    {
        errorMsg("[MPC_diffDrive_fblin.initialize] MPC and feedback linearization sampling times should be multiple");
        return false;
    }



    if (_n_obs > 0)
    {
        clustering(final_ref_pos);
    }

    std::ostringstream oss;
    oss << "[";

    for (int i = 0; i < _labels.size(); ++i) {
        oss << _labels[i];
        if (i + 1 < _labels.size())
            oss << ", ";
    }

    oss << "]";
    
    oss << _n_obs;

    oss << " ";

    oss << _n_groups;

    std::string s = oss.str();

    debugMsg(s);

    const std::vector<double> lb(2, -100.0);
    const std::vector<double> ub(2, +100.0);

    std::size_t totalSize = 2*_N + _n_groups + 2;
    _lowerBound.assign(totalSize, 0.0);
    _upperBound.assign(totalSize, 0.0);

    for (std::size_t i = 0; i < _N; ++i) {
        _lowerBound[2*i]     = lb[0];
        _lowerBound[2*i + 1] = lb[1];

        _upperBound[2*i]     = ub[0];
        _upperBound[2*i + 1] = ub[1];
    }

    for (std::size_t i = 0; i < _n_groups; ++i) {
        std::size_t idx = 2*_N + i;
        _lowerBound[idx] = 0.0;
        _upperBound[idx] = 1.0;
    }

    _lowerBound[2*_N + _n_groups]     = 0.0;
    _upperBound[2*_N + _n_groups]     = _dsl_a;
    _lowerBound[2*_N + _n_groups + 1] = 0.0;
    _upperBound[2*_N + _n_groups + 1] = _dsl_a;

    _H = Eigen::MatrixXd::Zero(2*_N + _n_groups +2, 2*_N + _n_groups + 2);
    _f = Eigen::VectorXd::Zero(2*_N + _n_groups +2);


    if (!_controllerInitialized)
    {
        infoMsg("[MPC_diffDrive_fblin] Initializing MPC controller");
        debugMsg("[MPC_diffDrive_fblin] Initializing plant matrices");
        // Initialize plant matrices
        _plant_A = Eigen::MatrixXd::Identity(2,2);
        _plant_B = _MPC_Ts*Eigen::MatrixXd::Identity(2,2);

        // Initialize MPC controller parameters
        _k = -1.0/(2.0*_MPC_Ts);
        _p = (_q+pow(_k, 2.0)*_r)/(1.0-pow(1.0+_MPC_Ts*_k, 2.0));

        // Initialize MPC controller matrices
        compute_AcalMatrix();
        compute_BcalMatrix();
        compute_QcalMatrix();
        compute_RcalMatrix();

        _Ain_vel = Eigen::MatrixXd::Zero(2*2*_N, 2*_N);
        _Bin_vel = Eigen::VectorXd::Zero(2*2*_N);

        // Initialize actual and reference data
        _actX = _actY = _actYaw = 0.0;
        _actXP = _actYP = 0.0;

        // Initialize the linearization controller
        _fblinController = new fblin_unicycle(_Pdist);
        _fblinSimController = new fblin_unicycle(_Pdist);

        // Set the initialization flag
        _controllerInitialized = true;

    }else
    {
        _solver->removeConstraint(totConstrain);

        totConstrain.clear();

        if (_solver)
        {
            delete _solver;
            _solver = NULL;
        }
    }

    _optimVect = Eigen::VectorXd::Zero(2*_N + _n_groups + 2);

    _solver = new GUROBIsolver(GUROBI_LICENSEID, GUROBI_USERNAME);
    if (!_solver->initProblem(2*_N + _n_groups + 2, _lowerBound, _upperBound))
    {
        errorMsg("[MPC_diffDrive_fblin.initialize] Error initializing the solver");
        return false;
    }

    return true;
}

void MPC_diffDrive_fblin::prediction()
{
    _XP_predicted = Eigen::VectorXd::Zero(_N);
    _YP_predicted = Eigen::VectorXd::Zero(_N);
    _v_long = Eigen::VectorXd::Zero(_N);
    _predictRobotState = Eigen::VectorXd::Zero(3*_N);

    Eigen::VectorXd theta = Eigen::VectorXd::Zero(_N);

    theta(0) = _actYaw;


    int j=2;
    for (int i = 1; i<_N; i++)
    {
        double omega;
        omega=(_optimVect_no_sl(j+1)*cos(theta(i-1))-_optimVect_no_sl(j)*sin(theta(i-1)))/_Pdist;
        theta(i)=theta(i-1)+omega*_MPC_Ts;
        _v_long(i-1)=_optimVect_no_sl(j)*cos(theta(i-1))+_optimVect_no_sl(j+1)*sin(theta(i-1));
        j=j+2;
    }

    _v_long(_N-1)=_v_long(_N-2);

    Eigen::VectorXd xp_est = Eigen::VectorXd::Zero(_N);
    Eigen::VectorXd yp_est = Eigen::VectorXd::Zero(_N);
    Eigen::VectorXd predictRobotState_P = Eigen::VectorXd::Zero(2*_N);


    predictRobotState_P = _Acal*Eigen::Vector2d(_prevXP, _prevYP) + _Bcal*_optimVect_no_sl;
    for (int i = 0; i<_N; i++)
    {
        xp_est(i) = predictRobotState_P(2*i);
        yp_est(i) = predictRobotState_P(2*i+1);
        _XP_predicted(i) = xp_est(i);
        _YP_predicted(i) = yp_est(i);
    }

    xp_est(0) = _actXP;
    yp_est(0) = _actYP;

    _XP_predicted(0) = xp_est(0);
    _YP_predicted(0) = yp_est(0);

    Eigen::VectorXd x = Eigen::VectorXd::Zero(_N);
    Eigen::VectorXd y = Eigen::VectorXd::Zero(_N);

    x(0) = _actX;
    y(0) = _actY;

    for (int i = 1; i<_N; i++)
    {
        x(i) = xp_est(i) - _Pdist*std::cos(theta(i));
        y(i) = yp_est(i) - _Pdist*std::sin(theta(i));
    }

    for (int i = 0; i<_N; i++)
    {
        _predictRobotState(3*i) = x(i);
        _predictRobotState(3*i+1) = y(i);
        _predictRobotState(3*i+2) = theta(i);
    }

}

pred_step MPC_diffDrive_fblin::prediction_step(double prevXP, double prevYP, double actX, double actY, double actYaw, double actXP, double actYP, Eigen::VectorXd optimVect_no_sl_samples, int k)
{
    Eigen::VectorXd XP_predicted = Eigen::VectorXd::Zero(_N);
    Eigen::VectorXd YP_predicted = Eigen::VectorXd::Zero(_N);
    Eigen::VectorXd v_long = Eigen::VectorXd::Zero(_N);

    pred_step prediction = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    Eigen::VectorXd theta = Eigen::VectorXd::Zero(_N);

    theta(0) = actYaw;


    int j=2;
    for (int i = 1; i<_N; i++)
    {
        double omega;
        omega=(optimVect_no_sl_samples(j+1)*cos(theta(i-1))-optimVect_no_sl_samples(j)*sin(theta(i-1)))/_Pdist;
        theta(i)=theta(i-1)+omega*_MPC_Ts;
        v_long(i-1)=optimVect_no_sl_samples(j)*cos(theta(i-1))+optimVect_no_sl_samples(j+1)*sin(theta(i-1));
        j=j+2;
    }

    v_long(_N-1)=v_long(_N-2);

    Eigen::VectorXd xp_est = Eigen::VectorXd::Zero(_N);
    Eigen::VectorXd yp_est = Eigen::VectorXd::Zero(_N);
    Eigen::VectorXd predictRobotState_P = Eigen::VectorXd::Zero(2*_N);


    predictRobotState_P = _Acal*Eigen::Vector2d(prevXP, prevYP) + _Bcal*optimVect_no_sl_samples;
    for (int i = 0; i<_N; i++)
    {
        xp_est(i) = predictRobotState_P(2*i);
        yp_est(i) = predictRobotState_P(2*i+1);
        XP_predicted(i) = xp_est(i);
        YP_predicted(i) = yp_est(i);
    }

    xp_est(0) = actXP;
    yp_est(0) = actYP;

    XP_predicted(0) = xp_est(0);
    YP_predicted(0) = yp_est(0);

    Eigen::VectorXd x = Eigen::VectorXd::Zero(_N);
    Eigen::VectorXd y = Eigen::VectorXd::Zero(_N);

    x(0) = actX;
    y(0) = actY;

    for (int i = 1; i<_N; i++)
    {
        x(i) = xp_est(i) - _Pdist*std::cos(theta(i));
        y(i) = yp_est(i) - _Pdist*std::sin(theta(i));
    }

    prediction.XP_pred = XP_predicted(k-1);
    prediction.YP_pred = YP_predicted(k-1);
    if (k < 20){
    prediction.vxP_pred = optimVect_no_sl_samples(2*k);
    prediction.vyP_pred = optimVect_no_sl_samples(2*k+1);
    }else{
    prediction.vxP_pred = optimVect_no_sl_samples(2*19);
    prediction.vyP_pred = optimVect_no_sl_samples(2*19+1);    
    }
    prediction.x_pred = x(k-1);
    prediction.y_pred = y(k-1);
    prediction.theta_pred = theta(k-1);
    prediction.vx_pred = v_long(k-1)*std::cos(theta(k-1));
    prediction.vy_pred = v_long(k-1)*std::sin(theta(k-1));

    return prediction;

}

void MPC_diffDrive_fblin::compute_wheelAccelerationConstraint()
{
    _A_dvel = Eigen::MatrixXd::Zero(4*_N, 2*_N + _n_groups + 2);
    _b_dvel = Eigen::VectorXd::Zero(4*_N);

    double deltaV_max = _a_max*_MPC_Ts;
    Eigen::VectorXd deltaV = deltaV_max*Eigen::VectorXd::Ones(4*_N);

    Eigen::MatrixXd A_var = Eigen::MatrixXd::Zero(4*_N, 2*_N);
    A_var.block(0, 0, 2*_N, 2*_N) = Eigen::MatrixXd::Identity(2*_N, 2*_N);
    A_var.block(2*_N, 0, 2*_N, 2*_N) = - Eigen::MatrixXd::Identity(2*_N, 2*_N);

    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(2*_N, 2*_N);

    for (int i = 0; i < 2*_N; i += 2)
    {
        V.block(i, i, 2, 2) = Eigen::MatrixXd::Identity(2, 2);

        _A_dvel.block(i, 2*_N + _n_groups, 2, 2) = - Eigen::MatrixXd::Identity(2, 2);
        _A_dvel.block(2*_N + i, 2*_N + _n_groups, 2, 2) = - Eigen::MatrixXd::Identity(2, 2);

        if (i < 2*_N - 3)
        {
            V.block(i + 2, i, 2, 2) = - Eigen::MatrixXd::Identity(2, 2);
        }
    }

    Eigen::VectorXd v0 = Eigen::VectorXd::Zero(2*_N);
    v0(0) = _optimVect_no_sl(0);
    v0(1) = _optimVect_no_sl(1);

    _A_dvel.block(0, 0, 4*_N, 2*_N) = A_var*V;

    _b_dvel = deltaV + A_var*v0;
}

personal_space MPC_diffDrive_fblin::personalSpaceFunction_mod(Eigen::RowVector2d v_ped, Eigen::Vector2d relative_velocity, Eigen::Vector2d robot_pos, Eigen::Vector2d ped_pos, double r_c)
{
    Eigen::Vector2d relative_velocity_wrt_robot = - relative_velocity;

    double c = - (ped_pos - robot_pos).dot(relative_velocity_wrt_robot);

    personal_space personalSpace = {Eigen::Vector2d::Zero()};

    if (c > 0.0)
    {
        if (v_ped.norm() != 0.0)
        {
            personalSpace.PSf(0) = std::min(1.5, std::max(r_c, 3*relative_velocity.norm()));
            personalSpace.PSf(1) = 1.5*personalSpace.PSf(0);
            personalSpace.scenario = 0;
        }else
        {
            personalSpace.PSf(0) = std::min(1.5, std::max(r_c, 3*relative_velocity.norm()));
            personalSpace.PSf(1) = personalSpace.PSf(0);
            personalSpace.scenario = 1;
        }
    }else
    {
        if (v_ped.norm() != 0.0)
        {
            personalSpace.PSf(0) = std::min(1.5, std::max(r_c, 3*v_ped.norm()));
            personalSpace.PSf(1) = 1.5*personalSpace.PSf(0);
            personalSpace.scenario = 2;
        }else
        {
            personalSpace.PSf(0) = r_c;
            personalSpace.PSf(1) = personalSpace.PSf(0);
            personalSpace.scenario = 3;
        }
    }

    personalSpace.theta_heading = std::atan2(v_ped(1), v_ped(0));

    return personalSpace;
}

tangent_local MPC_diffDrive_fblin::findTangentPoints_local_mod(Eigen::Vector2d robot_pos_local, double a, double b)
{
    tangent_local findTangentPoints_local = {Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero()};

    double x0 = robot_pos_local(0);
    double y0 = robot_pos_local(1);

    double R =std::hypot(x0/a, y0/b);

    if (R <= 1.0)
    {
        errorMsg("Punto non esterno all''ellisse.");
        return findTangentPoints_local;
    }

    double phi = std::atan2(y0/b, x0/a);
    double delta = std::acos(1/R);

    double t1 = phi + delta;
    double t2 = phi - delta;

    findTangentPoints_local.Q1 << a*std::cos(t1), b*std::sin(t1);
    findTangentPoints_local.Q2 << a*std::cos(t2), b*std::sin(t2);

    return findTangentPoints_local;

}

tangent MPC_diffDrive_fblin::findTangentPoints_mod(const Eigen::Vector2d& robot_pos, const Eigen::Vector2d& ped_pos, Eigen::Vector2d PSf, double theta_heading, double r_red, double r_c)
{
    tangent tangent_points = {Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero(), false};

    Eigen::Matrix2d R;
    R << std::cos(theta_heading), std::sin(theta_heading),
        -std::sin(theta_heading), std::cos(theta_heading);

    Eigen::Matrix2d R_inverse;
    R_inverse << std::cos(theta_heading), -std::sin(theta_heading),
                 std::sin(theta_heading), std::cos(theta_heading);

    Eigen::Vector2d robot_pos_local = R*(robot_pos - ped_pos);

    double delta = (r_c - r_red);

    double a = PSf(1) + delta;
    double b = PSf(0) + delta;

    if ((robot_pos_local(0) >= 0.0 && (std::pow(robot_pos_local(0), 2)/std::pow(a, 2) + std::pow(robot_pos_local(1), 2)/std::pow(b, 2) > 1.0)) || (robot_pos_local(0) < 0.0 && (std::pow(robot_pos_local(0), 2)/std::pow(b, 2) + std::pow(robot_pos_local(1), 2)/std::pow(b, 2) > 1.0)))
    {
        tangent_local tangents_circle_local = findTangentPoints_local_mod(robot_pos_local, b, b);

        if (tangents_circle_local.Q1(0) < 0.0 && tangents_circle_local.Q2(0) < 0.0)
        {
            tangent_points.Q1 = R_inverse*tangents_circle_local.Q1 + ped_pos;
            tangent_points.Q2 = R_inverse*tangents_circle_local.Q2 + ped_pos;
        }else
        {
            tangent_local tangents_ellips_local = findTangentPoints_local_mod(robot_pos_local, a, b);
            if  (tangents_circle_local.Q1(0) < 0.0 && tangents_circle_local.Q2(0) >=0.0)
            {
                tangent_points.Q1 = R_inverse*tangents_circle_local.Q1 + ped_pos;
                tangent_points.Q2 = R_inverse*tangents_ellips_local.Q2 + ped_pos;
            }else if (tangents_circle_local.Q1(0) >= 0.0 && tangents_circle_local.Q2(0) < 0.0)
            {
                tangent_points.Q1 = R_inverse*tangents_ellips_local.Q1 + ped_pos;
                tangent_points.Q2 = R_inverse*tangents_circle_local.Q2 + ped_pos;
            }else
            {
                tangent_points.Q1 = R_inverse*tangents_ellips_local.Q1 + ped_pos;
                tangent_points.Q2 = R_inverse*tangents_ellips_local.Q2 + ped_pos;
            }
        }
    }else{

        tangent_points.robot_inside_obst = true;

    }

    return tangent_points;
}

tangent MPC_diffDrive_fblin::findTangentPoints_mod_objects(const Eigen::Vector2d& robot_pos, double r_red, double r_c, const obstacle& obstacle)
{
    tangent tangent_points = {Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero(), false};

    double delta = (r_c - r_red);

    std::vector<point> polygon = obstacle.polygon;

    double signedArea = 0.0;
    for (size_t i = 0; i < polygon.size(); ++i)
    {
        const auto& p1 = polygon[i];
        const auto& p2 = polygon[(i+1) % polygon.size()];
        signedArea += (p1.x * p2.y - p2.x * p1.y);
    }
    signedArea *= 0.5;

    if (signedArea < 0)
    {
        std::reverse(polygon.begin(), polygon.end());
    }

    Eigen::VectorXd cross_products = Eigen::VectorXd::Zero(static_cast<int>(polygon.size()));
    Eigen::MatrixXd polygon_enlarged = Eigen::MatrixXd::Zero(static_cast<int>(polygon.size()), 2);
    Eigen::MatrixXd vect_Q_tot = Eigen::MatrixXd::Zero(static_cast<int>(polygon.size()), 2);
    Eigen::VectorXd vector_cross = Eigen::VectorXd::Zero(static_cast<int>(polygon.size()));

    for (int i = 0; i < polygon.size(); ++i)
    {
        const auto& p = polygon[i];
        double v_norm = std::hypot(p.x - obstacle.x, p.y - obstacle.y) + 1e-6;
        double vx = (p.x - obstacle.x)/v_norm;
        double vy = (p.y - obstacle.y)/v_norm;
        polygon_enlarged(i, 0) = p.x + vx*delta;
        polygon_enlarged(i, 1) = p.y + vy*delta;

    }

    for (int i = 0; i < polygon.size(); ++i)
    {
        const auto& p1 = polygon_enlarged.row(i);
        const auto& p2 =  polygon_enlarged.row((i+1) % static_cast<int>(polygon.size()));
        cross_products(i) = (p2(0) - p1(0))*(robot_pos(1) - p1(1)) - (p2(1) - p1(1))*(robot_pos(0) - p1(0));
    }

    if ((cross_products.array() >= 0).all())
    {

        tangent_points.robot_inside_obst = true;

    }else{

      vect_Q_tot = polygon_enlarged - robot_pos.transpose().replicate(static_cast<int>(polygon.size()), 1);

        tangent_points.robot_inside_obst = false;
        for (int k = 0; k < polygon.size(); k++)
        {
            vector_cross(k) = std::atan2(vect_Q_tot(0, 0)*vect_Q_tot(k, 1) - vect_Q_tot(0, 1)*vect_Q_tot(k, 0), vect_Q_tot(0, 0)*vect_Q_tot(k, 0) + vect_Q_tot(0, 1)*vect_Q_tot(k, 1));
        }

        int I_min, I_max;
        vector_cross.maxCoeff(&I_max);
        vector_cross.minCoeff(&I_min);

        tangent_points.Q1 = vect_Q_tot.block(I_min, 0, 1, 2).transpose() + robot_pos;
        tangent_points.Q2 = vect_Q_tot.block(I_max, 0, 1, 2).transpose() + robot_pos;

    }

    return tangent_points;
}

double MPC_diffDrive_fblin::compute_s_objects(const Eigen::RowVector2d& n, double r_red, double r_c, const obstacle& obstacle, const Eigen::Vector2d& pos)
{
    double delta = (r_c - r_red);

    std::vector<point> polygon = obstacle.polygon;

    Eigen::MatrixXd polygon_enlarged = Eigen::MatrixXd::Zero(static_cast<int>(polygon.size()), 2);
    Eigen::VectorXd s = Eigen::VectorXd::Zero(static_cast<int>(polygon.size()));
    double s_tot = 0.0;

    auto sign_nz = [](double x){ return (x >= 0.0) ? 1.0 : -1.0; };

    for (int i = 0; i < polygon.size(); ++i)
    {
        const auto& p = polygon[i];
        double vx;
        double vy;
        if (!obstacle.is_segment){
            double v_norm = std::hypot(p.x - obstacle.x, p.y - obstacle.y) + 1e-6;
            vx = (p.x - obstacle.x)/v_norm;
            vy = (p.y - obstacle.y)/v_norm;
            polygon_enlarged(i, 0) = p.x + vx*delta;
            polygon_enlarged(i, 1) = p.y + vy*delta;
        }else{
            Eigen::Vector2d normal = Eigen::Vector2d::Zero();
            normal << (obstacle.y2 - obstacle.y1), -(obstacle.x2 - obstacle.x1);
            Eigen::Vector2d longitudinal = Eigen::Vector2d::Zero();
            longitudinal << (obstacle.x2 - obstacle.x1), (obstacle.y2 - obstacle.y1);
            double normal_norm = normal.norm() + 1e-6;
            double longitudinal_norm = longitudinal.norm() + 1e-6;

            Eigen::Vector2d normal_normalized = normal/normal_norm;
            Eigen::Vector2d longitudinal_normalized = longitudinal/longitudinal_norm;

            Eigen::Vector2d delta_point;
            delta_point << (p.x - obstacle.x), (p.y - obstacle.y);

            polygon_enlarged(i, 0) = p.x + sign_nz(delta_point.dot(normal))*normal_normalized(0)*delta;
            polygon_enlarged(i, 1) = p.y + sign_nz(delta_point.dot(normal))*normal_normalized(1)*delta;
            polygon_enlarged(i, 0) = polygon_enlarged(i, 0) + sign_nz(delta_point.dot(longitudinal))*longitudinal_normalized(0)*delta;
            polygon_enlarged(i, 1) = polygon_enlarged(i, 1) + sign_nz(delta_point.dot(longitudinal))*longitudinal_normalized(1)*delta;
        }

    }

    for (int i = 0; i < polygon.size(); ++i)
    {
        double vx = (polygon_enlarged(i, 0) - obstacle.x);
        double vy = (polygon_enlarged(i, 1) - obstacle.y);
        Eigen::Vector2d v = Eigen::Vector2d::Zero();
        v << vx, vy;
        s(i) = n*v;
    }
    s_tot = n*pos + s.maxCoeff();

    return s_tot;
}

tangent MPC_diffDrive_fblin::findTangentPoints_mod_convex_hull(Eigen::Vector2d robot_pos, Eigen::MatrixXd ped_pos, Eigen::MatrixXd v_ped, Eigen::MatrixXd relative_velocity, double r_c, double r_red, Eigen::VectorXi idx, int i_group, int Na, Eigen::VectorX<obstacle> obstacles, Eigen::Vector2d P_pos)
{
    double tol = 1e-8;

    long n_obs_group = (idx.array() == i_group).count();

    Eigen::MatrixXd ped_pos_cluster = Eigen::MatrixXd::Zero(n_obs_group, 2);
    Eigen::MatrixXd v_ped_cluster = Eigen::MatrixXd::Zero(n_obs_group, 2);
    Eigen::MatrixXd relative_velocity_cluster = Eigen::MatrixXd::Zero(n_obs_group, 2);
    Eigen::VectorX<obstacle> obstacles_cluster;
    obstacles_cluster.resize(n_obs_group);
    obstacle zero = {};
    obstacles_cluster.setConstant(zero);

    int x = 0;

    for (int i = 0; i < _n_obs; i++)
    {
        if (idx(i) == i_group)
        {
            ped_pos_cluster.row(x) = ped_pos.row(i);
            v_ped_cluster.row(x) = v_ped.row(i);
            relative_velocity_cluster.row(x) = relative_velocity.row(i);
            obstacles_cluster(x) = obstacles(i);
            x++;
        }
    }

    tangent tangent_points = {Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero()};

    Eigen::MatrixXd PSf_transpose_cluster = Eigen::MatrixXd::Zero(2, n_obs_group);
    Eigen::VectorXd theta_heading_cluster = Eigen::VectorXd::Zero(n_obs_group);

    std::vector<Eigen::Matrix2d> D;
    D.resize(n_obs_group, Eigen::Matrix2d::Zero());

    std::vector<Eigen::Matrix2d> B;
    B.resize(n_obs_group, Eigen::Matrix2d::Zero());

    Eigen::MatrixXd Q1_tot = Eigen::MatrixXd::Zero( 3, n_obs_group);
    Eigen::MatrixXd Q2_tot = Eigen::MatrixXd::Zero( 3, n_obs_group);
    Eigen::Matrix<bool, Eigen::Dynamic, 1> robot_inside_obst(n_obs_group);
    robot_inside_obst.setConstant(false);

    Eigen::VectorXd vector_cross = Eigen::VectorXd::Zero(2*n_obs_group);

    Eigen::VectorXd theta = Eigen::VectorXd::LinSpaced(Na+1, 0, 2*M_PI).head(Na).eval();

    Eigen::VectorXd n_1 = theta.array().cos().matrix();
    Eigen::VectorXd n_2 = theta.array().sin().matrix();

    Eigen::MatrixXd n = Eigen::MatrixXd::Zero(Na, 2);
    n.col(0) = n_1;
    n.col(1) = n_2;

    Eigen::MatrixXd s = Eigen::MatrixXd::Zero(n_obs_group, Na);

    double delta = (r_c - r_red);

    for (int i = 0; i < n_obs_group; i++)
    {
        tangent tangent_ponits_temp;

        if (obstacles_cluster(i).is_person)
        {
            personal_space personalSpace = personalSpaceFunction_mod(v_ped_cluster.row(i), relative_velocity_cluster.row(i).transpose(), robot_pos, ped_pos_cluster.row(i).transpose(), r_c);
            tangent_ponits_temp = findTangentPoints_mod(robot_pos, ped_pos_cluster.row(i).transpose(), personalSpace.PSf, personalSpace.theta_heading, r_red, r_c);
        }else
        {
            tangent_ponits_temp = findTangentPoints_mod_objects(robot_pos, r_red, r_c, obstacles_cluster(i));
        }
        Q1_tot.block(0,i,2,1) = tangent_ponits_temp.Q1;
        Q2_tot.block(0,i,2,1) = tangent_ponits_temp.Q2;
        robot_inside_obst(i) = tangent_ponits_temp.robot_inside_obst;
    }

    for (int i = 0; i < n_obs_group; i++)
    {
        if (obstacles_cluster(i).is_person)
        {
            personal_space personalSpace = personalSpaceFunction_mod(v_ped_cluster.row(i), relative_velocity_cluster.row(i).transpose(), robot_pos, ped_pos_cluster.row(i).transpose(), r_c);
            PSf_transpose_cluster.col(i) = personalSpace.PSf;
            theta_heading_cluster(i) = personalSpace.theta_heading;

            double a = PSf_transpose_cluster(1, i);// + delta;
            double b = PSf_transpose_cluster(0, i);// + delta;

            D[i] << a, 0.0,
                    0.0, b;

            Eigen::Matrix2d R;
            R << std::cos(theta_heading_cluster(i)), -std::sin(theta_heading_cluster(i)),
                 std::sin(theta_heading_cluster(i)), std::cos(theta_heading_cluster(i));

            B[i] = R*D[i];

            Eigen::RowVector2d v;
            v << std::cos(theta_heading_cluster(i)), std::sin(theta_heading_cluster(i));

            for (int j = 0; j < n.rows(); j++)
            {
                if (n.row(j).dot(v) >= 0.0)
                {
                    s(i,j) = n.row(j)*ped_pos_cluster.row(i).transpose() + (B[i].transpose()*n.row(j).transpose()).norm();
                }else
                {
                    s(i,j) = n.row(j)*ped_pos_cluster.row(i).transpose() + b;
                }

            }
        }else
        {
            for (int j = 0; j < n.rows(); j++)
            {
                    s(i,j) = compute_s_objects(n.row(j), r_red, r_red, obstacles_cluster(i), ped_pos_cluster.row(i).transpose());
            }
        }

    }

    Eigen::VectorXd s_convex_hull = s.colwise().maxCoeff().transpose();
    Eigen::VectorXd dist = s_convex_hull - n*P_pos;

    bool isInsideGroup = (dist.array() >= -tol).all();

    Eigen::Vector3d robot_pos_w_0;
    robot_pos_w_0 << robot_pos, 0.0;
    Eigen::MatrixXd vect_Q_tot(3, 2*n_obs_group);
    vect_Q_tot << Q1_tot - robot_pos_w_0.replicate(1, n_obs_group), Q2_tot - robot_pos_w_0.replicate(1, n_obs_group);

    if (isInsideGroup)
    {
        tangent_points.robot_inside_obst = true;
    }else
    {
        tangent_points.robot_inside_obst = false;
        for (int k = 0; k < 2*n_obs_group; k++)
        {
            vector_cross(k) = std::atan2(vect_Q_tot(0, 0)*vect_Q_tot(1, k) - vect_Q_tot(1, 0)*vect_Q_tot(0, k), vect_Q_tot(0, 0)*vect_Q_tot(0, k) + vect_Q_tot(1, 0)*vect_Q_tot(1, k));
        }

        int I_min, I_max;
        vector_cross.maxCoeff(&I_max);
        vector_cross.minCoeff(&I_min);

        tangent_points.Q1 = vect_Q_tot.block(0, I_min, 2, 1) + robot_pos;
        tangent_points.Q2 = vect_Q_tot.block(0, I_max, 2, 1) + robot_pos;
    }

    return tangent_points;
}

double MPC_diffDrive_fblin::compute_m_objects(const Eigen::RowVector2d& n, double r_red, double r_c, const obstacle& obstacle, const Eigen::Vector2d& pos)
{
    double delta = (r_c - r_red);

    std::vector<point> polygon = obstacle.polygon;

    Eigen::MatrixXd polygon_enlarged = Eigen::MatrixXd::Zero(static_cast<int>(polygon.size()), 2);
    Eigen::VectorXd m = Eigen::VectorXd::Zero(static_cast<int>(polygon.size()));
    double m_tot = 0.0;

    auto sign_nz = [](double x){ return (x >= 0.0) ? 1.0 : -1.0; };

    for (int i = 0; i < polygon.size(); ++i)
    {
        const auto& p = polygon[i];
        double vx;
        double vy;
        if (!obstacle.is_segment){
            double v_norm = std::hypot(p.x - obstacle.x, p.y - obstacle.y) + 1e-6;
            vx = (p.x - obstacle.x)/v_norm;
            vy = (p.y - obstacle.y)/v_norm;
            polygon_enlarged(i, 0) = p.x + vx*delta;
            polygon_enlarged(i, 1) = p.y + vy*delta;
        }else{
            Eigen::Vector2d normal = Eigen::Vector2d::Zero();
            normal << (obstacle.y2 - obstacle.y1), -(obstacle.x2 - obstacle.x1);
            Eigen::Vector2d longitudinal = Eigen::Vector2d::Zero();
            longitudinal << (obstacle.x2 - obstacle.x1), (obstacle.y2 - obstacle.y1);
            double normal_norm = normal.norm() + 1e-6;
            double longitudinal_norm = longitudinal.norm() + 1e-6;

            Eigen::Vector2d normal_normalized = normal/normal_norm;
            Eigen::Vector2d longitudinal_normalized = longitudinal/longitudinal_norm;

            Eigen::Vector2d delta_point;
            delta_point << (p.x - obstacle.x), (p.y - obstacle.y);

            polygon_enlarged(i, 0) = p.x + sign_nz(delta_point.dot(normal))*normal_normalized(0)*delta;
            polygon_enlarged(i, 1) = p.y + sign_nz(delta_point.dot(normal))*normal_normalized(1)*delta;
            polygon_enlarged(i, 0) = polygon_enlarged(i, 0) + sign_nz(delta_point.dot(longitudinal))*longitudinal_normalized(0)*delta;
            polygon_enlarged(i, 1) = polygon_enlarged(i, 1) + sign_nz(delta_point.dot(longitudinal))*longitudinal_normalized(1)*delta;
        }

    }

    for (int i = 0; i < polygon.size(); ++i)
    {
        double vx = (polygon_enlarged(i, 0) - obstacle.x);
        double vy = (polygon_enlarged(i, 1) - obstacle.y);
        Eigen::Vector2d v = Eigen::Vector2d::Zero();
        v << vx, vy;
        m(i) = n*v;
    }
    m_tot = n*pos + m.minCoeff();

    return m_tot;
}

constraint MPC_diffDrive_fblin::constraintCoefficients_mod_mod_convex_hull(Eigen::Vector2d robot_pos, Eigen::MatrixXd ped_pos, Eigen::MatrixXd v_ped, Eigen::MatrixXd relative_velocity, double r_c, double r_red, Eigen::VectorXi idx, int i_group, Eigen::VectorX<obstacle> obstacles, int Na)
{
    double tol = 1e-8;

    constraint constraintCoefficients = {0.0, 0.0, 0.0, 0.0};

    long n_obs_group = (idx.array() == i_group).count();

    Eigen::MatrixXd ped_pos_cluster = Eigen::MatrixXd::Zero(n_obs_group, 2);
    Eigen::MatrixXd v_ped_cluster = Eigen::MatrixXd::Zero(n_obs_group, 2);
    Eigen::MatrixXd relative_velocity_cluster = Eigen::MatrixXd::Zero(n_obs_group, 2);
    Eigen::VectorX<obstacle> obstacles_cluster;
    obstacles_cluster.resize(n_obs_group);
    obstacle zero = {};
    obstacles_cluster.setConstant(zero);

    int x = 0;

    for (int i = 0; i < _n_obs; i++)
    {
        if (idx(i) == i_group)
        {
            ped_pos_cluster.row(x) = ped_pos.row(i);
            v_ped_cluster.row(x) = v_ped.row(i);
            relative_velocity_cluster.row(x) = relative_velocity.row(i);
            obstacles_cluster(x) = obstacles(i);
            x++;
        }
    }

    Eigen::MatrixXd PSf_transpose_cluster = Eigen::MatrixXd::Zero(2, n_obs_group);
    Eigen::VectorXd theta_heading_cluster = Eigen::VectorXd::Zero(n_obs_group);

    std::vector<Eigen::Matrix2d> D;
    D.resize(n_obs_group, Eigen::Matrix2d::Zero());

    std::vector<Eigen::Matrix2d> B;
    B.resize(n_obs_group, Eigen::Matrix2d::Zero());

    Eigen::VectorXd m = Eigen::VectorXd::Zero(n_obs_group);

    Eigen::RowVectorXd coefficients = Eigen::RowVectorXd::Zero(3);

    Eigen::VectorXd hx_tot_cluster = Eigen::VectorXd::Zero(n_obs_group);

    Eigen::VectorXd hy_tot_cluster = Eigen::VectorXd::Zero(n_obs_group);

    Eigen::VectorXd l_tot_cluster = Eigen::VectorXd::Zero(n_obs_group);

    Eigen::VectorXd dist_tot_cluster = Eigen::VectorXd::Zero(n_obs_group);

    Eigen::VectorXd theta = Eigen::VectorXd::LinSpaced(Na+1, 0, 2*M_PI).head(Na).eval();

    Eigen::VectorXd n_1 = theta.array().cos().matrix();
    Eigen::VectorXd n_2 = theta.array().sin().matrix();

    Eigen::MatrixXd n = Eigen::MatrixXd::Zero(Na, 2);
    n.col(0) = n_1;
    n.col(1) = n_2;

    double delta = (r_c - r_red);

    double x_centroid = ped_pos_cluster.col(0).mean();
    double y_centroid = ped_pos_cluster.col(1).mean();

    double hx_best = 0.0; 
    double hy_best = 0.0;
    double l_best = 0.0;
    double dist_best = 0.0;
    double d_min_best = 0.0;


    for (int j = 0; j < Na; j++){

        double d_min_temp = 0.0;
        double dist_temp = 0.0;

    coefficients(0) = - n(j,0);//(x_centroid - robot_pos(0))/(std::hypot(x_centroid - robot_pos(0), y_centroid - robot_pos(1)) + 1e-6);
    coefficients(1) = - n(j,1);//(y_centroid - robot_pos(1))/(std::hypot(x_centroid - robot_pos(0), y_centroid - robot_pos(1)) + 1e-6);

    coefficients(2) = coefficients(0)*x_centroid + coefficients(1)*y_centroid;

    /*std::vector<ClosestOnPolygonResult> p_near;

    for (int i = 0; i < n_obs_group; i++)
    {
        std::vector<Eigen::Vector2d> poly;
        for (int j = 0; j < obstacles_cluster(i).polygon.size(); j++)
        {
            Eigen::Vector2d v;
            v << obstacles_cluster(i).polygon[j].x - obstacles_cluster(i).x + ped_pos_cluster(i, 0), obstacles_cluster(i).polygon[j].y - obstacles_cluster(i).y + ped_pos_cluster(i, 1);
            poly.push_back(v);
        }
        p_near.push_back(closestPointOnPolygonBoundary(robot_pos, poly));
    }

    int best_idx = -1;
    double best_dist2 = std::numeric_limits<double>::infinity();

    for (int i = 0; i < static_cast<int>(p_near.size()); ++i)
    {
    if (p_near[i].dist2 < best_dist2)
    {
        best_dist2 = p_near[i].dist2;
        best_idx = i;
    }
    }

// best_idx = indice del minimo, oppure -1 se p_near e' vuoto

    double x_centroid = p_near[best_idx].closest_point(0);
    double y_centroid = p_near[best_idx].closest_point(1);

    coefficients(0) = (x_centroid - robot_pos(0))/(std::hypot(x_centroid - robot_pos(0), y_centroid - robot_pos(1)) + 1e-6);
    coefficients(1) = (y_centroid - robot_pos(1))/(std::hypot(x_centroid - robot_pos(0), y_centroid - robot_pos(1)) + 1e-6);

    coefficients(2) = coefficients(0)*x_centroid + coefficients(1)*y_centroid;*/


    /*for (int i = 0; i < n_obs_group; i++)
    {
    double x_centroid = ped_pos_cluster(i, 0);
    double y_centroid = ped_pos_cluster(i, 1);

    hx_tot_cluster(i) = -(x_centroid - robot_pos(0))/(std::hypot(x_centroid - robot_pos(0), y_centroid - robot_pos(1)) + 1e-6);
    hy_tot_cluster(i) = -(y_centroid - robot_pos(1))/(std::hypot(x_centroid - robot_pos(0), y_centroid - robot_pos(1)) + 1e-6);

    l_tot_cluster(i) = hx_tot_cluster(i)*x_centroid + hy_tot_cluster(i)*y_centroid;
    }

    dist_tot_cluster = hx_tot_cluster*robot_pos(0) + hy_tot_cluster*robot_pos(1) - l_tot_cluster;

    int I_coeff;
    dist_tot_cluster.minCoeff(&I_coeff);

    coefficients(0) = - hx_tot_cluster(I_coeff);
    coefficients(1) = - hy_tot_cluster(I_coeff);
    coefficients(2) = - l_tot_cluster(I_coeff);*/


    for (int i = 0; i < n_obs_group; i++)
    {
        if (obstacles_cluster(i).is_person)
        {
            personal_space personalSpace = personalSpaceFunction_mod(v_ped_cluster.row(i), relative_velocity_cluster.row(i).transpose(), robot_pos, ped_pos_cluster.row(i).transpose(), r_c);
            PSf_transpose_cluster.col(i) = personalSpace.PSf;
            theta_heading_cluster(i) = personalSpace.theta_heading;

            double a = PSf_transpose_cluster(1, i);// + delta;
            double b = PSf_transpose_cluster(0, i);// + delta;

            D[i] << a, 0.0,
                    0.0, b;

            Eigen::Matrix2d R;
            R << std::cos(theta_heading_cluster(i)), -std::sin(theta_heading_cluster(i)),
                 std::sin(theta_heading_cluster(i)), std::cos(theta_heading_cluster(i));

            B[i] = R*D[i];

            Eigen::RowVector2d v;
            v << std::cos(theta_heading_cluster(i)), std::sin(theta_heading_cluster(i));

            if (coefficients.segment(0, 2).dot(v) >= 0.0)
            {
                m(i) = coefficients.segment(0, 2)*ped_pos_cluster.row(i).transpose() - coefficients(2) - b;
            }else
            {
                m(i) = coefficients.segment(0, 2)*ped_pos_cluster.row(i).transpose() - coefficients(2) - (B[i].transpose()*coefficients.segment(0, 2).transpose()).norm();
            }
        }else
        {
            m(i) = compute_m_objects(coefficients.segment(0, 2), r_red, r_c, obstacles_cluster(i), ped_pos_cluster.row(i).transpose()) - coefficients(2);
        }


    }

    bool isHullTangentOriginal = (m.array() >= -tol).all();

    int I_m;
    m.minCoeff(&I_m);

    Eigen::RowVector2d h;
    h << coefficients(0), coefficients(1);

    if (isHullTangentOriginal)
    {
        d_min_temp = h*(ped_pos_cluster.row(I_m).transpose()) - coefficients(2) - r_c;
    }else
    {
        coefficients(2) =  coefficients(2) + m(I_m);
        if (obstacles_cluster(I_m).is_person)
        {
        d_min_temp = h*(ped_pos_cluster.row(I_m).transpose()) - coefficients(2) - r_c;
        }else{
            d_min_temp = delta;//- m(I_m) + compute_m_objects(coefficients.segment(0, 2), r_red, r_red, obstacles_cluster(I_m), ped_pos_cluster.row(I_m).transpose()) - coefficients(2);
        }
    }

    dist_temp = coefficients(0)*robot_pos(0) + coefficients(1)*robot_pos(1) - coefficients(2);

    if (dist_temp < dist_best || j == 0){
        hx_best = coefficients(0);
        hy_best = coefficients(1);
        l_best = coefficients(2);
        d_min_best = d_min_temp;
        dist_best = dist_temp;
    }
}

constraintCoefficients.hx = hx_best;
constraintCoefficients.hy = hy_best;
constraintCoefficients.l = l_best;
constraintCoefficients.d_min = d_min_best;


    return constraintCoefficients;
}

double MPC_diffDrive_fblin::safetyDistance_mod_convex_hull(double r_red, double Ts, std::vector<Eigen::MatrixXd> v_ped_samples, std::vector<Eigen::MatrixXd> p_ped_samples, Eigen::MatrixXd p_robot_samples, Eigen::MatrixXd v_robot_samples, int i, double r_c, double epsilon, Eigen::VectorXi idx, int i_group, Eigen::MatrixX<obstacle> obstacles_samples, int Na, Eigen::MatrixXd p_P_robot_samples, Eigen::VectorXd yaw_robot, Eigen::MatrixXd optimVect_no_sl_samples)
{
    long N_samples = p_ped_samples.front().rows();

    long N = p_ped_samples.front().rows() - 20;

    double dl_sl = 0.0;

    std::vector<Eigen::MatrixXd> p_ped_pred;
    p_ped_pred.resize(p_ped_samples.size(), Eigen::MatrixXd::Zero(N, 2));
    std::vector<Eigen::MatrixXd> v_ped_pred;
    v_ped_pred.resize(p_ped_samples.size(), Eigen::MatrixXd::Zero(N, 2));
    Eigen::MatrixX<obstacle> obstacles_pred;
    obstacles_pred.resize(N, static_cast<int>(p_ped_samples.size()));
    obstacle zero = {};
    obstacles_pred.setConstant(zero);
    Eigen::MatrixXd p_robot_pred = Eigen::MatrixXd::Zero(N, 2);
    Eigen::MatrixXd v_robot_pred = Eigen::MatrixXd::Zero(N, 2);
    Eigen::MatrixXd p_P_robot_pred = Eigen::MatrixXd::Zero(N, 2);

    for (int x = 0; x < p_ped_samples.size(); x++)
    {
        for (int j = 0; j < N; j++)
        {
            p_ped_pred[x](j,0) = p_ped_samples[x](N_samples-20+i-j-1,0);
            p_ped_pred[x](j,1) = p_ped_samples[x](N_samples-20+i-j-1,1);
            v_ped_pred[x](j,0) = v_ped_samples[x](N_samples-20+i-j-1,0);
            v_ped_pred[x](j,1) = v_ped_samples[x](N_samples-20+i-j-1,1);
            for (int k = 0; k < i; k++)
            {
                p_ped_pred[x](j,0) = p_ped_pred[x](j,0)+Ts*v_ped_pred[x](j,0);
                p_ped_pred[x](j,1) = p_ped_pred[x](j,1)+Ts*v_ped_pred[x](j,1);
                v_ped_pred[x](j,0) = v_ped_pred[x](j,0);
                v_ped_pred[x](j,1) = v_ped_pred[x](j,1);
            }
        }

    }

    obstacles_pred = obstacles_samples.block(i, 0, N, static_cast<int>(p_ped_samples.size()));

    for (int j = 0; j < N; j++)
    {
        pred_step prediction = prediction_step(p_P_robot_samples(N_samples-20+i-j-1,0), p_P_robot_samples(N_samples-20+i-j-1,1), p_robot_samples(N_samples-20+i-j-2,0), p_robot_samples(N_samples-20+i-j-2,1), yaw_robot(N_samples-20+i-j-2), p_P_robot_samples(N_samples-20+i-j-2,0), p_P_robot_samples(N_samples-20+i-j-2,1), optimVect_no_sl_samples.row(N_samples-20+i-j-1).transpose(), i);
        p_robot_pred(j,0) = prediction.x_pred;
        p_robot_pred(j,1) = prediction.y_pred;
        v_robot_pred(j,0) = prediction.vx_pred;
        v_robot_pred(j,1) = prediction.vy_pred;
        p_P_robot_pred(j,0) = prediction.XP_pred;
        p_P_robot_pred(j,1) = prediction.YP_pred;
    }

    Eigen::VectorXd hx = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd hy = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd l = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd hx_pred = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd hy_pred = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd l_pred = Eigen::VectorXd::Zero(N);

    for (int q = 0; q < N; q++)
    {
        Eigen::MatrixXd ped_pos(2, p_ped_pred.size());
        Eigen::MatrixXd v_ped(2, p_ped_pred.size());
        Eigen::MatrixXd relative_velocity(p_ped_pred.size(), 2);

        for (int z = 0; z < p_ped_pred.size(); z++)
        {
            ped_pos.block(0, z, 2, 1) = p_ped_pred[z].row(q).transpose();
            v_ped.block(0, z, 2, 1) = v_ped_pred[z].row(q).transpose();
        }

        relative_velocity = v_robot_pred.row(q).replicate(static_cast<long>(p_ped_pred.size()), 1) - v_ped.transpose();

        constraint constraintCoefficients_pred = constraintCoefficients_mod_mod_convex_hull(p_robot_pred.row(q).transpose(), ped_pos.transpose(), v_ped.transpose(), relative_velocity, r_c, r_red, idx, i_group, obstacles_pred.row((N-q)-1).transpose(), Na);

        hx_pred(q) = constraintCoefficients_pred.hx;
        hy_pred(q) = constraintCoefficients_pred.hy;
        l_pred(q) = constraintCoefficients_pred.l;

        Eigen::MatrixXd ped_pos_samples(2, p_ped_samples.size());
        Eigen::MatrixXd v_ped_samples_tot(2, p_ped_samples.size());
        Eigen::MatrixXd relative_velocity_samples(p_ped_samples.size(), 2);

        for (int z = 0; z < p_ped_pred.size(); z++)
        {
            ped_pos_samples.block(0, z, 2, 1) = p_ped_samples[z].row((N-q)-1).transpose();
            v_ped_samples_tot.block(0, z, 2, 1) = v_ped_samples[z].row((N-q)-1).transpose();
        }

        relative_velocity_samples = v_robot_samples.row((N-q)-1).replicate(static_cast<long>(p_ped_pred.size()), 1) - v_ped_samples_tot.transpose();

        constraint constraintCoefficients = constraintCoefficients_mod_mod_convex_hull(p_robot_samples.row((N-q)-1).transpose(), ped_pos_samples.transpose(), v_ped_samples_tot.transpose(), relative_velocity_samples, r_c, r_red, idx, i_group, obstacles_samples.row((N-q)-1).transpose(), Na);

        hx(q) = constraintCoefficients.hx;
        hy(q) = constraintCoefficients.hy;
        l(q) = constraintCoefficients.l;
    }

    Eigen::VectorXd Delta_l = l - l_pred;

    double mu_Delta_l = Delta_l.mean();
    double Var_Delta_l =  ((Delta_l.array() - mu_Delta_l).square().sum())/static_cast<double>(N - 1);

    Eigen::VectorXd Delta_hx = hx.array()*p_P_robot_samples.col(0).segment(0, N).reverse().array() - hx_pred.array()*p_P_robot_pred.col(0).array();//*p_robot_samples.col(0).segment(0, N).reverse().array();

    double mu_Delta_hx = Delta_hx.mean();
    double Var_Delta_hx = ((Delta_hx.array() - mu_Delta_hx).square().sum())/static_cast<double>(N - 1);

    Eigen::VectorXd Delta_hy = hy.array()*p_P_robot_samples.col(1).segment(0, N).reverse().array() - hy_pred.array()*p_P_robot_pred.col(1).array();//*p_robot_samples.col(1).segment(0, N).reverse().array();

    double mu_Delta_hy = Delta_hy.mean();
    double Var_Delta_hy = ((Delta_hy.array() - mu_Delta_hy).square().sum())/static_cast<double>(N - 1);

    double Cov_Delta_l_Delta_hx = ((Delta_l.array() - mu_Delta_l).array()*(Delta_hx.array() - mu_Delta_hx).array()).sum()/static_cast<double>(N - 1);

    double Cov_Delta_l_Delta_hy = ((Delta_l.array() - mu_Delta_l).array()*(Delta_hy.array() - mu_Delta_hy).array()).sum()/static_cast<double>(N - 1);

    double Cov_Delta_hx_Delta_hy = ((Delta_hx.array() - mu_Delta_hx).array()*(Delta_hy.array() - mu_Delta_hy).array()).sum()/static_cast<double>(N - 1);

    double mu_tot = mu_Delta_l - mu_Delta_hx - mu_Delta_hy;
    double Var_tot = Var_Delta_l + Var_Delta_hx + Var_Delta_hy - 2*Cov_Delta_l_Delta_hx - 2*Cov_Delta_l_Delta_hy + 2*Cov_Delta_hx_Delta_hy;

    dl_sl = - mu_tot + std::sqrt(std::max(0.0, Var_tot))*std::sqrt((1-epsilon)/epsilon);

    return dl_sl;
}


void MPC_diffDrive_fblin::compute_obstacle_avoidance_constraints()
{
    std::vector<Eigen::MatrixXd> H_obs;
    H_obs.resize(_n_groups, Eigen::MatrixXd::Zero(_N, 2*_N));

    Eigen::MatrixXd L_obs = Eigen::MatrixXd::Zero(_N, _n_groups);
    Eigen::MatrixXd dl = Eigen::MatrixXd::Zero(_N, _n_groups);
    Eigen::MatrixXd dl_min = Eigen::MatrixXd::Zero(_N, _n_groups);

    _A_obs = Eigen::MatrixXd::Zero(_n_groups*_N, 2*_N + _n_groups + 2);
    _b_obs = Eigen::VectorXd::Zero(_n_groups*_N);

    Eigen::MatrixXd x_ped = Eigen::MatrixXd::Zero(_n_obs, _N);
    Eigen::MatrixXd y_ped = Eigen::MatrixXd::Zero(_n_obs, _N);

    x_ped.col(0) = _x_ped_samples.row(0).transpose();
    y_ped.col(0) = _y_ped_samples.row(0).transpose();

    for (int i = 1; i < _N; i++)
    {
        x_ped.col(i) = x_ped.col(i-1) + _MPC_Ts*_vx_ped_samples.row(0).transpose();
        y_ped.col(i) = y_ped.col(i-1) + _MPC_Ts*_vy_ped_samples.row(0).transpose();
    }

    int Na = 200;

    Eigen::MatrixXd p_robot_samples(_N_samples, 2);
    p_robot_samples << _x_robot_samples, _y_robot_samples;

    Eigen::MatrixXd v_robot_samples(_N_samples, 2);
    v_robot_samples << _vx_robot_samples, _vy_robot_samples;

    Eigen::MatrixXd p_P_robot_samples(_N_samples, 2);
    p_P_robot_samples << _xP_robot_samples, _yP_robot_samples;

    std::vector<Eigen::MatrixXd> p_ped_samples;
    p_ped_samples.resize(_n_obs, Eigen::MatrixXd::Zero(_N_samples, 2));

    std::vector<Eigen::MatrixXd> v_ped_samples;
    v_ped_samples.resize(_n_obs, Eigen::MatrixXd::Zero(_N_samples, 2));

    Eigen::VectorX<obstacle> obstacles = _obstacles_samples.row(0).transpose();

    for (int i = 0; i < _n_obs; i++)
    {
        p_ped_samples[i] << _x_ped_samples.col(i), _y_ped_samples.col(i);
        v_ped_samples[i] << _vx_ped_samples.col(i), _vy_ped_samples.col(i);
    }

      std::vector<ConstraintLine> lines;

    for (int i = 0; i < _n_groups; i++)
    {
        Eigen::MatrixXd E_sl = Eigen::MatrixXd::Zero(_N, _n_groups);

        int j = 0;
        if (std::hypot(_x_near(i) - _x_robot_samples(0), _y_near(i) - _y_robot_samples(0)) <= _sensorRange)
        {
            for (int k = 0; k < _N; k++)
            {
                Eigen::MatrixXd relative_velocity = Eigen::MatrixXd::Zero(_n_obs, 2);
                Eigen::RowVectorXd relative_velocity_sum = Eigen::RowVectorXd::Zero(2);
                Eigen::RowVectorXd relative_velocity_cluster = Eigen::RowVectorXd::Zero(2);

                double r_c;

                if (_v_long(k) >= 0.0){
                    r_c = 0.8;
                }else{
                    r_c = 1.3;
                }

                relative_velocity.col(0) = _v_long(k)*cos(_predictRobotState(3*k+2))*Eigen::VectorXd::Ones(_n_obs) - _vx_ped_samples.row(0).transpose();
                relative_velocity.col(1) = _v_long(k)*sin(_predictRobotState(3*k+2))*Eigen::VectorXd::Ones(_n_obs) - _vy_ped_samples.row(0).transpose();

                int x = 0;

                for (int l = 0; l < _n_obs; l++)
                {
                    if (_labels(l) == i)
                    {
                        relative_velocity_sum += relative_velocity.row(l);
                        x++;
                    }
                }

                if (x > 0) {
                    relative_velocity_cluster(0) = relative_velocity_sum(0)/x;
                    relative_velocity_cluster(1) = relative_velocity_sum(1)/x;
                }else
                {
                    relative_velocity_cluster(0) = 0.0;
                    relative_velocity_cluster(1) = 0.0;
                }

                Eigen::Vector2d robot_pos = Eigen::Vector2d::Zero();
                robot_pos(0) = _predictRobotState(3*k);
                robot_pos(1) = _predictRobotState(3*k+1);

                Eigen::Vector2d P_pos = Eigen::Vector2d::Zero();
                P_pos(0) = _XP_predicted(k);
                P_pos(1) = _YP_predicted(k);

                Eigen::MatrixXd ped_pos = Eigen::MatrixXd::Zero(_n_obs, 2);
                ped_pos.col(0) = x_ped.col(k);
                ped_pos.col(1) = y_ped.col(k);

                Eigen::MatrixXd v_ped(_n_obs, 2);
                v_ped << _vx_ped_samples.row(0).transpose(), _vy_ped_samples.row(0).transpose();

                tangent tangent_points = findTangentPoints_mod_convex_hull(robot_pos, ped_pos, v_ped, relative_velocity, r_c, _r_red, _labels, i, Na, obstacles, P_pos);

                Eigen::Vector2d v1;
                v1 << tangent_points.Q1(0) - robot_pos(0), tangent_points.Q1(1) - robot_pos(1);

                Eigen::Vector2d v2;
                v2 << tangent_points.Q2(0) - robot_pos(0), tangent_points.Q2(1) - robot_pos(1);

                Eigen::Matrix2d V_12;
                V_12 << v1, v2;

                Eigen::Vector2d alfa_12;

                if (!tangent_points.robot_inside_obst)
                {
                    alfa_12 = V_12.colPivHouseholderQr().solve(relative_velocity_cluster.transpose());
                }else
                {
                    alfa_12 << 0.0, 0.0;
                }

                if (k==0){

                    _vx_robot = _v_long(k)*cos(_predictRobotState(3*k+2));
                    _vy_robot = _v_long(k)*sin(_predictRobotState(3*k+2));
                    _v_long_robot = _v_long(k);

                }

                if (tangent_points.robot_inside_obst || (alfa_12(0) > 0.0 && alfa_12(1) > 0.0))
                {
                    constraint constraintCoefficients = constraintCoefficients_mod_mod_convex_hull(robot_pos, ped_pos, v_ped, relative_velocity, r_c, _r_red, _labels, i, obstacles, Na);

                    H_obs[i].block(k, j, 1, 2) << constraintCoefficients.hx, constraintCoefficients.hy;

                    L_obs(k, i) = constraintCoefficients.l;

                    if (k==0){

                    lines.push_back({constraintCoefficients.hx, constraintCoefficients.hy, constraintCoefficients.l, i, k});

                    }

                    double dl_sl = safetyDistance_mod_convex_hull(_r_red, _MPC_Ts, v_ped_samples, p_ped_samples, p_robot_samples, v_robot_samples, k + 1, r_c, _e, _labels, i, _obstacles_samples, Na, p_P_robot_samples, _theta_robot_samples, _optimVect_no_sl_samples);
                    dl(k,i) = dl_sl;
                    dl_min(k,i) = constraintCoefficients.d_min;
                    debugMsg(std::to_string(dl_sl));
                    //debugMsg(std::to_string(dl_min(k,i)));

                    if (k==0){

                    lines.push_back({constraintCoefficients.hx, constraintCoefficients.hy, constraintCoefficients.l - dl(k,i), i, k});

                    }
                }

                E_sl(k,i)=-dl(k,i)-dl_min(k,i);

                j=j+2;
            }

        }

        Eigen::MatrixXd A_obsj(_N, 2*_N + _n_groups + 2);
        A_obsj << H_obs[i]*_Bcal, E_sl, Eigen::MatrixXd::Zero(_N,2);
        Eigen::VectorXd b_obsj(_N);
        Eigen::Vector2d csi;
        csi << _actXP, _actYP;
        b_obsj = L_obs.col(i) - dl.col(i) - H_obs[i]*_Acal*csi;

        _A_obs.block((i)*_N, 0, _N, 2*_N + _n_groups +2) = A_obsj;
        _b_obs.segment((i)*_N, _N) = b_obsj;

    }

    if (_constraints_cb) {
    _constraints_cb(lines);
  }

}

void MPC_diffDrive_fblin::compute_slack_variables_constraints()
{
    Eigen::MatrixXd Asl_a = Eigen::MatrixXd::Zero(4,2*_N + _n_groups + 2);
    Eigen::Vector4d bsl_a;
    bsl_a << _dsl_a, _dsl_a, 0.0, 0.0;

    Asl_a.block(0, 2*_N + _n_groups, 2, 2) = Eigen::MatrixXd::Identity(2,2);
    Asl_a.block(2, 2*_N + _n_groups, 2, 2) = - Eigen::MatrixXd::Identity(2,2);

    _Asl = Eigen::MatrixXd::Zero(2*_n_groups + 4,2*_N + _n_groups + 2);
    _bsl = Eigen::VectorXd::Zero(2*_n_groups + 4);

    if (_n_groups > 0)
    {
        Eigen::MatrixXd Asl_p = Eigen::MatrixXd::Zero(2*_n_groups,2*_N + _n_groups + 2);
        Eigen::VectorXd bsl_p = Eigen::VectorXd::Zero(2*_n_groups);

        Asl_p.block(0, 2*_N, _n_groups, _n_groups) = Eigen::MatrixXd::Identity(_n_groups,_n_groups);
        Asl_p.block(_n_groups, 2*_N, _n_groups, _n_groups) = - Eigen::MatrixXd::Identity(_n_groups,_n_groups);
        bsl_p.segment(0, _n_groups) = Eigen::VectorXd::Ones(_n_groups);
        _Asl << Asl_p, Asl_a;
        _bsl << bsl_p, bsl_a;
    }else
    {
        _Asl = Asl_a;
        _bsl = bsl_a;
    }



}

void MPC_diffDrive_fblin::extractTraj(Eigen::MatrixXd ref, int count, bool trajectory_flag, bool replanTrigger)
{

    _refMPCstate = Eigen::VectorXd::Zero(2*(_N));
    Eigen::VectorXd x_pr = ref.col(0);
    Eigen::VectorXd y_pr = ref.col(1);

    Eigen::Vector2d v;
    v << x_pr(x_pr.size() -1), y_pr(y_pr.size() -1);

    if (!trajectory_flag)
    {
        _refMPCstate = v.replicate(_N, 1);
    }else
    {
        if (count == 0 || replanTrigger)
        {
            _kt = 1;
        }

        if (!_is_blocked_inside_hull){

        _refMPCstate = v.replicate(_N, 1);

        Eigen::VectorXd Csi_r_temp_x = x_pr.segment(std::min(static_cast<long>(_kt),x_pr.size() - 1), std::min(static_cast<long>(_kt+_N),x_pr.size()) - std::min(static_cast<long>(_kt+1),x_pr.size()) + 1);
        Eigen::VectorXd Csi_r_temp_y = y_pr.segment(std::min(static_cast<long>(_kt),y_pr.size() - 1), std::min(static_cast<long>(_kt+_N),y_pr.size()) - std::min(static_cast<long>(_kt+1),y_pr.size()) + 1);

        int j = 0;
        for (int i=0; i<std::min(static_cast<long>(2*_N), 2*(x_pr.size() - _kt)); i += 2)
        {
            _refMPCstate(i) = Csi_r_temp_x(j);
            _refMPCstate(i+1) = Csi_r_temp_y(j);
            if (j < Csi_r_temp_x.size() - 1)
            {
                j++;
            }
        }
    }else{
        _refMPCstate = _ref_temp.replicate(_N, 1);
    }

        _kt = _kt + 1;
    }

}

bool MPC_diffDrive_fblin::executeMPCcontroller() {
    /** Preliminary checks */
    if (!_controllerInitialized)
    {
        errorMsg("[MPC_diffDrive_fblin.executeMPCcontroller] Call initialize() before calling executeMPCcontroller()");
        return false;
    }

    prediction();

    // Compute constraint matrices
    compute_wheelVelocityConstraint();
    compute_wheelAccelerationConstraint();

    if (_n_groups > 0)
    {
        compute_obstacle_avoidance_constraints();
    }

    compute_slack_variables_constraints();

    _A_tot = Eigen::MatrixXd::Zero(8*_N + (_N + 2)*_n_groups + 4, 2*_N + _n_groups + 2);
    _B_tot = Eigen::VectorXd::Zero(8*_N + (_N + 2)*_n_groups + 4);

    if (_n_groups > 0)
    {
        _A_tot << _Ain_vel,
                  _A_dvel,
                  _A_obs,
                  _Asl;
        _B_tot << _Bin_vel,
                  _b_dvel,
                  _b_obs,
                  _bsl;
    }else
    {
        _A_tot << _Ain_vel,
                  _A_dvel,
                  _Asl;
        _B_tot << _Bin_vel,
                  _b_dvel,
                  _bsl;
    }

    if (totConstrain.size()==0) {
        if (!_solver->addConstraint(_A_tot, _B_tot, totConstrain))
        {
            errorMsg("[MPC_diffDrive_fblin.executeMPCcontroller] Error setting the wheel velocity constraint");
            return false;
        }
    }
    else {
        if (!_solver->addConstraint(_A_tot, _B_tot, totConstrain))
        {
            errorMsg("[MPC_diffDrive_fblin.executeMPCcontroller] Error setting the wheel velocity constraint");
            return false;
        }
    }

    // Compute cost function matrices
    compute_objectiveMatrix();
    if (!_solver->setObjective(_H, _f))
    {
        errorMsg("[MPC_diffDrive_fblin.executeMPCcontroller] Error setting the MPC objective");
        return false;
    }

    /*** ADD THIS BLOCK — QP DUMP FOR DEBUG ***/
    static int dump_id = 0;
    // (optional) only dump sometimes: if (dump_id % 10 == 0) { ... }
    char fname[128];
    std::snprintf(fname, sizeof(fname), "/home/alessandro/MPC_ros2/nav2_mpc_fblin_controller_0925/development/cmake-build-debug/test_MPC_Cpp/mpc_step_%04d", dump_id++);
    _solver->writeProblem(fname);   // writes /tmp/mpc_step_XXXX.lp (and .mps)
    /*** END ADDED BLOCK ***/

    // Solve optimization problem
    int optimizerStatus;
    double objectiveValue;
    if (!_solver->solveProblem(_optimVect, objectiveValue, optimizerStatus))
    {
        _optimVect.setZero();
        if (!_is_blocked){
        _n_fail++;
        }

        if (_n_fail > 5){
            _is_blocked = true;
            _n_fail = 50;
            _sp = 1.0e3;
        }


        switch (optimizerStatus) {
            case GUROBIsolver::INFEASIBLE:
                errorMsg("[MPC_diffDrive_fblin.executeMPCcontroller] Error solving the optimization problem (infeasible problem)");
                break;
            case GUROBIsolver::OTHER:
                errorMsg("[MPC_diffDrive_fblin.executeMPCcontroller] Error solving the optimization problem (not optimal but not infeasible)");
                break;
        }
        return false;
    }else{
        _n_fail = std::max(_n_fail - 1, 0);
        if (_n_fail == 0){
            _is_blocked = false;
            _sp = 1.0e9;
        }
    }

    return true;
}

bool MPC_diffDrive_fblin::executeLinearizationController() {
    /** Preliminary checks */
    if (!_controllerInitialized) {
        errorMsg("[MPC_diffDrive_fblin.executeLinearizationController] Call initialize() before calling executeLinearizationController()");
        return false;
    }

    // Execute the feedback linearization law
    _fblinController->control_transformation(_optimVect(0), _optimVect(1), _linearVelocity, _angularVelocity);

    return true;
}

void MPC_diffDrive_fblin::set_actualRobotState(double x, double y, double yaw) {
    _actX = x;
    _actY = y;
    _actYaw = yaw;

    // Update the linearizing controller state
    _fblinController->set_unicycleState(_actX, _actY, _actYaw);

    if (_init_P)
    {
        _prevXP = _actXP;
        _prevYP = _actYP;
    }

    // Update the MPC state
    _fblinController->ouput_transformation(_actXP, _actYP);

    if (!_init_P)
    {
        _prevXP = _actXP;
        _prevYP = _actYP;
        _init_P = true;
    }

}

void MPC_diffDrive_fblin::sample_Pedestrian_and_robot_states(const Eigen::VectorX<obstacle>& obstacles, double x_robot, double y_robot, double theta_robot)//, double sim_time, Eigen::MatrixXd u_prev, int count)
{
    if (!_is_started)
    {
        _x_ped_samples = Eigen::MatrixXd::Zero(_N_samples, _n_obs);
        _y_ped_samples = Eigen::MatrixXd::Zero(_N_samples, _n_obs);
        _vx_ped_samples = Eigen::MatrixXd::Zero(_N_samples, _n_obs);
        _vy_ped_samples = Eigen::MatrixXd::Zero(_N_samples, _n_obs);
        _x_robot_samples = Eigen::VectorXd::Zero(_N_samples);
        _y_robot_samples = Eigen::VectorXd::Zero(_N_samples);
        _theta_robot_samples = Eigen::VectorXd::Zero(_N_samples);
        _vx_robot_samples = Eigen::VectorXd::Zero(_N_samples);
        _vy_robot_samples = Eigen::VectorXd::Zero(_N_samples);
        _xP_robot_samples = Eigen::VectorXd::Zero(_N_samples);
        _yP_robot_samples = Eigen::VectorXd::Zero(_N_samples);
        _optimVect_no_sl_samples = Eigen::MatrixXd::Zero(_N_samples, 2*_N);
        _id_samples = Eigen::VectorXi::Zero(_n_obs);

        for (int i = 0; i < _n_obs; i++)
        {
            _x_ped_samples.col(i).setConstant(obstacles(i).x);
            _y_ped_samples.col(i).setConstant(obstacles(i).y);
            _vx_ped_samples.col(i).setConstant(0.0);
            _vy_ped_samples.col(i).setConstant(0.0);
            _vx_ped_samples(0, i) = obstacles(i).vx;
            _vy_ped_samples(0, i) = obstacles(i).vy;
            _id_samples(i) = obstacles(i).id;
        }

        _x_robot_samples.setConstant(x_robot);
        _y_robot_samples.setConstant(y_robot);
        _theta_robot_samples.setConstant(theta_robot);
        _xP_robot_samples.setConstant(x_robot + _Pdist*std::cos(theta_robot));
        _yP_robot_samples.setConstant(y_robot + _Pdist*std::sin(theta_robot));
        _vx_robot_samples.setConstant(0.0);
        _vy_robot_samples.setConstant(0.0);

        _is_started = true;
    }else
    {
        _optimVect_no_sl = _optimVect.segment(0, 2*_N); // u_prev.col(count-1);

        _optimVect_no_sl_samples.block(1, 0, _N_samples-1, 2*_N) = _optimVect_no_sl_samples.block(0, 0, _N_samples-1, 2*_N).eval();
        _optimVect_no_sl_samples.row(0) = _optimVect_no_sl.transpose();

        Eigen::MatrixXd _x_ped_samples_temp = Eigen::MatrixXd::Zero(_N_samples, _n_obs);
        Eigen::MatrixXd _y_ped_samples_temp = Eigen::MatrixXd::Zero(_N_samples, _n_obs);
        Eigen::MatrixXd _vx_ped_samples_temp = Eigen::MatrixXd::Zero(_N_samples, _n_obs);
        Eigen::MatrixXd _vy_ped_samples_temp = Eigen::MatrixXd::Zero(_N_samples, _n_obs);
        Eigen::VectorXi _id_samples_temp = Eigen::VectorXi::Zero(_n_obs);
        Eigen::VectorX<bool> _is_static_temp(_n_obs);
        _is_static_temp.setConstant(false);

        Eigen::VectorX<bool> is_marked_sample(_id_samples_obstacles.size());
        is_marked_sample.setConstant(false);

        for (int i = 0; i < _n_obs; i++)
        {
            bool is_memorized = false;
            for (int j = 0; j < _id_samples.size(); j++)
            {
                if (!is_marked_sample(j) && (obstacles(i).id == _id_samples(j)  || std::hypot(obstacles(i).x - (_x_ped_samples(0, j) + _vx_ped_samples(0, j)*_MPC_Ts), obstacles(i).y - (_y_ped_samples(0, j) + _vy_ped_samples(0, j)*_MPC_Ts)) < _r_red*2 + std::hypot(_vx_ped_samples(0, j)*_MPC_Ts, _vy_ped_samples(0, j)*_MPC_Ts)))
                {
                    _x_ped_samples_temp.col(i).segment(1, _N_samples-1) = _x_ped_samples.col(j).segment(0, _N_samples-1).eval();
                    _x_ped_samples_temp(0, i) = obstacles(i).x;
                    _y_ped_samples_temp.col(i).segment(1, _N_samples-1) = _y_ped_samples.col(j).segment(0, _N_samples-1).eval();
                    _y_ped_samples_temp(0, i) = obstacles(i).y;
                    _vx_ped_samples_temp.col(i).segment(1, _N_samples-1) = _vx_ped_samples.col(j).segment(0, _N_samples-1).eval();
                    _vx_ped_samples_temp(0, i) = obstacles(i).vx;
                    _vy_ped_samples_temp.col(i).segment(1, _N_samples-1) = _vy_ped_samples.col(j).segment(0, _N_samples-1).eval();
                    _vy_ped_samples_temp(0, i) = obstacles(i).vy;
                    is_memorized = true;
                    is_marked_sample(j) = true;
                    break;
                }

            }
            if (!is_memorized)
            {
                _x_ped_samples_temp.col(i).setConstant(obstacles(i).x);
                _y_ped_samples_temp.col(i).setConstant(obstacles(i).y);
                _vx_ped_samples_temp.col(i).setConstant(0.0);
                _vy_ped_samples_temp.col(i).setConstant(0.0);
                _vx_ped_samples_temp(0, i) = obstacles(i).vx;
                _vy_ped_samples_temp(0, i) = obstacles(i).vy;
            }
            _id_samples_temp(i) = obstacles(i).id;
        }

        _x_ped_samples = _x_ped_samples_temp;
        _y_ped_samples = _y_ped_samples_temp;
        _vx_ped_samples = _vx_ped_samples_temp;
        _vy_ped_samples = _vy_ped_samples_temp;
        _id_samples = _id_samples_temp;

        _x_robot_samples.segment(1, _N_samples-1) = _x_robot_samples.segment(0, _N_samples-1).eval();
        _x_robot_samples(0) = x_robot;
        _y_robot_samples.segment(1, _N_samples-1) = _y_robot_samples.segment(0, _N_samples-1).eval();
        _y_robot_samples(0) = y_robot;
        _theta_robot_samples.segment(1, _N_samples-1) = _theta_robot_samples.segment(0, _N_samples-1).eval();
        _theta_robot_samples(0) = theta_robot;
        _xP_robot_samples.segment(1, _N_samples-1) = _xP_robot_samples.segment(0, _N_samples-1).eval();
        _xP_robot_samples(0) = x_robot + _Pdist*std::cos(theta_robot);
        _yP_robot_samples.segment(1, _N_samples-1) = _yP_robot_samples.segment(0, _N_samples-1).eval();
        _yP_robot_samples(0) = y_robot + _Pdist*std::sin(theta_robot);

        double curr_v_long = _optimVect_no_sl(0)*std::cos(theta_robot) + _optimVect_no_sl(1)*std::sin(theta_robot);

        _vx_robot_samples.segment(1, _N_samples-1) = _vx_robot_samples.segment(0, _N_samples-1).eval();
        _vx_robot_samples(0) = curr_v_long*std::cos(theta_robot);
        _vy_robot_samples.segment(1, _N_samples-1) = _vy_robot_samples.segment(0, _N_samples-1).eval();
        _vy_robot_samples(0) = curr_v_long*std::sin(theta_robot);
    }
/*if (std::fabs(sim_time - 60.0) < 1e-6)
{
    _is_started = false;
}*/
}

void MPC_diffDrive_fblin::sample_Pedestrian_and_robot_states_complete(Eigen::VectorX<obstacle>& obstacles)//, double sim_time)
{
    obstacle zero = {};

    if (!_is_started_complete)
    {
        _obstacles_samples.resize(_N_samples, _n_obs);
        _obstacles_samples.setConstant(zero);
        _id_samples_obstacles = Eigen::VectorXi::Zero(_n_obs);

        for (int i = 0; i < _n_obs; i++)
        {
            _obstacles_samples.col(i).setConstant(obstacles(i));

            _id_samples_obstacles(i) = obstacles(i).id;
        }

        _is_started_complete = true;
    }else
    {

        Eigen::MatrixX<obstacle> _obstacles_samples_temp;
        _obstacles_samples_temp.resize(_N_samples, _n_obs);
        _obstacles_samples_temp.setConstant(zero);

        Eigen::VectorXi _id_samples_obstacles_temp = Eigen::VectorXi::Zero(_n_obs);
        Eigen::VectorX<bool> _is_static_obstacles_temp(_n_obs);
        _is_static_obstacles_temp.setConstant(false);

        Eigen::VectorX<bool> is_marked(_id_samples_obstacles.size());
        is_marked.setConstant(false);

        for (int i = 0; i < _n_obs; i++)
        {
            bool is_memorized = false;
            for (int j = 0; j < _id_samples_obstacles.size(); j++)
            {
                if (!is_marked(j) && (obstacles(i).id == _id_samples_obstacles(j) || std::hypot(obstacles(i).x - (_obstacles_samples(0, j).x + _obstacles_samples(0, j).vx*_MPC_Ts), obstacles(i).y - (_obstacles_samples(0, j).y + _obstacles_samples(0, j).vy*_MPC_Ts)) < _r_red*2 + std::hypot(_obstacles_samples(0, j).vx*_MPC_Ts, _obstacles_samples(0, j).vy*_MPC_Ts)))
                {
                    _obstacles_samples_temp.col(i).segment(1, _N_samples-1) = _obstacles_samples.col(j).segment(0, _N_samples-1).eval();
                    _obstacles_samples_temp(0, i) = obstacles(i);
                    is_marked(j) = true;

                    if (_obstacles_samples_temp(1, i).is_person && _obstacles_samples_temp(2, i).is_person && _obstacles_samples_temp(3, i).is_person){

                        _obstacles_samples_temp(0, i).is_person = true;
                        obstacles(i).is_person = true;
                        
                    }

                    is_memorized = true;
                    break;
                }

            }
            if (!is_memorized)
            {
                _obstacles_samples_temp.col(i).setConstant(obstacles(i));

            }
            _id_samples_obstacles_temp(i) = obstacles(i).id;
        }

        _obstacles_samples = _obstacles_samples_temp;

        _id_samples_obstacles = _id_samples_obstacles_temp;

    }
/*if (std::fabs(sim_time - 60.0) < 1e-6)
{
    _is_started = false;
}*/
}

/*void MPC_diffDrive_fblin::clustering()
{
    Eigen::VectorXd x_ped = _x_ped_samples.row(0).transpose();
    Eigen::VectorXd y_ped = _y_ped_samples.row(0).transpose();
    double x_robot = _x_robot_samples(0);
    double y_robot = _y_robot_samples(0);

    arma::uword n_obs = x_ped.size();

    arma::rowvec x_ped_arma(n_obs);
    arma::rowvec y_ped_arma(n_obs);

    for (arma::uword i = 0; i < n_obs; ++i) {
        x_ped_arma(i) = x_ped(static_cast<Eigen::Index>(i));
        y_ped_arma(i) = y_ped(static_cast<Eigen::Index>(i));
    }

    arma::mat ped_pos(2, n_obs, arma::fill::zeros);
    ped_pos.row(0) = x_ped_arma;
    ped_pos.row(1) = y_ped_arma;

    mlpack::dbscan::DBSCAN<> dbscan(2.7, 1);
    arma::Row<std::size_t> labels;

    std::size_t n_groups = dbscan.Cluster(ped_pos, labels);

    _labels = Eigen::VectorXi::Zero(static_cast<int>(n_obs));

    for (arma::uword i = 0; i < n_obs; ++i) {
        _labels(static_cast<int>(i)) = static_cast<int>(labels(i));
    }
    _n_groups = static_cast<int>(n_groups);

    _x_near = Eigen::VectorXd::Zero(static_cast<int>(n_groups));
    _y_near = Eigen::VectorXd::Zero(static_cast<int>(n_groups));

    for (std::size_t i = 0; i < n_groups; ++i)
    {
        std::vector<double> X_group;
        std::vector<double> Y_group;
        std::vector<double> dist;
        for (arma::uword j = 0; j < n_obs; ++j)
        {
            if (labels(j) == i)
            {
                const double x = x_ped(static_cast<Eigen::Index>(j));
                const double y = y_ped(static_cast<Eigen::Index>(j));

                 X_group.push_back(x);
                 Y_group.push_back(y);
                dist.push_back(std::hypot((x_robot - x), (y_robot - y)));

            }


        }

        if (!dist.empty())
        {
            auto it = std::min_element(dist.begin(), dist.end());
            long I = std::distance(dist.begin(), it);
            _x_near(static_cast<Eigen::Index>(i)) = X_group[I];
            _y_near(static_cast<Eigen::Index>(i)) = Y_group[I];
        }
    }
}*/

bool MPC_diffDrive_fblin::is_inside_convex_hull(Eigen::Vector2d robot_pos, Eigen::MatrixXd ped_pos, Eigen::MatrixXd v_ped, Eigen::MatrixXd relative_velocity, double r_c, double r_red, Eigen::VectorXi idx, int i_group, int Na, Eigen::VectorX<obstacle> obstacles, Eigen::Vector2d final_ref_pos, const double& tau)//, bool& is_blocked)
{
    double tol = 1e-8;

    long n_obs_group = (idx.array() == i_group).count();

    Eigen::MatrixXd ped_pos_cluster = Eigen::MatrixXd::Zero(n_obs_group, 2);
    Eigen::MatrixXd v_ped_cluster = Eigen::MatrixXd::Zero(n_obs_group, 2);
    Eigen::MatrixXd relative_velocity_cluster = Eigen::MatrixXd::Zero(n_obs_group, 2);
    Eigen::VectorX<obstacle> obstacles_cluster;
    obstacles_cluster.resize(n_obs_group);
    obstacle zero = {};
    obstacles_cluster.setConstant(zero);

    int x = 0;

    for (int i = 0; i < _n_obs; i++)
    {
        if (idx(i) == i_group)
        {
            ped_pos_cluster.row(x) = ped_pos.row(i);
            v_ped_cluster.row(x) = v_ped.row(i);
            relative_velocity_cluster.row(x) = relative_velocity.row(i);
            obstacles_cluster(x) = obstacles(i);
            x++;
        }
    }

    

    Eigen::MatrixXd PSf_transpose_cluster = Eigen::MatrixXd::Zero(2, n_obs_group);
    Eigen::VectorXd theta_heading_cluster = Eigen::VectorXd::Zero(n_obs_group);

    std::vector<Eigen::Matrix2d> D;
    D.resize(n_obs_group, Eigen::Matrix2d::Zero());

    std::vector<Eigen::Matrix2d> B;
    B.resize(n_obs_group, Eigen::Matrix2d::Zero());

    Eigen::VectorXd vector_cross = Eigen::VectorXd::Zero(2*n_obs_group);

    Eigen::VectorXd theta = Eigen::VectorXd::LinSpaced(Na+1, 0, 2*M_PI).head(Na).eval();

    Eigen::VectorXd n_1 = theta.array().cos().matrix();
    Eigen::VectorXd n_2 = theta.array().sin().matrix();

    Eigen::MatrixXd n = Eigen::MatrixXd::Zero(Na, 2);
    n.col(0) = n_1;
    n.col(1) = n_2;

    Eigen::MatrixXd s = Eigen::MatrixXd::Zero(n_obs_group, Na);

    double delta = (r_c - r_red);

    bool is_inside = false;


    for (int i = 0; i < n_obs_group; i++)
    {
        if (obstacles_cluster(i).is_person)
        {
            personal_space personalSpace = personalSpaceFunction_mod(v_ped_cluster.row(i), relative_velocity_cluster.row(i).transpose(), robot_pos, ped_pos_cluster.row(i).transpose(), r_c);
            PSf_transpose_cluster.col(i) = personalSpace.PSf;
            theta_heading_cluster(i) = personalSpace.theta_heading;

            double a = PSf_transpose_cluster(1, i);// + delta;
            double b = PSf_transpose_cluster(0, i);// + delta;

            D[i] << a, 0.0,
                    0.0, b;

            Eigen::Matrix2d R;
            R << std::cos(theta_heading_cluster(i)), -std::sin(theta_heading_cluster(i)),
                 std::sin(theta_heading_cluster(i)), std::cos(theta_heading_cluster(i));

            B[i] = R*D[i];

            Eigen::RowVector2d v;
            v << std::cos(theta_heading_cluster(i)), std::sin(theta_heading_cluster(i));

            for (int j = 0; j < n.rows(); j++)
            {
                if (n.row(j).dot(v) >= 0.0)
                {
                    s(i,j) = n.row(j)*ped_pos_cluster.row(i).transpose() + (B[i].transpose()*n.row(j).transpose()).norm();
                }else
                {
                    s(i,j) = n.row(j)*ped_pos_cluster.row(i).transpose() + b;
                }

            }
        }else
        {
            for (int j = 0; j < n.rows(); j++)
            {
                    s(i,j) = compute_s_objects(n.row(j), r_red, r_c, obstacles_cluster(i), ped_pos_cluster.row(i).transpose());
            }
        }

    }

    Eigen::Vector2d csi;
    csi << _actXP, _actYP;

    Eigen::VectorXd s_convex_hull = s.colwise().maxCoeff().transpose();
    Eigen::VectorXd dist = s_convex_hull - n*final_ref_pos;
    Eigen::VectorXd dist_robot_pos = s_convex_hull - n*csi;

    bool isInsideGroup = (dist.array() >= -tol).all();
    bool isInsideGroup_robot_pos = (dist_robot_pos.array() >= -tol).all();

    if (isInsideGroup || isInsideGroup_robot_pos || _is_blocked)
    {
        is_inside = true;
    }else
    {
        is_inside = false;
    }

    //is_blocked = false;

    /*if (!isInsideGroup && (isInsideGroup_robot_pos && _is_blocked))
    {
        
        is_blocked = true;

    }*/

    /*if (isInsideGroup_robot_pos || _is_blocked)
    {
        
        is_blocked = true;

    }*/

    return is_inside;
}

Eigen::Vector2d MPC_diffDrive_fblin::compute_new_ref(Eigen::Vector2d robot_pos, Eigen::MatrixXd ped_pos, Eigen::MatrixXd v_ped, Eigen::MatrixXd relative_velocity, double r_c, double r_red, Eigen::VectorXi idx, int i_group, Eigen::VectorX<obstacle> obstacles, Eigen::Vector2d final_ref_pos, bool& is_closed)
{
    double tol = 1e-8;

    Eigen::Vector2d new_ref = Eigen::Vector2d::Zero();

    long n_obs_group = (idx.array() == i_group).count();

    Eigen::MatrixXd ped_pos_cluster = Eigen::MatrixXd::Zero(n_obs_group, 2);
    Eigen::MatrixXd v_ped_cluster = Eigen::MatrixXd::Zero(n_obs_group, 2);
    Eigen::MatrixXd relative_velocity_cluster = Eigen::MatrixXd::Zero(n_obs_group, 2);
    Eigen::VectorX<obstacle> obstacles_cluster;
    obstacles_cluster.resize(n_obs_group);
    obstacle zero = {};
    obstacles_cluster.setConstant(zero);

    int x = 0;

    for (int i = 0; i < _n_obs; i++)
    {
        if (idx(i) == i_group)
        {
            ped_pos_cluster.row(x) = ped_pos.row(i);
            v_ped_cluster.row(x) = v_ped.row(i);
            relative_velocity_cluster.row(x) = relative_velocity.row(i);
            obstacles_cluster(x) = obstacles(i);
            x++;
        }
    }

    

    Eigen::MatrixXd PSf_transpose_cluster = Eigen::MatrixXd::Zero(2, n_obs_group);
    Eigen::VectorXd theta_heading_cluster = Eigen::VectorXd::Zero(n_obs_group);

    std::vector<Eigen::Matrix2d> D;
    D.resize(n_obs_group, Eigen::Matrix2d::Zero());

    std::vector<Eigen::Matrix2d> B;
    B.resize(n_obs_group, Eigen::Matrix2d::Zero());

    Eigen::VectorXd m = Eigen::VectorXd::Zero(n_obs_group);

    Eigen::RowVectorXd coefficients = Eigen::RowVectorXd::Zero(3);

    Eigen::VectorXd hx_tot_cluster = Eigen::VectorXd::Zero(n_obs_group);

    Eigen::VectorXd hy_tot_cluster = Eigen::VectorXd::Zero(n_obs_group);

    Eigen::VectorXd l_tot_cluster = Eigen::VectorXd::Zero(n_obs_group);

    Eigen::VectorXd dist_tot_cluster = Eigen::VectorXd::Zero(n_obs_group);

    double delta = (r_c - r_red);

    double x_centroid = ped_pos_cluster.col(0).mean();
    double y_centroid = ped_pos_cluster.col(1).mean();

    std::vector<P> pts;

    for (int i = 0; i < n_obs_group; i++)
    {
        P pts_temp = {0.0, 0.0, 0.0};
        
        if (!obstacles_cluster(i).is_segment){
        
        pts_temp.x = obstacles_cluster(i).x;
        pts_temp.y = obstacles_cluster(i).y;
        pts_temp.radius = obstacles_cluster(i).radius;
        pts.push_back(pts_temp);
        }else{
            double v_norm = std::hypot(obstacles_cluster(i).x2 - obstacles_cluster(i).x1, obstacles_cluster(i).y2 - obstacles_cluster(i).y1) + 1e-6;
            int N_segments = static_cast<int>(std::ceil(v_norm/(2*delta*std::sqrt(2))));
            double vx = (obstacles_cluster(i).x2 - obstacles_cluster(i).x1)/v_norm;
            double vy = (obstacles_cluster(i).y2 - obstacles_cluster(i).y1)/v_norm;
            double new_delta = v_norm/N_segments;

            for (int j=0; j <= N_segments; ++j){
                pts_temp.x = obstacles_cluster(i).x1 + j*vx*new_delta;
                pts_temp.y = obstacles_cluster(i).y1 + j*vy*new_delta;
                pts_temp.radius = delta*std::sqrt(2);
                pts.push_back(pts_temp);   
            }

        }    
    }

    dedupKeepMaxParam(pts);

    auto H = convexHull(pts);
    auto gaps = gapsOnHullAboveDelta(H, delta);

    /*double ax = 0.0, ay = 0.0, a_radius = 0.0;
    double bx = 0.0, by = 0.0, b_radius = 0.0;
    double clearance = 0.0;

    if (edgeRes.ia >= 0) {
    ax = edgeRes.a.x;
    ay = edgeRes.a.y;
    a_radius = edgeRes.a.radius;

    bx = edgeRes.b.x;
    by = edgeRes.b.y;
    b_radius = edgeRes.b.radius;

    clearance = edgeRes.clearance; // già = dist - ra - rb
} else {
    is_closed = true;
    return robot_pos;
}*/

    if (!gaps.empty()){
        is_closed = false;

        Eigen::Vector2d n_best = Eigen::Vector2d::Zero();
        LongestEdgeResult gap_best = {};
        double dist_best = std::numeric_limits<double>::infinity();

        bool found = false;

        for (int i = 0; i < gaps.size(); i++){
        
        double dist_temp = 0.0;    

        Eigen::Vector2d e(gaps[i].b.x - gaps[i].a.x, gaps[i].b.y - gaps[i].a.y);
        Eigen::Vector2d n(e(1), -e(0));
        Eigen::Vector2d a_point;
        a_point << gaps[i].a.x, gaps[i].a.y;

        Eigen::Vector2d centroid;
        centroid << x_centroid, y_centroid;


        double s = (centroid - a_point).dot(n);
        if (s < 0){
            n = -n;
        }

        dist_temp = std::hypot((gaps[i].b.x + gaps[i].a.x)/2 - _actXP, (gaps[i].b.y + gaps[i].a.y)/2 - _actYP);

        if (dist_temp < dist_best){ //&& (final_ref_pos - Eigen::Vector2d((gaps[i].b.x + gaps[i].a.x)/2, (gaps[i].b.y + gaps[i].a.y)/2)).dot(n) >= 0){
            n_best = n;
            gap_best = gaps[i];
            dist_best = dist_temp;
            found = true;
        }
    }
    
    if (n_best.norm() < 1e-12 || !found) { is_closed = true; return robot_pos; }
    
    double n_norm = n_best.norm();

    coefficients(0) = n_best(0)/n_norm;
    coefficients(1) = n_best(1)/n_norm;

    coefficients(2) = coefficients(0)*gap_best.a.x + coefficients(1)*gap_best.a.y;


    for (int i = 0; i < n_obs_group; i++)
    {
        if (obstacles_cluster(i).is_person)
        {
            personal_space personalSpace = personalSpaceFunction_mod(v_ped_cluster.row(i), relative_velocity_cluster.row(i).transpose(), robot_pos, ped_pos_cluster.row(i).transpose(), r_c);
            PSf_transpose_cluster.col(i) = personalSpace.PSf;
            theta_heading_cluster(i) = personalSpace.theta_heading;

            double a = PSf_transpose_cluster(1, i) + delta;
            double b = PSf_transpose_cluster(0, i) + delta;

            D[i] << a, 0.0,
                    0.0, b;

            Eigen::Matrix2d R;
            R << std::cos(theta_heading_cluster(i)), -std::sin(theta_heading_cluster(i)),
                 std::sin(theta_heading_cluster(i)), std::cos(theta_heading_cluster(i));

            B[i] = R*D[i];

            Eigen::RowVector2d v;
            v << std::cos(theta_heading_cluster(i)), std::sin(theta_heading_cluster(i));

            if (coefficients.segment(0, 2).dot(v) >= 0.0)
            {
                m(i) = coefficients.segment(0, 2)*ped_pos_cluster.row(i).transpose() - coefficients(2) - b;
            }else
            {
                m(i) = coefficients.segment(0, 2)*ped_pos_cluster.row(i).transpose() - coefficients(2) - (B[i].transpose()*coefficients.segment(0, 2).transpose()).norm();
            }
        }else
        {
            m(i) = compute_m_objects(coefficients.segment(0, 2), r_red, r_c, obstacles_cluster(i), ped_pos_cluster.row(i).transpose()) - coefficients(2);
        }


    }



    bool isHullTangentOriginal = (m.array() >= -tol).all();

    double hx = coefficients(0);
    double hy = coefficients(1);
    double l;

    Eigen::Vector2d h;
    h << hx, hy;

    int I_m;
    m.minCoeff(&I_m);

    if (isHullTangentOriginal)
    {
        l = coefficients(2) - 2.0;
    }else
    {
        l =  coefficients(2) + m(I_m) - 2.0;
    }
    double h_norm = h.norm();
    new_ref(0) = -(hx/h_norm)*std::abs(coefficients(2) - l) + (gap_best.a.x + gap_best.b.x)/2;
    new_ref(1) = -(hy/h_norm)*std::abs(coefficients(2) - l) + (gap_best.a.y + gap_best.b.y)/2;
}else{
    is_closed = true;
    return robot_pos;
}   

    return new_ref;
}

void MPC_diffDrive_fblin::clustering(const Eigen::Vector2d& final_ref_pos)
{
    Eigen::VectorX<obstacle> obstacles = _obstacles_samples.row(0).transpose();
    double x_robot = _x_robot_samples(0);
    double y_robot = _y_robot_samples(0);

    int n_obs = obstacles.size();

    int n_circles = 0;
    int n_segments = 0;

    for (int i = 0; i < n_obs; ++i) {
        if (!obstacles(i).is_segment){
        n_circles++;
        }else{
        n_segments++;  
        }
    }
    
    std::vector<Circle> circles;
    circles.reserve(n_circles);

    std::vector<Segment> segments;
    segments.reserve(n_segments);

    double tau = 2.0;

    for (int i = 0; i < n_obs; ++i) {
        if (!obstacles(i).is_segment){
        Circle circle = {0.0, 0.0, 0.0};
        circle.x = obstacles(i).x;
        circle.y = obstacles(i).y;
        circle.r = obstacles(i).radius;
        circles.push_back(circle);
        }else{
        Segment segment = {0.0, 0.0, 0.0, 0.0};
        segment.x1 = obstacles(i).x1;
        segment.y1 = obstacles(i).y1;
        segment.x2 = obstacles(i).x2;
        segment.y2 = obstacles(i).y2;
        segments.push_back(segment);   
        }
    }

    auto clusters = cluster_circles_and_segments_long_separate(circles, segments, tau, std::numeric_limits<int>::max(), 2.0);

    _labels = Eigen::VectorXi::Zero(n_obs);

    for (int cid = 0; cid < (int)clusters.size(); ++cid) {
        for (int idx : clusters[cid])
        {
            _labels(idx) = cid;
        }
    }

    _n_groups = (int)clusters.size();

    Eigen::Vector2d robot_pos = Eigen::Vector2d::Zero();
    robot_pos << x_robot, y_robot;

    Eigen::MatrixXd ped_pos = Eigen::MatrixXd::Zero(_n_obs, 2);
    ped_pos.col(0) = _x_ped_samples.row(0).transpose();
    ped_pos.col(1) = _y_ped_samples.row(0).transpose();

    Eigen::MatrixXd v_ped(_n_obs, 2);
    v_ped << _vx_ped_samples.row(0).transpose(), _vy_ped_samples.row(0).transpose();

    Eigen::MatrixXd relative_velocity = Eigen::MatrixXd::Zero(_n_obs, 2);
    relative_velocity.col(0) = _vx_robot_samples(0)*Eigen::VectorXd::Ones(_n_obs) - _vx_ped_samples.row(0).transpose();
    relative_velocity.col(1) = _vy_robot_samples(0)*Eigen::VectorXd::Ones(_n_obs) - _vy_ped_samples.row(0).transpose();
    
    int Na = 50;

    Eigen::Vector2d P_position;
    P_position << _actXP, _actYP;

    bool is_inside_flag = true;
    bool is_inside_flag_first = true;
    bool is_blocked_inside_hull_internal = false;

    while (is_inside_flag){

        is_inside_flag = false;

    Eigen::VectorXi labels_base = _labels;

    //debugMsg(std::to_string(_ref_temp(0)));
    //debugMsg(std::to_string(_ref_temp(1)));
    //debugMsg(std::to_string(_is_blocked));

    for (int i = 0; i < _n_groups; ++i)
    {
        //bool is_blocked_inside_hull;
        bool is_inside = is_inside_convex_hull(robot_pos, ped_pos, v_ped, relative_velocity, _r_red, _r_red, labels_base, i, Na, obstacles, final_ref_pos, tau);//, is_blocked_inside_hull);

        //debugMsg(std::to_string(is_blocked_inside_hull));
        
        /*if (is_blocked_inside_hull && is_inside_flag_first){
                //_is_blocked_inside_hull = true;
                //_ref_temp = compute_new_ref(robot_pos, ped_pos, v_ped, relative_velocity, _r_c, _r_red, labels_base, i, obstacles, final_ref_pos, _is_closed);
                is_blocked_inside_hull_internal = true;
                is_inside_flag_first = false;
        }*/
        
        if (is_inside){

            std::vector<Circle> circles_temp;
            std::vector<Segment> segments_temp;
            std::vector<int> global_idx;

            int x = 0;

            for (int j = 0; j < n_circles; j++)
            {
                if (_labels(j) == i)
                {
                circles_temp.push_back(circles[j]);
                global_idx.push_back(j);
                x++;
                }
            }

            for (int j = 0; j < n_segments; j++)
            {
                if (_labels(j + n_circles) == i)
                {
                segments_temp.push_back(segments[j]);
                global_idx.push_back(j + n_circles);
                x++;
                }
            }

            auto clusters_spezzati = cluster_circles_and_segments_long_separate(circles_temp, segments_temp, tau, (x+1)/2, 2.0);

            if (clusters_spezzati.size() <= 1) {
                continue;
            }

            is_inside_flag = true;

            bool is_first = true;

            int new_id = _labels.maxCoeff() + 1;

            for (int cid = 0; cid < (int)clusters_spezzati.size(); ++cid) {
                if (is_first){
                        is_first = false;
                }else{
                    for (int idx : clusters_spezzati[cid]){
                        _labels(global_idx[idx]) = new_id;
                    }
                    new_id++;
                }
            }
        }
    }
    
    _n_groups = _labels.maxCoeff() + 1;
}
/*if (!is_blocked_inside_hull_internal && _is_blocked){
_is_blocked = false;
_n_fail = 0;
_sp = 1.0e9;
//_is_blocked_inside_hull = false;
}*/


    _x_near = Eigen::VectorXd::Zero(_n_groups);
    _y_near = Eigen::VectorXd::Zero(_n_groups);

    for (int i = 0; i < _n_groups; ++i)
    {
        std::vector<double> X_group;
        std::vector<double> Y_group;
        std::vector<double> dist;
        for (int j = 0; j < n_obs; ++j)
        {
            if (_labels(j) == i)
            {
                const double x = obstacles(j).x;
                const double y = obstacles(j).y;

                 X_group.push_back(x);
                 Y_group.push_back(y);
                dist.push_back(std::max(0.0, std::hypot((x_robot - x), (y_robot - y)) - obstacles(j).radius));

            }


        }

        if (!dist.empty())
        {
            auto it = std::min_element(dist.begin(), dist.end());
            long I = std::distance(dist.begin(), it);
            _x_near(i) = X_group[I];
            _y_near(i) = Y_group[I];
        }
    }
}

void MPC_diffDrive_fblin::set_referenceRobotState(double x, double y, double yaw)
{
    // Transform from robot state to point P state
    double xP_ref, yP_ref;
    _fblinController->reference_transformation(x, y, yaw, xP_ref, yP_ref);

    // Fill in the reference vector along the prediction horizon with a constant reference
    for (auto i=0; i<_N+1; i++) {
        _refMPCstate.block(i*2, 0, 2, 1) = Eigen::Vector2d(xP_ref, yP_ref);
    }
}

void MPC_diffDrive_fblin::set_referenceRobotState(const Eigen::VectorXd& refRobotState)
{
    if (refRobotState.size()!=3*(_N+1))
    {
        errorMsg("[MPC_diffDrive_fblin.set_referenceRobotState] The refState variable has the wrong size");
        _refMPCstate = refRobotState.segment(0, 3*(_N+1));
    }
    else
    {
        // Transform from robot state to point P state and store in reference vector
        for (auto i=0; i<_N+1; i++) {
            _fblinController->reference_transformation(refRobotState(3*i), refRobotState(3*i+1),
                                                       refRobotState(3*i+2), _refMPCstate(2*i), _refMPCstate(2*i+1));
        }
    }
}

void MPC_diffDrive_fblin::get_actualMPCControl(double& vPx, double& vPy)
{
    vPx = _optimVect(0);
    vPy = _optimVect(1);
}

void MPC_diffDrive_fblin::get_actualMPCstate(double& xP, double& yP)
{
    xP = _actXP;
    yP = _actYP;
}

void MPC_diffDrive_fblin::get_referenceMPCstate(Eigen::VectorXd& refState)
{
    refState = _refMPCstate;
}

void MPC_diffDrive_fblin::get_actualControl(double &linVelocity, double &angVelocity)
{
    linVelocity = _linearVelocity;
    angVelocity = _angularVelocity;
}

void MPC_diffDrive_fblin::get_v_robot(double &vx, double &vy, double &v_long)
{
    vx = _vx_robot;
    vy = _vy_robot;
    v_long = _v_long_robot;
}

/** Private methods */
void MPC_diffDrive_fblin::compute_AcalMatrix() {
    // Initialize Acal matrix
    _Acal = Eigen::MatrixXd::Zero(_N*2, 2);

    // Compute Acal matrix
    _Acal.block(0, 0, 2, 2) = Eigen::MatrixXd::Identity(2, 2);

    for (int k=0; k<_N; k++) {
        _Acal.block(2*k, 0, 2, 2) = _plant_A.pow(k+1);
    }

    // Check matrix
    // saveMatrixToFile("Acal_matrix.csv", _Acal);
}

void MPC_diffDrive_fblin::compute_BcalMatrix() {
    // Initialize Bcal matrix
    _Bcal = Eigen::MatrixXd::Zero(_N*2, _N*2);

    // Compute Bcal matrix
    for (int i=1; i<=_N; i++)
        for (int j=1; j<=i; j++)
            _Bcal.block(2*(i-1), 2*(j-1), 2, 2) = _plant_A.pow(i-j)*_plant_B;

    // Check matrix
    //  saveMatrixToFile("Bcal_matrix.csv", _Bcal);
}

void MPC_diffDrive_fblin::compute_QcalMatrix() {
    // Initialize Qcal matrix
    _Qcal = Eigen::MatrixXd::Zero((_N)*2, (_N)*2);

    // Compute Qcal matrix
    _Qcal.diagonal().segment(0, _N*2-2) = Eigen::VectorXd::Ones(_N*2-2)*_q;
    _Qcal.diagonal().segment(_N*2-2, 2) = Eigen::VectorXd::Ones(2)*_p;

    // Check matrix
//    saveMatrixToFile("Qcal_matrix.csv", _Qcal);
}

void MPC_diffDrive_fblin::compute_RcalMatrix() {
    // Initialize Rcal matrix
    _Rcal = Eigen::MatrixXd::Zero(_N*2, _N*2);

    // Compute Rcal matrix
    _Rcal.diagonal() = Eigen::VectorXd::Ones(_N*2)*_r;

    // Check matrix
//    saveMatrixToFile("Rcal_matrix.csv", _Rcal);
}

void MPC_diffDrive_fblin::compute_objectiveMatrix() {
    // Initialize H and f matrices
    _H_no_sl = Eigen::MatrixXd::Zero(2*_N, 2*_N);
    _f_no_sl = Eigen::VectorXd::Zero(_H.rows());

    // Compute H and f matrices
    _H_no_sl = _Bcal.transpose()*_Qcal.block(0, 0, 2*_N, 2*_N)*_Bcal+_Rcal;
    _f_no_sl = (_Acal*Eigen::Vector2d(_actXP, _actYP)-_refMPCstate).transpose()*_Qcal.block(0, 0, 2*_N, 2*_N)*_Bcal;

    _H.block(0, 0, 2*_N, 2*_N) = _H_no_sl;
    _H.block(2*_N, 2*_N, _n_groups, _n_groups) = Eigen::MatrixXd::Identity(_n_groups, _n_groups)*_sp;
    _H.block(2*_N + _n_groups, 2*_N + _n_groups, 2, 2) = Eigen::MatrixXd::Identity(2, 2)*_sa;

    _f.segment(0, 2*_N) = _f_no_sl;
    _f.segment(2*_N, _n_groups + 2) = Eigen::VectorXd::Zero(_n_groups + 2);

    // Check matrix
//    saveMatrixToFile("H_matrix.csv", _H);
//    saveMatrixToFile("f_matrix.csv", _f);
}

void MPC_diffDrive_fblin::compute_wheelVelocityConstraint() {
    // Initialize Ain_vel and Bin_vel matrices
    _Ain_vel = Eigen::MatrixXd::Zero(2*2*_N, 2*_N + _n_groups + 2);
    _Bin_vel = Eigen::VectorXd::Zero(2*2*_N);

    // Compute constraint matrices
    for (auto k=0; k<_N; k++) {
        Eigen::Matrix2d Abar = (1/(2.0*_wheelRadius))*Eigen::Matrix2d({{2.0*std::cos(_predictRobotState(3*k+2))-_track/_Pdist*std::sin(_predictRobotState(3*k+2)),
                               2.0*std::sin(_predictRobotState(3*k+2))+_track/_Pdist*std::cos(_predictRobotState(3*k+2))},
                              {2.0*std::cos(_predictRobotState(3*k+2))+_track/_Pdist*std::sin(_predictRobotState(3*k+2)),
                               2.0*std::sin(_predictRobotState(3*k+2))-_track/_Pdist*std::cos(_predictRobotState(3*k+2))}});
        _Ain_vel.block(2*k, 2*k, 2, 2) = Abar;
        _Ain_vel.block(2*(k+_N), 2*k, 2, 2) = -Abar;

        _Bin_vel(2*k) = _wheelVelMax;
        _Bin_vel(2*k+1) = _wheelVelMax;
        _Bin_vel(2*(k+_N)) = -_wheelVelMin;
        _Bin_vel(2*(k+_N)+1) = -_wheelVelMin;
    }

    // Check matrix
//    saveMatrixToFile("Ain_vel_matrix.csv", _Ain_vel);
//    saveMatrixToFile("Bin_vel_matrix.csv", _Bin_vel);
}

void MPC_diffDrive_fblin::saveMatrixToFile(std::string fileName, Eigen::MatrixXd matrix) {
    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");

    std::ofstream file(fileName);
    if (file.is_open())
    {
        file << matrix.format(CSVFormat);
        file.close();
    }
}

void MPC_diffDrive_fblin::set_ErrorMsgCallback(ErrorMsgCallback callback) {
    _error = callback;
}

void MPC_diffDrive_fblin::set_InfoMsgCallback(InfoMsgCallback callback) {
    _info = callback;
}

void MPC_diffDrive_fblin::set_DebugMsgCallback(DebugMsgCallback callback) {
    _debug = callback;
}

void MPC_diffDrive_fblin::errorMsg(const std::string& message) {
    if (_error) {
        _error(message);
    }
}

void MPC_diffDrive_fblin::infoMsg(const std::string& message) {
    if (_info) {
        _info(message);
    }
}

void MPC_diffDrive_fblin::debugMsg(const std::string& message) {
    if (_debug) {
        _debug(message);
    }
}