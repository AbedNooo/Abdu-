#ifndef MPC_DIFFDRIVE_FBLIN_H
#define MPC_DIFFDRIVE_FBLIN_H

#include <boost/function.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <algorithm>

#include "solver/GUROBIsolver.h"
#include "controller/fblin_unicycle.h"
#include "multitarget_tracker/HungarianAlg.h"
#include <mlpack/methods/dbscan/dbscan.hpp>
#include <functional>
#include <utility>
#include <limits>

#define GUROBI_LICENSEID 2745424
#define GUROBI_USERNAME  "alessandro"

// Message handlers
typedef boost::function<void(const std::string &)> DebugMsgCallback;
typedef boost::function<void(const std::string &)> InfoMsgCallback;
typedef boost::function<void(const std::string &)> ErrorMsgCallback;

struct tangent
{
    Eigen::Vector2d Q1;
    Eigen::Vector2d Q2;
    bool robot_inside_obst;
};

struct personal_space
{
    Eigen::Vector2d PSf;
    double theta_heading;
    int scenario;
};

struct tangent_local
{
    Eigen::Vector2d Q1;
    Eigen::Vector2d Q2;
};

struct constraint
{
    double hx;
    double hy;
    double l;
    double d_min;
};

struct point
{
    double x;
    double y;
    double z = 0;
};

struct obstacle
{
    double x;
    double y;
    double x1;
    double y1;
    double x2;
    double y2;
    double vx;
    double vy;
    double radius;
    std::vector<point> polygon;
    bool is_person;
    bool is_segment;
    int id;
};

struct Track
{
  obstacle state;         // stato “pubblicabile”
  int missed = 0;         // frames senza match
};

struct AssocParams
{
  double dt = 0.2;

  double gate_dist = 0.6;     // [m] soglia gating (rifiuta match lontani)
  int max_missed = 3;         // dopo quanti frame “muore” un track

  double alpha_pos = 0.8;     // [0..1] quanto segui la misura per x,y
  double vel_decay = 0.7;     // [0..1] decadimento verso 0 (vx,vy)
};

struct P {
  double x, y;
  double radius;
};

struct LongestEdgeResult {
  P a, b;
  double clearance, length;
  int ia, ib;
};

struct DSUCap {
    std::vector<int> parent, rank, sz;

    DSUCap(int n) : parent(n), rank(n, 0), sz(n, 1) {
        for (int i = 0; i < n; ++i) parent[i] = i;
    }

    int find(int x) {
        if (parent[x] != x) parent[x] = find(parent[x]);
        return parent[x];
    }

    // prova a unire rispettando il massimo numero di elementi per cluster
    bool unite(int x, int y, int max_size) {
        int rx = find(x), ry = find(y);
        if (rx == ry) return false;

        if (sz[rx] + sz[ry] > max_size) return false; // <-- vincolo

        // union by rank
        if (rank[rx] < rank[ry]) std::swap(rx, ry);
        parent[ry] = rx;
        sz[rx] += sz[ry];
        if (rank[rx] == rank[ry]) rank[rx]++;

        return true;
    }
};

struct Circle {
    double x, y, r;
};

struct Segment {
    double x1, y1, x2, y2;
};

static inline double dist_centri(const Circle& a, const Circle& b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return std::sqrt(dx*dx + dy*dy);
}


struct Edge {
    int u, v;
    double w;
};

struct pred_step{
    double XP_pred;
    double YP_pred;
    double vxP_pred;
    double vyP_pred;
    double x_pred;
    double y_pred;
    double theta_pred;
    double vx_pred;
    double vy_pred;
};

std::vector<std::vector<int>> cluster_cerchi_metricaB_max(const std::vector<Circle>& circles, double tau, int max_cluster_size);

static double cross2(const point& O, const point& A, const point& B);

std::vector<point> convexHull(std::vector<point> P);

Eigen::Vector2d polygonAreaCentroid(const std::vector<point>& poly);

static void dedupKeepMaxParam(std::vector<P>& pts);

static std::vector<P> convexHull(std::vector<P> pts);

//static LongestEdgeResult bestGapOnHull(const std::vector<P>& H);

static std::vector<LongestEdgeResult> gapsOnHullAboveDelta(const std::vector<P>& H, double delta);

static inline double norm2(double x, double y);

static double recomputeRadiusFromCentroid(const obstacle& o);

static void predictTracks(std::vector<Track>& tracks, const AssocParams& p);

static distMatrix_t buildCostMatrix(const std::vector<Track>& tracks, const std::vector<obstacle>& dets);

static assignments_t solveHungarian(const distMatrix_t& Cost, size_t N, size_t M);

static void applyGating(assignments_t& assignment, const distMatrix_t& Cost, size_t N, double gate_dist);

static void updateAssignedTracks(std::vector<Track>& tracks, const std::vector<obstacle>& dets, const assignments_t& assignment, const AssocParams& p);

static void markAndDeleteLost(std::vector<Track>& tracks, const assignments_t& assignment, int max_missed);

static std::vector<bool> computeUsedDetections(const assignments_t& assignment, size_t M);

static void spawnNewTracks(std::vector<Track>& tracks, const std::vector<obstacle>& dets, const std::vector<bool>& used, int& next_id);

void associateAndTrack(std::vector<Track>& tracks, int& next_id, const std::vector<obstacle>& detections, const AssocParams& p);

//std::vector<std::vector<int>> split_cluster_branches_mst(const std::vector<Circle>& circles, const std::vector<int>& cluster_idx, double tau, bool keep_junctions_as_own_cluster, int min_branch_size, double theta_cut_deg);

//std::vector<std::vector<int>> recluster_split_branches(const std::vector<Circle>& circles, const std::vector<std::vector<int>>& dsu_clusters, double tau, bool keep_junctions_as_own_cluster, int min_branch_size, double theta_cut_deg);

struct ClosestOnPolygonResult {
  Eigen::Vector2d closest_point; // Q
  int edge_index;                // i = edge (vi -> v(i+1))
  double t;                      // parametro su segmento [0,1]
  double dist2;                  // distanza al quadrato
};

/**
 * @class MPC_diffDrive_fblin
 * @brief MPC regulator/trajectory tracking controller with feedback linearization for differential drive robots
 */
class MPC_diffDrive_fblin {

public:
    /**
     * @brief Constructor for MPC_diffDrive_fblin
     */
    MPC_diffDrive_fblin();

    /**
     * @brief Destructor for MPC_diffDrive_fblin
     */
    ~MPC_diffDrive_fblin();

    /**
     * @brief Set robot model parameters
     * @param wheelVelMax maximum wheel angular velocity
     * @param wheelVelMin minimum wheel angular velocity
     * @param wheelRadius wheel nominal radius
     * @param track distance between wheel contact points
     */


     struct ConstraintLine {
  double hx, hy, l;
  int group; // quale cluster/ostacolo
  int k;     // quale step di MPC
};

using ConstraintsDebugCallback =
  std::function<void(const std::vector<ConstraintLine>&)>;


void setConstraintsDebugCallback(ConstraintsDebugCallback cb) {
  _constraints_cb = std::move(cb);
}




    void set_robotParams(double wheelVelMax, double wheelVelMin, double wheelRadius, double track);

    /**
     * @brief Set MPC controller parameters
     * @param samplingTime controller sampling time (must be an integer multiple of feedback linearization sampling time)
     * @param predictionHorizon length of the prediction horizon (number of samples)
     * @param q diagonal element of Q matrix
     * @param r diagonal element of R matrix
     * @param variableLB lower bounds of optimization variables
     * @param variableUB upper bounds of optimization variables
     */
    void set_MPCparams(double samplingTime, int predictionHorizon, double q, double r, double a_max, double sp, double sa, double r_red, double T_gp, double sensorRange, double e, int N_samples);
    //void set_MPCparams(double samplingTime, int predictionHorizon, double q, double r, double a_max, double sp, int n_obs, double sa, double r_c, double r_red, double T_gp, double sensorRange, double e, int N_samples);

    void set_n_obs(int n_obs);
    void set_r_c(const double& r_c);

    /**
     * @brief Set feedback linearization parameters
     * @param samplingTime feedback linearization sampling time
     * @param pointPdistance distance of point P with respect to the wheel axis
     */
    void set_FBLINparams(double samplingTime, double pointPdistance);

    /**
     * @brief Initialize MPC controller
     * @return          false in case of an error, true otherwise
     */
    bool initialize(const Eigen::Vector2d& final_ref_pos);

    /**
     * @brief Executes MPC controller computing control values (x/y velocity of point P)
     * @return          false in case of an error, true otherwise
     */
    bool executeMPCcontroller();

    /**
     * @brief Executes the feedback linearizing controller computing robot velocities from x/y velocity of point P
     * @return          false in case of an error, true otherwise
     */
    bool executeLinearizationController();

    /**
     * @brief Set the actual robot pose
     * @param x actual robot x-position
     * @param y actual robot y-position
     * @param yaw actual robot orientation
     */
    void set_actualRobotState(double x, double y, double yaw);

    /**
     *  Sampling of pedestrians velocities and positions*/

    void sample_Pedestrian_and_robot_states(const Eigen::VectorX<obstacle>& obstacles, double x_robot, double y_robot, double theta_robot);//, double sim_time, Eigen::MatrixXd u_prev, int count);

    void sample_Pedestrian_and_robot_states_complete(Eigen::VectorX<obstacle>& obstacles);//, double sim_time);

    void extractTraj(Eigen::MatrixXd ref, int count, bool trajectory_flag, bool replanTrigger);
    /**
     * @brief Set the reference robot state for a regulation problem (constant reference along the prediction horizon)
     * @param x reference robot x-position
     * @param y reference robot y-position
     * @param yaw reference robot orientation
     */
    void set_referenceRobotState(double x, double y, double yaw);

    /**
     * @brief Set the reference robot state along the prediction horizon for a trajectory tracking problem
     * @param refRobotState vector of reference robot states [x(k), y(k), yaw(k), ..., x(k+N), y(k+N), yaw(k+N)]
     * where N is the length of the prediction horizon
     */
    void set_referenceRobotState(const Eigen::VectorXd &refRobotState);

    /**
     * @brief Retrieve the actual MPC state (point P position)
     * @param xP point P x-position
     * @param yP point P y-position
     */
    void get_actualMPCstate(double& xP, double& yP);

    /**
     * @brief Retrieve the reference MPC state (point P position) along the prediction horizon
     * @param refState vector of reference states (N+1 elements, where N is the length of the prediction horizon)
     */
    void get_referenceMPCstate(Eigen::VectorXd &refState);

    /**
     * @brief Retrieve the actual MPC control vector (velocity of point P)
     * @param vPx x-component of point P velocity
     * @param vPy y-component of point P velocity
     */
    void get_actualMPCControl(double& vPx, double& vPy);

    /**
     * @brief Retrieve the actual robot control vector (linear and angular robot velocity)
     * @param linVelocity linear velocity of the robot
     * @param angVelocity angular velocity of the robot
     */
    void get_actualControl(double &linVelocity, double &angVelocity);

    /**
     * @brief Callbacks used by the class to print error/info/debug messages
     * The callback functions should have the following prototype
     * void debugMsgCallback(const std::string &message)
     * void infoMsgCallback(const std::string &message)
     * void errorMsgCallback(const std::string &message)
     * @param callback a pointer to the callback function that prints the message
     */

    void get_v_robot(double &vx, double &vy, double &v_long);

    void set_ErrorMsgCallback(ErrorMsgCallback callback);
    void set_InfoMsgCallback(InfoMsgCallback callback);
    void set_DebugMsgCallback(DebugMsgCallback callback);

private:
ConstraintsDebugCallback _constraints_cb;

    // MPC parameters
    int _N, _n_obs, _N_samples, _n_groups;
    double _MPC_Ts;
    double _q, _r, _k, _p, _s, _a_max, _sp, _sa, _r_c, _r_red, _T_gp, _sensorRange, _e, _dsl_a;
    bool _MPCparamsInitialized;
    bool _is_started, _is_started_complete;
    int _n_fail;
    bool _is_closed;
    bool _is_blocked_inside_hull = false;
    bool _is_blocked = false;

    Eigen::MatrixXd _plant_A, _plant_B;
    std::vector<double> _lowerBound, _upperBound;

    Eigen::MatrixXd _Acal, _Bcal;
    Eigen::MatrixXd _Qcal, _Rcal;
    Eigen::MatrixXd _H, _H_no_sl;
    Eigen::VectorXd _f, _f_no_sl;
    Eigen::MatrixXd _Ain_vel;
    Eigen::VectorXd _Bin_vel;
    Eigen::VectorXd _x_robot_samples, _y_robot_samples, _vx_robot_samples, _vy_robot_samples, _x_near, _y_near, _xP_robot_samples, _yP_robot_samples, _theta_robot_samples;
    Eigen::MatrixXd _x_ped_samples, _y_ped_samples, _vx_ped_samples, _vy_ped_samples;
    Eigen::VectorXi _labels, _id_samples, _id_samples_obstacles;
    Eigen::VectorXd _XP_predicted, _YP_predicted, _v_long;
    Eigen::MatrixXd _A_dvel;
    Eigen::VectorXd _b_dvel;
    Eigen::MatrixXd _A_obs;
    Eigen::VectorXd _b_obs;
    Eigen::MatrixXd _Asl;
    Eigen::VectorXd _bsl;
    Eigen::MatrixXd _A_tot;
    Eigen::VectorXd _B_tot;
    Eigen::MatrixX<obstacle> _obstacles_samples;
    Eigen::Vector2d _ref_temp;
    Eigen::MatrixXd _optimVect_no_sl_samples;

    int _kt;

    double _actXP, _actYP, _actX, _actY, _actYaw, _prevXP, _prevYP;
    Eigen::VectorXd _predictRobotState, _refMPCstate, _optimVect, _optimVect_no_sl;
    bool _init_P;

    GUROBIsolver *_solver;

    // Feedback linearization parameters
    double _fblin_Ts, _Pdist;

    fblin_unicycle *_fblinController, *_fblinSimController;
    bool _FBLINparamsInitialized;

    // Robot parameters
    double _wheelVelMax, _wheelVelMin;
    double _track, _wheelRadius;
    bool _robotParamsInitialized;

    // Controller variables
    bool _controllerInitialized;
    double _linearVelocity, _angularVelocity;
    double _vx_robot = 0.0, _vy_robot = 0.0, _v_long_robot = 0.0;

    // Controller constraints
    std::vector<GRBConstr> totConstrain;

    // Message handler function pointers
    DebugMsgCallback _debug;
    InfoMsgCallback _info;
    ErrorMsgCallback _error;

    // Private member functions
    void compute_AcalMatrix();
    void compute_BcalMatrix();
    void compute_QcalMatrix();
    void compute_RcalMatrix();
    void compute_objectiveMatrix();
    void compute_wheelVelocityConstraint();
    void prediction();
    pred_step prediction_step(double prevXP, double prevYP, double actX, double actY, double actYaw, double actXP, double actYP, Eigen::VectorXd optimVect_no_sl_samples, int k);
    void compute_wheelAccelerationConstraint();
    void compute_obstacle_avoidance_constraints();
    tangent findTangentPoints_mod_convex_hull(Eigen::Vector2d robot_pos, Eigen::MatrixXd ped_pos, Eigen::MatrixXd v_ped, Eigen::MatrixXd relative_velocity, double r_c, double r_red, Eigen::VectorXi idx, int i_group, int Na, Eigen::VectorX<obstacle> obstacles, Eigen::Vector2d P_pos);
    personal_space personalSpaceFunction_mod(Eigen::RowVector2d v_ped, Eigen::Vector2d relative_velocity, Eigen::Vector2d robot_pos, Eigen::Vector2d ped_pos, double r_c);
    tangent findTangentPoints_mod(const Eigen::Vector2d& robot_pos, const Eigen::Vector2d& ped_pos, Eigen::Vector2d PSf, double theta_heading, double r_red, double r_c);
    tangent_local findTangentPoints_local_mod(Eigen::Vector2d robot_pos_local, double a, double b);
    constraint constraintCoefficients_mod_mod_convex_hull(Eigen::Vector2d robot_pos, Eigen::MatrixXd ped_pos, Eigen::MatrixXd v_ped, Eigen::MatrixXd relative_velocity, double r_c, double r_red, Eigen::VectorXi idx, int i_group, Eigen::VectorX<obstacle> obstacles, int Na);
    double safetyDistance_mod_convex_hull(double r_red, double Ts, std::vector<Eigen::MatrixXd> v_ped_samples, std::vector<Eigen::MatrixXd> p_ped_samples, Eigen::MatrixXd p_robot_samples, Eigen::MatrixXd v_robot_samples, int i, double r_c, double epsilon, Eigen::VectorXi idx, int i_group, Eigen::MatrixX<obstacle> obstacles_samples, int Na, Eigen::MatrixXd p_P_robot_samples, Eigen::VectorXd yaw_robot, Eigen::MatrixXd optimVect_no_sl_samples);
    tangent findTangentPoints_mod_objects(const Eigen::Vector2d& robot_pos, double r_red, double r_c, const obstacle& obstacle);
    double compute_s_objects(const Eigen::RowVector2d& n, double r_red, double r_c, const obstacle& obstacle, const Eigen::Vector2d& pos);
    double compute_m_objects(const Eigen::RowVector2d& n, double r_red, double r_c, const obstacle& obstacle, const Eigen::Vector2d& pos);
    bool is_inside_convex_hull(Eigen::Vector2d robot_pos, Eigen::MatrixXd ped_pos, Eigen::MatrixXd v_ped, Eigen::MatrixXd relative_velocity, double r_c, double r_red, Eigen::VectorXi idx, int i_group, int Na, Eigen::VectorX<obstacle> obstacles, Eigen::Vector2d final_ref_pos, const double& tau);//, bool& is_blocked);
    Eigen::Vector2d compute_new_ref(Eigen::Vector2d robot_pos, Eigen::MatrixXd ped_pos, Eigen::MatrixXd v_ped, Eigen::MatrixXd relative_velocity, double r_c, double r_red, Eigen::VectorXi idx, int i_group, Eigen::VectorX<obstacle> obstacles, Eigen::Vector2d final_ref_pos, bool& is_closed);
    void compute_slack_variables_constraints();
    void clustering(const Eigen::Vector2d& final_ref_pos);

    void saveMatrixToFile(std::string fileName, Eigen::MatrixXd matrix);

    void errorMsg(const std::string &message);
    void infoMsg(const std::string &message);
    void debugMsg(const std::string &message);
};

#endif //MPC_DIFFDRIVE_FBLIN_H
