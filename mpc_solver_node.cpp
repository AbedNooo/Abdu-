#include "rclcpp/rclcpp.hpp"
#include "nav_msgs/msg/odometry.hpp"  // Changed from pose_stamped to odometry
#include "nav_msgs/msg/path.hpp"
#include "nav2_mpc_fblin_controller/msg/control_inputs.hpp"
#include "controller/MPC_diffDrive_fblin.h"
#include "tf2/utils.h"
#include "tf2_geometry_msgs/tf2_geometry_msgs.hpp"
#include "obstacle_detector/msg/obstacles.hpp"
#include "std_msgs/msg/string.hpp"
#include <Eigen/Dense>
#include <tf2_ros/buffer.h>
#include <tf2_ros/transform_listener.h>
#include "geometry_msgs/msg/vector3_stamped.hpp"
#include "visualization_msgs/msg/marker.hpp"
#include "geometry_msgs/msg/pose_array.hpp"
#include <visualization_msgs/msg/marker_array.hpp>

class MPCSolverNode : public rclcpp::Node {
public:
    MPCSolverNode() : Node("mpc_solver_node"), tf_buffer_(this->get_clock()), tf_listener_(tf_buffer_) {
        /*// Parameters declaration
        declare_parameter("prediction_horizon", 5);
        declare_parameter("q", 2.0);
        declare_parameter("r", 0.5);
        declare_parameter("p_dist", 0.25);
        declare_parameter("w_min", -0.5);
        declare_parameter("w_max", 0.5);
        declare_parameter("wheel_radius", 0.175);
        declare_parameter("track_width", 0.53);
        declare_parameter("mpc_frequency", 5.0);
        
        RCLCPP_INFO(this->get_logger(), "MPC Solver Node initialized with frequency: %f", 
                   get_parameter("mpc_frequency").as_double());
        
        prediction_horizon_ = get_parameter("prediction_horizon").as_int();*/

        marker_pub_ = this->create_publisher<visualization_msgs::msg::Marker>(
  "mpc_constraints", 10);

        
        // Changed to Odometry subscriber
        odom_sub_ = create_subscription<nav_msgs::msg::Odometry>(
            "/diff_drive_controller/odom", 10,  // Increased queue size
            std::bind(&MPCSolverNode::odomCallback, this, std::placeholders::_1));
            
        ref_traj_sub_ = create_subscription<nav_msgs::msg::Path>(
            "reference_path", 1,
            std::bind(&MPCSolverNode::trajCallback, this, std::placeholders::_1));

        control_pub_ = create_publisher<nav2_mpc_fblin_controller::msg::ControlInputs>(
            "mpc_control_outputs", 1);

        subscription_ = this->create_subscription<obstacle_detector::msg::Obstacles>(
      "obstacles",
      10,
      [this](obstacle_detector::msg::Obstacles::SharedPtr msg) {

        std::lock_guard<std::mutex> lock(mutex_);

        n_obs = msg->circles.size() + msg->segments.size();

        obstacle zero ={};

        obstacles.resize(n_obs);

        for (int k = 0; k < n_obs; ++k){ 
        obstacles(k) = zero;
        }

        int i = 0;

    for (const auto& circle_index : msg->circles) {
    obstacle o;

    o.x = circle_index.center.x;
    o.y = circle_index.center.y;

    o.x1 = circle_index.center.x;
    o.y1 = circle_index.center.y;

    o.x2 = circle_index.center.x;
    o.y2 = circle_index.center.y;

    o.vx = circle_index.velocity.x;
    o.vy = circle_index.velocity.y;

    o.radius = circle_index.radius;

    o.is_person = false;
    
    o.is_segment = false;

    o.id = circle_index.uid;

    std::vector<point> P;

    P.reserve(_k);

    for (int j = 0; j < _k; ++j){
      
      point p;

      double theta = ((2*M_PI)/_k)*j;

      p.x = circle_index.center.x + circle_index.radius*std::cos(theta);
      p.y = circle_index.center.y + circle_index.radius*std::sin(theta);
      p.z = 0.0;

      P.push_back(p);

    }

    o.polygon = P;

    obstacles(i) = o;

    i++;
  }

  for (const auto& segment_index : msg->segments) {
    obstacle o;

    o.x = (segment_index.first_point.x + segment_index.last_point.x)/2;
    o.y = (segment_index.first_point.y + segment_index.last_point.y)/2;

    o.x1 = segment_index.first_point.x;
    o.y1 = segment_index.first_point.y;

    o.x2 = segment_index.last_point.x;
    o.y2 = segment_index.last_point.y;

    o.vx = (segment_index.first_velocity.x + segment_index.last_velocity.x)/2;
    o.vy = (segment_index.first_velocity.y + segment_index.last_velocity.y)/2;

    o.radius = std::hypot((segment_index.first_point.x - segment_index.last_point.x), (segment_index.first_point.y - segment_index.last_point.y))/2;

    o.is_person = false;
    
    o.is_segment = true;

    o.id = segment_index.uid;

    std::vector<point> P;

    P.reserve(4);

    Eigen::Vector2d n = Eigen::Vector2d::Zero();
    n << (segment_index.last_point.y - segment_index.first_point.y), -(segment_index.last_point.x - segment_index.first_point.x);
    Eigen::Vector2d l = Eigen::Vector2d::Zero();
    l << (segment_index.last_point.x - segment_index.first_point.x), (segment_index.last_point.y - segment_index.first_point.y);
    double n_norm = n.norm() + 1e-6;
    double l_norm = l.norm() + 1e-6;

    point p;

    double delta = 3.0e-2;

    p.x = segment_index.first_point.x + (n(0)/n_norm)*delta;
    p.y = segment_index.first_point.y + (n(1)/n_norm)*delta;
    p.x = p.x - (l(0)/l_norm)*delta;
    p.y = p.y - (l(1)/l_norm)*delta;
    p.z = 0.0;

    P.push_back(p);

    p.x = segment_index.last_point.x + (n(0)/n_norm)*delta;
    p.y = segment_index.last_point.y + (n(1)/n_norm)*delta;
    p.x = p.x + (l(0)/l_norm)*delta;
    p.y = p.y + (l(1)/l_norm)*delta;
    p.z = 0.0;

    P.push_back(p);

    p.x = segment_index.last_point.x - (n(0)/n_norm)*delta;
    p.y = segment_index.last_point.y - (n(1)/n_norm)*delta;
    p.x = p.x + (l(0)/l_norm)*delta;
    p.y = p.y + (l(1)/l_norm)*delta;
    p.z = 0.0;

    P.push_back(p);

    p.x = segment_index.first_point.x - (n(0)/n_norm)*delta;
    p.y = segment_index.first_point.y - (n(1)/n_norm)*delta;
    p.x = p.x - (l(0)/l_norm)*delta;
    p.y = p.y - (l(1)/l_norm)*delta;
    p.z = 0.0;

    P.push_back(p);    

    o.polygon = convexHull(P);

    obstacles(i) = o;

    i++;
  }

  });

  lidar_poses_ = this->create_subscription<geometry_msgs::msg::PoseArray>(
  "/detected_object_pose",
  10,
  [this](const geometry_msgs::msg::PoseArray::SharedPtr msg)
  {
    std::lock_guard<std::mutex> lock(mutex_);

    // msg->poses è un std::vector<geometry_msgs::msg::Pose>
    // msg->header.frame_id = lidar_frame (nel nodo che hai)

    msg_lidar_ = *msg;
    lidar_received_ = true;


  }
);

ellipse_pub_ = this->create_publisher<visualization_msgs::msg::MarkerArray>(
  "mpc_personal_space", 10);

        // Initialize MPC solver
        mpc_solver_ = std::make_unique<MPC_diffDrive_fblin>();

        mpc_solver_->setConstraintsDebugCallback(
  [this](const std::vector<ConstraintLine>& lines) {
    last_lines_ = lines;
    lines_ready_ = true;
  }
);

        
        // Set parameters
        /*mpc_solver_->set_robotParams(
            get_parameter("w_max").as_double(),
            get_parameter("w_min").as_double(),
            get_parameter("wheel_radius").as_double(),
            get_parameter("track_width").as_double());
            
        mpc_solver_->set_MPCparams(
            1.0/get_parameter("mpc_frequency").as_double(),
            get_parameter("prediction_horizon").as_int(),
            get_parameter("q").as_double(),
            get_parameter("r").as_double());
            
        mpc_solver_->set_FBLINparams(
            0.01,
            get_parameter("p_dist").as_double());*/

        mpc_solver_->set_ErrorMsgCallback([this](const std::string& message) {
  this->errorMsgCallback(message);
});

mpc_solver_->set_DebugMsgCallback([this](const std::string& message) {
  this->debugMsgCallback(message);
});

mpc_solver_->set_InfoMsgCallback([this](const std::string& message) {
  this->infoMsgCallback(message);
});

        mpc_solver_->set_MPCparams(Ts_MPC, N, q, r, a_max, sp, sa, r_red, T_gp, sensorRange, e, N_samples);
        mpc_solver_->set_FBLINparams(Ts_MPC, p_dist);
        mpc_solver_->set_robotParams(wMax, wMin, R, d);
            
        /*if (!mpc_solver_->initialize()) {
            RCLCPP_ERROR(get_logger(), "MPC solver initialization failed!");
            rclcpp::shutdown();
        }*/

        timer_ = create_wall_timer(
            std::chrono::milliseconds(static_cast<int>(Ts_MPC * 1000.0)),
            std::bind(&MPCSolverNode::solveMPC, this));
    }

private:
    using ConstraintLine = MPC_diffDrive_fblin::ConstraintLine;


nav_msgs::msg::Odometry toMapPoseAndTwist(
  const nav_msgs::msg::Odometry & odom_msg,
  tf2_ros::Buffer & tf_buffer,
  rclcpp::Logger logger)
{
  nav_msgs::msg::Odometry out = odom_msg;

  // Pose: odom.header.frame_id -> map
  geometry_msgs::msg::PoseStamped pose_in, pose_out;
  pose_in.header = odom_msg.header;  // tipico: frame_id="odom"
  pose_in.pose   = odom_msg.pose.pose;

  // Twist: 2 vettori nel child_frame_id (tipico: "base_link")
  geometry_msgs::msg::Vector3Stamped lin_in, lin_out;
  lin_in.header.stamp = odom_msg.header.stamp;
  lin_in.header.frame_id = odom_msg.child_frame_id;
  lin_in.vector = odom_msg.twist.twist.linear;

  geometry_msgs::msg::Vector3Stamped ang_in, ang_out;
  ang_in.header.stamp = odom_msg.header.stamp;
  ang_in.header.frame_id = odom_msg.child_frame_id;
  ang_in.vector = odom_msg.twist.twist.angular;

  try {
    // Pose in map
    pose_out = tf_buffer.transform(pose_in, "map", tf2::durationFromSec(0.1));

    // Vettori velocità in map
    lin_out = tf_buffer.transform(lin_in, "map", tf2::durationFromSec(0.1));
    ang_out = tf_buffer.transform(ang_in, "map", tf2::durationFromSec(0.1));

    out.header.frame_id = "map";
    out.header.stamp    = pose_out.header.stamp;
    out.pose.pose       = pose_out.pose;

    // Ora la twist è espressa in "map"
    out.twist.twist.linear  = lin_out.vector;
    out.twist.twist.angular = ang_out.vector;

  } catch (const tf2::TransformException & ex) {
    RCLCPP_WARN(logger, "TF odom/twist -> map failed: %s", ex.what());
  }

  return out;
}

Eigen::VectorX<obstacle> removeIndex(const Eigen::VectorX<obstacle>& v, int k)
{
  const int n = v.size();
  if (k < 0 || k >= n) return v;            // oppure throw
  Eigen::VectorX<obstacle> out(n - 1);

  out.head(k) = v.head(k);
  out.tail(n - k - 1) = v.tail(n - k - 1);
  return out;
}

void publishPersonalSpaces(const rclcpp::Time& stamp)
{
  visualization_msgs::msg::MarkerArray arr;

  // opzionale: pulisce i marker vecchi
  visualization_msgs::msg::Marker clear;
  clear.action = visualization_msgs::msg::Marker::DELETEALL;
  arr.markers.push_back(clear);

  const double r_red = 0.3;   // usa i TUOI valori reali
  double r_c_int;
  //const double delta = (r_c - r_red);

  if (v_long_robot_ >= 0.0){
    r_c_int = 0.8;
  }else{
    r_c_int = 1.3;
  }

  const double rx = current_odom_.pose.pose.position.x;
  const double ry = current_odom_.pose.pose.position.y;

  const double vrx = vx_robot_;//current_odom_.twist.twist.linear.x; // in map (tu già trasformi odom->map)
  const double vry = vy_robot_;//current_odom_.twist.twist.linear.y;

  const int N = 60;
  const double TWO_PI = 6.283185307179586;

  int fallback_id = 0;

  for (int i = 0; i < obstacles_vis.size(); ++i)
  {
    if (!obstacles_vis(i).is_person) continue;

    const double px = obstacles_vis(i).x;
    const double py = obstacles_vis(i).y;

    const double vpx = obstacles_vis(i).vx;
    const double vpy = obstacles_vis(i).vy;

    // replica personalSpaceFunction_mod()
    const double relx = vrx - vpx;
    const double rely = vry - vpy;

    const double c = -((px - rx) * (-relx) + (py - ry) * (-rely));

    const double vped_norm = std::hypot(vpx, vpy);
    const double rel_norm  = std::hypot(relx, rely);

    double b0 = r_c_int;
    double a0 = r_c_int;

    if (c > 0.0) {
      if (vped_norm > 1e-6) { b0 = std::max(r_c_int, 3.0 * rel_norm); a0 = 1.5 * b0; }
      else                  { b0 = std::max(r_c_int, 3.0 * rel_norm); a0 = b0; }
    } else {
      if (vped_norm > 1e-6) { b0 = std::max(r_c_int, 3.0 * vped_norm); a0 = 1.5 * b0; }
      else                  { b0 = r_c_int; a0 = b0; }
    }

    // se vuoi visualizzare l’ellisse “enlarged” come fai nel resto del codice:
    const double a = a0;// + delta;
    const double b = b0;// + delta;

    const double yaw = std::atan2(vpy, vpx);

    visualization_msgs::msg::Marker mk;
    mk.header.frame_id = "map";   // oppure il frame che usi in RViz
    mk.header.stamp = stamp;
    mk.ns = "personal_space";
    mk.id = (obstacles_vis(i).id >= 0) ? obstacles_vis(i).id : (100000 + fallback_id++);
    mk.type = visualization_msgs::msg::Marker::LINE_STRIP;
    mk.action = visualization_msgs::msg::Marker::ADD;

    mk.scale.x = 0.03;     // spessore linea
    mk.color.r = 1.0;
    mk.color.g = 0.2;
    mk.color.b = 0.2;
    mk.color.a = 0.9;

    mk.lifetime = rclcpp::Duration::from_seconds(0.25);

    mk.points.reserve(N + 1);
    for (int k = 0; k <= N; ++k) {
        const double t  = TWO_PI * double(k) / double(N);
        const double ct = std::cos(t);
        const double st = std::sin(t);

        // metà "front" ellisse, metà "back" cerchio
        double lx, ly;
        if (ct >= 0.0) {
            // front: ellisse
            lx = a * ct;
            ly = b * st;
        } else {
            // back: cerchio (r = b)
            lx = b * ct;
            ly = b * st;
        }

    geometry_msgs::msg::Point p;
    p.x = px + lx * std::cos(yaw) - ly * std::sin(yaw);
    p.y = py + lx * std::sin(yaw) + ly * std::cos(yaw);
    p.z = 0.02;
    mk.points.push_back(p);
}

    arr.markers.push_back(std::move(mk));
  }

  ellipse_pub_->publish(arr);
}

    // Changed to odom callback
    void odomCallback(const nav_msgs::msg::Odometry::SharedPtr msg) {
        std::lock_guard<std::mutex> lock(mutex_);
        //current_odom_ = *msg;
        current_odom_ = toMapPoseAndTwist(*msg, tf_buffer_, this->get_logger());
        odom_received_ = true;

        if (msg->twist.twist.linear.x >= 0){
            r_c = 0.8;
        }else{
            r_c = 1.3;
        }
        
        // Debug output
        RCLCPP_DEBUG(get_logger(), "Received odom - x: %.2f, y: %.2f, theta: %.2f",
                    msg->pose.pose.position.x,
                    msg->pose.pose.position.y,
                    tf2::getYaw(msg->pose.pose.orientation));
    }

    void trajCallback(const nav_msgs::msg::Path::SharedPtr msg) {
        std::lock_guard<std::mutex> lock(mutex_);
        reference_trajectory_ = *msg;

        ref = Eigen::MatrixXd::Zero(reference_trajectory_.poses.size(), 2);

        for (size_t k = 0; k < reference_trajectory_.poses.size(); k++) {
           
            ref(k, 0) = reference_trajectory_.poses[k].pose.position.x;
            ref(k, 1) = reference_trajectory_.poses[k].pose.position.y;
        }

        traj_received_ = true;
    }

    void solveMPC() {
        std::lock_guard<std::mutex> lock(mutex_);
    
        if (!odom_received_) {
            RCLCPP_WARN_THROTTLE(get_logger(), *get_clock(), 2000, 
                               "No odometry data received yet.");
            return;
        }
    
        if (reference_trajectory_.poses.empty()) {
            RCLCPP_WARN_THROTTLE(get_logger(), *get_clock(), 2000, 
                               "Reference trajectory is empty.");
            return;
        }
    
        /*// Use odometry pose data
        mpc_solver_->set_actualRobotState(
            current_odom_.pose.pose.position.x,
            current_odom_.pose.pose.position.y,
            tf2::getYaw(current_odom_.pose.pose.orientation));
    
        Eigen::VectorXd ref_state(3 * (prediction_horizon_ + 1));
        for (int i = 0; i < prediction_horizon_ + 1; ++i) {
            size_t idx = std::min(static_cast<size_t>(i), reference_trajectory_.poses.size() - 1);
            ref_state(3 * i) = reference_trajectory_.poses[idx].pose.position.x;
            ref_state(3 * i + 1) = reference_trajectory_.poses[idx].pose.position.y;
            ref_state(3 * i + 2) = tf2::getYaw(reference_trajectory_.poses[idx].pose.orientation);
        }
    
        mpc_solver_->set_referenceRobotState(ref_state);*/

        if (lidar_received_){
        geometry_msgs::msg::PoseArray pose_array_out;

        pose_array_out.header = msg_lidar_.header;
        pose_array_out.header.frame_id = "map";
        pose_array_out.poses.reserve(msg_lidar_.poses.size());

        for (const auto& pose : msg_lidar_.poses){
            geometry_msgs::msg::PoseStamped pose_in, pose_out;
            pose_in.header = msg_lidar_.header;
            pose_in.pose   = pose;

        try {

            pose_out = tf_buffer_.transform(pose_in, "map", tf2::durationFromSec(0.1));

            pose_array_out.poses.push_back(pose_out.pose);
    
        } catch (const tf2::TransformException & ex) {
            RCLCPP_WARN(get_logger(), "TF lidar_3D_frame -> map failed: %s", ex.what());
        }

        }

        if (obstacles.size() > 0){
        for (const auto& pose_map : pose_array_out.poses){

        for (int i = 0; i < obstacles.size(); i++){

            if (std::hypot(obstacles(i).x - pose_map.position.x, obstacles(i).y - pose_map.position.y) < r_red*2){
                obstacles(i).is_person = true;
            }

            }
        }

        for (int i = 0; i < obstacles.size(); i++){

            if (obstacles(i).is_person == true){
                for (int j = obstacles.size() - 1; j > i; j--){
                    //if (obstacles(j).is_person == true && std::hypot(obstacles(j).x - obstacles(i).x, obstacles(j).y - obstacles(i).y) < r_red){
                    if (std::hypot(obstacles(j).x - obstacles(i).x, obstacles(j).y - obstacles(i).y) < r_red*2){    
                        obstacles = removeIndex(obstacles, j);
                    }
                }            
            }

            for (int j = i - 1; j >= 0; --j){
                    //if (obstacles(j).is_person == true && std::hypot(obstacles(j).x - obstacles(i).x, obstacles(j).y - obstacles(i).y) < r_red){
                    if (std::hypot(obstacles(j).x - obstacles(i).x, obstacles(j).y - obstacles(i).y) < r_red*2){    
                        obstacles = removeIndex(obstacles, j);
                        --i;
                    }
                }

            }
        }
    lidar_received_ = false;
    }

        mpc_solver_->set_n_obs(obstacles.size());

        mpc_solver_->set_r_c(r_c);

        mpc_solver_->sample_Pedestrian_and_robot_states(obstacles, current_odom_.pose.pose.position.x, current_odom_.pose.pose.position.y, tf2::getYaw(current_odom_.pose.pose.orientation));//, time, u_prev, MPC_count);

        mpc_solver_->sample_Pedestrian_and_robot_states_complete(obstacles);//, double sim_time);//CAMBIARE!!!

        obstacles_vis = obstacles;
        
        if (traj_received_){
            mpc_solver_->extractTraj(ref, MPC_count, true, true);
            traj_received_ = false;
        }else{
            mpc_solver_->extractTraj(ref, MPC_count, true, false);
        } //CAMBIARE!!!

        mpc_solver_->initialize(ref.row(ref.rows() - 1).transpose());

        mpc_solver_->set_actualRobotState(current_odom_.pose.pose.position.x, current_odom_.pose.pose.position.y, tf2::getYaw(current_odom_.pose.pose.orientation));

        mpc_solver_->executeMPCcontroller();
            nav2_mpc_fblin_controller::msg::ControlInputs control_msg;
            mpc_solver_->get_actualMPCControl(control_msg.vpx, control_msg.vpy);
            control_pub_->publish(control_msg);
            
            // Debug output
            RCLCPP_DEBUG(get_logger(), "Published control - vpx: %.2f, vpy: %.2f",
                        control_msg.vpx, control_msg.vpy);
       
        //mpc_solver_->executeMPCcontroller();

        //mpc_solver_->get_actualMPCControl(vPx_act, vPy_act);

        //mpc_solver_->executeLinearizationController();

        mpc_solver_->get_v_robot(vx_robot_, vy_robot_, v_long_robot_);

        MPC_count++;

        auto stamp = this->now();

        if (lines_ready_) {
            publishConstraintLines(last_lines_, stamp);
            lines_ready_ = false;
        }

        publishPersonalSpaces(stamp);
    }

    inline void debugMsgCallback(const std::string &message) {
    std::cout << message << std::endl;
    }

    inline void infoMsgCallback(const std::string &message) {
    std::cout << message << std::endl;
    }

    inline void errorMsgCallback(const std::string &message) {
    std::cout << message << std::endl;
    }

    void publishConstraintLines(const std::vector<ConstraintLine>& lines, const rclcpp::Time& stamp)
{
  if (!marker_pub_) return;

  visualization_msgs::msg::Marker m;
  m.header.frame_id = "map";        // deve essere coerente con hx/hy/l
  m.header.stamp = stamp;//this->now();
  m.ns = "mpc_constraints";
  m.id = 0;
  m.type = visualization_msgs::msg::Marker::LINE_LIST;
  m.action = visualization_msgs::msg::Marker::ADD;
  m.pose.orientation.w = 1.0;

  m.scale.x = 0.03;
  m.color.a = 1.0;
  m.color.r = 1.0;

  m.lifetime = rclcpp::Duration::from_seconds(0.25);

  const double span = 5.0;

  for (const auto& c : lines) {
    const double n2 = c.hx*c.hx + c.hy*c.hy;
    if (n2 < 1e-12) continue;

    const double x0 = c.l * c.hx / n2;
    const double y0 = c.l * c.hy / n2;

    const double n = std::sqrt(n2);
    const double tx = -c.hy / n;
    const double ty =  c.hx / n;

    geometry_msgs::msg::Point p1, p2;
    p1.x = x0 + span*tx;  p1.y = y0 + span*ty;  p1.z = 0.05;
    p2.x = x0 - span*tx;  p2.y = y0 - span*ty;  p2.z = 0.05;

    m.points.push_back(p1);
    m.points.push_back(p2);
  }

  marker_pub_->publish(m);
}



    /*****************************************************/

    /************* Controller initialization *************/
    // MPC parameters
    const int N = 20;
    const double Ts_MPC = 0.2;
    const double q = 1.0;
    const double r = 1.0;
    const double a_max = 0.20;
    const double sp = 1e9;
    int n_obs = 0;
    const double sa = 1e3;
    const double r_red = 0.3;
    double r_c = 0.8;
    const double T_gp = 4.0;
    const double sensorRange = 10.0;
    const double e = 0.3;
    const int N_samples = 40;
    const int maxInfeasibleSolution = 2;

    // Feedback linearization parameters
    const double p_dist = 0.5;
    const double Ts_fblin = 0.01;

    // Robot parameters
    const double v_max = 0.55;
    const double R = 0.175;
    const double d = 0.565;
    const double wMax = v_max/R;
    const double wMin = -wMax;

    // Changed members to use odometry
    std::mutex mutex_;
    nav_msgs::msg::Odometry current_odom_;
    nav_msgs::msg::Path reference_trajectory_;
    int prediction_horizon_;
    bool odom_received_ = false;
    bool traj_received_ = false;
    bool lidar_received_ = false;

    double vPx_act = 0.0, vPy_act = 0.0, v_act = 0.0, w_act = 0.0;

    int MPC_count = 0;
    
    rclcpp::Subscription<nav_msgs::msg::Odometry>::SharedPtr odom_sub_;
    rclcpp::Subscription<nav_msgs::msg::Path>::SharedPtr ref_traj_sub_;
    rclcpp::Publisher<nav2_mpc_fblin_controller::msg::ControlInputs>::SharedPtr control_pub_;
    rclcpp::TimerBase::SharedPtr timer_;
    std::unique_ptr<MPC_diffDrive_fblin> mpc_solver_;
    rclcpp::Subscription<obstacle_detector::msg::Obstacles>::SharedPtr subscription_;
    Eigen::VectorX<obstacle> obstacles;
    Eigen::VectorX<obstacle> obstacles_vis;
    Eigen::MatrixXd ref;
    rclcpp::Subscription<geometry_msgs::msg::PoseArray>::SharedPtr lidar_poses_;

    tf2_ros::Buffer tf_buffer_;
    tf2_ros::TransformListener tf_listener_;

    rclcpp::Publisher<visualization_msgs::msg::Marker>::SharedPtr marker_pub_;
    rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr ellipse_pub_;

    geometry_msgs::msg::PoseArray msg_lidar_;

    std::vector<ConstraintLine> last_lines_;
    bool lines_ready_ = false;

    double vx_robot_; 
    double vy_robot_;
    double v_long_robot_;

    int _k = 10;



};

int main(int argc, char** argv) {
    rclcpp::init(argc, argv);
    rclcpp::spin(std::make_shared<MPCSolverNode>());
    rclcpp::shutdown();
    return 0;
}
