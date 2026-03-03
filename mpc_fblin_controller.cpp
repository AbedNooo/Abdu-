#include <string>
#include <vector>
#include <array>
#include <mutex>
#include <algorithm>

#include "nav2_mpc_fblin_controller/mpc_fblin_controller.hpp"
#include "nav2_core/exceptions.hpp"
#include "nav2_util/node_utils.hpp"
#include "nav2_util/geometry_utils.hpp"
#include "alglib/interpolation.h"
#include "nav2_costmap_2d/costmap_filters/filter_values.hpp"
#include "nav2_mpc_fblin_controller/msg/control_inputs.hpp"
#include "visualization_msgs/msg/marker.hpp"
#include <tf2/LinearMath/Quaternion.h>
#include <tf2_geometry_msgs/tf2_geometry_msgs.hpp>

using std::hypot;
using std::min;
using std::max;
using std::abs;
using nav2_util::declare_parameter_if_not_declared;
using nav2_util::geometry_utils::euclidean_distance;
using namespace nav2_costmap_2d;
using rcl_interfaces::msg::ParameterType;

namespace nav2_mpc_fblin_controller {

void handleDebugMessages(const std::string& message) {
    RCLCPP_DEBUG(rclcpp::get_logger("MPC_fblin_diffDrive"), "%s", message.c_str());
}

void handleInfoMessages(const std::string& message) {
    RCLCPP_INFO(rclcpp::get_logger("MPC_fblin_diffDrive"), "%s", message.c_str());
}

void handleErrorMessages(const std::string& message) {
    RCLCPP_ERROR(rclcpp::get_logger("MPC_fblin_diffDrive"), "%s", message.c_str());
}

MPC_fblin_diffDrive::MPC_fblin_diffDrive()
: node_(),
  tf_(nullptr),
  costmap_ros_(nullptr),
  costmap_(nullptr),
  logger_(rclcpp::get_logger("MPC_fblin_diffDrive")),
  clock_(nullptr),
  plugin_name_(""),
  p_dist_(0.5),
  fblin_frequency_(50.0),
  w_min_(-3.143),
  w_max_(3.143),
  wheel_radius_(0.175),
  track_width_(0.565),
  use_tracking_controller_(true),
  path_sampling_threshold_(1.0),
  next_goal_threshold_(0.1),
  max_linear_velocity_(0.55),
  max_angular_velocity_(1.0),
  MPC_frequency_(5.0),
  q_(1.0),
  r_(1.0),
  acc_lin_max_(0.2),
  prediction_horizon_(20),
  max_infeasible_sol_(5),
  path_idx_(0),
  path_time_(0),
  path_duration_(0.0),
  path_length_(0.0),
  mpc_sub_(nullptr),
  collision_checker_(nullptr),
  dyn_params_handler_(nullptr),
  next_goal_pub_(nullptr),
  reference_path_pub_(nullptr)
{
}

MPC_fblin_diffDrive::~MPC_fblin_diffDrive()
{
  // Clean up resources
}

void MPC_fblin_diffDrive::configure(
    const rclcpp_lifecycle::LifecycleNode::WeakPtr &parent,
    std::string name, std::shared_ptr<tf2_ros::Buffer> tf,
    std::shared_ptr<nav2_costmap_2d::Costmap2DROS> costmap_ros)
{
    auto node = parent.lock();
    node_ = parent;
    if (!node) {
        throw nav2_core::PlannerException("Unable to lock node!");
    }

    costmap_ros_ = costmap_ros;
    costmap_ = costmap_ros_->getCostmap();
    tf_ = tf;
    plugin_name_ = name;
    logger_ = node->get_logger();
    clock_ = node->get_clock();

    // Controller parameters
    declare_parameter_if_not_declared(node, "controller_frequency", rclcpp::ParameterValue(100.0));
    node->get_parameter("controller_frequency", fblin_frequency_);

    declare_parameter_if_not_declared(node, plugin_name_ + ".p_dist", rclcpp::ParameterValue(0.25));
    node->get_parameter(plugin_name_ + ".p_dist", p_dist_);

    declare_parameter_if_not_declared(node, plugin_name_ + ".w_min", rclcpp::ParameterValue(-0.1));
    node->get_parameter(plugin_name_ + ".w_min", w_min_);

    declare_parameter_if_not_declared(node, plugin_name_ + ".w_max", rclcpp::ParameterValue(0.1));
    node->get_parameter(plugin_name_ + ".w_max", w_max_);

    declare_parameter_if_not_declared(node, plugin_name_ + ".wheel_radius", rclcpp::ParameterValue(0.1));
    node->get_parameter(plugin_name_ + ".wheel_radius", wheel_radius_);

    declare_parameter_if_not_declared(node, plugin_name_ + ".track_width", rclcpp::ParameterValue(0.5));
    node->get_parameter(plugin_name_ + ".track_width", track_width_);

    declare_parameter_if_not_declared(node, plugin_name_ + ".use_tracking_controller", rclcpp::ParameterValue(false));
    node->get_parameter(plugin_name_ + ".use_tracking_controller", use_tracking_controller_);

    declare_parameter_if_not_declared(node, plugin_name_ + ".path_sampling_threshold", rclcpp::ParameterValue(1.0));
    node->get_parameter(plugin_name_ + ".path_sampling_threshold", path_sampling_threshold_);

    declare_parameter_if_not_declared(node, plugin_name_ + ".next_goal_threshold", rclcpp::ParameterValue(0.1));
    node->get_parameter(plugin_name_ + ".next_goal_threshold", next_goal_threshold_);

    declare_parameter_if_not_declared(node, plugin_name_ + ".max_linear_velocity", rclcpp::ParameterValue(0.5));
    node->get_parameter(plugin_name_ + ".max_linear_velocity", max_linear_velocity_);

    declare_parameter_if_not_declared(node, plugin_name_ + ".max_angular_velocity", rclcpp::ParameterValue(1.0));
    node->get_parameter(plugin_name_ + ".max_angular_velocity", max_angular_velocity_);

    // MPC parameters
    declare_parameter_if_not_declared(node, plugin_name_ + ".MPC_frequency", rclcpp::ParameterValue(10.0));
    node->get_parameter(plugin_name_ + ".MPC_frequency", MPC_frequency_);

    declare_parameter_if_not_declared(node, plugin_name_ + ".q", rclcpp::ParameterValue(1.0));
    node->get_parameter(plugin_name_ + ".q", q_);

    declare_parameter_if_not_declared(node, plugin_name_ + ".r", rclcpp::ParameterValue(0.1));
    node->get_parameter(plugin_name_ + ".r", r_);

    declare_parameter_if_not_declared(node, plugin_name_ + ".acc_lin_max", rclcpp::ParameterValue(0.5));
    node->get_parameter(plugin_name_ + ".acc_lin_max", acc_lin_max_);

    declare_parameter_if_not_declared(node, plugin_name_ + ".prediction_horizon", rclcpp::ParameterValue(10));
    node->get_parameter(plugin_name_ + ".prediction_horizon", prediction_horizon_);

    declare_parameter_if_not_declared(node, plugin_name_ + ".max_infeasible_sol", rclcpp::ParameterValue(5));
    node->get_parameter(plugin_name_ + ".max_infeasible_sol", max_infeasible_sol_);

    // Create publishers
    next_goal_pub_ = node->create_publisher<geometry_msgs::msg::PoseStamped>("next_goal_point", 1);
    reference_path_pub_ = node->create_publisher<nav_msgs::msg::Path>("reference_path", 1);

    // Create subscriber for MPC outputs
    mpc_sub_ = node->create_subscription<nav2_mpc_fblin_controller::msg::ControlInputs>(
        "mpc_control_outputs",
        rclcpp::SystemDefaultsQoS(),
        [this](const nav2_mpc_fblin_controller::msg::ControlInputs::SharedPtr msg) {
            std::lock_guard<std::mutex> lock(mpc_mutex_);
            current_mpc_output_ = {msg->vpx, msg->vpy};
        });

    // Initialize collision checker
    collision_checker_ = std::make_unique<nav2_costmap_2d::FootprintCollisionChecker<nav2_costmap_2d::Costmap2D*>>(costmap_);
    collision_checker_->setCostmap(costmap_);

    path_duration_ = path_length_ = 0.0;
    path_idx_ = path_time_ = 0;
}

void MPC_fblin_diffDrive::cleanup() {
    RCLCPP_INFO(logger_, "Cleaning up controller: %s", plugin_name_.c_str());
    next_goal_pub_.reset();
    reference_path_pub_.reset();
}

void MPC_fblin_diffDrive::activate() {
    RCLCPP_INFO(logger_, "Activating controller: %s", plugin_name_.c_str());
    next_goal_pub_->on_activate();
    reference_path_pub_->on_activate();
    auto node = node_.lock();
    if (node) {
        dyn_params_handler_ = node->add_on_set_parameters_callback(
            std::bind(&MPC_fblin_diffDrive::dynamicParametersCallback, this, std::placeholders::_1));
    }
}

void MPC_fblin_diffDrive::deactivate() {
    RCLCPP_INFO(logger_, "Deactivating controller: %s", plugin_name_.c_str());
    next_goal_pub_->on_deactivate();
    reference_path_pub_->on_deactivate();
    dyn_params_handler_.reset();
}

geometry_msgs::msg::TwistStamped MPC_fblin_diffDrive::computeVelocityCommands(
     const geometry_msgs::msg::PoseStamped &pose,
     const geometry_msgs::msg::Twist& /*speed*/, nav2_core::GoalChecker* /*goal_checker*/)
{
     std::lock_guard<std::mutex> lock_reinit(mutex_);
     //nav2_costmap_2d::Costmap2D *costmap = costmap_ros_->getCostmap();
     //std::unique_lock<nav2_costmap_2d::Costmap2D::mutex_t> lock(*(costmap->getMutex()));

     // Get latest MPC outputs
     double vpx, vpy;
     {
         std::lock_guard<std::mutex> lock(mpc_mutex_);
         vpx = current_mpc_output_[0];
         vpy = current_mpc_output_[1];
     }

     geometry_msgs::msg::PoseStamped pose_map;
    try {
  
        tf_->transform(pose, pose_map, "map");
    } catch (const tf2::TransformException & ex) {
        RCLCPP_WARN(logger_, "Cannot transform pose to map: %s", ex.what());
        return geometry_msgs::msg::TwistStamped();  // o cmd_vel a zero
    }

    double theta = tf2::getYaw(pose_map.pose.orientation);

     // TODO: Implement your control logic here using vpx and vpy
     double linear_vel = vpx*std::cos(theta) + vpy*std::sin(theta);
     double angular_vel = (vpy*std::cos(theta) - vpx*std::sin(theta))/p_dist_;

     // Apply speed limits
     linear_vel = std::clamp(linear_vel, -max_linear_velocity_, max_linear_velocity_);
     angular_vel = std::clamp(angular_vel, -max_angular_velocity_, max_angular_velocity_);

     // Populate and return message
     geometry_msgs::msg::TwistStamped cmd_vel;
     cmd_vel.header = pose.header;
     cmd_vel.header.stamp = clock_->now();
     cmd_vel.twist.linear.x = linear_vel;
     cmd_vel.twist.angular.z = angular_vel;
     return cmd_vel;
}

/*geometry_msgs::msg::TwistStamped MPC_fblin_diffDrive::computeVelocityCommands(
    const geometry_msgs::msg::PoseStamped &pose,
    const geometry_msgs::msg::Twist &velocity,
    nav2_core::GoalChecker *goal_checker)
{
    (void)velocity;
    (void)goal_checker;
    geometry_msgs::msg::TwistStamped cmd_vel;

    if (plan_.empty()) {
        RCLCPP_WARN(logger_, "No plan received.");
        return cmd_vel;
    }

    auto node = node_.lock();
    if (!node) {
        throw std::runtime_error("Failed to lock node");
    }

    // Get current robot pose
    double x = pose.pose.position.x;
    double y = pose.pose.position.y;
    double theta = tf2::getYaw(pose.pose.orientation);

    // Track next waypoint
    const auto &target_pose = plan_[path_idx_];
    double x_ref = target_pose.pose.position.x;
    double y_ref = target_pose.pose.position.y;

    // Error in world frame
    double ex = x_ref - x;
    double ey = y_ref - y;

    // Rotate error into body frame
    double ex_r = std::cos(theta) * ex + std::sin(theta) * ey;
    double ey_r = -std::sin(theta) * ex + std::cos(theta) * ey;

    // Feedback linearization
    double k1 = 1.0, k2 = 2.0;
    double v = k1 * ex_r;
    double omega = k2 * ey_r;

    // Set command
    cmd_vel.twist.linear.x = v;
    cmd_vel.twist.angular.z = omega;
    cmd_vel.header.stamp = node->now();
    cmd_vel.header.frame_id = "base_link";

    // Optional: publish goal marker
    geometry_msgs::msg::PoseStamped next_goal = target_pose;
    next_goal.header.stamp = node->now();
    next_goal_pub_->publish(next_goal);

    // Move to next waypoint if close enough
    double distance = std::hypot(ex, ey);
    if (distance < 0.1 && path_idx_ < plan_.size() - 1) {
        path_idx_++;
    }

    return cmd_vel;
}*/


void MPC_fblin_diffDrive::setPlan(const nav_msgs::msg::Path &path) {
    if (path.poses.empty()) {
        RCLCPP_WARN(logger_, "Received empty path");
        return;
    }

    plan_ = path.poses;
    path_idx_ = 0;

    std::lock_guard<std::mutex> lock(mutex_);

    nav_msgs::msg::Path path_P;
    
    if (use_tracking_controller_) {
        // Natural coordinate computation
        std::vector<double> path_x(path.poses.size(), 0.0);
        std::vector<double> path_y(path.poses.size(), 0.0);
        std::vector<double> path_yaw(path.poses.size(), 0.0);
        std::vector<double> path_s(path.poses.size(), 0.0);

        path_x[0] = path.poses[0].pose.position.x;
        path_y[0] = path.poses[0].pose.position.y;
        path_yaw[0] = tf2::getYaw(path.poses[0].pose.orientation);
        
        for (size_t k = 1; k < path.poses.size(); k++) {
            path_s[k] = path_s[k-1] + 
                       std::hypot(path.poses[k].pose.position.x - path.poses[k-1].pose.position.x,
                                  path.poses[k].pose.position.y - path.poses[k-1].pose.position.y);
            path_x[k] = path.poses[k].pose.position.x;
            path_y[k] = path.poses[k].pose.position.y;
            path_yaw[k] = tf2::getYaw(path.poses[k].pose.orientation);
        }

        alglib::spline1dinterpolant path_sp_x_, path_sp_y_, path_sp_yaw_;

        // Interpolation of the path components
        try {
            alglib::real_1d_array s, x, y, yaw;
            s.setcontent(path_s.size(), path_s.data());
            x.setcontent(path_x.size(), path_x.data());
            y.setcontent(path_y.size(), path_y.data());
            yaw.setcontent(path_yaw.size(), path_yaw.data());

            spline1dbuildcubic(s, x, path_sp_x_);
            spline1dbuildcubic(s, y, path_sp_y_);
            spline1dbuildcubic(s, yaw, path_sp_yaw_);
        } catch(alglib::ap_error &e) {
            RCLCPP_ERROR(logger_, "ALGLIB exception: %s", e.msg.c_str());
        }

        path_length_ = path_s.back();
        path_duration_ = std::max(1.5 * path_length_ / (w_max_ * wheel_radius_), 
                          std::sqrt(6.0 * path_length_ / acc_lin_max_));
        path_time_ = 0;
        double step = 1/MPC_frequency_;
        int i = 0;
        for (auto time = 0.0; time <= path_duration_; time += step) {
           geometry_msgs::msg::Pose pose_P;
           geometry_msgs::msg::PoseStamped pose_stamped;
           double x_interp;
           double y_interp;
           double yaw_interp;
           double yaw_ref;
           double x_interp_step;
           double y_interp_step;
           double yaw_interp_step;
           interpolateTrajectory(time, path_sp_x_, path_sp_y_, path_sp_yaw_, x_interp, y_interp, yaw_interp);

           if (time < path_duration_){
            interpolateTrajectory(time + step, path_sp_x_, path_sp_y_, path_sp_yaw_, x_interp_step, y_interp_step, yaw_interp_step);
            yaw_ref = std::atan2(y_interp_step - y_interp, x_interp_step - x_interp);
           }else{
            interpolateTrajectory(time - step, path_sp_x_, path_sp_y_, path_sp_yaw_, x_interp_step, y_interp_step, yaw_interp_step);
            yaw_ref = std::atan2(y_interp - y_interp_step, x_interp - x_interp_step);
           }
           pose_P.position.x = x_interp;// + p_dist_*std::cos(yaw_ref);
           pose_P.position.y = y_interp;// + p_dist_*std::sin(yaw_ref);
           pose_P.position.z = 0.0;
           double roll  = 0.0;
           double pitch = 0.0;
           double yaw   = yaw_interp;
            
           tf2::Quaternion q;
           q.setRPY(roll, pitch, yaw);
           
           geometry_msgs::msg::Quaternion orientation = tf2::toMsg(q);

           pose_P.orientation = orientation;

           pose_stamped.pose = pose_P;

           path_P.poses.push_back(pose_stamped);
        }
    } else {
        std::vector<double> path_x;
        std::vector<double> path_y;
        std::vector<double> path_yaw;
        std::vector<double> path_s;

        path_x.push_back(path.poses[0].pose.position.x);
        path_y.push_back(path.poses[0].pose.position.y);
        path_yaw.push_back(tf2::getYaw(path.poses[0].pose.orientation));
        path_s.push_back(0.0);

        for (size_t k = 1; k < path.poses.size(); k++) {
            double dist = std::hypot(path.poses[k].pose.position.x - path_x.back(),
                                    path.poses[k].pose.position.y - path_y.back());
            if (dist >= path_sampling_threshold_ || k == path.poses.size() - 1) {
                path_x.push_back(path.poses[k].pose.position.x);
                path_y.push_back(path.poses[k].pose.position.y);
                path_yaw.push_back(tf2::getYaw(path.poses[k].pose.orientation));
                path_s.push_back(dist + path_s.back());
            }
        }
        path_idx_ = 0;
        alglib::spline1dinterpolant path_sp_x_, path_sp_y_, path_sp_yaw_;

        // Interpolation of the path components
        try {
            alglib::real_1d_array s, x, y, yaw;
            s.setcontent(path_s.size(), path_s.data());
            x.setcontent(path_x.size(), path_x.data());
            y.setcontent(path_y.size(), path_y.data());
            yaw.setcontent(path_yaw.size(), path_yaw.data());

            spline1dbuildcubic(s, x, path_sp_x_);
            spline1dbuildcubic(s, y, path_sp_y_);
            spline1dbuildcubic(s, yaw, path_sp_yaw_);
        } catch(alglib::ap_error &e) {
            RCLCPP_ERROR(logger_, "ALGLIB exception: %s", e.msg.c_str());
        }

        path_length_ = path_s.back();
        path_duration_ = std::max(1.5 * path_length_ / (w_max_ * wheel_radius_), 
                          std::sqrt(6.0 * path_length_ / acc_lin_max_));
        path_time_ = 0;
        double step = 1/MPC_frequency_;
        int i = 0;
        for (auto time = 0.0; time <= path_duration_; time += step) {
           geometry_msgs::msg::Pose pose_P;
           geometry_msgs::msg::PoseStamped pose_stamped;
           double x_interp;
           double y_interp;
           double yaw_interp;
           interpolateTrajectory(time, path_sp_x_, path_sp_y_, path_sp_yaw_, x_interp, y_interp, yaw_interp);
           pose_P.position.x = x_interp;// + p_dist_*std::cos(yaw_interp);
           pose_P.position.y = y_interp;// + p_dist_*std::sin(yaw_interp);
           pose_P.position.z = 0.0;
           double roll  = 0.0;
           double pitch = 0.0;
           double yaw   = yaw_interp;
            
           tf2::Quaternion q;
           q.setRPY(roll, pitch, yaw);
           
           geometry_msgs::msg::Quaternion orientation = tf2::toMsg(q);

           pose_P.orientation = orientation;

           pose_stamped.pose = pose_P;

           path_P.poses.push_back(pose_stamped);
        }
    }

    path_P.header.stamp = clock_->now();
    path_P.header.frame_id = "map";

    // Publish reference path
    if (reference_path_pub_) {
        nav_msgs::msg::Path ref_path;
        ref_path.header = path_P.header;
        ref_path.poses = path_P.poses;
        reference_path_pub_->publish(ref_path);
    }
}

void MPC_fblin_diffDrive::setSpeedLimit(const double &speed_limit, const bool &percentage) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (speed_limit == nav2_costmap_2d::NO_SPEED_LIMIT) {
        // Restore default value
        max_linear_velocity_ = w_max_ * wheel_radius_;
        max_angular_velocity_ = (w_max_ - w_min_) / track_width_ * wheel_radius_;
    } else {
        if (percentage) {
            max_linear_velocity_ = w_max_ * wheel_radius_ * speed_limit / 100.0;
            max_angular_velocity_ = (max_linear_velocity_ / wheel_radius_ - w_min_) / track_width_ * wheel_radius_;
        } else {
            max_linear_velocity_ = speed_limit;
            max_angular_velocity_ = (max_linear_velocity_ / wheel_radius_ - w_min_) / track_width_ * wheel_radius_;
        }
    }
}

rcl_interfaces::msg::SetParametersResult MPC_fblin_diffDrive::dynamicParametersCallback(
    std::vector<rclcpp::Parameter> parameters) 
{
    rcl_interfaces::msg::SetParametersResult result;
    result.successful = true;
    std::lock_guard<std::mutex> lock_reinit(mutex_);

    for (auto parameter : parameters) {
        const auto &type = parameter.get_type();
        const auto &name = parameter.get_name();

        if (type == ParameterType::PARAMETER_DOUBLE) {
            if (name == plugin_name_ + ".MPC_frequency") {
                MPC_frequency_ = parameter.as_double();
            } else if (name == "controller_frequency") {
                fblin_frequency_ = parameter.as_double();
            } else if (name == plugin_name_ + ".q") {
                q_ = parameter.as_double();
            } else if (name == plugin_name_ + ".r") {
                r_ = parameter.as_double();
            } else if (name == plugin_name_ + ".p_dist") {
                p_dist_ = parameter.as_double();
            } else if (name == plugin_name_ + ".w_min") {
                w_min_ = parameter.as_double();
            } else if (name == plugin_name_ + ".w_max") {
                w_max_ = parameter.as_double();
            } else if (name == plugin_name_ + ".acc_lin_max") {
                acc_lin_max_ = parameter.as_double();
            } else if (name == plugin_name_ + ".wheel_radius") {
                wheel_radius_ = parameter.as_double();
            } else if (name == plugin_name_ + ".track_width") {
                track_width_ = parameter.as_double();
            } else if (name == plugin_name_ + ".path_sampling_threshold") {
                path_sampling_threshold_ = parameter.as_double();
            } else if (name == plugin_name_ + ".next_goal_threshold") {
                next_goal_threshold_ = parameter.as_double();
            }
        } else if (type == ParameterType::PARAMETER_INTEGER) {
            if (name == plugin_name_ + ".prediction_horizon") {
                prediction_horizon_ = parameter.as_int();
            } else if (name == plugin_name_ + ".max_infeasible_sol") {
                max_infeasible_sol_ = parameter.as_int();
            }
        } else if (type == ParameterType::PARAMETER_BOOL) {
            if (name == plugin_name_ + ".use_tracking_controller") {
                use_tracking_controller_ = parameter.as_bool();
            }
        }
    }

    return result;
}

void MPC_fblin_diffDrive::interpolateTrajectory(const double& t, const alglib::spline1dinterpolant& path_sp_x_, const alglib::spline1dinterpolant& path_sp_y_, const alglib::spline1dinterpolant& path_sp_yaw_, double& x, double& y, double& yaw) {
    double s = path_length_ * (3.0 * std::pow(std::min(t/path_duration_, 1.0), 2.0) - 
                               2.0 * std::pow(std::min(t/path_duration_, 1.0), 3.0));

    x = spline1dcalc(path_sp_x_, s);
    y = spline1dcalc(path_sp_y_, s);
    yaw = spline1dcalc(path_sp_yaw_, s);
}

}  // namespace nav2_mpc_fblin_controller

PLUGINLIB_EXPORT_CLASS(
    nav2_mpc_fblin_controller::MPC_fblin_diffDrive,
    nav2_core::Controller)