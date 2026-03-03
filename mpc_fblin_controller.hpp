#ifndef NAV2_MPC_FBLIN_CONTROLLER__MPC_FBLIN_DIFFDRIVE_CONTROLLER_HPP_
#define NAV2_MPC_FBLIN_CONTROLLER__MPC_FBLIN_DIFFDRIVE_CONTROLLER_HPP_

#include <string>
#include <vector>
#include <mutex>
#include <array>
#include "alglib/interpolation.h"

#include "nav2_costmap_2d/footprint_collision_checker.hpp"
#include "nav2_core/controller.hpp"
#include "rclcpp/rclcpp.hpp"
#include "pluginlib/class_list_macros.hpp"
#include "nav2_mpc_fblin_controller/msg/control_inputs.hpp"

namespace nav2_mpc_fblin_controller
{

class MPC_fblin_diffDrive : public nav2_core::Controller
{
public:
  MPC_fblin_diffDrive();
  ~MPC_fblin_diffDrive() override;

  void configure(
    const rclcpp_lifecycle::LifecycleNode::WeakPtr & parent,
    std::string name, std::shared_ptr<tf2_ros::Buffer> tf,
    std::shared_ptr<nav2_costmap_2d::Costmap2DROS> costmap_ros) override;

  void cleanup() override;
  void activate() override;
  void deactivate() override;

  geometry_msgs::msg::TwistStamped computeVelocityCommands(
    const geometry_msgs::msg::PoseStamped & pose,
    const geometry_msgs::msg::Twist & velocity,
    nav2_core::GoalChecker * goal_checker) override;

  void setPlan(const nav_msgs::msg::Path & path) override;
  //void setPlan(const nav_msgs::msg::Path &plan) override;

  void setSpeedLimit(const double & speed_limit, const bool & percentage) override;

protected:
  rcl_interfaces::msg::SetParametersResult dynamicParametersCallback(std::vector<rclcpp::Parameter> parameters);
  void interpolateTrajectory(const double& t, const alglib::spline1dinterpolant& path_sp_x_, const alglib::spline1dinterpolant& path_sp_y_, const alglib::spline1dinterpolant& path_sp_yaw_, double& x, double& y, double& yaw);

  // ROS 2 infrastructure
  rclcpp_lifecycle::LifecycleNode::WeakPtr node_;
  std::shared_ptr<tf2_ros::Buffer> tf_;
  std::shared_ptr<nav2_costmap_2d::Costmap2DROS> costmap_ros_;
  nav2_costmap_2d::Costmap2D * costmap_;
  rclcpp::Logger logger_{rclcpp::get_logger("MPC_fblin_diffDrive")};
  rclcpp::Clock::SharedPtr clock_;
  std::string plugin_name_;

  // Configuration parameters
  double p_dist_;
  double fblin_frequency_;
  double w_min_, w_max_, wheel_radius_, track_width_;
  bool use_tracking_controller_;
  double path_sampling_threshold_;
  double next_goal_threshold_;
  double max_linear_velocity_;
  double max_angular_velocity_;
  double MPC_frequency_;
  double q_;
  double r_;
  double acc_lin_max_;
  int prediction_horizon_;
  int max_infeasible_sol_;

  // Plan bookkeeping (order matters)
  std::vector<geometry_msgs::msg::PoseStamped> plan_;
  size_t path_idx_ {0};
  size_t path_time_ {0};
  double path_duration_ {0.0};
  double path_length_ {0.0};



  // Trajectory tracking
  // double path_duration_, path_length_;
  alglib::spline1dinterpolant path_sp_x_, path_sp_y_, path_sp_yaw_;
  std::vector<double> path_x_, path_y_, path_yaw_;
  // long unsigned int path_idx_, path_time_;

  // MPC outputs subscription
  rclcpp::Subscription<nav2_mpc_fblin_controller::msg::ControlInputs>::SharedPtr mpc_sub_;
  std::array<double, 2> current_mpc_output_{0,0};
  std::mutex mpc_mutex_;

  // Collision checking
  std::unique_ptr<nav2_costmap_2d::FootprintCollisionChecker<nav2_costmap_2d::Costmap2D*>> collision_checker_;

  // Dynamic parameters
  std::mutex mutex_;
  rclcpp::node_interfaces::OnSetParametersCallbackHandle::SharedPtr dyn_params_handler_;

  // Publishers
  std::shared_ptr<rclcpp_lifecycle::LifecyclePublisher<geometry_msgs::msg::PoseStamped>> next_goal_pub_;
  std::shared_ptr<rclcpp_lifecycle::LifecyclePublisher<nav_msgs::msg::Path>> reference_path_pub_;
};

}  // namespace nav2_mpc_fblin_controller

#endif  // NAV2_MPC_FBLIN_CONTROLLER__MPC_FBLIN_DIFFDRIVE_CONTROLLER_HPP_
