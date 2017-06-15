#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;
// ditto
using namespace std;

namespace
{
  // Debug
  template<typename T>
  ostream& operator<<(ostream&o, const vector<T> &v) {
    o << "[ ";
    for(auto &x : v)
      o << x << ", ";
    o << "]";
    return o;
  }
  
  // For converting back and forth between radians and degrees.
  constexpr double pi() { return M_PI; }
  double deg2rad(double x) { return x * pi() / 180; }
  double rad2deg(double x) { return x * 180 / pi(); }

  // Checks if the SocketIO event has JSON data.
  // If there is data the JSON object in string format will be returned,
  // else the empty string "" will be returned.
  string hasData(string s) {
    auto found_null = s.find("null");
    auto b1 = s.find_first_of("[");
    auto b2 = s.rfind("}]");
    if (found_null != string::npos) {
      return "";
    } else if (b1 != string::npos && b2 != string::npos) {
      return s.substr(b1, b2 - b1 + 2);
    }
    return "";
  }

  // Evaluate a polynomial.
  double polyeval(Eigen::VectorXd coeffs, double x) {
    int i = coeffs.size() - 1;
    double result = coeffs[i];
    while (i > 0)
      result = result * x + coeffs[--i];
    return result;
  }
  
  // Fit a polynomial.
  Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                          int order) {
    assert(xvals.size() == yvals.size());
    assert(order >= 1 && order <= xvals.size() - 1);
    Eigen::MatrixXd A = Eigen::MatrixXd::Ones(xvals.size(), order + 1);
    for (int j = 1; j < order + 1; ++j)
      A.col(j) = A.col(j - 1).cwiseProduct(xvals);
    return A.householderQr().solve(yvals);
  }
  

  // Autonomous Driver using MPC

  // Refactored outside the lambda for clarity
  struct AutonomousDriver
  {
    struct Inputs
    {
      double steer;
      double throttle;
    };
    
    MPC mpc;
    
    // last inputs
    Inputs inputs;
    // current state
    Eigen::VectorXd state;
    
    // Enable/disable the generation of the reference, modeled curve
    bool visualizations{true};
    // Display the MPC predicted trajectory
    vector<double> mpc_x_vals;
    vector<double> mpc_y_vals;
    // Display the waypoints/reference line
    vector<double> next_x_vals;
    vector<double> next_y_vals;
    
    
    AutonomousDriver(size_t N, double dt)
    : mpc(N, dt), state(6)
    {}
    
    template <typename T>
    Inputs Update(const T& telemetry) {
      // Waypoints in global coordinates
      vector<double> ptsx = telemetry["ptsx"];
      vector<double> ptsy = telemetry["ptsy"];
      // Position in global coordinates
      double px = telemetry["x"];
      double py = telemetry["y"];
      // Orientation in radians (math coord, x right hand)
      double psi = telemetry["psi"];
      // Orientation in radians (navigation coord, y left hand)
      //double psi_nav = telemetry["psi_unity"];
      // Current speed in mph
      double v = telemetry["speed"];
      
      // Current steering angle in radians
      //double steering_angle = telemetry["steering_angle"];
      // Current throttle value normalized [-1, 1]
      //double throttle = telemetry["throttle"];

      assert(ptsx.size() == 6);
      assert(ptsy.size() == 6);
      
    #if VERBOSE_DEBUG
      cout << "psi = " << psi
           << "; px = " << px
           << "; py = " << py << endl;
      cout << "ptsx = " << ptsx << endl;
      cout << "ptsy = " << ptsy << endl;
    #endif
      
      // Convert waypoints to car coordinates
      Eigen::VectorXd ptsx_v(6);
      Eigen::VectorXd ptsy_v(6);
      
      for(int i = 0 ; i < 6; i++) {
        const double
          dx = ptsx[i] - px,
          dy = ptsy[i] - py,
          // The compiler should put these guys outside the loop.
          cos_psi = cos(psi),
          sin_psi = sin(psi);
        ptsx_v[i] =   dx * cos_psi + dy * sin_psi;
        ptsy_v[i] = - dx * sin_psi + dy * cos_psi;
      }

      // We are now in the car reference:  psi = px = py = 0.

      // Fit the waypoints to a 3-order polynomial
      auto coeffs = polyfit(ptsx_v, ptsy_v, 3);
      
      // cte: shortest distance from the origin to the polynomial curve.
      // For now we approximate it with the horizontal distance.
      double cte = polyeval(coeffs, 0);
      
      // epsi in car coordinates
      double epsi = -atan(coeffs[1]);
      
    #if VERBOSE_DEBUG
      cout << "cte = " << cte << "; epsi = " << epsi << endl;
    #endif
      
      // Update the state
      state << /*px*/0, /*py*/0, /*psi*/0, v, cte, epsi;
      
      // Solve the MPC for this state
      auto vars = mpc.Solve(state, coeffs);
      
      // vars = x, y, psi, v, cte, epsi, delta, a
      inputs.steer = vars[0];
      inputs.throttle = vars[1];

      // Handle visualizations
      mpc_x_vals.clear();
      mpc_y_vals.clear();
      
      next_x_vals.clear();
      next_y_vals.clear();
      
      if (visualizations) {
        //for(int i = 2; i < vars.size(); i++) {
        //}
        
        // Reference line
        for(int i = 1 ; i < 25; i++) {
          double x = 2.5 * i;
          double y = polyeval(coeffs, x);
          next_x_vals.push_back(x);
          next_y_vals.push_back(y);
        }
      }
      
      // Return the inputs to drive the car another step
      return inputs;
    }
  };
  
} // (anonymous)



int main(int argc, char **argv) {
  uWS::Hub h;
  
  // timestamp configuration
  size_t N = 10;
  double dt = 0.1;
  
  // Parse command line arguments
  for (int i = 1; i < argc; i++) {
    if (!strncmp(argv[i], "-N=", 3)) {
      N = atoi(argv[i]+3);
      if (N <= 0 || N > 100) {
        cerr << "Invalid N parameter." << endl;
        return -1;
      }
    } else if (!strncmp(argv[i], "-dt=", 4)) {
      dt = atof(argv[i]+4);
      if (dt <= 0 || dt > 1) {
        cerr << "Invalid dt parameter." << endl;
        return -1;
      }
    } else {
      cout << "Usage: mpc [options]" << endl
           << "Options:" << endl
           << "  -N=<INT>       Timestamp length (default " << N << ")." << endl
           << "  -dt=<FLOAT>    Timestamp duration (default " << dt << ")." << endl
           << endl;
      if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "-help"))
        return 0;
      return -1;
    }
  }
  
  cout << "Parameters:" << endl
       << "  N  = " << N << endl
       << "  dt = " << dt << endl;
  
  // MPC is initialized here!
  AutonomousDriver driver(N, dt);
  
  h.onMessage([&driver](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
          // Compute the inputs
          auto inputs = driver.Update(j[1]);
          
          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25)] instead of [-1, 1].
          msgJson["steering_angle"] = inputs.steer / deg2rad(25);
          msgJson["throttle"] = inputs.throttle;

          msgJson["mpc_x"] = driver.mpc_x_vals;
          msgJson["mpc_y"] = driver.mpc_y_vals;

          msgJson["next_x"] = driver.next_x_vals;
          msgJson["next_y"] = driver.next_y_vals;

          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          //cout << msg.substr(2) << endl;

          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
