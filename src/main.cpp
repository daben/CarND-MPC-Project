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
#include "utils.h"

// for convenience
using json = nlohmann::json;

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


// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                       int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A = Eigen::MatrixXd::Ones(xvals.size(), order + 1);
  for (int j = 1; j < order + 1; ++j)
    A.col(j) = A.col(j - 1).cwiseProduct(xvals);
  return A.householderQr().solve(yvals);
}

int main() {
  uWS::Hub h;

  // Latency in seconds. We are simulating a latency in the actuator
  // of 100 ms + an extra to account for the time involved in computing.
  double latency = 0.15;
  {
    // allow the user to specify the latency
    char *value = getenv("latency");
    if (value) latency = atof(value);
    printf("Latency=%.2f; ", latency);
  }

  // MPC is initialized here!
  MPC mpc;

  // For the waypoints/reference line
  vector<double> next_x_vals(25);
  vector<double> next_y_vals(25);
  for (int i = 0; i < next_x_vals.size(); i++) {
    // precomputed x
    next_x_vals[i] = 2.5 * i;
    // y will be computed later
    next_y_vals[i] = 0.0;
  }

  h.onMessage([&](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                  uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object

          // waypoints, I assume they come in meters...
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          // car position in global coordinates
          double px = j[1]["x"];
          double py = j[1]["y"];
          // orientation in radians
          double psi = j[1]["psi"];
          // speed in MPH
          double speed = j[1]["speed"];

          // Current steering angle in radians
          double steer_value = j[1]["steering_angle"];
          // Current throttle value normalized [-1, 1]
          double throttle_value = j[1]["throttle"];

          // The speed comes in MPH
          double v = mph2mps(speed);
          double accel = throttle_value;

          /*
          * Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */
          benchmark clock;
          clock.reset();

          // Handle the latency predicting the state of the system
          // from the current state.

          // Kinematic model
          if (latency > 1e-3) {
            // from MPC.cpp...
            const double Lf = 2.67;
            const double
              cos_psi = cos(psi),
              sin_psi = sin(psi);

            px  += v * cos_psi * latency + 0.5 * accel * cos_psi * pow(latency, 2);
            py  += v * sin_psi * latency + 0.5 * accel * sin_psi * pow(latency, 2);
            v   += accel * latency;
            // Note the minus sign to apply correct meaning of the steer_value
            psi -= v * steer_value / Lf * latency;
          }

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
          // We could compute the actual distance
          // double cte = polycte(coeffs, 3);
          //
          // But approximating it with the horizontal distance at the origin
          // it's quite good already.
          //
          // This should be `cte = polyeval(coeffs, x) - y`, but x and y
          // are zero in car coordinates:
          double cte = coeffs[0];

          // cout << coeffs.transpose() << " ⊢ " << cte << endl;

          // epsi in car coordinates
          //  `epsi = psi - atan(c1 + 2 c2 x + 3 c3 x^2)`
          // with psi = x = 0:
          double epsi = -atan(coeffs[1]);

          // Radius of curvature 3m ahead, computed over a polynomial
          // of second degree
          double radius = polyroc(polyfit(ptsx_v, ptsy_v, 2), 2, 3);

          // Update the state
          Eigen::VectorXd state(6);
          state << /*px*/0, /*py*/0, /*psi*/0, v, cte, epsi;

          // Solve the MPC for this state
          MPC::Result result;

          bool ok = mpc.Solve(state, coeffs, radius, result);

          const auto mpc_elapsed = clock.elapsed();

          if (ok) {
            steer_value = -result.delta / deg2rad(25);
            throttle_value = result.a;

            fprintf(stdout,
                    "Cost = %8.2e; δ = %6.3f; a = %6.3f; "
                    "cte = %7.3f; eψ = %11.8f; "
                    "x = %7.2f; y = %7.2f; speed = %6.2f; "
                    "Radius = %7.3f; "
                    "Solve(ms) = %6.3f\n",
                    result.cost, steer_value, throttle_value,
                    cte, epsi, px, py, speed,
                    radius,
                    mpc_elapsed * 1e3);

          } else {
            // We should do something better here...
            cerr << "*** WARNING ***" << endl;
          }

          // Waypoints curve
          for (int i = 0; i < next_x_vals.size(); i++) {
            next_y_vals[i] = polyeval(coeffs, next_x_vals[i]);
          }

          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          // Display the MPC predicted trajectory
          msgJson["mpc_x"] = result.x;
          msgJson["mpc_y"] = result.y;

          // Display the waypoints polynomial curve
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          //std::cout << msg << std::endl;
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
    //ws.close();
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
