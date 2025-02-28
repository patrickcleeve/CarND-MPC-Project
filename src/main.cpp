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
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];


          // Flags for debugging
          bool debug_mode = false;
          bool controller_data_mode = false;
          bool equation_mode = true;
          
          // Vehicle constant
          const double Lf = 2.67;

          // Convert points from map to car coordinates, shift reference angle
          for (int i = 0; i < ptsx.size(); i++ ) {

            // shift car reference angle to 90 degrees
            double shift_x = ptsx[i] - px;
            double shift_y = ptsy[i] - py;

            ptsx[i] = (shift_x * cos(0-psi) - shift_y * sin(0-psi));
            ptsy[i] = (shift_x * sin(0-psi) + shift_y * cos(0-psi));

          }

          double* ptrx = &ptsx[0];
          Eigen::Map<Eigen::VectorXd> ptsx_transform(ptrx, 6);

          double* ptry = &ptsy[0];
          Eigen::Map<Eigen::VectorXd> ptsy_transform(ptry, 6);

          auto coeffs = polyfit(ptsx_transform, ptsy_transform, 3);

          // Calculate cte and epsi
          double cte = polyeval(coeffs, 0);
          double epsi = -atan(coeffs[1]);

          
          // steer value needs to be inverted (map coords are LH, car is RH)
          double steer_value = j[1]["steering_angle"];
          double throttle_value =  j[1]["throttle"];


          // dt in seconds, latency in ms
          double dt = 0.1;
          int latency = dt * 1000;

          // Use equations of motion to predict where vehicle should be 
          // after time delay (latency)
          // Readings we receive are from 100ms in the past, 
          // therefore to predict the current location to best accuracy
          // use equations of motion.

          // Because of change in coordiante system
          // System now in car coordinates. centre of car > 0, 0, psi0 = 0
          // car orientated straight ahead in car coordinate system
          // x_t = 0
          // y_t = 0
          // psi_t = 0
          // v_t = v
          // cte_t = cte
          // epsi_t = epsi = -atan(coeffs[1])
          // a_t = throttle_value
          // delta_t = steer_angle

          // x_t+1    = x_t + v * cos(psi) * dt
          // y_t+1    = y_t + v * sin(psi) * dt
          // psi_t+1  = psi_t - v_t / Lf * delta_t * dt
          // v_t+1    = v_t + a_t * dt
          // cte_t+1  = f(x_t) - y_t + v_t * sin(epsi_t) * dt
          // epsi_t+1 = psi_t - psi_t_des + v_t / Lf * delta_t * dt

          // Latency predicted equations of motion
          // These can be simplified

          double x_dt     = 0 + v * cos(0) * dt;
          double y_dt     = 0 + v * sin(0) * dt;
          double psi_dt   = 0 - v / Lf * steer_value * dt;
          double v_dt     = v + throttle_value * dt;
          double cte_dt   = cte - 0 + v * sin(epsi) * dt;
          double epsi_dt  = 0 - epsi + v / Lf * steer_value * dt;
          
          if (equation_mode) {

            cout << "Latency Update Equations:"                                                 << endl;
            cout << "       t     t+dt         diff"                                            << endl;
            cout << "X    : " << "0"  << "  ||  "  << x_dt    << "  || " << abs(0 - x_dt)       << endl;
            cout << "Y    : " << "0"  << "  ||  "  << y_dt    << "  || " << abs(0 - y_dt)       << endl;
            cout << "PSI  : " << "0"  << "  ||  "  << psi_dt  << "  || " << abs(0 - psi_dt)     << endl;
            cout << "V    : " << v    << "  ||  "  << v_dt    << "  || " << abs(v - v_dt)       << endl;
            cout << "CTE  : " << cte  << "  ||  "  << cte_dt  << "  || " << abs(cte - cte_dt)   << endl;
            cout << "EPSI : " << epsi << "  ||  "  << epsi_dt << "  || " << abs(epsi - epsi_dt) << endl;

          } 

          // Update state equation
          Eigen::VectorXd state(6);
          state << x_dt, y_dt, psi_dt, v_dt, cte_dt, epsi_dt;

          if (debug_mode) { cout << "SOLVE START " << endl; };
          
          // Solve MPC
          auto vars = mpc.Solve(state, coeffs);

          if (debug_mode) { cout << "SOLVE FINISHED" << endl; };


          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = vars[0] / (deg2rad(25) * Lf);
          msgJson["throttle"] = vars[1];

          // Display the MPC predicted trajectory (Green)
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line

          for (int i = 2; i < vars.size(); i++) {

            if(i % 2 == 0) {
              
              mpc_x_vals.push_back(vars[i]);
            
            } else {

              mpc_y_vals.push_back(vars[i]);

            }
          }

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;


          // Display the waypoints/reference line (Yellow)
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line

          double poly_inc = 2.5;
          int num_points = 25;
          for (int i = 1; i < num_points; i++) {

            next_x_vals.push_back(poly_inc * i);
            next_y_vals.push_back(polyeval(coeffs, poly_inc * i));

          }

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          if (controller_data_mode) { std::cout << msg << std::endl; };

          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.

          this_thread::sleep_for(chrono::milliseconds(latency));
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
