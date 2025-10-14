# 🧮 Simple Pendulum Simulation

This module simulates the motion of a simple pendulum using numerical methods. It serves as a foundational example of classical mechanics and demonstrates how Lagrangian dynamics can be applied to derive and solve equations of motion.

## 🎓 Why the Simple Pendulum?

The simple pendulum is one of the few mechanical systems that every physicist should be able to solve analytically. It’s a classic example of a system governed by a second-order differential equation, and it introduces key concepts like:

- Harmonic motion
- Energy conservation
- Small-angle approximation
- Nonlinear dynamics (for larger angles)

In *Mechanics I*, we studied the Lagrangian formulation:

$$
L = T - U = \\frac{1}{2} m l^2 \\dot{\\theta}^2 - mgl(1 - \\cos\\theta)
$$

From this, we derived the equation of motion using the Euler-Lagrange equation:

$$
\\frac{d}{dt} \\left( \\frac{\\partial L}{\\partial \\dot{\\theta}} \\right) - \\frac{\\partial L}{\\partial \\theta} = 0
$$

This leads to the nonlinear differential equation:

$$
\\ddot{\\theta} + \\frac{g}{l} \\sin\\theta = 0
$$

## 🎞️ Animation

To better visualize the pendulum dynamics, you can watch the following MP4 animations:

- [Pendulum_simulation.mp4](Pendulum_simulation.mp4):  
  A clean animation of the simple pendulum’s motion, showing angle evolution over time.

- [Combination_simulation.mp4](Combination_simulation.mp4):  
  A composite animation that compares different initial conditions or behaviors.

> 💡 Tip: Click the links to download or view the videos in your browser.

## 📦 Contents

This folder includes the following files:

| File                             | Description                                                                 |
|----------------------------------|-----------------------------------------------------------------------------|
| `simulacion_pendulo_simple.m`    | MATLAB script for simulating the simple pendulum                           |
| `combinacion_pendulos.m`         | MATLAB script combining multiple pendulum behaviors                        |
| `Description_analisis_pendulo_simple.m` | MATLAB file with descriptive analysis of the simple pendulum (ode45, implicit Euler and RK4)        |
| `Pendulum_simulation.mp4`        | MP4 animation showing the pendulum’s motion                                |
| `Combination_simulation.mp4`     | MP4 animation comparing different pendulum behaviors                       |


> All simulations are self-contained. You can modify initial conditions and parameters directly in the scripts.


## 🌀 Alternative Uses of This Calculus

The same mathematical tools used here—Lagrangian mechanics and numerical integration—are widely applicable:

- **Engineering**: Modeling suspension systems, robotic arms, or oscillating components.
- **Quantum Mechanics**: The classical pendulum is a stepping stone to understanding the quantum harmonic oscillator.
- **Cybersecurity**: While less direct, the ability to model dynamic systems and solve differential equations is foundational for cryptographic algorithms and signal analysis.
- **Freelance Projects**: This simulation can be adapted for educational content, interactive visualizations, or physics-based animations in web apps.

## 🚀 How to Run

1. Clone the repository:
   ```bash
   git clone https://github.com/Jorbanejah/Pendulos.git
   cd Pendulos/simple_pendulum
