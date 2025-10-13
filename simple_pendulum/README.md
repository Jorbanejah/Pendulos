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

## 🧰 What This Folder Contains

- `simple_pendulum.py`: Python script that numerically solves the pendulum equation using the Runge-Kutta method.
- `energy_plot.png`: Visualizes total, kinetic, and potential energy over time.
- `trajectory.gif`: Animation of the pendulum’s motion.
- `README.md`: You’re reading it!

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
