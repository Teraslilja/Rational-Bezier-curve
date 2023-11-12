# Rational-Bezier-curve
Calculate rational Bezier curves (2D or 3D)

See https://en.wikipedia.org/wiki/B%C3%A9zier_curve#Rational_B%C3%A9zier_curves

Purpose of this repository is study, how to implement rational Bezier curve calculation with C++20.

The classes and structures
* curve::bezier::rational::Rational manages control points.
* curve::bezier::rational::ValidateRational validates the control points.
* curve::bezier::rational::internal::CalculateRational implements the calculation related to the curve defined by control points.
* curve::bezier::rational::DelegationInterface defines the delagation interface between Rational and ValidateRational.

The class Rational can be replaced with other implementations written by e.g. Python or C# as it only manages control points and passes control points as a span to class ValidateRational.
