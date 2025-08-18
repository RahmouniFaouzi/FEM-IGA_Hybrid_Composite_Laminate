## 📄 Publication

This repository contains the MATLAB implementation for the research article:

**Isogeometric analysis based on non-uniform rational B-splines technology of stress and failure strength in inter-ply hybrid laminated composite**  
*Journal of Composite Materials, 2022, Vol. 56(18), 2921–2932*

🔗 [Read the full article](https://journals.sagepub.com/doi/abs/10.1177/00219983221105313)


##  Overview

This project implements **Isogeometric Analysis (IGA)** using **Non-Uniform Rational B-Splines (NURBS)** to predict:
* Stress distribution in hybrid laminated composites
* First ply failure strength under uniaxial tensile and compressive loads
* Failure envelopes using different failure criteria:
  - Maximum Stress
  - Maximum Strain
  - Tsai-Wu
  - Tsai-Hill
  - Hoffman

The developed IGA integrates **CAD-based NURBS representation** with numerical analysis for **accurate geometry modeling and efficient simulations**, avoiding the meshing limitations of FEM.

## Features

* NURBS-based isogeometric discretization  
* Supports unidirectional and angle-ply hybrid laminates  
* Predicts failure strength for different fiber orientations  
* Implements multiple failure theories  
* Evaluates effect of:

    - Polynomial order
    - Mesh density
    - Ply thickness  
*  Compares results with **FEM (ANSYS)** and **Autodesk Helius Composite**

## Test
This project includes a C implementation of the function NURBS2Dders for computing NURBS basis function derivatives efficiently. To use this version, you need to compile the MEX file.
This code is developed under Matlab® implementation, it should be able to run on any operating system.

## License
Distributed under the GNU LGPL v.3.0.

## Authors:
Rahmouni Faouzi: 📫 **rahmounifaouzi01@gmail.com**  <br />
Prof. Elajrami Mohamed: 📫 **eladjrami_mohamed@yahoo.fr** <br />
Prof. Kouider Madani: 📫 **koumad10@yahoo.fr** <br />
Prof. Raul Campilho <br />
