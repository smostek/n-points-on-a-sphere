## Problem Overview
Let $C$ be a set of $n$ unit vectors. Each point exerts a Coulomb force on all its neighbors such that the net force on each point $\hat r_{io}$ is given by

$$\vec F_i = \sum_{k\neq i}\frac{\vec r_{ik}}{r_{ik}^3}$$

where $r_{ik}=||\hat r_{io}-\hat r_{ko}||$ The component of this force that is tangent to the unit sphere is $\vec F_i^{\perp}=\left(\mathbb{U}-\hat r_{io}\hat r_{io}^T\right)\vec F_i$.

If a distribution of points is equilibrium, then

$$CS(C) = \sum_{i=1}^{n} F_i^{\perp} = 0$$

where $CS$ denotes the 'Cross Sum,' so called because I originally formulated it using a cross-product.

Most 'normal' shapes can be easily arranged so as to be in equilibrium: prisms, pyramids, regular polygons, and of course the Platonic solids. What one quickly notices about this solution set is that they are all well-described as a set of parallel planar polygons. For instance, the icosahedron has one point both poles and two pentagonal bands, as shown:

![icosahedron](https://github.com/smostek/n-points-on-a-sphere/assets/162070478/55374875-46d0-4747-a048-f7b04c5e3825)

This observation, while not strictly true in general, is at the core of the naming scheme used to describe and generate shapes: the icosahedron is given the name (1 5 5' 1), where the prime (') indicates that that second band has been rotated.

## Code Overview
The central file of this project is `shape.m`. Pass it a name as a character vector - basically any sequence of numbers; backticks are used to denote rotated bands - and it will make the corresponding shape. It will automatically adjust the positions of each of the polygonal bands to minimize the Cross Sum, which uses my custom gradient descent algorithm, `gradMin`. Unfortunately, most of the bands in most shapes are not regular, and thus cannot be accessed directly. That's why there's the `forceMin` method: run this and the computer will cautiously move each of the points a small distance in the direction of their force vectors until it either finds the minimum or gets stuck in a valley. One of the principle areas of focus has long been to prevent `forceMin` from getting stuck.

Once you have your minimized distribution, call the `see` method to display it. I have developed a number of different visualizations, and you can pick which one you want to use by passing its name as a char vector to the `see` method. For instance, you can view the polyhedron in its circumsphere by specifying `'solid'`, or you can have the computer highligt the shape's symmetries by passing `'symmetry'`. 

That's just one shape, though, and this project is all about finding as many shapes as possible. Run `TheAlgo` with your favorite value of $n$, and the computer will attempt to find every single equilibrium distribution with that many points. Be warned however, that the algorithm is not fast and I wouldn't recommend giving it a value of $n$ greater than 16 (the number of starting positions increases exponentially with $n$).

This only scratches the surface of the complexity this problem, and there remains much progress to be made, but this overview should be enough to give you a grasp of what's going on (in spite of my pitiful dearth of code comments). If you're really interested or have a question for me, feel free to reach out at smostek0407@gmail.com. Thanks!
