## Problem Overview
Let $C$ be a set of $n$ points on a unit sphere with its center at point $o$. Each point exerts a Coulomb force on all its neighbors such that the net force vector associated with each point $i$ is given by

$$\vec F_i = \sum_{k\neq i}\frac{\vec r_{ik}}{r_{ik}^3}$$

where a vector without a hat (e.g. $r_{ik}$) indicates that vector's magnitude.

The component of this force that is tangent to the unit sphere is 

$$\vec F_i^{\perp}=\left(\mathbb{1}-\hat r_{io}\hat r_{io}^T\right)\vec F_i$$

If a point is in equilibrium, then this perpendicular force component will have zero magnitude. So, we define the 'Cross Sum' (so called because I orignially formulated the idea using a cross-product) as

$$CS(C) = \sum_{i=1}^{n} F_i^{\perp}$$

The goal is to find all distributions $C$ for which $CS(C)=0$. Phrased that way, it doesn't sound very hard. Phrased as a search for an unknown number of predominantly unstable zeros of a highly non-linear function of $2n$ variables, the difficulty becomes clearer.

Most 'normal' shapes can be easily arranged so as to be in equilibrium: prisms, pyramids, regular polygons, and of course the Platonic solids. What one quickly notices about this solution set is that they are all well-described as a set of parallel planar polygons. For instance, the icosahedron has one point on both poles and two pentagonal bands, as shown:

![icosahedron](https://github.com/smostek/n-points-on-a-sphere/assets/162070478/55374875-46d0-4747-a048-f7b04c5e3825)

This observation, while not strictly true in general, has proven very helpful for finding and categorizing solutions and it remains the means by which shapes in the program are named: the icosahedron has the name (1 5 5' 1) where the prime (') indecates that the lower pentagon has been rotated. Nevertheless, I am transitioning away from this towards the more rigorous language of Group Theory, though there remain several kinks to work out.

## Code Overview
The central file of this project is `shape.m`. Pass it a name as a character vector - basically any sequence of numbers; backticks are used as primes - and it will make the corresponding shape as per the naming convention just discussed. It will automatically adjust the positions of each of the polygonal bands to minimize the Cross Sum, using my custom gradient descent algorithm, `gradMin`. Unfortunately, most of the bands in most shapes are not regular, and thus cannot be accessed directly. That's why there's the `forceMin` method: run this and the computer will cautiously move each of the points a small distance in the direction of their force vectors until it either finds the minimum or gets stuck in a valley. One of the principle areas of focus has long been to prevent `forceMin` from getting stuck.

Once you have your minimized distribution, call the `see` method to display it. I have developed a number of different visualizations, and you can pick which one you want to use by passing its name as a char vector to the `see` method. For instance, you can view the polyhedron formed by the points by specifying `'solid'`, or you can have the computer highligt the shape's symmetries by passing `'symmetry'`. 

That's just one shape, though, and this project is all about finding as many shapes as possible. Run `TheAlgo` with your favorite value of $n$, and the computer will attempt to find every single equilibrium distribution with that many points. Be warned however, that the algorithm is not fast and I wouldn't recommend giving it a value of $n$ greater than 16 (the number of starting positions increases exponentially with $n$).

This only scratches the surface of the complexity this problem, and there remains much progress to be made, but this overview should be enough to give you a grasp of what's going on (in spite of my pitiful dearth of code comments). If you're really interested or have a question for me, feel free to reach out at smostek0407@gmail.com. Thanks!
