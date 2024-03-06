## Problem Overview
Let $C$ be a set of $n$ unit vectors. Each point exerts a Coulomb force on all its neighbors such that the net force on each point $\vec r_i$ is given by 
$$\vec F_i = \sum_{k\neq i}\frac{\vec r_i-\vec r_k}{||\vec r_i-\vec r_k||^3}$$
The component of this force that is tangent to the unit sphere is $\vec T_i=\left(I-\vec r_i\vec r_i^T\right)\vec F_i$. If a distribution of points is equilibrium, then
$$CS(C) = \sum_{i=1}^{n}||\vec T_i|| = 0$$
where $CS$ denotes the 'Cross Sum,' so called because I originally formulated it using a cross-product.

Most 'normal' shapes can be easily arranged so as to be in equilibrium: prisms, pyramids, regular polygons, and of course the Platonic solids. What one quickly notices about this solution set is that they are all well-described as a set of parallel planar polygons. For instance, the icosahedron has one point both poles and two pentagonal bands, as shown:

![icosahedron](https://github.com/smostek/n-points-on-a-sphere/assets/162070478/55374875-46d0-4747-a048-f7b04c5e3825)

This observation, while not strictly true in general, is at the core of the naming scheme used to describe and generate shapes: the icosahedron is given the name (1 5 5' 1), where the prime (') indicates that that second band has been rotated.
## Code Overview
The central file of this project is `shape.m`. Pass it a name as a character vector - basically any sequence of numbers; backticks are used to denote rotated bands - and it will make the corresponding shape. It will automatically adjust the positions of each of the polygonal bands to minimize the Cross Sum, which uses my custom gradient descent algorithm, `gradMin`. Unfortunately, most shapes are not symmetric enough for this to work, and that's why there's the `forceMin` method: run this and the computer will cautiously move each of the points a small distance in the direction of their force vectors until it either finds the minimum or gets stuck in a valley. One of the principle areas of focus has long been to prevent `forceMin` from getting stuck.

Once you have your minimized distribution, call the `rename` method to analyze the shape. This method does a lot, but the main ultimate effect is that it finds every other valid name for the distribution. What does that mean? Call the `see` method and pass it `'band'` to find out: use the dropdown to see your shape from a new perspective. And there are plenty of other visualizations I've developed: pass it `'solid'` to see the polyhedron defined by the distribution's coordinates; give it `'highlight'` to have it group the polyhedron's similar faces and edges; or try `'contourPlot'` to see how the strength of the Coulomb potential varies over the surface.

That's just one shape, though, and this project is all about finding as many shapes as possible. Run `TheAlgo` with your favorite value of $n$, though be warned that problems start cropping in above $n=10$ so try not to push the program too hard. That might sound underwhelming, but when I first tried to find the solutions for $n=12$ a year ago it took over half and hour to calculate, and now it takes less than 10 seconds, so there's been a lot of progress.

This only scratches the surface of the complexity this problem, and there remains much progress to be made, but this overview should be enough to give you a grasp of what's going on (in spite of my pitiful dirth of code comments). If you're really interested or have a question for me, feel free to reach out at smostek0407@gmail.com. Thanks!
