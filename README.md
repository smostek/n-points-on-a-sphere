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
The central file of this project is `shape.m`. 
