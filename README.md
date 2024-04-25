# Apollonian Gaskets

## Tangent Circles

Start by considering 3 circles that are pairwise tangent at 3 distinct points. Denote the signed curvatures of these circles by $\kappa_i \in \mathbb{R}$ and the centres by $z_i \in \mathbb{C}$, for $i = 1,2,3$. The curvature is just the reciprocal of the radius and the sign is positive if we are taking the perspective of being outside the circle and negative if we are taking the perspective of being inside the circle.

Descarte's theorem states that there exists exactly 2 circles that are tangent to all three of these circles. Furthermore, the signed curvature $\kappa_4$ of these circles satisfy
$$\left( \sum_{i=1}^4 \kappa_i \right)^2 = 2 \sum_{i=1}^4 \kappa_i^2$$
More recently, it was shown that the centres satisfy
$$\left( \sum_{i=1}^4 \kappa_i z_i \right)^2 = 2 \sum_{i=1}^4 (\kappa_i z_i)^2$$
Note that each equation gives 2 solutions, giving 4 in total. However, only 2 of these circles will be tangent to all 3 of the original circles.

## Constructing a Gasket

If we start with 3 pairwise tangent circles and use this result to generate 2 more, we now have 5 circles. From these five circles we can then take triplets of pairwise tangent circles and generate even more circles. After $n$ steps, there will be $3^n + 2$ circles. 

Lets start with 3 circles with centres
$$(-1,0) \qquad (+1, 0) \qquad (0, +\sqrt{3})$$
and all having signed curvature equal to $+1$.

First write a function `is_tangent` to check if any 2 given circles are tangent to each other. It should return `True` for each pair of our circles.

Next, write a function `get_circles` that takes 3 pairwise tangent circles and returns 2 new circles that are each tangent to all 3 of the input circles using the formulae above. You may use `is_tangent` to check that any circles you find are valid.

You can now repeatedly apply `get_circles` to generate the gasket. You will need to set some form of stopping condition to make the program run n a finite amount of time.

The gasket can be displayed using tools from `matplotlib`.

## Extensions

Try generating gaskets for different triplets of initial circles. To generate more initial triplets, use [Soddy circles of a triangle](https://en.wikipedia.org/wiki/Soddy_circles_of_a_triangle).

Think about how to optimise this procedure by modifying either the iterative process or the stopping condition. For example, some of the circles we find will be so small that they are barely visible in the rendered gasket so they could possibly be excluded.

This problem can also be extended to higher dimensions (and even different geometries) but the centres can no longer be found by just solving a quadratic equation. For further details see [here](https://arxiv.org/abs/math/0101066).
