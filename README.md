# Algorithms
A set of four algorithms to find the shortest path between two points of a bidimensional space under a potential, this four methods are:
- Method 1: The first approximation of the algorithm is to look at each of the neighbours given a point P, moving to the one that reduces the weight function and continue this way until the end is reached.
- Method 2: This method can be conceived as an approximation to Bellman-Ford algorithm reducing the number of times the relaxation function is executed for each pair points, here we only relax the vertices that verify the condition dist(Pj, Pf) < dist(P, Pf) where Pj is one of the four neighbours of P.
- Dijkstra's shortest path algoithm.
- Bellman-Ford's shortest path algoithm.

There are some parameters to play with:
- `h`: the step size, the lower its value the smoother the path but for this will cause some algorithms to crash under a `Max-Depth Limit Excedeed` exception due to the large number of vertex the graph will have.
- `ALPHA`: the correction coeficient on the distance, this makes the distance be more sensitive in changes of potential `> 1` or less sensitive `< 1`. For cases near an infinite divergence in the function, values of ALPHA greater than 1 give better results than the default one.
- Weight function: this is used to calculate the weight of each edge, there are two weight function implemented that give better results in some cases and worst in other so there is no clear winner in the case.

A complete explanation of this methods and a discussion on them can be found in this [article](/discussion.pdf).

# License
    Copyright 2017 labay11

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
