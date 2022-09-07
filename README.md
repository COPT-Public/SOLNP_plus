#  SOLNP README
## Basic Introduction
This is C implementation of SOLNP algorithm, which is proposed by Yinyu Ye (1989) and originally implemented in Matlab. Various improvements have been made to increase the robustness and reduce number of function evaluations compared with the original version. 
    
SOLNP solves the general derivative-free nonlinear constrained problems of the form:

$$
\begin{aligned}
    \min\_x\ &f(x) \\ 
    \text{s.t.} \ & g(x) = 0, \\
                  & l\_h\leq h(x)\leq u\_h, \\
                  & l\_x \leq x \leq u\_x,
\end{aligned}
$$

where $f(x),g(x),h(x)$ are smooth functions.
## Getting Started
Currently, SOLNP only provides Matlab interface. The user can learn the details of installing and using in `SOLNP USER GUIDE.html`. The user can choose OSQP solver to solve the QP subproblem.
## Developing Team
- Jiyuan Tan Email: JiyuanTan19@gmail.com
- Tianhao Liu Email:liu.tianhao@163.sufe.edu.cn
- Jinsong Liu Email: liujinsong@163.sufe.edu.cn
## Reference
- Ye, Y. (1988). Interior algorithms for linear, quadratic, and linearly constrained convex programming. Stanford University.
- [Original source](https://web.stanford.edu/~yyye/matlab/) of SOLNP by Professor Ye.
- Stellato, B., Banjac, G., Goulart, P., Bemporad, A., & Boyd, S. (2020). OSQP: An operator splitting solver for quadratic programs. Mathematical Programming Computation, 12(4), 637-672.
