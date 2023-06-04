# Practical 3: Smoothing with basis expansions and penalties
Smoothing and function estimation play important parts in applied statistics and data science. One approach
combines basis expansions and penalized regression, both techniques with much wider application. In this practical
you (working individually) will write R functions for smoothing x, y data. The idea is that you have a model:
$$y_i = f(x_i) + ϵ_i , i = 1, . . . , n$$
where $x_i$ and $y_i$ are observed, f is an unknown smooth function, and $ϵ_i$ a zero mean error term with variance $σ^2$.
