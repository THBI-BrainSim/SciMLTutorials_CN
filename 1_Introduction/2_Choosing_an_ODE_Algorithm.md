---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.10.0
  kernelspec:
    display_name: Julia 1.5.3
    language: julia
    name: julia-1.5
---

# Choosing an ODE Algorithm

虽然默认算法以及`alg_hints = [:stiff]`在大多数情况下就足够了，但有时您可能需要施加更多的控制。本教程的这一部分的目的是向您介绍一些最广泛使用的算法选择以及应该在什么时候使用它们。文档中相应的页面是[ODE Solvers](https://docs.sciml.ai/dev/solvers/ode_solve/)页面，该页面更加深入。


## Diagnosing Stiffness

关于算法选择，需要知道的一个关键问题是你要求解的问题是否是刚性的。让我们以Van Der Pol方程为例：

```julia
using DifferentialEquations, ParameterizedFunctions
van! = @ode_def VanDerPol begin
  dy = μ*((1-x^2)*y - x)
  dx = 1*y
end μ

prob = ODEProblem(van!,[0.0,2.0],(0.0,6.3),1e6)
```

一个提示该模型可能是刚性的因素是参数`1e6`。参数大通常意味着刚性模型。如果我们尝试使用默认方法来求解此问题：

```julia
sol = solve(prob,Tsit5())
```

在这里，显示以达到最大迭代次数。另一个可能发生的情况是求解器返回的solution是不稳定的（爆炸到无穷大或`dt`变得太小）。如果发生了这些情况，首先要做的是检查您的模型是否正确。很有可能是您犯了一个导致模型不稳定的错误！

如果是模型的问题，则刚性可能就是原因。因此，我们可以提示求解器使用适当的方法：

```julia
sol = solve(prob,alg_hints = [:stiff])
```

理解刚度的另一种方法是查看solution。

```julia
using Plots; gr()
sol = solve(prob,alg_hints = [:stiff],reltol=1e-6)
plot(sol,denseplot=false)
```

让我们放大y轴，看看发生了什么：

```julia
plot(sol,ylims = (-10.0,10.0))
```

注意一些极端的垂直移动是如何发生的。这些垂直位移是导数项非常大的地方，这表示刚度。这是一个突出这种行为的极端例子，但是可以将这一总体思想带到您的问题中。如有疑问，只需尝试同时使用刚性求解器和非刚性求解器进行计时，看看哪种效率更高。

为了验证这一点，让我们使用BenchmarkTools，这是一个使我们能相对可靠地对代码块进行计时的软件包。

```julia
function lorenz!(du,u,p,t)
    σ,ρ,β = p
    du[1] = σ*(u[2]-u[1])
    du[2] = u[1]*(ρ-u[3]) - u[2]
    du[3] = u[1]*u[2] - β*u[3]
end
u0 = [1.0,0.0,0.0]
p = (10,28,8/3)
tspan = (0.0,100.0)
prob = ODEProblem(lorenz!,u0,tspan,p)
```

现在，让我们使用基准工具中的`@btime`宏来比较此问题上使用非刚性和刚性求解器的情况。

```julia
using BenchmarkTools
@btime solve(prob);
```

```julia
@btime solve(prob,alg_hints = [:stiff]);
```

在这种特殊情况下，我们可以看到非刚性求解器使我们更快地找到了solution。


## The Recommended Methods

在选择一种方法时，一般规则如下：

* 较高的阶数在较低的公差范围内效率更高，较低的阶数在较高的公差范围内效率更高
* 在大多数现实情况下，适应性都是必不可少的
* Runge-Kutta方法适用于非刚性方程，Rosenbrock方法适用于较小的刚性方程，BDF方法适用于较大的刚性方程

虽然规则总有例外，但这些都是很好的指导原则。基于这些，选择求解算法的一个简单方法是：

* 默认是`Tsit5()`，这是5阶的非刚性Runge-Kutta方法
* 如果使用低公差（`1e-8`），请尝试使用`Vern7()`或`Vern9()`
* 如果使用高公差，请尝试`BS3()`
* 如果问题是刚性的，请尝试`Rosenbrock23()`、`Rodas5()`或`CVODE_BDF()`
* 如果不确定，请使用`AutoTsit5(Rosenbrock23())`或`AutoVern9(Rodas5())`

```julia

```
