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

# An Intro to DifferentialEquations.jl


## Basic Introduction Via Ordinary Differential Equations


本笔记将通过向您介绍用于求解微分方程（ODE）的功能，来帮助您开始使用DifferentialEquations.jl。相应的英文文档页面是：[ODE tutorial](https://docs.sciml.ai/dev/tutorials/ode_example/)。尽管某些语法可能与其他类型的方程不同，但在每种情况下都遵循相同的一般原则。我们的目标是提供一种温和而透彻的介绍，以突出这些原则，从而帮助您概括所学的内容。


### Background


如果您不熟悉微分方程，那么快速阅读一下[常微分方程的定义](https://en.wikipedia.org/wiki/Ordinary_differential_equation)可能会有所帮助。我们将常微分方程定义为描述变量$u$变化方式的方程，即：

$$u' = f(u,p,t)$$

其中，$p$是模型的参数，$t$是时间变量，$f$是$u$如何变化的非线性模型。另外还有初始值的问题，它包含了关于起始值的信息：

$$u(t_0) = u_0$$

总之，如果您知道起始值并且知道该值是如何随时间变化的，那么您就会知道该值在未来的任意时刻将是什么。这是微分方程的直观定义。


### First Model: Exponential Growth

<!-- #region -->
我们的第一个模型将是典型的指数增长模型。这个模型中变化率与当前值成比例，就像这样：

$$u' = au$$

其中，我们有一个起始值$u(0)=u_0$。假设我们把1美元投入比特币，比特币以每年98%的速度增长。现在调用$t=0$，并以年为单位测量时间，我们的模型是：

$$u' = 0.98u$$

并且$u(0) = 1.0$，我们注意到，在这个设置中我们将其编码为Julia代码的形式与常规表示形式相匹配。

```julia
f(u,p,t) = 0.98u
```

其中$u_0 = 1.0$。如果想要在`t=0.0`到`t=1.0`的时间范围内求解此模型，则可以通过指定这个函数`f`、这个初始条件`u0`和这个时间跨度来定义一个`ODEProblem`，如下所示：
<!-- #endregion -->

```julia
using DifferentialEquations
f(u,p,t) = 0.98u
u0 = 1.0
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
```

为了解这个ODEProblem，我们使用solve命令。

```julia
sol = solve(prob)
```

就这样，我们成功地解决了我们的第一个ODE！


#### Analyzing the Solution


当然，solution类型（这里的sol变量）本身并不有趣。我们想要深入了解solution！详细说明分析solution的功能的文档页面为[Solution Handling](https://docs.sciml.ai/dev/basics/solution/)。这里我们将介绍一些基础知识。您可以使用Plots.jl提供的绘图方法绘制solution：

```julia
using Plots; gr()
plot(sol)
```

从图中我们可以看到，solution是一条与我们的直觉相符合的指数曲线。作为绘图方法，我们可以使用任何Plots.jl的属性来注释结果。例如：

```julia
plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
     xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!") # legend=false
```

使用可变的plot!命令，我们可以添加其他部分到我们的图中。对于此ODE，我们知道真正的解是$u(t) = u_0 exp(at)$，因此让我们在图中添加这一真实解：

```julia
plot!(sol.t, t->1.0*exp(0.98t),lw=3,ls=:dash,label="True Solution!")
```

在上一个命令中，我演示了`sol.t`，其获取并保存了solution中的时间点数组：

```julia
sol.t
```

我们可以使用`sol.u`获取solution中的值数组：

```julia
sol.u
```

`sol.u[i]`是时间`sol.t[i]`时解的值。我们可以使用标准的推导表达式来计算solution值的函数数组，例如：

```julia
[t+u for (u,t) in tuples(sol)]
```

然而，一个有趣的特性是，默认情况下，solution是一个连续函数。如果我们再次检查打印输出：

```julia
sol
```

你可以看到，solution有一个顺序变化的插值。默认算法会自动在方法之间切换，以处理所有类型的问题。对于非刚性方程（如我们正在求解的方程），它是4阶精度的连续函数。我们可以将解称为时间`sol(t)`的函数。例如，要获得`t=0.45`处的值，可以使用以下命令：

```julia
sol(0.45)
```

#### Controlling the Solver


DifferentialEquations.jl在其算法中有一组通用的求解器控件，可以在[Common Solver Options](https://docs.sciml.ai/dev/basics/common_solver_opts/)页面上找到。我们将详细介绍一些最广泛使用的选项。

最有用的选项是公差`abstol`和`reltol`。这两个参数告诉内部自适应时间步进引擎您想要的solution的精度。一般来说，`reltol`是相对精度，而`abstol`是`u`接近零时的精度。这些公差是局部公差，因此不做全局保证。但是，一个好的经验法则是，总的solution方案精度比相对公差小1-2位。因此，对于默认值`abstol=1e-6`和`reltol=1e-3`，您可以期望大约1-2位的全局精度。如果我们想获得6位左右的精度，可以使用以下命令：

```julia
sol = solve(prob,abstol=1e-8,reltol=1e-8)
```

现在，我们看不出与真实解有明显的区别：

```julia
plot(sol)
plot!(sol.t, t->1.0*exp(0.98t),lw=3,ls=:dash,label="True Solution!")
```

请注意，通过减小公差，求解器必须执行的步骤数为9，而不是之前的5。在精度和速度之间进行权衡，由您来确定什么是解决问题的正确平衡。

另一个常见的选择是使用`saveat`来使求解器在特定时间点保存。例如，如果我们希望在整数`k`的`t=0.1k`的均匀网格处求解，则可以使用以下命令：

```julia
sol = solve(prob,saveat=0.1)
```

注意，当使用`saveat`时，不再保存连续输出变量，因此插值的`sol(t)`仅为一阶。我们可以通过向`saveat`传递一组值再不均匀的网格点中进行保存，例如：

```julia
sol = solve(prob,saveat=[0.2,0.7,0.9])
```

如果需要减少存储，我们还可以通过`dense=false`直接关闭连续输出：

```julia
sol = solve(prob,dense=false)
```

要关闭所有中间保存，我们可以使用`save_everystep=false`：

```julia
sol = solve(prob,save_everystep=false)
```

如果我们想求解并仅保存最终值，我们甚至可以设置`save_start=false`。

```julia
sol = solve(prob,save_everystep=false,save_start=false)
```

请注意，类似地，在另一端也有类似的`save_end=false`。

可以通过[Callback Library](https://docs.sciml.ai/dev/features/callback_library/#saving_callback-1)中的`SavingCallback`处理更高级的保存行为，例如保存solution的函数，这将在本教程的后面部分讨论。


#### Choosing Solver Algorithms


对于微分方程，没有最好的数值求解算法。当您调用`solve(prob)`时，DifferentialEquations.jl会根据您要求的属性（公差，要保存的信息等）来为您的问题找到一个好的算法。但是，在许多情况下，您可能需要更直接的控制。稍后的笔记将帮助介绍DifferentialEquations.jl中的各种*算法*，但现在让我们先介绍一下语法。

选择数值方法的最关键决定因素是模型的刚度。刚度大致由特征值是一个具有大特征值的雅可比矩阵`f`。这是相当数学化的，我们可以更直观地想象：如果您在`f`中有很大的数字（例如`1e5`阶的参数），那么他可能是刚性的。或者，正如MATLAB ODE套件的常见这Lawrence Shampine喜欢定义的那样，如果标准算法很慢，那么它就是刚性的。我们将在后面的教程中更深入地了解对刚性的判断，但是现在请注意，如果您认为模型可能是刚性的，则可以通过`alg_hints = [:stiff]`提示算法选择器。

```julia
sol = solve(prob,alg_hints=[:stiff])
```

刚性算法必须在每个步骤中求解隐式方程和线性系统，因此仅应在需要时使用它们。

如果我们要直接选择一种算法，则可以在问题参数prob之后传递算法类型，就像`solve(prob,alg)`。例如，让我们使用`Tsit5()`算法来求解此问题，为了便于演示，我们同时将相对公差更改为`1e-6`。

```julia
sol = solve(prob,Tsit5(),reltol=1e-6)
```

### Systems of ODEs: The Lorenz Equation


现在，让我们转到ODE系统。洛仑兹方程是催生混沌理论的著名“蝴蝶吸引子”。它由ODE系统定义：

$$
\begin{align}
\frac{dx}{dt} &= \sigma (y - x)\\
\frac{dy}{dt} &= x (\rho - z) -y\\
\frac{dz}{dt} &= xy - \beta z
\end{align}
$$

为了在DifferentialEquations.jl中定义一个微分方程组，我们把`f`定义为具有向量初始条件的向量函数。因此，对于向量`u = [x,y,z]`，我们有导数函数：

```julia
function lorenz!(du,u,p,t)
    σ,ρ,β = p
    du[1] = σ*(u[2]-u[1])
    du[2] = u[1]*(ρ-u[3]) - u[2]
    du[3] = u[1]*u[2] - β*u[3]
end
```

注意，这里我们使用了in-place格式，该格式将输出写入到预先分配的向量`du`中。对于方程组，in-place格式更快。我们使用初始条件$u_0 = [1.0,0.0,0.0]$，如下所示：

```julia
u0 = [1.0,0.0,0.0]
```

最后，对于该模型，我们使用了参数p，我们也需要在ODEProblem中设置这个值。对于我们的模型，我们想使用参数$\sigma = 10$、$\rho = 28$和$\beta = 8/3$进行求解，因此我们建立了参数集合：

```julia
p = (10,28,8/3) # we could also make this an array, or any other type!
```

现在我们生成ODEProblem类型。在本例中，因为我们有参数，所以我们将参数值添加到构造函数调用的末尾。让我们在`t=0`到`t=100`的时间范围内解这个问题：

```julia
tspan = (0.0,100.0)
prob = ODEProblem(lorenz!,u0,tspan,p)
```

现在，和之前一样，我们solve这个问题：

```julia
sol = solve(prob)
```

同样的solution处理特性也适用于这种情况。因此，`sol.t`存储时间点，而`sol.u`是一个数组，存储相应时间点的solution。

然而，在处理方程组时，还有一些额外的特性是需要知道的。首先，`sol`还充当数组的角色，`sol[i]`返回第`i`个时间点的解。

```julia
sol.t[10],sol[10]
```

另外，solution类似于矩阵，其中`sol[j,i]`是第`j`个变量在第`i`个时刻的值：

```julia
sol[2,10]
```

我们可以通过执行转换来获得一个真正的矩阵：

```julia
A = Array(sol)
```

这与sol相同，即`sol[i,j] = A[i,j]`，但现在它是一个真正的矩阵。默认情况下，绘图将显示每个变量的时间序列：

```julia
plot(sol)
```

如果我们想要绘制相互对照的值，我们可以使用`vars`命令。让我们将变量`1`相对变量`2`相对变量`3`进行绘制：

```julia
plot(sol,vars=(1,2,3))
```

这是经典的洛伦兹吸引子图，其中`x`轴为`u[1]`、`y`轴为`u[2]`、`z`轴为`u[3]`。注意，默认情况下绘图方法会使用插值，但是我们可以将其关闭：

```julia
plot(sol,vars=(1,2,3),denseplot=false)
```

哎呀！这显示了通过仅计算稀疏解并填充值来计算连续解可以节省大量计算工作！请注意，在vars中，`0=time`，因此我们可以绘制单个组件的时间序列，如下所示：

```julia
plot(sol,vars=(0,2))
```

## Internal Types


要研究的最后一个基本用户界面特性是类型的选择。DifferentialEquations.jl根据输入类型来确定所使用的内部类型，因此，由于在前面的例子中，当我们使用`Float64`值作为初始条件时，这意味着内部值将使用`Float64`来求解。我们确保时间是通过`Float64`值指定的，这意味着时间步也将使用64位浮点数。但是，通过简单地改变这些类型，我们可以改变内部使用的东西。

举个简单的例子，假设我们要解决由矩阵定义的ODE。为此，我们可以简单地使用矩阵作为输入。

```julia
A  = [1. 0  0 -5
      4 -2  4 -3
     -4  0  0  1
      5 -2  2  3]
u0 = rand(4,2)
tspan = (0.0,1.0)
f(u,p,t) = A*u
prob = ODEProblem(f,u0,tspan)
sol = solve(prob)
```

这与我们之前所作的没有什么区别，但是在这个例子中，`u0`是一个`4x2`的矩阵。因此，每个时间点的solution都是矩阵：

```julia
sol[3]
```

在DifferentialEquations.jl中，您可以使用定义了`+`、`-`、`*`、`/`，并具有适当范数的任何类型。例如，如果我们要使用任意精度的浮点数，则可以将输入更改为`BigFloat`的矩阵：

```julia
big_u0 = big.(u0)
```

我们可以使用该初始条件来求解具有任意精度的`ODEProblem`。

```julia
prob = ODEProblem(f,big_u0,tspan)
sol = solve(prob)
```

```julia
sol[1,3]
```

要真正利用这一点，我们要把`abstol`和`reltol`变小！请注意，“时间”的类型与因变量的类型不同，这样可以通过保持多个精度来优化算法。通过使用`BigFloat`变量定义时间范围，我们也可以将时间转换为任意精度：

```julia
prob = ODEProblem(f,big_u0,big.(tspan))
sol = solve(prob)
```

最后，让我们展示一个更复杂的类型用法。对于小型数组，通过包[StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl)对静态数组进行操作通常会更快。语法类似于普通数组，但是对于这些特殊的数组，我们使用`@SMatrix`宏来表示我们想要创建一个静态数组。

```julia
using StaticArrays
A  = @SMatrix [ 1.0  0.0 0.0 -5.0
                4.0 -2.0 4.0 -3.0
               -4.0  0.0 0.0  1.0
                5.0 -2.0 2.0  3.0]
u0 = @SMatrix rand(4,2)
tspan = (0.0,1.0)
f(u,p,t) = A*u
prob = ODEProblem(f,u0,tspan)
sol = solve(prob)
```

```julia
sol[3]
```

## Conclusion


这些是DifferentialEquations.jl中的基本控件。所有方程都是通过一个problem类型来定义的，`solve`命令与算法选择（或使用默认算法）一起使用以获取solution。每个solution的行为都相同，例如带有`sol.t[i]`的数组`sol[i]`，以及带有绘图命令`plot(sol)`的连续函数`sol(t)`。通用求解器选项可用于控制任意方程类型的求解。最后，数值求解中使用的类型由输入类型确定，并且可以用于任意精度求解并添加其他优化（例如，可以用于通过GPU求解！）。虽然这在ODE上得到了证明，但这些技术也可以推广到其他类型的方程上。
