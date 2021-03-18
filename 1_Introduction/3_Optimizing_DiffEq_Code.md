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

# Optimizing DiffEq Code

在本Notebook中，我们将逐步介绍一些用于优化代码的主要工具，以便有效地求解DifferentialEquations.jl。用户端优化很重要，因为对于足够困难的问题，大部分时间将花费在您要解决的函数`f`内。“高效的”积分器是那些减少所需的`f`调用次数以达到容错能力的积分器。优化DiffEq代码或任何Julia函数的主要思路如下：

* 不做内存分配
* 对于小型数组使用StaticArrays
* 使用广播融合
* 使其类型稳定
* 减少冗余计算
* 利用BLAS调用
* 优化算法选择

我们将在小型和大型系统的背景下讨论这些策略。让我们从小型系统开始。


## Optimizing Small Systems (<100 DEs)

让我们以之前的经典Lorenz系统为例。让我们从天真地以不合适的形式编写系统开始：

```julia
function lorenz(u,p,t)
 dx = 10.0*(u[2]-u[1])
 dy = u[1]*(28.0-u[3]) - u[2]
 dz = u[1]*u[2] - (8/3)*u[3]
 [dx,dy,dz]
end
```

这里，lorenz返回一个对象`[dx,dy,dz]`，它是在lorenz的内部创建的。

这是高级语言中常见的代码模式，如MATLAB、SciPy或R的deSolve。然而，这种形式的问题在于它在每个步骤都分配了一个向量`[dx,dy,dz]`。让我们用这个函数作为求解过程的测试基准：

```julia
using DifferentialEquations, BenchmarkTools
u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz,u0,tspan)
@benchmark solve(prob,Tsit5())
```

BenchmarkTools程序包的`@benchmark`通过多次运行该代码以获取准确的测量结果。最短时间是指您的操作系统和其他后台进程没有进行妨碍时所需要的时间。注意，在这种情况下，大约需要5ms来求解并分配11MB内存。然而，如果我们要在真正的用户代码中使用它，我们将看到大量时间花费在了垃圾收集（GC）上，以清理我们创建的所有数组。即使关闭保存，我们也有这些内存分配。

```julia
@benchmark solve(prob,Tsit5(),save_everystep=false)
```

当然，问题在于每次调用派生函数时都会创建数组。此功能每步被多次调用，因此是内存使用的主要来源。为了解决这个问题，我们可以使用in-place的形式来**使我们的代码不分配内存**：

```julia
function lorenz!(du,u,p,t)
 du[1] = 10.0*(u[2]-u[1])
 du[2] = u[1]*(28.0-u[3]) - u[2]
 du[3] = u[1]*u[2] - (8/3)*u[3]
end
```

这里，我们不是每次都创建一个数组，而是利用了缓存数组`du`。当使用inplace形式时，DifferentialEquations.jl采用不同的内部路由，该路由也使内部分配最小化。当我们对该函数进行基准测试时，我们将看到很大的不同。

```julia
u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz!,u0,tspan)
@benchmark solve(prob,Tsit5())
```

```julia
@benchmark solve(prob,Tsit5(),save_everystep=false)
```

仅此变化就有数倍的时差！请注意，仍有一些内存分配，这是由于集成缓存的构造所致。但这并不会挥着问题的规模而扩展：

```julia
tspan = (0.0,500.0) # 5x longer than before
prob = ODEProblem(lorenz!,u0,tspan)
@benchmark solve(prob,Tsit5(),save_everystep=false)
```

因为那只是设置时的分配。


#### But if the system is small we can optimize even more

只有当分配是“堆分配”时，分配才是昂贵的。对于更深入的堆分配定义，[有很多在线资源](http://net-informations.com/faq/net/stack-heap.htm)。但是，一个很好的工作定义是，堆分配是必须由指针指向的可变大小的内存块，这种指针间接操作会花费时间。此外，必须对堆进行管理，垃圾控制器必须主动跟踪堆中的内容。

然而，有一种替代堆分配的方法，称为栈分配。栈是静态大小的（在编译时已知），因此它的访问速度很快。此外，编译器提前知道确切的内存块，因此重用内存很便宜。这意味着在栈上进行分配基本上没有成本！

数组必须在堆中分配，因为它们的大小（以及它们占用的内存数量）是在运行时确定的。但是在Julia中有一些结构是栈分配的。例如，`struct`是栈分配的值类型。`Tuple`是一个栈分配的集合。对于DiifEq而言，最有用的数据结构是来自软件包[StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl)中的`StaticArray`。这些数据的长度在编译时确定。它们是使用附加到普通数组表达式的宏创建的，例如：

```julia
using StaticArrays
A = @SVector [2.0,3.0,5.0]
```

注意，`SVector`后面的`3`给出了`SVector`的大小，无法更改。另外，`SVector`是不可变的，因此我们必须创建一个新的`SVector`来更改值。但是请记住，我们不必担心分配，因为此数据结构是栈分配的。`SArray`还有很多额外的优化功能：它们有快速的矩阵乘法，快速的QR分解等。它们直接利用了有关数组大小的信息。因此，应尽可能使用它们。

不幸的是，静态数组只能用于足够小的数组。在达到一定大小后，它们将在执行某些指令及其时限之后被迫堆分配。因此，如果系统具有超过100个变量，则不应使用静态数组。此外，只有本机的Julia算法才能充分利用静态数组。

让我们**使用静态数组来优化lorenz**。注意，在这种情况下，我们想使用out-of-place分配形式，但是这次我们要输出一个静态数组：

```julia
function lorenz_static(u,p,t)
 dx = 10.0*(u[2]-u[1])
 dy = u[1]*(28.0-u[3]) - u[2]
 dz = u[1]*u[2] - (8/3)*u[3]
 @SVector [dx,dy,dz]
end
```

为了让求解器在内部使用静态数组，我们只需要给它一个静态数组作为初始条件：

```julia
u0 = @SVector [1.0,0.0,0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz_static,u0,tspan)
@benchmark solve(prob,Tsit5())
```

```julia
@benchmark solve(prob,Tsit5(),save_everystep=false)
```

差不多就这些了。对于静态数组，您不必担心分配问题，因此使用像`*`这样的操作，也不必担心融合操作（将在下一节中讨论）。使用R/MATLAB/Python的“向量化代码”，或是直接使用数字/值，您的代码在这种情况下会很快。

```julia

```

```julia

```

```julia

```

```julia

```

```julia

```
