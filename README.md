# 18.330-Problem-set-8
18.330 Problem set 8

**Download Link:https://programming.engineering/product/18-330-problem-set-8/**


Description
5/5 – (2 votes)
Exercise 1 : Eigenvalue solvers for a special matrix

In this problem we will consider a very special symmetric matrix. Recall that the second-order finite difference scheme for a function ( ) is given by

″

=

+1 − 2 + −1

ℎ0 .

= ( 0

+ )

ℎ =

2

−

where

and

Define the vector f such that f = . Then the corresponding derivative vector using finite differences can be written as f″ = f . Some care must be taken about the boundary conditions; for simplicity we will assume that 0 = = 0. Under these assumptions the matrix is a tridiagonal matrix:

′

−2

1

1

1′

2′

1−21

2

3

1

−2 1

3

′ ⋮

=

⋱ ⋱ ⋱

1

⋮

′

1

−2

−2

−2

1

−1

−2 −1

You can construct the matrix in Julia using the diagm function. For simplicity keep everything dense in this problem, although there clearly a lot of structure that could be exploited.

Consider the matrix . Using the ansatz (i.e. hypothesis) v = sin( / ), where 1 ≤ ≤ − 1, show that v is an eigenvector of . What is the corresponding eigenvalue? Remember that the eigenvalues should be real!

In class we discussed several ways to find eigenvalues of matrices. The simplest algorithm for finding eigenvalues is the power method. Supppse that a symmetric matrix has the eigendecomposition = Λ −1. Recall that for a symmetric matrix the eigenvectors associated with distinct eigenvalues are orthogonal. Starting with a random intitial vector x which can be written as a sum over the eigenvectors,

x = ∑ v


show that

= ∑

x v

where we order the eigenvalues by their magnitude, | 1| > | 2| > ⋯ > | | and hence

x = 1 1 (v1 + ( 2 v2)) .

1

This shows why the power method converges to the leading eigenvector.

̃

3. The power method gives an approximation for the leading eigenvector v1.

Since it is only an approximation it is not necessarily true that v ̃ =

1

̃

1v1 exactly.

Instead we will say that the best approximation to 1 is given by the that satifies

1 =

min

v

v

v

v

2

||

1̃−

1̃||2 = ∑ (∑ (

1̃)−(

1̃)) .

v

By differentiating this

expression with respect to

, show that

v

v

v

v

v

v⊤ v

1̃⋅ ( 1̃)

1̃ 1̃

1 ≈

=

⊤

.

1̃⋅

1̃

1̃

1̃

This is called the Rayleigh quotient.

Implement the power method power_method(A, N, x) for a symmetric matrix which iterates for a fixed number of iterations on the intial

vector x and returns ṽand approximated by the Rayleigh quotient

1 1

above.

Run this on the matrix 10 which is of size (9 × 9). Use the true largest magnitude vector from part 1. Plot the relative error between your result and the true value as a function of to get a smooth plot; use the same initial vector each time.

Remember to normalize the vector at each iteration to avoid overflow! This initial random vector should be complex. (Extra credit: show that the relationship is what you would expect analytically!)


A more advanced method that we discussed to find all the eigenvalues was the QR algorithm. We define (0) = and then

() ()= ()

(+1)= () ()

This will normally converge to an upper-triangular matrix. How does this work? We call two matrices and similar if = ⊤ where is orthogonal. Show that if has eigenvalue pairs ( , v ) then has eigenpairs ( , v ).

Show that in the QR algorithm ( +1) is similar to ( ) and hence . Therefore if the algorithm converges to an upper-triangular matrix we can read off the eigenvalues from the diagonal and the eigenvectors will be given by the columns of (1) (2) (3) ⋯ ( ).

Implement the QR algorithm in a function qr_algorithm(A, N) for a ma-trix with iterations. It should return the resulting matrix ( ) and the eigenvector matrix. Run the algorithm on a random symmetric matrix of size 10 × 10 for 100 iterations. How close to upper-triangular is the resulting matrix?

Run the QR algorithm on the matrix 10 of size (9 × 9) for a suitable number of iterations to get convergence. Compare the results to the the-oretical one. Can you find all the eigenvalues and eigenvectors?

Using the fact that the eigendecomposition of a symmetric matrix gives orthogonal matrices (which are easy to invert) propose a method to solve the linear system

x = b

Solve the system for b = 0.01^2*sin.(2π*(0.01:0.01:0.99)) and

as a 99 × 99 matrix. Plot the resulting vector x. Plot on the same axis sin.(2π*(0.01:0.01:0.99))/4π^2. Is it similar to x? We have just solved the boundary value problem ″ = sin(2 ) with (0) = (2 ) = 0. In the next problem set we will see how to do this even quicker using Fourier analysis.

Exercise 2: Low-rank approximation

In this problem we will use the SVD to produce low-rank approximations for matrices. We will then use this to compress images and write fast multiplication algorithms.


Plot the singular values as a function of . What do they tell us? What is a suitable rank- approximation to take? What does it look like?

Now let’s apply these ideas to a real image. Use the Images.jl package to load a test image using the following code. (Remember that you may need to install the relevant packages):

using Images, TestImages

img = testimage(“mandrill”)

imgmat = Array(channelview(Float64.(Gray.(img))))

heatmap(imgmat, c=ColorGradient(:grays), showaxis=false, clims=(0,1), yflip=true)

Here imgmat is a standard Julia matrix that contains greyscale pixel val-ues.

Plot the rank-50 approximation to imgmat (using heatmap). How does it compare to the original?

Plot the singular values as a function of . You should see an “elbow” shape. Approximately where does it occur?

Create an interactive visualization that shows the low-rank approximation for different values of . What do you observe? After which are you happy with the quality of the image. Is it related to where the elbow is?

Exercise 3: Dynamic mode decomposition

In this problem we will use the SVD to predict the future!

Suppose that we have some data that satisfies an ODE,

ẋ= f(x)

We can solve the ODE and take time snapshots of the result which we can stack into a matrix as follows:

,

= [

x

|

x

|

x

|

( )

( +1) ⋯

( )]

|

|

|

We will assume that the dynamics is linear. This means that we expect that two consecutive snapshots should satisfy

x +1 = x

for some matrix .


If we have snapshots then we can write this as a matrix equation,

2, = 1,( −1)

Suppose that we have an eigen-decomposition for , i.e. = Λ −1, and that we can write the initial time snapshot as a sum over the eigenbasis of ,

x0 = ∑ w

where w is the th eigenvector of . Show that future time snapshots can be written as

x = ∑ w

where is the th eigenvalue of .

We are now going to try and find the eigendecomposition of without ever constructing it! Start by calculating the SVD of 1,( −1) = Σ ⊤. Find an expression for ⊤ in terms of 2, , Σ, , . We can then calculate the eigenspectrum of A by taking the eigenvalue decomposition of ⊤ since they are similar and using use the result in [1.5].

Write a function that calculates the eigenspectrum of given 1, . In-stead of using the full SVD, use a truncated SVD keeping only the first singular values and singular vectors; i.e. use the matrices Vr = V[:, 1:r], Ur = U[:, 1:r] and Σr = Σ[1:r, 1:r] in the expression above for , , Σ. Your function should be of the form aspectrum(X, r). The reason for truncating is in case A is not full rank, in which case some terms in Σ might be nearly 0 and dividing by them would be a bad idea.

Test your function on a random 10 × 10 matrix generating some data for of size (10×11) from it starting from a random 0. Compare the re-sults of eigen(A) with aspectrum(X, 10). Remember that eigenvectors are the same upto a multiplicative constant.

We are now going to apply this to some dynamical data. Use the ODE integrator code below to generate some data for coupled oscillators with 1(0) = 0.5 and all the others (0) = 0 for N = 10 from = 0 to = 10. The coupled system can be written as

1̈

−2

1

1

1

2̈

1

−2

2

3

=

1

−2 1

3

̈

⋮

⋱ ⋱ ⋱

⋮

̈

1

−2


or as a system of first order equations,

1

1̇

0

1

1

1̇

−2

0

1

1

2̇

0

0

0

2

2

1

0

−2

0

1

2

̇

3

3

=

0

0

0

0

1

3̇

1

3

̇

1

0

−2

0

⋮

⋱⋱⋱⋱⋱

⋮

0

0

0

0

1

The output

̇

̇

1

0

−2

0

1, 1̇ , 2,

2̇ , ⋯

is a data matrix with rows

Generate a plot of 1( ) to check that everything went according to plan.

Split the data into two parts 1 and 2, the first half from = 0 to = 5 = 5 to = 10. Calculate the spefctrum of withandthesecondhalf

= 10 using 1.

Use the first column of 2 as the initial condition. Use the spectrum you found to predict the future dynamics. [Hint: use the initial condition to find the s, which is a matrix solve. Then use the equations in part 1.1 to calculate the prediction.]

Plot the prediction for the 10 springs on the same axis as the true solution. What happens?

Repeat [3.6–3.7] for = 15 and = 20. What do you observe?

ODE Code:

struct RKMethod{T}

c::Vector{T}

b::Vector{T}

a::Matrix{T}

s::Int

# Make sure that the matrices and vectors have the correct size

function RKMethod(c::Vector{T}, b::Vector{T}, a::Matrix{T}, s::Int) where T lc = length(c); lb = length(b); sa1, sa2 = size(a)

if lc == sa1 == sa2 == s-1 && lb == s

new{T}(c, b, a, s)

else

error(“Sizes should be (s = $s) : \n length(c) = s-1 ($lc) \n length(b) = s

end

end

end


function (method::RKMethod)(f, x, t, h)

# Extract the parameters

a, b, c, s = method.a, method.b, method.c, method.s

Vector to hold the k terms k = [f(t, x)]

for i in 1:s-1

tn = t + c[i]*h

xn = x + h*sum(a[i, j] * k[j] for j in 1:i)

push!(k, f(tn, xn))

end

return x + h * sum(b.*k)

end

function integrate(method, f, x0, t0, tf, h)

Calculate the number of time steps N = ceil(Int, (tf – t0) / h)

hf = tf – (N – 2)*h

#initiate tracking vectors

xs = [copy(x0)]

ts = [t0]

#step

for i in 1:N-1

push!(xs, method(f, xs[i], ts[i], h))

push!(ts, ts[i] + h)

end

# Special last step

push!(xs, method(f, xs[N-1], ts[N-1], hf))

push!(ts, ts[N-1] + hf)

return ts, xs

end

c = [1/2, 1/2, 1]

b = (1/6) .* [1, 2, 2, 1]

a = [1/2 0 0; 0 1/2 0; 0 0 1]

rk4 = RKMethod(c, b, a, 4)

function build_springmat(N)

springmat = zeros(2N,2N)

for i = 1:2N


(i < 2N) && (springmat[i, i+1] = 1.0)

if iseven(i)

springmat[i, i-1] = –2.0

(i > 3) && (springmat[i, i-3] = 1.0)

end

end

springmat

end

N=10

const spmat = build_springmat(N)

spf(t, x) = spmat*x

x0 = zeros(2N); x0[1] = 0.5

ts, xs = integrate(rk4, spf, x0, 0.0, 20.0, 0.005); X = hcat(xs…)

plot(ts, X[1:2:2N, :]’, leg=:outertopright, box=:on)

Exercise 4: Fourier integrals

In lectures we saw that we could write periodic functions in the form

̂

( ) = ∑

=−

where ̂= 1 ∫2 ( ) − .

2 0

1. Consider the saw-tooth function,

( ) =

mod

0 ≤ <

{2 −

≤ <2

{

( , 2 ))

(

else

Calculate the Fourier series coefficients analytically.

Write a function fourier_coefficients(f::Vector, n::Int) that takes in a vector of samples of uniformly distributed over [0, 2 ] and returns an approximation to ̂by calculating the integral using the trapezoidal rule.

Now calculate ̂using your trapezoidal code using 100 points and = −3, −2, … , 3. How do they compare to the theoretical results?

Fix to be 1. Plot the relative error between the theoretical result and the result from using fourier_coefficients for calculating 1̂using a number of points between 10 and 1000. What does the convergence look like? Does it agree with what we discussed in class?


Now consider the smooth periodic function exp(cos( )). Repeat [4.3– 4.4] for this function. The analytical result is given by ̂ = | |(1). You can calculate this using the besseli(abs(n), 1) in SpecialFunc-tions.jl.

Plot the magnitude of ̂as a function of for the two functions. How do the coefficients decay? Is this what you expected?

10

