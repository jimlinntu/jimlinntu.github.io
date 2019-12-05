---
layout: post
title: Dive into Fast Fourier Transform
---

Recently, I saw this codeforce problem <https://codeforces.com/problemset/problem/1257/G>, 
so I start to figure out how FFT works.

# Iterative FFT

Let $$ a $$ be the original sequence and $$ b $$ is the sequence after divide and conquer,

Then we can observe that:

$$ b_i = a_{reverse(i)} \iff b_{reverse(i)} = a_i $$

Let's initialize $$ j = 0 $$ index for $$ a $$ and $$ i = 0 $$ index for $$ b $$,
note that we can accelerate the procedure by incremental and in-place fashion, that is

$$ b_{i = 0} = a_{reverse(i) == 0} $$


Given the fact that:

$$  b_{i + 1} = a_{reverse(i+1)} $$

, we can actually exploit the previous step $$ j == reverse(i) $$ result.

Let us rewrite it to: $$ reverse(j) == i $$, then we can plug $$ j $$ into the update step.

$$ b_{i + 1} = a_{reverse(reverse(j) + 1)} $$

Though $$ reverse(reverse(j) + 1)) $$ seems crazy, it can actually be implemented fast by bitwise operations in `reverse_add` below.

Horay! Let's see the code

~~~c++
int reverse_add(int x, int bit_length)
{
    // try to find the most first bit that is 0
    for(int l = 1 << bit_length; (x ^= l) < l; l >>= 1);
    return x;
}
~~~

~~~c++
void bit_reverse(int n, complex_t *x)
{
    int bit_length = (int) log2(n);
    for(int i = 0, j = 0; i != n; ++i)
    {
        if(i > j) swap(x[i], x[j]);
        j = reverse_add(j, bit_length);
    }
}
~~~

Let us try this example:

$$
    w^{0}_{8}, w^{1}_{8}, w^{2}_{8}, w^{3}_{8}, w^{4}_{8}, w^{5}_{8}, w^{6}_{8}, w^{7}_{8}
$$

where $$ w^{k}_{n} = e^{\frac{2 \pi i k}{n}} $$.

And a polynomial ( I use the vector notation introduced by MIT open course):

$$
    A(x) = a_0 x^0 + a_1 x^1 + a_2 x^2 + .... a_7 x^7 = < a_0, a_1, ..., a_7 >(x)
$$

Let's see the computational graph:

![fft](/img/blog/fft.jpg)


Code:


~~~c++
void transform(int n, complex_t *x, complex_t *w)
{
    // Prepare the leaf nodes ( you can see this as a bottom-up approach )
    bit_reverse(n, x);
    for(int stride = 2; stride <= n; stride <<= 1)
    {
 
        int half_stride = stride >> 1;
        // the start index of a stride region
        for(int start = 0; start < n; start += stride)
        {
            // this loop will finish up [start, start + stride] region
            for(int k = 0; k < half_stride; ++k)
            {
                complex_t z = x[start + half_stride + k] * w[n / stride * k];
                // Use the phase property (see below)
                x[start + half_stride + k] = x[start + k] - z;
                x[start + k] = x[start + k] + z;
            }
        }
    }
}
~~~

Note, in the code,

$$
    z = ... \cdot w^{k}_{stride} =  ... \cdot  w^{k * \frac{n}{stride}}_{n}
$$

and the phase property:

$$
    w^{k + \frac{n}{2}}_{n} = e^{\frac{2 \pi i}{n} \cdot \frac{n}{2}} \cdot w^{k}_{n} = - w^{k}_{n}
$$



# References
* <http://blog.miskcoo.com/2015/04/polynomial-multiplication-and-fast-fourier-transform>
* [MIT 3. Divide & Conquer: FFT](https://www.youtube.com/watch?v=iTMn0Kt18tg)
* [Number Theoretic Transform](https://ccrma.stanford.edu/~jos/mdft/Number_Theoretic_Transform.html)
