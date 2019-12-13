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

# Recursively FFT
Here is the pseudocode:
~~~
fft(<a0, a1, a2, a3 ... , an-1>){
    result = [None] * n
    even = fft(<a0, a2, ... an-2>) // even[i] == A(w^{i}_{n/2})
    assert(len(even) == n/2)
    odd = fft(<a1, a3, ..., an-1>)
    assert(len(odd) == n/2)
    
    for(int i = 0; i < n; i++){
        // A(w^i_n) = A_even(w^{2i}_n) + w^i_n * A_odd(w^{2i}_n)
        // See the fact below.
        result[i] = even[i % (n/2)] + w^i_n * odd[i % (n/2)]
    }

    return result
}
~~~

$$
    (w^{i}_n)^2 = w^{2i}_n = w^{i}_{\frac{n}{2}} = w^{i \mod \frac{n}{2}}_{\frac{n}{2}}
$$

This is the code written by miskoo (See reference):

~~~c++
void fft(int n, complex<double>* buffer, int offset, int step, complex<double>* epsilon)
{
    if(n == 1) return;
    int m = n >> 1;
    fft(m, buffer, offset, step << 1, epsilon);
    fft(m, buffer, offset + step, step << 1, epsilon);
    for(int k = 0; k != m; ++k)
    {
        int pos = 2 * step * k;
        temp[k] = buffer[pos + offset] + epsilon[k * step] * buffer[pos + offset + step];
        temp[k + m] = buffer[pos + offset] - epsilon[k * step] * buffer[pos + offset + step];
    }
 
    for(int i = 0; i != n; ++i)
        buffer[i * step + offset] = temp[i];
}
~~~

For some people, this code is hard to understand at the first glance.
So let us visualize it!:

![recurseive-fft](/img/blog/fft-recursive.jpg)

Note that for simplicity, I denote:

$$
    A_{i:j:step}(x) = <a_i, a_{i+step}, ... >(x) \text{(Not include $j$)}
$$

# Number Theoretic Transform
TODO




# References
* <http://blog.miskcoo.com/2015/04/polynomial-multiplication-and-fast-fourier-transform>
* [MIT 3. Divide & Conquer: FFT](https://www.youtube.com/watch?v=iTMn0Kt18tg)
* [Number Theoretic Transform](https://ccrma.stanford.edu/~jos/mdft/Number_Theoretic_Transform.html)
