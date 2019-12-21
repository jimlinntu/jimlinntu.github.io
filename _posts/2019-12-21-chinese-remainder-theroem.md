---
layout: post
title: Chinese Remainder Theorem
---
# Chinese Remainder Theorem (CRT)

In the last post, I discuss how to compute NTT under a finite field.

But as [从多项式乘法到快速傅里叶变换](http://blog.miskcoo.com/2015/04/polynomial-multiplication-and-fast-fourier-transform) post mentioned,
sometimes you will be asked to compute modulus under a number that might not be prime.

In this case, you must use [Chinese Remainder Theorem](https://www.wikiwand.com/en/Chinese_remainder_theorem),
that is, use multiple prime numbers $$ p_1, p_2, ..., p_k $$ such that $$ p_1 \cdot p_2 \cdot ... p_k $$ is larger than the required number.


## Theorem

$$
    x \equiv a_1 \mod p_1 \\
    x \equiv a_2 \mod p_2 \\
           ...\\
    x \equiv a_{k-1} \mod p_{k-1} \\
    x \equiv a_{k} \mod p_{k}
$$

Let

$$
    P = p_1 \cdot ... p_{k-1} \cdot p_k
$$

There exists a unique $$ x \in [0, P) $$,

$$
    x = \sum_{i=1}^{k} a_i M_i^{-1} N_i^ \mod P
$$

where  

$$
    \begin{align}
        M_i^{-1} = N_i^{-1} \mod p_i
    \end{align}
$$

and

$$
    N_i = \frac{P}{p_i}
$$

## Proof

Because 

$$
    N_j \mid p_i  \text{ if } i \neq j
$$

Therefore, the only thing we need to verify is:

$$
    a_i M_i^{-1} N_i \mod p_i = a_i \mod p_i
$$

## Practical Suggestion

In practice, you can just choose two $$ p_1, p_2 $$ to perform NTT and then you CRT to restore the original coefficient!

# References
* [b451. 圖片匹配 by Morris](http://morris821028.github.io/2015/07/21/zj-b451/#NTT-FNT)

