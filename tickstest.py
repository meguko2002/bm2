import numpy as np
import math

def tx(N):
    delta = 0
    mod = N
    for i in (5, 6, 7, 8):  # グラフを等分する個数
        x = N / i
        k = int(math.log10(x))
        if x<1: k-=1
        for j in (1,2,5):  # 2,5,10の倍数で等分する
            tmpdelta = j*10**k
            tmpmod = N - tmpdelta * i
            if abs(mod) > abs(tmpmod):  # あまりが小さいのものを選ぶ
                mod = tmpmod
                delta = tmpdelta

    return delta, N // delta

input =0.069
print(input , tx(input))
