import numpy as np
import test4
print("testing without memory mod:")
a = np.array([1,2,3.0])
b = test4.square_vec(a)
print("a:", a)
print("b:", b)

print("testing with mem mod:")

a = np.array([1,2,3.0])
b = test4.square_vec_noreturn(a)
print("a:", a)
print("b:", b)

print("testing with mem mod:")

