#!/usr/bin/env sage

var('J L v0 vp x T1 T2 A D')

# No constant velocity, no constant acceleration
sol = solve([
    L == 0.5*(v0 + vp) * x,
    vp == v0 + 0.25*J*x^2,
], v0, x)

P = (v0 - sol[0][0]) * (v0 - sol[1][0]) * (v0 - sol[2][0])
Pc = P.right().coefficients(v0, sparse=False)

cs = ["d", "c", "b", "a"]

print("No constant velocity, no constant acceleration case (T1=T3, T2=T4=T5=T6=T7=0)")
print("v0:")
for i in range(3, -1, -1):
    print(f"auto {cs[i]} = {Pc[i].simplify()};")




# No constant velocity, constant acceleration
sol = solve([
    L == 0.5*(v0 + vp) * x,
    vp == v0 - A^2/J + A*x,
    T1 == A/J, T2 == x - 2*T1
], v0, x, T1, T2)

P = (v0 - sol[0][0]) * (v0 - sol[1][0])

Pc = P.right().coefficients(v0, sparse=False)

cs = ["c", "b", "a"]

print()
print("No constant velocity, constant acceleration case (T1=T3, T2 > 0, T4=T5=T6=T7=0)")
print("v0:")
for i in range(2, -1, -1):
    print(f"auto {cs[i]} = {Pc[i].simplify()};")