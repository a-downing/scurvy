#!/usr/bin/env sage

from sage.symbolic.ring import SymbolicRing

var('J L v0 vp x T1 T2 A D')

sol = solve([
    L == 0.5*(v0 + vp) * x,
    vp == v0 + 0.25*J*x^2,
], v0, x)

print((-1/72*4^(1/3)*3^(1/3)*((I*4^(1/3)*3^(5/6) - 4^(1/3)*3^(1/3))*J^2*(-(9*L - sqrt(3)*sqrt((27*J*L^2 - 32*vp^3)/J))/J)^(4/3) + 2*4^(2/3)*3^(2/3)*J*vp*(-(9*L - sqrt(3)*sqrt((27*J*L^2 - 32*vp^3)/J))/J)^(2/3) + 48*vp^2*(-I*sqrt(3) - 1))/(J*(-(9*L - sqrt(3)*sqrt((27*J*L^2 - 32*vp^3)/J))/J)^(2/3))).simplify_full())

print(f"{len(sol)} solutions")

for i in range(len(sol)):
    print(f"Solution {i}")
    for eq in sol[i]:
        print(eq.simplify_full())

P = (x - sol[0][1]) * (x - sol[1][1]) * (x - sol[2][1])
Pc = P.right().coefficients(x, sparse=False)

cs = ["d", "c", "b", "a"]

print("No constant velocity, no constant acceleration case (T1=T3, T2=T4=T5=T6=T7=0)")
print("v0:")
for i in range(3, -1, -1):
    print(f"auto {cs[i]} = {Pc[i].simplify_full() if isinstance(Pc[i].parent(), SymbolicRing) else Pc[i]};")

print()
# No constant velocity, constant acceleration
sol = solve([
    L == 0.5*(v0 + vp) * x,
    vp == v0 - A^2/J + A*x,
    T1 == A/J, T2 == x - 2*T1
], v0, x, T1, T2)

print(f"{len(sol)} solutions")

for i in range(len(sol)):
    print(f"Solution {i}")
    for eq in sol[i]:
        print(eq.simplify_full())

P = (x - sol[0][1]) * (x - sol[1][1])

Pc = P.right().coefficients(x, sparse=False)

cs = ["c", "b", "a"]

print()
print("No constant velocity, constant acceleration case (T1=T3, T2 > 0, T4=T5=T6=T7=0)")
print("v0:")
for i in range(2, -1, -1):
    print(f"auto {cs[i]} = {Pc[i].simplify()};")

#auto a = 1;
#auto b = -(A^2 + 2*J*vp)/(A*J);
#auto c = 2*L/A;

print((-4/3*sqrt(1/3)*sqrt(27*L**2 - 32*vp**3/J)/J + 4*L/J - 128/9*vp**3/(J**3*(sqrt(1/3)*sqrt(27*L**2 - 32*vp**3/J)/J - 3*L/J))).simplify_full())
print((-1/18*sqrt(3)*sqrt(27*L**2 - 32*vp**3/J)*J*L + 1/2*J*L**2 - 19/27*vp**3 + 256/9*vp**6/(J**3*(sqrt(3)*sqrt(27*L**2 - 32*vp**3/J)/J - 9*L/J)**2)).simplify_full())