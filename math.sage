#!/usr/bin/env sage

var('J L v0 vp x T1 T2 A D')

R = RealField(1000)

#J = R(83.263479748143922)/R(60);
#L = R(-68.675129777373016);
#v0 = R(-0.0034879858661581052)/R(60);

# No constant velocity, no constant acceleration
sol = solve([
    L == R(0.5)*(v0 + vp) * x,
    vp == v0 + R(0.25)*J*x^2
], vp, x)

#print(sol)

P = (vp - sol[0][0]) * (vp - sol[1][0]) * (vp - sol[2][0])
Q = (x - sol[0][1]) * (x - sol[1][1]) * (x - sol[2][1])

Pc = P.right().coefficients(vp, sparse=False)

cs = ["d", "c", "b", "a"]

print("No constant velocity, no constant acceleration case (T1=T3, T2=T4=T5=T6=T7=0)")
print("vp:")
for i in range(3, -1, -1):
    print(f"auto {cs[i]} = {Pc[i].simplify()};")

Qc = Q.right().coefficients(x, sparse=False)

print("x:")
for i in range(3, -1, -1):
    print(f"auto {cs[i]} = {Qc[i]};")




# No constant velocity, constant acceleration
sol = solve([L == 0.5*(v0 + vp) * x, vp == v0 - A^2/J + A*x, T1 == A/J, T2 == x - 2*T1], vp, x, T1, T2)

#print(sol)

P = (vp - sol[0][0]) * (vp - sol[1][0])
Q = (x - sol[0][1]) * (x - sol[1][1])

Pc = P.right().coefficients(vp, sparse=False)

cs = ["c", "b", "a"]

print()
print("No constant velocity, constant acceleration case (T1=T3, T2 > 0, T4=T5=T6=T7=0)")
print("vp:")
for i in range(2, -1, -1):
    print(f"auto {cs[i]} = {Pc[i].simplify()};")

print()

Qc = Q.right().coefficients(x, sparse=False)

print("x:")
for i in range(2, -1, -1):
    print(f"auto {cs[i]} = {Qc[i]};")