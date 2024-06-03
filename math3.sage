#!/usr/bin/env sage

var('V A D J L v0 vp vf x x_hat x_bar C')

# No constant velocity, no constant acceleration
sol = solve([
    L == 0.5*(v0 + vp)*x + 0.5*(vp + vf)*x_bar + vp*x_hat,
    vp == V,
    V == v0 - A^2/J + A*x,
    V == vf - D^2/J + D*x_bar,
    x == (V - v0)/A + A/J,
    x_bar == (V - vf)/D + D/J,
    x_hat == (2*L - (v0 + V)*x - (V + vf)*x_hat)/2*V,
], v0, x, x_hat, x_bar, vp)

print(f"{len(sol)} solutions")

P = (v0 - sol[0][0]) * (v0 - sol[1][0])
Pc = P.right().coefficients(v0, sparse=False)

cs = ["c", "b", "a"]

for i in range(2, -1, -1):
    print(f"auto {cs[i]} = {Pc[i].simplify_full().simplify_full()};")
