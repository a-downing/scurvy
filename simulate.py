import matplotlib.pyplot as plt
import sys

def simulate(periods, v_0, J, plot):
    if v_0 >= 0:
        J_acc = J
        J_dec = -J
        v = v_0
    else:
        J_acc = -J
        J_dec = J
        v = -v_0

    T1, T2, T3, T4, T5, T6, T7 = periods
    t = 0
    a = 0
    d = 0
    dt = 0.0001

    xs = []
    y1s = []
    y2s = []
    y3s = []

    while t < sum(periods):
        _d = d + v*dt + a*(dt**2*0.5)

        if t > T1+T2+T3+T4+T5+T6:
            _a = a + J_acc * dt
        elif t > T1+T2+T3+T4+T5:
            None
        elif t > T1+T2+T3+T4:
            _a = a + J_dec * dt
        elif t > T1+T2+T3:
            None
        elif t > T1+T2:
            _a = a + J_dec * dt
        elif t > T1:
            None
        else:
            _a = a + J_acc * dt

        v = v + (a + _a) * dt*0.5
        d = _d
        a = _a

        t += dt

        xs.append(t)
        y1s.append(v * 60)
        y2s.append(d)
        y3s.append(a * 60)

    d = y2s[-1]
    v = y1s[-1]

    if plot:
        plt.plot(xs, y1s, label="Velocity")
        plt.plot(xs, y2s, label="Distance")
        plt.plot(xs, y3s, label="Acceleration")
        plt.legend()
        plt.show()

#simulate.py T1 T2 T3 T4 T5 T6 T7 J v_0
periods = tuple(map(lambda x: float(x), sys.argv[1:8]))
J = float(sys.argv[8])
v_0 = float(sys.argv[9])
simulate(periods, v_0, J, True)