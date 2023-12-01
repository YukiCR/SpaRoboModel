function next_state = runge_kutta_step(func, t, state, u, dt)
    k1 = func(state, t, u);
    k2 = func(state + 0.5 * dt * k1, t + 0.5 * dt, u);
    k3 = func(state + 0.5 * dt * k2, t + 0.5 * dt, u);
    k4 = func(state + dt * k3, t + dt, u);

    next_state = state + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end