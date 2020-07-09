function tau = lagrange(L,q,dq)
    syms q1 q2 q3
    syms dq1 dq2 dq3 ddq1 ddq2 ddq3 real;
    syms th1(t) th2(t) th3(t)
    
    L_ddq = diff(L, dq);
    L_t = subs(L_ddq,[q1 q2 q3 dq1 dq2 dq3], [th1, th2, th3, ...
        diff(th1(t),t), diff(th2(t),t), diff(th3(t), t)]);
    L_dt = diff(L_t, t);
    
    L1 = subs(L_dt,[th1 th2 th3 diff(th1(t),t) diff(th2(t),t), ...
        diff(th3(t),t) diff(th1(t),t,t) diff(th2(t),t,t), ...
        diff(th3(t), t,t)],[q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);
    L2 = diff(L, q);
    
    tau = L1 - L2;
end