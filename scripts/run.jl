using diode_sim

params=diode_sim.Parameters(
    t_steps=1000,
    dt=1e-22,
    n=256,
    Gr=0.0,
    dx=1e-6,
    N_A=2.1e18,
    C_Rr=1.1e-3,
    N_D=2.1e18,
    v_start=-2,
)
#s1=diode_sim.run_sim(params,-2.0,0.0);
#diode_sim.plotstate(params,s1);
s2=diode_sim.run_sim(params,0.0,0.0);
diode_sim.plotstate(params,s2);
#v,currents=diode_sim.iv_curve(params);
#Plot semilog in y
