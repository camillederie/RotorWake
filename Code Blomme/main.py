from create_geometry import BladeGeometry
from lifting_line_solver import LiftingLineSolver

# single rotor
prop_geo = BladeGeometry(radius=50.0, tsr=8, v_inf=10.0, n_blades=3, n_span=20, n_theta=200,
                spacing="equal",a=0.33,xshift=0,yshift=0,zshift=0,phase_diff=0,double_rotor=False)
solver = LiftingLineSolver(geo=prop_geo, r_rotor=50, tol=1e-5, n_iter=100, weight=0.5)

solver.run_solver()

# double rotor
prop_geo = BladeGeometry(radius=50.0, tsr=8, v_inf=10.0, n_blades=3, n_span=20, n_theta=200,
                spacing="equal",a=0.33,xshift=0,yshift=200,zshift=0,phase_diff=0,double_rotor=True)
solver = LiftingLineSolver(geo=prop_geo, r_rotor=50, tol=1e-5, n_iter=100, weight=0.5)

solver.run_solver()