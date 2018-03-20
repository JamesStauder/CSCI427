from support.expressions import *
from support.momentum import *
from support.fenics_optimizations import *
import fenics as fc
import dolfin as df
from constants import *

'''
Class: IceCube
Argument list: fileName with all of the data
Purpose: This is the black box that runs the model. It runs the model and stores the output in various arrays.
This will run the model all at once instead of updating a GUI every tick.

LEGEND:
BB, HH, TD, TB, TX, TY, TZ, us, ub






Return types, values: None
Dependencies: dolfin, fenics, support functions(written by Dr. Jesse Johnson), h5 file with the initial datasets
Creator: James Stauder
Date created:2/23/18
Last edited: 2/23/18
'''


# Unit Test: Calculating outside of range
class IceCube:
    def __init__(self, fileName, timeEnd, timeStep):

        self.times = []
        self.BB = []
        self.HH = []
        self.TD = []
        self.TB = []
        self.TX = []
        self.TY = []
        self.TZ = []
        self.us = []
        self.ub = []

        ##########################################################
        ################           MESH          #################
        ##########################################################
        # TODO: Probably do not have to save then open mesh
        self.mesh = df.Mesh()
        self.inFile = fc.HDF5File(self.mesh.mpi_comm(), fileName, "r")
        self.inFile.read(self.mesh, "/mesh", False)

        #########################################################
        #################  FUNCTION SPACES  #####################
        #########################################################
        self.E_Q = df.FiniteElement("CG", self.mesh.ufl_cell(), 1)
        self.Q = df.FunctionSpace(self.mesh, self.E_Q)
        self.E_V = df.MixedElement(self.E_Q, self.E_Q, self.E_Q)
        self.V = df.FunctionSpace(self.mesh, self.E_V)

        self.assigner_inv = fc.FunctionAssigner([self.Q, self.Q, self.Q], self.V)
        self.assigner = fc.FunctionAssigner(self.V, [self.Q, self.Q, self.Q])

        self.U = df.Function(self.V)
        self.dU = df.TrialFunction(self.V)
        self.Phi = df.TestFunction(self.V)
        self.u, self.u2, self.H = df.split(self.U)
        self.phi, self.phi1, self.xsi = df.split(self.Phi)

        self.un = df.Function(self.Q)
        self.u2n = df.Function(self.Q)

        self.zero_sol = df.Function(self.Q)

        self.S0 = df.Function(self.Q)
        self.B = df.Function(self.Q)
        self.H0 = df.Function(self.Q)
        self.A = df.Function(self.Q)

        self.inFile.read(self.S0.vector(), "/surface", True)
        self.inFile.read(self.B.vector(), "/bed", True)
        self.inFile.read(self.A.vector(), "/smb", True)

        self.H0.assign(self.S0 - self.B)

        self.Hmid = theta * self.H + (1 - theta) * self.H0

        self.S = self.B + self.Hmid

        self.width = df.interpolate(Width(degree=2), self.Q)

        self.strs = Stresses(self.U, self.Hmid, self.H0, self.H, self.width, self.B, self.S, self.Phi)

        self.R = -(self.strs.tau_xx + self.strs.tau_xz + self.strs.tau_b + self.strs.tau_d + self.strs.tau_xy) * df.dx

        #############################################################################
        ########################  MASS CONSERVATION  ################################
        #############################################################################
        self.h = df.CellSize(self.mesh)
        self.D = self.h * abs(self.U[0]) / 2.
        self.area = self.Hmid * self.width

        self.mesh_min = self.mesh.coordinates().min()
        self.mesh_max = self.mesh.coordinates().max()

        # Define boundaries
        self.ocean = df.FacetFunctionSizet(self.mesh, 0)
        self.ds = fc.ds(subdomain_data=self.ocean)  # THIS DS IS FROM FENICS! border integral

        for f in df.facets(self.mesh):
            if df.near(f.midpoint().x(), self.mesh_max):
                self.ocean[f] = 1
            if df.near(f.midpoint().x(), self.mesh_min):
                self.ocean[f] = 2

        self.R += ((self.H - self.H0) / dt * self.xsi
                   - self.xsi.dx(0) * self.U[0] * self.Hmid
                   + self.D * self.xsi.dx(0) * self.Hmid.dx(0)
                   - (self.A - self.U[0] * self.H / self.width * self.width.dx(0))
                   * self.xsi) * df.dx + self.U[0] * self.area * self.xsi * self.ds(1) \
                  - self.U[0] * self.area * self.xsi * self.ds(0)

        #####################################################################
        #########################  SOLVER SETUP   ###########################
        #####################################################################

        # Bounds
        self.l_thick_bound = df.project(Constant(thklim), self.Q)
        self.u_thick_bound = df.project(Constant(1e4), self.Q)

        self.l_v_bound = df.project(-10000.0, self.Q)
        self.u_v_bound = df.project(10000.0, self.Q)

        self.l_bound = df.Function(self.V)
        self.u_bound = df.Function(self.V)

        self.assigner.assign(self.l_bound, [self.l_v_bound] * 2 + [self.l_thick_bound])
        self.assigner.assign(self.u_bound, [self.u_v_bound] * 2 + [self.u_thick_bound])

        # This should set the velocity at the divide (left) to zero
        self.dbc0 = df.DirichletBC(self.V.sub(0), 0, lambda x, o: df.near(x[0], self.mesh_min) and o)
        # Set the velocity on the right terminus to zero
        self.dbc1 = df.DirichletBC(self.V.sub(0), 0, lambda x, o: df.near(x[0], self.mesh_max) and o)
        # overkill?
        self.dbc2 = df.DirichletBC(self.V.sub(1), 0, lambda x, o: df.near(x[0], self.mesh_max) and o)
        # set the thickness on the right edge to thklim
        self.dbc3 = df.DirichletBC(self.V.sub(2), thklim, lambda x, o: df.near(x[0], self.mesh_max) and o)

        # Define variational solver for the mass-momentum coupled problem
        self.J = df.derivative(self.R, self.U, self.dU)

        self.coupled_problem = df.NonlinearVariationalProblem(self.R, self.U, bcs=[self.dbc0, self.dbc1, self.dbc3], \
                                                              J=self.J)

        self.coupled_problem.set_bounds(self.l_bound, self.u_bound)

        self.coupled_solver = df.NonlinearVariationalSolver(self.coupled_problem)

        # Acquire the optimizations in fenics optimizations
        set_solver_options(self.coupled_solver)

        self.t = 0
        self.timeEnd = float(timeEnd)
        self.dtFloat = float(timeStep)

        self.inFile.close()

    def runAllSteps(self):
        while self.t < self.timeEnd:
            self.runNextStep()

    def runNextStep(self):
        self.coupled_problem = df.NonlinearVariationalProblem(self.R, self.U, bcs=[self.dbc0, self.dbc1, self.dbc3],
                                                              J=self.J)
        self.coupled_problem.set_bounds(self.l_bound, self.u_bound)
        self.coupled_solver = df.NonlinearVariationalSolver(self.coupled_problem)

        # Optimizations in fenics optimizations
        set_solver_options(self.coupled_solver)

        try:
            self.coupled_solver.solve(set_solver_options())
        except:
            self.coupled_solver.parameters['snes_solver']['error_on_nonconvergence'] = False
            self.assigner.assign(self.U, [self.zero_sol, self.zero_sol, self.H0])
            self.coupled_solver.solve()
            self.coupled_solver.parameters['snes_solver']['error_on_nonconvergence'] = True

        self.assigner_inv.assign([self.un, self.u2n, self.H0], self.U)
        self.times.append(self.t)

        self.BB.append(self.strs.B.compute_vertex_values())
        self.HH.append(self.strs.H0.compute_vertex_values())
        self.TD.append(df.project(self.strs.tau_d_plot).compute_vertex_values())
        self.TB.append(df.project(self.strs.tau_b_plot).compute_vertex_values())
        self.TX.append(df.project(self.strs.tau_xx_plot).compute_vertex_values())
        self.TY.append(df.project(self.strs.tau_xy_plot).compute_vertex_values())
        self.TZ.append(df.project(self.strs.tau_xz_plot).compute_vertex_values())
        self.us.append(df.project(self.strs.u(0)).compute_vertex_values())
        self.ub.append(df.project(self.strs.u(1)).compute_vertex_values())

        self.t += self.dtFloat

        return self.BB[-1], self.HH[-1], self.TD[-1], self.TB[-1], self.TX[-1], self.TY[-1], self.TZ[-1], self.us[-1], \
               self.ub[-1]
