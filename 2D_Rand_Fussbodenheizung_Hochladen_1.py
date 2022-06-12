"""Wenn Simulation fertig, dann ist das Gitter weg. Metrik: L in mm T in Sec """
#### Netzgenerierung und Mesh Laden ###########################################
#   Erzeuge in einem CAD-Programm eine 2D-Form. In SolidWorks (nur für Windows)
#    kann unter "Gußwerkzeuge" einer Skizze eine Fläche zugewiesen werden. 
#   Speichere dieses Objekt als STEP-Datei ab, damit es in GMSH(Linux-Version) 
#   eingelesen und vernetzt werden kann. Exportiere das Netz als msh-Datei.
#   Nutze nun das Zusatzpaket: "Meshio", um diese msh-Datei in Phyton einzulesen.
#   Speichere diese eingelesene Datei als xml-Datei ab.

#import meshio    # -->Infos:  https://pypi.org/project/meshio/2.3.5/
#msh = meshio.read("2D_Rand_msh.msh") 
#meshio.write("2D_Rand.xml",msh)# !!! Deshalb gibt es eine kleine Warnung in der Konsole.

#   Lade nun mit Dolfin diese xml-Datei in Phyton ein.
import dolfin as dlfn
mesh = dlfn.Mesh("2D_Rand.xml")
#   "et voila", Netz ist da. :) 
#### andere nützliche Pakete laden ############################################
from dolfin import dot, inner, outer, grad, div #....extra definieren, damit etwaige Gleichung übersichtlich bleiben.
from vedo.dolfin import * # ....damit die Simulation sofort in einem Plot sichtbar ist.
import matplotlib.pyplot as plt

#### Startbedingungen für die Temperatursimulation definieren #################
Start_T_Sockel = 0 # Feld-Temperatur zum Zeitpunkt t=0
Vorlauf_T = 45  # Vorlauftemperatur  ...ist eine Dirichlet-RB
Nachlauf_T = 40 # Rücklauftemperatur ...ist eine Dirichlet-RB
Halbrohr_T = 44 # Diese Temperatur hängt von der Verlegeform ab. Schnecke versus Alternative
K = dlfn.Constant(0.994)  # Temperaturleitfähigkeit für Beton: https://de.wikipedia.org/wiki/Temperaturleitf%C3%A4higkeit
T_end = 6600 # Gesamt-Simulations-Zeit
dt = 100 # Zeitschritt in Sekunden

#### Experiment: Neumann-RBs als Relation modellieren in Anlehnung an DIN 1264 Teil 1 bis 5 
R_1 = 0.75 # minimal zulässiger  Wärmeleitwiderstand durch untere Wände zur Umgebung (anderer Raum)
R_4 = 0.75*2 # Annahme: minimal zulässiger  Wärmeleitwiderstand durch Wände zur Umgebung (Draußen)
R_3 = 0.15 # maximal zulässiger  Wärmeleitwiderstand der tragenden Fussbodenschicht (Drinne)
Skalierer_Neumann = 0.050 # Sakliere diesen Wert solange bis die Fussboden Temperatur die maximal zulässigen 29°C erreicht
#   Equation coefficients für die schwache Form 
g_1 = dlfn.Constant(-1/R_1*Skalierer_Neumann)  # Neumann heat flux unten
g_4 = dlfn.Constant(-1/R_4*Skalierer_Neumann)  # Neumann heat flux links zur Wand hin
g_3 = dlfn.Constant(-1/R_3*Skalierer_Neumann)  # Neumann heat flux oben aus dem Fussboden raus
# ---> Simuliere solange, bis sich ein Gleichgewicht einstellt.<---Schaue bitte auf die Colorbar



#### Gebiete für Dirichlet- und Neumann-RBs ansteuern. #######################
#### Die folgenden Zeilen sind in Anlehnung an Sebastians Rund-Mail entstanden. 
#    ....Funktionieren aber nicht so wirklich einfach, wie gedacht.
space_dim = mesh.geometry().dim()
mvc = dlfn.MeshValueCollection("size_t", mesh, space_dim - 1)
facet_markers = dlfn.cpp.mesh.MeshFunctionSizet(mesh, mvc)
#facet_markers = dlfn.cpp.mesh.DomainBoundary(mesh)
#facet_markers = dlfn.cpp.mesh.SubDomain(map_tol=1e-1)
dsN_1 = dlfn.Measure("ds", subdomain_id=1, subdomain_data=facet_markers)
dsN_2 = dlfn.Measure("ds", subdomain_id=2, subdomain_data=facet_markers)
dsN_3 = dlfn.Measure("ds", subdomain_id=3, subdomain_data=facet_markers)
dsN_4 = dlfn.Measure("ds", subdomain_id=4, subdomain_data=facet_markers)
dsN_5 = dlfn.Measure("ds", subdomain_id=5, subdomain_data=facet_markers)
dsN_6 = dlfn.Measure("ds", subdomain_id=6, subdomain_data=facet_markers)
dsN_7 = dlfn.Measure("ds", subdomain_id=7, subdomain_data=facet_markers)

####  Ok, jetzt anders nach Recherche im Internet: ...stundenlang suchen und ausprobieren 
tdim = mesh.topology().dim()
boundary_parts = dlfn.MeshFunction("size_t", mesh, tdim - 1)
#boundary_parts.set_all(2) # ...lustiger Befehl, um alle Ränder auf einmal ansteuern.
#   Um die genauen geometrischen Daten in den folgenden Zeilen einzugeben, lohnt sich ein Blig auf die technische Zeichnung des CAD-Models 
unten = dlfn.CompiledSubDomain("near(x[1], 0.0) && on_boundary")
rechts = dlfn.CompiledSubDomain("near(x[0], 132.82) && on_boundary")
oben = dlfn.CompiledSubDomain("near(x[1], 50.0) && on_boundary")
links = dlfn.CompiledSubDomain("near(x[0], 0.0) && on_boundary")
Kreis_Vorlauf = dlfn.CompiledSubDomain("near(pow(pow(x[0]-37.50,2)+pow(x[1]-12.50,2),0.5)> 7.6, 0) && on_boundary")
Kreis_Nachlauf = dlfn.CompiledSubDomain("near(pow(pow(x[0]-82.82,2)+pow(x[1]-12.50,2),0.5)> 7.6, 0) && on_boundary")
Halbkreis = dlfn.CompiledSubDomain("near(pow(pow(x[0]-132.82,2)+pow(x[1]-12.50,2),0.5)> 7.6, 0) && on_boundary")

unten.mark(boundary_parts, 1) # ...für Neumann-RB
rechts.mark(boundary_parts, 2)# ...für Neumann-RB ..geht noch nicht anwendbar anzusteuern, weil es wahrscheinlich aus  zwei Teilen besteht
oben.mark(boundary_parts, 3)  # ...für Neumann-RB
links.mark(boundary_parts, 4) # ...für Neumann-RB
# Define the circle boundary
Kreis_Vorlauf.mark(boundary_parts, 5)  # ...für Dirichlet-RB
Kreis_Nachlauf.mark(boundary_parts, 6) # ...für Dirichlet-RB   
Halbkreis.mark(boundary_parts, 7)      # ...für Dirichlet-RB
 
#### Zuweisung der RBs  und Operator-Definition ###############################
V = dlfn.FunctionSpace(mesh, "CG", 1)      
bc = (dlfn.DirichletBC(V,dlfn.Expression('Vorlauf_T',Vorlauf_T=Vorlauf_T,degree = 1), boundary_parts,5),
      dlfn.DirichletBC(V,dlfn.Expression('Nachlauf_T',Nachlauf_T=Nachlauf_T,degree = 1), boundary_parts,6),
      dlfn.DirichletBC(V,dlfn.Expression('Halbrohr_T',Halbrohr_T=Halbrohr_T,degree = 1), boundary_parts,7))
#   ...und noch die Neumann-RBs
dsN_1 = dlfn.Measure("ds", subdomain_id=1, subdomain_data=boundary_parts)
dsN_3 = dlfn.Measure("ds", subdomain_id=3, subdomain_data=boundary_parts)
dsN_4 = dlfn.Measure("ds", subdomain_id=4, subdomain_data=boundary_parts)
#   Define steady part of the equation
def operator(u, v):
    return ((K * inner(grad(u), grad(v))) * dlfn.dx 
           - K * g_1 * v * dsN_1
           - K * g_3 * v * dsN_3
           - K * g_4 * v * dsN_4)

#### Define trial and test function and solution at previous time-step  #######
u = dlfn.TrialFunction(V)
v = dlfn.TestFunction(V)
u0 = dlfn.Function(V)

#   Time-stepping parameters
theta = dlfn.Constant(1)  # Crank-Nicolson scheme (0 für exliziten, 1 für impliziten Euler)

#   Define time discretized equation
F = ((1.0 / dt) * inner(u - u0, v) * dlfn.dx
    + theta * operator(u, v)
    + (1.0 - theta) * operator(u0, v) )

#   Prepare solution function and solver
u = dlfn.Function(V) # ...jetzt wird u überschrieben
problem = dlfn.LinearVariationalProblem(dlfn.lhs(F), dlfn.rhs(F), u, bc)
solver = dlfn.LinearVariationalSolver(problem)

#   Initial condition 
ic = dlfn.Expression("""1*pow(x[0] - 0.25, 0) + 1*pow(x[1] - 0.25, 0) < 0.2*0.2 
                   ? -25.0 * ((pow(x[0] - 0.25, 2) + pow(x[1] - 0.25, 2)) - 0.2*0.2)
                   : Start_T_Sockel""",Start_T_Sockel=Start_T_Sockel, degree=1)

#   Prepare initial condition
u0.interpolate(ic)
u.interpolate(ic)


t = 0.0
while t < T_end:
        
    solver.solve()
     
    plot(u,
         wireframe = True,
         interactive=False,
         text=__doc__+"\nTemperatur nach Minute = %g" % float(t/60))
#    plot(u,
#         text=__doc__+"\nTemperature at t = %g" % t,
#         style=2,
#         axes=2,
#         lw=0, # no mesh edge lines
#         warpZfactor=0.05,
#         isolines={"n": 12, "lw":1, "c":'black', "alpha":0.1},
#         wireframe = False,
#         interactive=False)
        # Move to next time step
    u0.assign(u)
    t += dt


plot(u,ytitle="mm",xtitle="mm",
     text=__doc__+"\nTemperatur nach Minute = %g" % float(t/60),
     style=1,    
     lw=0,# no mesh edge lines
     title="Temperaturfeld") 

interactive()


