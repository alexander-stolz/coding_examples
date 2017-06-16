#!/usr/bin/env python

import vtk
import time
import random
from math import sin, cos

scale  = 1.e-34
G      = 6.67384e-11
AU     = 149597871000.
M      = 1.98892e30
ve     = 29780
me     = 5.972e24
dt     = 3600. / 10000.

ANZAHL     = 20
MAX_X      = 10.
MAX_V      = 1.
MAX_M      = 10.
RK         = True
LINES      = False
RESETCAM   = True
RECHNUNGEN = 10

MSONNE    = None
MPLANETEN = None
R         = None
E         = None
PHI       = None


class Stern(vtk.vtkActor):
    def __init__(self, parent, m, x=[0., 0., 0.], v=[0., 0., 0.], farbe=None):
        self.x = x
        self.v = v
        self.m = m

        self.source = vtk.vtkSphereSource()
        self.source.SetCenter(0., 0., 0.)
        self.source.SetRadius((scale * m)**.33)
        self.source.SetThetaResolution(20)
        self.source.SetPhiResolution(20)

        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInputConnection(self.source.GetOutputPort())

        self.SetMapper(self.mapper)
        self.SetPosition([xi / AU for xi in self.x])

        if farbe:
            self.GetProperty().SetColor(farbe)

        if LINES:   parent.Pfade.append(Path(parent, self))

    def setx(self, x):
        self.x = x
        self.SetPosition([xi / AU for xi in self.x])

    def killme(self):
        del self.source
        del self.mapper


class OSD(vtk.vtkTextActor):
    def __init__(self, parent, text="", farbe=(1., 1., 1.)):
        self.parent = parent
        self.SetInput(text)
        self.GetTextProperty().SetFontSize(12)
        self.GetTextProperty().SetColor(farbe)

    def setText(self, text):
        self.SetInput(text)

    def update(self):
        self.SetInput(
            "Anzahl Sterne    = {}\n"
            "Runge-Kutta      = {}\n"
            "Reset Camera  = {}\n"
            "dt                          = {} Tage"
            .format(len(self.parent.Sterne), RK, RESETCAM, dt / 3600. / 24.))


class Path(vtk.vtkActor):
    def __init__(self, parent, stern):
        self.parent  = parent
        self.pointID = 0

        self.points = vtk.vtkPoints()
        self.points.InsertPoint(0, stern.x)

        self.lines = vtk.vtkCellArray()
        self.lines.InsertNextCell(1)    # 1 Punkt
        self.lines.InsertCellPoint(0)

        self.track = vtk.vtkPolyData()
        self.track.SetPoints(self.points)
        self.track.SetLines(self.lines)

        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInput(self.track)

        self.SetMapper(self.mapper)
        self.GetProperty().SetColor(1., 1., 1.)
        self.GetProperty().SetOpacity(.1)

    def nextPoint(self, stern):
        if LINES:
            self.pointID += 1

            self.points.InsertPoint(self.pointID, stern.x)
            self.points.Modified()

            self.lines.InsertNextCell(2)
            self.lines.InsertCellPoint(self.pointID - 1)
            self.lines.InsertCellPoint(self.pointID)
            self.lines.Modified()

    def killme(self):
        del self.points
        del self.track
        del self.lines
        del self.mapper


###############
## Animation ##
###############

class Animation():
    def __init__(self):
        self.fs         = False
        self.render     = True
        self.rechnungen = 1.

        self.setSterne()

        self.ren = vtk.vtkRenderer()
        self.ren.SetBackground(.1*.75, .2*.75, .3*.75)

        for stern in self.Sterne:
            self.ren.AddActor(stern)

        if LINES:
            for pfad in self.Pfade:
                self.ren.AddActor(pfad)

        axisactor = vtk.vtkCubeAxesActor()
        axisactor.SetBounds((0, 1., 0, 1., 0, 1.))
        #axisactor.SetXTitle("AU")
        axisactor.YAxisVisibilityOn()
        axisactor.ZAxisVisibilityOn()
        axisactor.SetCamera(self.ren.GetActiveCamera())
        #self.ren.AddActor(axisactor)

        self.osd = OSD(self, text="")
        self.osd.update()
        self.ren.AddActor2D(self.osd)

        #self.ren.GetActiveCamera().SetPosition(0., 0., 10. * MAX_X)
        #self.ren.GetActiveCamera().SetFocalPoint(0., 0., 0.)
        #self.ren.GetActiveCamera().SetClippingRange(000., -2000.)

        self.ren.ResetCamera()
        self.ren.ResetCameraClippingRange()

        self.renwin = vtk.vtkRenderWindow()
        self.renwin.AddRenderer(self.ren)
        self.renwin.SetSize(800, 600)
        self.renwin.LineSmoothingOn()
        self.renwin.PolygonSmoothingOn()
        self.renwin.SetStereoTypeToAnaglyph()

        if self.fs:    self.renwin.FullScreenOn()

        inter = vtk.vtkRenderWindowInteractor()
        inter.SetRenderWindow(self.renwin)
        inter.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())

        inter.Initialize()
        self.renwin.Render()

        inter.AddObserver("TimerEvent", self.animate)
        inter.AddObserver("KeyPressEvent", self.keypress)
        self.timer = inter.CreateRepeatingTimer(20)

        self.renwin.SetWindowName("+/- = in/decrease dt,  \
            space = Runge-Kutta on/off,  \
            c = toggle cam mode,  \
            r = reset")

        inter.Start()

        inter.DestroyTimer(self.timer)
        time.sleep(.1)

        for stern in self.Sterne:
            stern.killme()
        for pfad in self.Pfade:
            pfad.killme()

        del inter
        del self.renwin
        del self.Sterne
        del self.Pfade
        del self.osd
        del axisactor
        del self.ren
        del self.timer

    def setSterne(self):
        self.Sterne = []
        self.Pfade  = []
        if ANZAHL == 0:
            #self.Sterne.append(Stern(self, m=M, farbe=(1, .5, 0.)))
            #self.Sterne.append(Stern(self, m=.171*M, x=(AU*cos(PHI), AU*sin(PHI), 0), v=(-ve*sin(PHI), ve*cos(PHI), 0), farbe=(0, .5, 1)))
            if not PHI:
                self.setSternSystem(
                    mSonne=1.98892e30,
                    mPlaneten=[4.785e26, 3.247e26],
                    r=[2.09e10, 3.366e10],
                    e=[0.6, .8],
                    phi=[180., 0.]
                )
            else:
                self.setSternSystem(
                    mSonne=MSONNE,
                    mPlaneten=MPLANETEN,
                    r=R,
                    e=E,
                    phi=PHI
                )

        random.seed()
        ran = random.uniform
        for i in xrange(ANZAHL):
            self.Sterne.append(Stern(
                self,
                m=MAX_M * M * ran(0, 10),
                x=(MAX_X * ran(-1, 1) * AU, MAX_X * ran(-1, 1) * AU, MAX_X * ran(-1, 1) * AU),
                v=(MAX_V * ran(-1, 1) * ve, MAX_V * ran(-1, 1) * ve, MAX_V * ran(-1, 1) * ve),
                farbe=(ran(0, 1), ran(0, 1), ran(0, 1)))
            )

    def setSternSystem(self, mSonne, mPlaneten, r, e, phi):
        self.Sterne.append(Stern(self, m=mSonne, farbe=(1, 1, 0.)))
        for i in xrange(len(mPlaneten)):
            _a   = (1. - e[i]) * r[i] / (1. - (e[i])**2)
            _v   = (G * (mSonne + mPlaneten[i]) * (2. / r[i] - 1. / _a))**.5
            _phi = phi[i] / 180. * 3.14159265359

            xPlanet = [r[i] * cos(_phi), r[i] * sin(_phi), 0.]
            vPlanet = [- _v * sin(_phi), _v   * cos(_phi), 0.]

            self.Sterne.append(Stern(
                self,
                m=mPlaneten[i],
                x=xPlanet,
                v=vPlanet,
                farbe=(0, .5, .5))
            )

        ### Transformation ins Schwerpunktsystem
        m_gesamt = 0.
        mx       = [0., 0., 0.]
        for stern in self.Sterne:
            m_gesamt += stern.m
            mx        = [mx[i] + stern.m * stern.x[i] for i in xrange(3)]
        x_mitte = [mx[i] / m_gesamt for i in xrange(3)]

        for stern in self.Sterne:
            stern.setx([stern.x[i] - x_mitte[i] for i in xrange(3)])

        # Anfangsgeschwindigkeit Sonne, damit SP nicht wandert.
        vm = [0., 0., 0.]
        for planet in self.Sterne[1:]:
            vm = [vm[i] + planet.v[i] * planet.m for i in xrange(3)]
        self.Sterne[0].v = [-vm[i] / mSonne for i in xrange(3)]

    def resetSterne(self):
        ran = random.uniform
        for i in xrange(ANZAHL):
            self.Sterne[i].v = MAX_V * ran(-1, 1) * ve, MAX_V * ran(-1, 1) * ve, MAX_V * ran(-1, 1) * ve
            self.Sterne[i].setx(x=(MAX_X * ran(-1, 1) * AU, MAX_X * ran(-1, 1) * AU, MAX_X * ran(-1, 1) * AU))

    def keypress(self, o, e):
        global dt, RK, RESETCAM
        if o.GetKeySym() == "space":
            RK = not RK
        if o.GetKeySym() == "plus":
            dt *= 1.1
        if o.GetKeySym() == "minus":
            dt *= .9
        if o.GetKeySym() == "r":
            self.resetSterne()
        if o.GetKeySym() == "c":
            RESETCAM = not RESETCAM
        self.osd.update()

    def animate(self, obj=None, event=None):
        if self.render:
            self.render = False

            t = gettime()
            for steps in xrange(int(self.rechnungen)):
                if RK:  x_temp = self.rk4_step()
                else:   x_temp = self.easy_step()

                for i in xrange(len(self.Sterne)):
                    self.Sterne[i].setx(x_temp[i])

            if LINES:
                for i in xrange(len(self.Sterne)):
                    self.Pfade[i].nextPoint(self.Sterne[i])
            t = gettime() - t

            self.rechnungen = self.rechnungen * 20. / (t * 1000.)
            if self.rechnungen < 1.: self.rechnungen = 1.

            if RESETCAM:    self.ren.ResetCamera()

            self.renwin.Render()
            self.render = True

#################
## Integration ##
#################

    def easy_step(self):
        x_temp = []
        for stern in self.Sterne:
            a = self.acceleration(stern, stern.x)
            stern.v = [stern.v[i] + a[i] * dt for i in xrange(3)]
            x_temp.append([stern.x[i] + stern.v[i] * dt for i in xrange(3)])
        return x_temp

    def rk4_step(self):
        x_temp = [None] * len(self.Sterne)
        s = -1
        for stern in self.Sterne:
            s += 1
            x_temp[s], v = self.rk4(stern, dt)
            stern.v = v
        return x_temp

    def rk4(self, stern, dt):
        x = stern.x
        v = stern.v

        x1 = x
        v1 = v
        a1 = self.acceleration(stern, x1)

        x2 = [x[i] + .5 * v1[i] * dt for i in xrange(3)]
        v2 = [v[i] + .5 * a1[i] * dt for i in xrange(3)]
        a2 = self.acceleration(stern, x2)

        x3 = [x[i] + .5 * v2[i] * dt for i in xrange(3)]
        v3 = [v[i] + .5 * a2[i] * dt for i in xrange(3)]
        a3 = self.acceleration(stern, x3)

        x4 = [x[i] + v3[i] * dt for i in xrange(3)]
        v4 = [v[i] + a3[i] * dt for i in xrange(3)]
        a4 = self.acceleration(stern, x4)

        x_next = [x[i] + (dt/6.) * (v1[i] + 2.*(v2[i] + v3[i]) + v4[i]) for i in xrange(3)]
        v_next = [v[i] + (dt/6.) * (a1[i] + 2.*(a2[i] + a3[i]) + a4[i]) for i in xrange(3)]

        return x_next, v_next

    def acceleration(self, stern, x):
        f0 = .000000001e-15
        e0 = 1
        a = [0., 0., 0.]
        for stern2 in self.Sterne:
            if stern2 is stern:
                continue
            r = ((x[0] - stern2.x[0]) ** 2 + (x[1] - stern2.x[1]) ** 2 + (x[2] - stern2.x[2]) ** 2) ** .5
            a = [a[i] - G * stern2.m * (x[i] - stern2.x[i]) / r ** 3 for i in xrange(3)]
        r0 = (x[0] ** 2 + x[1] ** 2 + x[2] ** 2) ** .5
        a = [a[0] - f0 * x[0] * r0 ** e0, a[1] - f0 * x[1] * r0 ** e0, a[2] - f0 * x[2] * r0 ** e0]
        return a


def timefunc():
    t = time.time()
    if time.time() - t > 0.:
        return time.time
    else:
        return time.clock

gettime = timefunc()


if __name__ == "__main__":
    _answer = raw_input("Zeichne Linien?\t[y / ENTER]\t> ")
    if _answer == "y":  LINES = True
    try:
        ANZAHL = int(input("Anzahl Planeten? [Zahl / ENTER]\t> "))
        _answer = raw_input("Eigene Planeten? [y / ENTER]\t> ")
        if _answer == "y":
            _mp  = []
            _r   = []
            _e   = []
            _phi = []
            MSONNE = float(raw_input("mSonne?\t[kg]\t> "))
            for i in xrange(ANZAHL):
                print "\nPlanet Nr. {}\n##################".format(i + 1)
                _mp.append(float(raw_input("Masse Planet?\t[kg]\t> ")))
                _r.append(float(raw_input("Abstand zur Sonne? [m]\t> ")))
                _e.append(float(raw_input("Exzentrizitaet?\t\t> ")))
                _phi.append(float(raw_input("Winkel?\t[grad]\t\t> ")))
            MPLANETEN = _mp
            R = _r
            E = _e
            PHI = _phi
            ANZAHL = 0
    except:
        ANZAHL = 0
        print "Beispiel wurde geladen."
    program = Animation()
