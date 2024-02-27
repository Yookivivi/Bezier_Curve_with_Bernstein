import numpy as np
from scipy.special import binom
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import customtkinter as ctk
import tkinter
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)

import time


ctk.set_appearance_mode("light")  # Modes: system (default), light, dark
ctk.set_default_color_theme("blue")  # Themes: blue (default), dark-blue, green

class BezierBuilder(object):
    """Bézier curve interactive builder.
    """
    def __init__(self, control_polygon, ax_bernstein):
        """Constructor.
        Receives the initial control polygon of the curve.
        """
        self.control_polygon = control_polygon
        self.points = list(control_polygon.get_data())
        # points initial with 2 blank point, I do not want it, so I pop twice
        self.points.pop()
        self.points.pop()
        self.xp = list(control_polygon.get_xdata())
        self.yp = list(control_polygon.get_ydata())
        self.canvas = control_polygon.figure.canvas
        self.ax_main = control_polygon.axes
        self.ax_bernstein = ax_bernstein

        # Event handler for mouse clicking
        self.connect("all")
        self.disconnect("all")
        self.press = None
        
        # Create Bézier curve
        line_bezier = Line2D([], [],
                             c=control_polygon.get_markeredgecolor())
        self.bezier_curve = self.ax_main.add_line(line_bezier)

        self.has_point=False
        self.t=0
        self.ann_list_bernstein= []
        self.ann_list_bezier= []
        self.ann_list_points= []
        self.scatter_bezier=[]
        self.scatter_bernstein=[]


    def get_has_point(self):
        return self.has_point
    
    # used for add new point
    def on_press(self, event):

        # Ignore clicks outside axes
        if event.inaxes != self.control_polygon.axes:
            return

        # Add point. I only get the int part of th point here.
        x=int(event.xdata)
        y=int(event.ydata)
        self.xp.append(x)
        self.yp.append(y)
        self.control_polygon.set_data(self.xp, self.yp)
        self.points.append([x,y])

        ann=self.ax_main.annotate(("{:.2f}".format(x),"{:.2f}".format(y)), (x, y))
        self.ann_list_points.append(ann)

        self.has_point=True

        # Rebuild Bézier curve and update canvas
        self.bezier_curve.set_data(*self._build_bezier())
        self._update_bernstein()
        self._update_bezier()
        if(self.t!=0):
            self._plot_bezier_with_t(self.t)
            self._plot_bernstein_with_t(self.t)



    # used in "moving point", checking which point the user clicks
    def on_short_press(self, event):
        if event.inaxes != self.control_polygon.axes:
            return
        point=None
        x=int(event.xdata)
        y=int(event.ydata)
        i=0
        for p in self.points:
            if len(p)>1 and abs(p[0]-x)<2 and abs(p[1]-y)<2:
                point=p
                break
            i+=1
        if(point==None):
            return
        self.press = point, i, (x, y) # i is the index of clicked point in self.points
    
    def on_motion(self, event):
        """Move the point if the mouse is over."""

        if event.inaxes != self.control_polygon.axes:
            return
        # ignore if the user don't moving the point in self.points
        if self.press==None:
            return
        
        (x0, y0), index, (xpress, ypress) = self.press
        dx = int(event.xdata) - xpress
        dy = int(event.ydata) - ypress

        self.points[index]=[x0+dx,y0+dy]
        self.xp[index]=x0+dx
        self.yp[index]=y0+dy

    def on_release(self, event):
        """Clear button press information."""
        if event.inaxes != self.control_polygon.axes:
            return
        self.press = None

        # Rebuild Bézier curve and update canvas
        self.control_polygon.set_data(self.xp, self.yp)
        self.bezier_curve.set_data(*self._build_bezier())

        
        self._update_bernstein()
        self._update_bezier()

        if self.ann_list_points!=[]:
            for i, a in enumerate(self.ann_list_points):
                a.remove()
            self.ann_list_points[:] = []    
        for x in list(self.xp):
            index=self.xp.index(x)
            ann=self.ax_main.annotate(("{:.2f}".format(x),"{:.2f}".format(self.yp[index])), (x, self.yp[index]))
            self.ann_list_points.append(ann)
        
        self._plot_bezier_with_t(self.t)
        self._plot_bernstein_with_t(self.t)

    def connect(self, event):
        if(event=="all"):
            self.cidpress = self.canvas.mpl_connect(
                'button_press_event', self.on_press)
            self.cidmotion = self.canvas.mpl_connect(
                'motion_notify_event', self.on_motion)
            self.cidrelease = self.canvas.mpl_connect(
                'button_release_event', self.on_release)
        elif(event=="press1"):
            self.cidpress = self.canvas.mpl_connect(
                'button_press_event', self.on_press)
        elif(event=="press2"):
            self.cidpress = self.canvas.mpl_connect(
                'button_press_event', self.on_short_press)
        elif(event=="move"):
            self.cidmotion = self.canvas.mpl_connect(
                'motion_notify_event', self.on_motion)
        else:
            self.cidrelease = self.canvas.mpl_connect(
                'button_release_event', self.on_release)
        
    def disconnect(self, event):
        if(event=="all"):
            self.canvas.mpl_disconnect(self.cidpress)
            self.canvas.mpl_disconnect(self.cidmotion)
            self.canvas.mpl_disconnect(self.cidrelease)
        elif(event=="press"):
            self.canvas.mpl_disconnect(self.cidpress)
        elif(event=="move"):
            self.canvas.mpl_disconnect(self.cidmotion)
        else:
            self.canvas.mpl_disconnect(self.cidrelease)

    def _build_bezier(self):
        x, y = Bezier(list(zip(self.xp, self.yp))).T
        return x, y

    def _update_bezier(self):
        self.canvas.draw()

    def _update_bernstein(self):
        N = len(self.xp) - 1
        t = np.linspace(0, 1, num=200)
        ax = self.ax_bernstein
        ax.clear()
        for kk in range(N + 1):
            ax.plot(t, Bernstein(N, kk)(t))
        ax.set_title("Bernstein basis, N = {}".format(N))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)


    def _plot_bezier_with_t(self,t):
        self.t=t
        if self.scatter_bezier!=[]:
            for i, a in enumerate(self.scatter_bezier):
                a.remove()
            self.scatter_bezier[:] = []    
        if self.ann_list_bezier!=[]:
            for i, a in enumerate(self.ann_list_bezier):
                a.remove()
            self.ann_list_bezier[:] = []
        x, y = Bezier_with_t(list(zip(self.xp, self.yp)), t).T
        self.scatter_bezier.append(self.ax_main.scatter(x,y, c="black", marker="D"))
        ann=self.ax_main.annotate(("{:.2f}".format(x),"{:.2f}".format(y)), (x, y))
        self.ann_list_bezier.append(ann)

        self.canvas.draw()

    def _plot_bernstein_with_t(self,t):
        self.t=t
        N = len(self.xp) - 1
        if self.scatter_bernstein!=[]:
            self.scatter_bernstein.remove()
            self.scatter_bernstein=[]
        if self.ann_list_bernstein!=[]:
            for i, a in enumerate(self.ann_list_bernstein):
                a.remove()
            self.ann_list_bernstein[:] = []    
        points_list=list()
        for kk in range(N + 1):
            points_list.append((t, Bernstein(N, kk)(t)))
        self.scatter_bernstein = plt.plot(*zip(*points_list), c="black", marker="D")[0]
        for i,a in enumerate(points_list):
            ann=plt.annotate("{:.2f}".format(points_list[i][1]), points_list[i])
            self.ann_list_bernstein.append(ann)

        self.canvas.draw()


def Bernstein(n, k):
    """Bernstein polynomial.
    """
    coeff = binom(n, k)

    def _bpoly(t):
        return coeff * t ** k * (1 - t) ** (n - k)

    return _bpoly


def Bezier(points, num=200):
    """Build Bézier curve from points.
    """
    N = len(points)
    t = np.linspace(0, 1, num=num)
    curve = np.zeros((num, 2))
    for ii in range(N):
        curve += np.outer(Bernstein(N - 1, ii)(t), points[ii])
    return curve

def Bezier_with_t(points,t, num=200):
    """Build Bézier curve from points.
    """
    N = len(points)
    curve = np.zeros((num, 2))
    for ii in range(N):
        curve += np.outer(Bernstein(N - 1, ii)(t), points[ii])
    return curve[0]


class App(ctk.CTk):
        
    def __init__(self):
        super().__init__()
        self.geometry(f"{1100}x{580}")
        self.title("Dynamic Bezier Curve")
        self.update()

        self.frame = ctk.CTkFrame(master=self,
                         height= 550,
                         width = 825,
                         fg_color="blue")
        self.frame.place(relx=0.33, rely=0.025)

        self.add_button = ctk.CTkButton(master=self, text="switch to add point", width=300,height=50,command=self.add_point)
        self.add_button.place(relx=0.15, rely=0.1, anchor=ctk.CENTER)

        self.move_button = ctk.CTkButton(master=self, text="switch to move point", width=300,height=50,command=self.move_point)
        self.move_button.place(relx=0.15, rely=0.2, anchor=ctk.CENTER)

        self.clear_button = ctk.CTkButton(master=self, text="clear the canvas", width=300,height=50,command=self.clear)
        self.clear_button.place(relx=0.15, rely=0.3, anchor=ctk.CENTER)


        t_value = tkinter.DoubleVar()
        t_value.set('0.000')
        self.slider = ctk.CTkSlider(master=self,
                                    variable = t_value,
                                    width=300,
                                    height=20,
                                    from_=0, to=1,
                                    command=self.slider_function)
        self.slider.set(0)
        self.slider.place(relx= 0.025,rely=0.55) 

        self.description_text = ctk.CTkLabel(master=self, text="Slide to change the t value")
        self.description_text.place(relx= 0.025,rely=0.45)

        self.val_text = ctk.CTkLabel(master=self, text="t =")
        self.val_text.place(relx= 0.025,rely=0.5)
        self.slider_val = ctk.CTkLabel(master=self, textvariable=t_value)
        self.slider_val.place(relx= 0.043,rely=0.5)
        
        self.draw()

        self.protocol("WM_DELETE_WINDOW", self.close_window)
        self.mainloop()

    def close_window(self):
        self.withdraw()
        self.quit()

    def clear(self):
        for widget in self.frame.winfo_children():
            widget.destroy()
        self.draw()

    def draw(self):
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12))

        # Empty line
        line = Line2D([], [], ls='--', c='#666666',
                      marker='x', mew=2, mec='#204a87')
        ax1.add_line(line)
        # Canvas limits
        ax1.set_xlim(0, 100)
        ax1.set_ylim(0, 100)
        ax1.set_title("Bézier curve")
        # Bernstein plot
        ax2.set_title("Bernstein basis")

        canvas = FigureCanvasTkAgg(fig, master=self.frame)  # A tk.DrawingArea.
        canvas.draw()
        canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)


        self.bezier=BezierBuilder(line, ax2)
    
    def slider_function(self,t):
        if self.bezier.get_has_point():
            self.bezier._plot_bernstein_with_t(t)
            self.bezier._plot_bezier_with_t(t)
    
    def add_point(self):
        print("add new points")
        self.bezier.disconnect("all")
        self.bezier.connect("press1")
        
    def move_point(self):
        print("move points")
        self.bezier.disconnect("all")
        self.bezier.connect("press2")
        self.bezier.connect("move")
        self.bezier.connect("release")


if __name__ == '__main__':
    app = App()