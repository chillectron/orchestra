from tkinter import *
from math import sqrt
import math
import time
import numpy as np
from scipy.special import comb
import datetime
from math import sin, cos, tan

master = Tk()
cheight = 855
cwidth = 1440
framerate = 75

c = Canvas(master, bg='#202020', bd=0, height=cheight, width=cwidth, highlightthickness=0, relief='ridge')
c.pack()

elements = ['res','nmos','wire','pmos','idc','vdc','vcvs','vccs','plot','plottrace','text']
dict_of_elements = dict([(y,x+1) for x,y in enumerate(elements)])

def eqn_interpolator(eqn_in_x, steps, step_size):

	x = 0
	steps = round(steps/step_size)
	xs = [0]*steps
	ys = [0]*steps
	for i in range(0, steps):
		ys[i] = eval(eqn_in_x)
		x+=step_size
		xs[i] = x

	return xs, [i * steps/(max(ys)-min(ys)) for i in ys]

def bernstein_poly(i, n, t):

    return comb(n, i) * ( t**(n-i) ) * (1 - t)**i


def bezier_curve(points, nTimes=1000):

	nPoints = len(points)
	xPoints = np.array([p[0] for p in points])
	yPoints = np.array([p[1] for p in points])

	t = np.linspace(0.0, 1.0, nTimes)

	polynomial_array = np.array([ bernstein_poly(i, nPoints-1, t) for i in range(0, nPoints)])

	xvals = np.dot(xPoints, polynomial_array)
	yvals = np.dot(yPoints, polynomial_array)

	return xvals, yvals

def bezier(t):
    return t**2 * (3.0 - 2.0 * t)

def ease(t):
    sqt = t**2
    return sqt / (2 * (sqt - t) + 1)

def math_rotate(oldx, oldy, centerx, centery, angle):

	out = complex(centerx, -1*centery) + complex(math.cos(math.radians(angle)), math.sin(math.radians(angle))) * (complex(oldx, -1*oldy) - complex(centerx, -1*centery))
	return [out.real, -1*out.imag]

def math_translate(oldx, oldy, deltax, deltay):
	out = complex(oldx, -1*oldy) + complex(deltax, -1*deltay)
	return [out.real, -1*out.imag]

def math_scale(oldx, oldy, centerx, centery, scale_x, scale_y):
	out = complex(oldx, -1*oldy) - complex(centerx,  -1*centery)
	out = complex(out.real*scale_x, out.imag*scale_y) + complex(oldx, -1*oldy)
	return [out.real, -1*out.imag]

def math_appear(oldx, oldy, centerx, centery, scale_x, scale_y):
	out = complex(oldx, -1*oldy) - complex(centerx,  -1*centery)
	out = out * complex(scale_x, 0) + complex(centerx, -1*centery)
	return [out.real, -1*out.imag]

def reCenter(obj):
	obj.rotate_around_x_old = getCenter(obj)[0]
	obj.rotate_around_y_old = getCenter(obj)[1]

	for i in range(0, obj.items):
		obj.drawing_coords[i] = (obj.drawing_coords[i][0] + obj.posx - obj.rotate_around_x_old, obj.drawing_coords[i][1] + obj.posy - obj.rotate_around_y_old, obj.drawing_coords[i][2] + obj.posx - obj.rotate_around_x_old, obj.drawing_coords[i][3] + obj.posy -obj.rotate_around_y_old)

	obj.rotate_around_x = getCenter(obj)[0]
	obj.rotate_around_y = getCenter(obj)[1]

def getCenter(obj):
	if obj.name == 'text':
		return [obj.drawing_coords[0][0], obj.drawing_coords[0][1]]
	if obj.name == 'plottrace':
		return [obj.drawing_coords[0][0], obj.drawing_coords[0][1]]

	for i in range(0, obj.items):
		obj.xmin_tmp[i] = min(obj.drawing_coords[i][0], obj.drawing_coords[i][2])
		obj.xmax_tmp[i] = max(obj.drawing_coords[i][0], obj.drawing_coords[i][2])
		obj.ymin_tmp[i] = min(obj.drawing_coords[i][1], obj.drawing_coords[i][3])
		obj.ymax_tmp[i] = max(obj.drawing_coords[i][1], obj.drawing_coords[i][3])
	return [(min(obj.xmin_tmp) + max(obj.xmax_tmp))/2, (min(obj.ymin_tmp) + max(obj.ymax_tmp))/2]

def draw(obj, position):

	if position==1:
		obj.drawing[0] = c.create_line(obj.drawing_coords[0], arrow='none', fill=obj.color, width=obj.width, capstyle=ROUND)
		obj.drawing[1] = c.create_line(obj.drawing_coords[1], arrow='none', fill=obj.color, width=obj.width, capstyle=ROUND)
		obj.drawing[2] = c.create_line(obj.drawing_coords[2], arrow='none', fill=obj.color, width=obj.width, capstyle=ROUND)
		obj.drawing[3] = c.create_line(obj.drawing_coords[3], arrow='none', fill=obj.color, width=obj.width, capstyle=ROUND)
		obj.drawing[4] = c.create_line(obj.drawing_coords[4], arrow='none', fill=obj.color, width=obj.width, capstyle=ROUND)
		obj.drawing[5] = c.create_line(obj.drawing_coords[5], arrow='none', fill=obj.color, width=obj.width, capstyle=ROUND)
		obj.drawing[6] = c.create_line(obj.drawing_coords[6], arrow='none', fill=obj.color, width=obj.width, capstyle=ROUND)
		obj.drawing[7] = c.create_line(obj.drawing_coords[7], arrow='none', fill=obj.color, width=obj.width, capstyle=ROUND)
		obj.drawing[8] = c.create_line(obj.drawing_coords[8], arrow='none', fill=obj.color, width=obj.width, capstyle=ROUND)
		obj.drawing[9] = c.create_line(obj.drawing_coords[9], arrow='none', fill=obj.color, width=obj.width, capstyle=ROUND)
		obj.drawing[10] = c.create_line(obj.drawing_coords[10], arrow='none', fill=obj.color, width=obj.width, capstyle=ROUND)
		obj.drawing[11] = c.create_line(obj.drawing_coords[11], arrow='none', fill=obj.color, width=obj.width, capstyle=ROUND)

	if position==2:
		obj.drawing[0]= c.create_line(obj.drawing_coords[0], arrow='none', fill=obj.color, width=obj.width, capstyle=PROJECTING)
		obj.drawing[1]= c.create_line(obj.drawing_coords[1], arrow='none', fill=obj.color, width=obj.width, capstyle=PROJECTING)
		obj.drawing[2]= c.create_line(obj.drawing_coords[2], arrow='none', fill=obj.color, width=obj.width, capstyle=PROJECTING)
		obj.drawing[3]= c.create_line(obj.drawing_coords[3], arrow='none', fill=obj.color, width=obj.width, capstyle=PROJECTING)
		obj.drawing[4]= c.create_line(obj.drawing_coords[4], arrow='none', fill=obj.color, width=obj.width, capstyle=PROJECTING)
		obj.drawing[5]= c.create_line(obj.drawing_coords[5], arrow='none', fill=obj.color, width=obj.width, capstyle=PROJECTING)
		obj.drawing[6]= c.create_line(obj.drawing_coords[6], arrow='none', fill=obj.color, width=obj.width, capstyle=PROJECTING)
		obj.drawing[7] = c.create_line(obj.drawing_coords[7], arrow='none', fill=obj.color, width=obj.width, capstyle=PROJECTING)
		obj.drawing[8] = c.create_line(obj.drawing_coords[8], arrow='none', fill=obj.color, width=obj.width, capstyle=PROJECTING)

	if position==3:
		obj.drawing[0] = c.create_line(obj.drawing_coords[0], arrow='none', fill=obj.color, width=obj.width, joinstyle='bevel')

	if position==4:
		obj.drawing[0]= c.create_line(obj.drawing_coords[0], arrow='none', fill=obj.color, width=obj.width, capstyle=PROJECTING)
		obj.drawing[1]= c.create_line(obj.drawing_coords[1], arrow='none', fill=obj.color, width=obj.width, capstyle=PROJECTING)
		obj.drawing[2]= c.create_line(obj.drawing_coords[2], arrow='none', fill=obj.color, width=obj.width, capstyle=PROJECTING)
		obj.drawing[3]= c.create_line(obj.drawing_coords[3], arrow='none', fill=obj.color, width=obj.width, capstyle=PROJECTING)
		obj.drawing[4]= c.create_line(obj.drawing_coords[4], arrow='none', fill=obj.color, width=obj.width, capstyle=PROJECTING)
		obj.drawing[5]= c.create_line(obj.drawing_coords[5], arrow='none', fill=obj.color, width=obj.width, capstyle=PROJECTING)
		obj.drawing[6]= c.create_line(obj.drawing_coords[6], arrow='none', fill=obj.color, width=obj.width, capstyle=PROJECTING)
		obj.drawing[7] = c.create_line(obj.drawing_coords[7], arrow='none', fill=obj.color, width=obj.width, capstyle=PROJECTING)
		obj.drawing[8] = c.create_line(obj.drawing_coords[8], arrow='none', fill=obj.color, width=obj.width, capstyle=PROJECTING)

	if position==5:
		obj.drawing[0] = c.create_line(obj.drawing_coords[0], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[1] = c.create_line(obj.drawing_coords[1], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[2] = c.create_oval(obj.drawing_coords[2], outline='#FFFFFF',fill=obj.color, width=obj.width)
		obj.drawing[3] = c.create_oval(obj.drawing_coords[3], outline='#FFFFFF',fill='#000000', width=obj.width)
		obj.drawing[4] = c.create_line(obj.drawing_coords[4], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[5] = c.create_line(obj.drawing_coords[5], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[6] = c.create_line(obj.drawing_coords[6], arrow='none', fill=obj.color, width=obj.width)

	if position==6:
		obj.drawing[0] = c.create_line(obj.drawing_coords[0], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[1] = c.create_line(obj.drawing_coords[1], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[2] = c.create_oval(obj.drawing_coords[2], outline='#FFFFFF',fill=obj.color, width=obj.width)
		obj.drawing[3] = c.create_oval(obj.drawing_coords[3], outline='#FFFFFF',fill='#000000', width=obj.width)
		obj.drawing[4] = c.create_line(obj.drawing_coords[4], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[5] = c.create_line(obj.drawing_coords[5], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[6] = c.create_line(obj.drawing_coords[6], arrow='none', fill=obj.color, width=obj.width)
	
	if position==7:
		obj.drawing[0] = c.create_line(obj.drawing_coords[0], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[1] = c.create_line(obj.drawing_coords[1], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[2] = c.create_line(obj.drawing_coords[2], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[3] = c.create_line(obj.drawing_coords[3], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[4] = c.create_line(obj.drawing_coords[4], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[5] = c.create_line(obj.drawing_coords[5], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[6] = c.create_line(obj.drawing_coords[6], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[7] = c.create_line(obj.drawing_coords[7], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[8] = c.create_line(obj.drawing_coords[8], arrow='none', fill=obj.color, width=obj.width)

	if position==8:
		obj.drawing[0] = c.create_line(obj.drawing_coords[0], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[1] = c.create_line(obj.drawing_coords[1], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[2] = c.create_line(obj.drawing_coords[2], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[3] = c.create_line(obj.drawing_coords[3], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[4] = c.create_line(obj.drawing_coords[4], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[5] = c.create_line(obj.drawing_coords[5], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[6] = c.create_line(obj.drawing_coords[6], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[7] = c.create_line(obj.drawing_coords[7], arrow='none', fill=obj.color, width=obj.width)
		obj.drawing[8] = c.create_line(obj.drawing_coords[8], arrow='none', fill=obj.color, width=obj.width)

	if position==9:
		for i in range(0,obj.items):
			obj.drawing[i] = c.create_line(obj.drawing_coords[i], arrow='none', fill=obj.color, width=obj.width)

	if position==10:
		for i in range(0,len(obj.xvals)-1):
			obj.drawing[i] = c.create_line(obj.drawing_coords[i], arrow='none', fill=obj.color, width=obj.width)

	if position==11:
		obj.drawing[0] = c.create_text(obj.posx, obj.posy, text=obj.text, fill=obj.color, font=("Arial", obj.fontsize), angle=obj.angle)


def drawings(obj, position, draw):

	if position==1:
		if draw:
			obj.drawing_coords[0] = (obj.posx, obj.posy, obj.posx+15, obj.posy)
			obj.drawing_coords[1] = (obj.posx+15, obj.posy, obj.posx+15+4, obj.posy+9)
			obj.drawing_coords[2] = (obj.posx+15+4, obj.posy+9, obj.posx+15+4+4, obj.posy-9)
			obj.drawing_coords[3] = (obj.posx+15+4+4, obj.posy-9, obj.posx+15+4+4+4, obj.posy+9)
			obj.drawing_coords[4] = (obj.posx+15+4+4+4, obj.posy+9, obj.posx+15+4+4+4+4, obj.posy-9)
			obj.drawing_coords[5] = (obj.posx+15+4+4+4+4, obj.posy-9, obj.posx+15+4+4+4+4+4, obj.posy+9)
			obj.drawing_coords[6] = (obj.posx+15+4+4+4+4+4, obj.posy+9, obj.posx+15+4+4+4+4+4+4, obj.posy-9)
			obj.drawing_coords[7] = (obj.posx+15+4+4+4+4+4+4, obj.posy-9, obj.posx+15+4+4+4+4+4+4+4, obj.posy+9)
			obj.drawing_coords[8] = (obj.posx+15+4+4+4+4+4+4+4, obj.posy+9, obj.posx+15+4+4+4+4+4+4+4+4, obj.posy-9)
			obj.drawing_coords[9] = (obj.posx+15+4+4+4+4+4+4+4+4, obj.posy-9, obj.posx+15+4+4+4+4+4+4+4+4+4, obj.posy+9)
			obj.drawing_coords[10] = (obj.posx+15+4+4+4+4+4+4+4+4+4, obj.posy+9, obj.posx+15+4+4+4+4+4+4+4+4+4+4, obj.posy)
			obj.drawing_coords[11] = (obj.posx+15+4+4+4+4+4+4+4+4+4+4, obj.posy, obj.posx+15+4+4+4+4+4+4+4+4+4+4+15, obj.posy)
		return 12

	if position==2:
		if draw:
			obj.drawing_coords[0] = (obj.posx, obj.posy, obj.posx, obj.posy+15)
			obj.drawing_coords[1] = (obj.posx, obj.posy+15, obj.posx-26*obj.flip, obj.posy+15)
			obj.drawing_coords[2] = (obj.posx-26*obj.flip, obj.posy+15, obj.posx-26*obj.flip, obj.posy+15+26)
			obj.drawing_coords[3] = (obj.posx-26*obj.flip, obj.posy+15+26, obj.posx, obj.posy+15+26)
			obj.drawing_coords[4] = (obj.posx, obj.posy+15+26, obj.posx, obj.posy+15+26+15)
			obj.drawing_coords[5] = (obj.posx-30*obj.flip, obj.posy+15, obj.posx-30*obj.flip, obj.posy+15+26)
			obj.drawing_coords[6] = (obj.posx-30*obj.flip, obj.posy+15+26/2, obj.posx-30*obj.flip-15*obj.flip, obj.posy+15+26/2)
			obj.drawing_coords[7] = (obj.posx-3*obj.flip, obj.posy+15+26, obj.posx-10*obj.flip, obj.posy+15+26-7)
			obj.drawing_coords[8] = (obj.posx-3*obj.flip, obj.posy+15+26, obj.posx-10*obj.flip, obj.posy+15+26+7)
		return 9

	if position==3:
		if draw:
			obj.drawing_coords[0] = (obj.startx, obj.starty, obj.endx, obj.endy)
		return 1

	if position==4:
		if draw:
			obj.drawing_coords[0] = (obj.posx, obj.posy, obj.posx, obj.posy+15)
			obj.drawing_coords[1] = (obj.posx, obj.posy+15, obj.posx-26*obj.flip, obj.posy+15)
			obj.drawing_coords[2] = (obj.posx-26*obj.flip, obj.posy+15, obj.posx-26*obj.flip, obj.posy+15+26)
			obj.drawing_coords[3] = (obj.posx-26*obj.flip, obj.posy+15+26, obj.posx, obj.posy+15+26)
			obj.drawing_coords[4] = (obj.posx, obj.posy+15+26, obj.posx, obj.posy+15+26+15)
			obj.drawing_coords[5] = (obj.posx-30*obj.flip, obj.posy+15, obj.posx-30*obj.flip, obj.posy+15+26)
			obj.drawing_coords[6] = (obj.posx-30*obj.flip, obj.posy+15+26/2, obj.posx-30*obj.flip-15*obj.flip, obj.posy+15+26/2)
			obj.drawing_coords[7] = (obj.posx-10*obj.flip, obj.posy+15, obj.posx-3*obj.flip, obj.posy+15-7)
			obj.drawing_coords[8] = (obj.posx-10*obj.flip, obj.posy+15, obj.posx-3*obj.flip, obj.posy+15+7)
		return 9

	if position==5:
		if draw:
			obj.drawing_coords[0] = (obj.posx, obj.posy, obj.posx, obj.posy+20)
			obj.drawing_coords[1] = (obj.posx, obj.posy+60, obj.posx, obj.posy+80)
			obj.drawing_coords[2] = (obj.posx-20, obj.posy+20, obj.posx+20, obj.posy+60)
			obj.drawing_coords[3] = (obj.posx-19, obj.posy+21, obj.posx+19, obj.posy+59)
			obj.drawing_coords[4] = (obj.posx, obj.posy+20+10, obj.posx, obj.posy+20+30)
			obj.drawing_coords[5] = (obj.posx, obj.posy+20+30, obj.posx-8, obj.posy+20+30-8)
			obj.drawing_coords[6] = (obj.posx, obj.posy+20+30, obj.posx+8, obj.posy+20+30-8)
		return 7

	if position==6:
		if draw:
			obj.drawing_coords[0] = (obj.posx, obj.posy, obj.posx, obj.posy+20)
			obj.drawing_coords[1] = (obj.posx, obj.posy+60, obj.posx, obj.posy+80)
			obj.drawing_coords[2] = (obj.posx-20, obj.posy+20, obj.posx+20, obj.posy+60)
			obj.drawing_coords[3] = (obj.posx-19, obj.posy+21, obj.posx+19, obj.posy+59)
			obj.drawing_coords[4] = (obj.posx, obj.posy+20+10, obj.posx, obj.posy+20+18)
			obj.drawing_coords[5] = (obj.posx-4, obj.posy+20+10+4, obj.posx+4, obj.posy+20+18-4)
			obj.drawing_coords[6] = (obj.posx-4, obj.posy+20+30, obj.posx+4, obj.posy+20+30)
		return 7

	if position==7:
		if draw:
			obj.drawing_coords[0] = (obj.posx, obj.posy, obj.posx, obj.posy+20)
			obj.drawing_coords[1] = (obj.posx, obj.posy+60, obj.posx, obj.posy+80)
			obj.drawing_coords[2] = (obj.posx, obj.posy+20, obj.posx-20, obj.posy+40)
			obj.drawing_coords[3] = (obj.posx, obj.posy+20, obj.posx+20, obj.posy+40)
			obj.drawing_coords[4] = (obj.posx-20, obj.posy+40, obj.posx, obj.posy+60)
			obj.drawing_coords[5] = (obj.posx+20, obj.posy+40, obj.posx, obj.posy+60)
			obj.drawing_coords[6] = (obj.posx, obj.posy+28, obj.posx, obj.posy+38)
			obj.drawing_coords[7] = (obj.posx-5, obj.posy+28+5, obj.posx+5, obj.posy+38-5)
			obj.drawing_coords[8] = (obj.posx-5, obj.posy+30+5+15, obj.posx+5, obj.posy+40-5+15)
		return 9

	if position==8:
		if draw:
			obj.drawing_coords[0] = (obj.posx, obj.posy, obj.posx, obj.posy+20)
			obj.drawing_coords[1] = (obj.posx, obj.posy+60, obj.posx, obj.posy+80)
			obj.drawing_coords[2] = (obj.posx, obj.posy+20, obj.posx-20, obj.posy+40)
			obj.drawing_coords[3] = (obj.posx, obj.posy+20, obj.posx+20, obj.posy+40)
			obj.drawing_coords[4] = (obj.posx-20, obj.posy+40, obj.posx, obj.posy+60)
			obj.drawing_coords[5] = (obj.posx+20, obj.posy+40, obj.posx, obj.posy+60)
			obj.drawing_coords[6] = (obj.posx, obj.posy+28, obj.posx, obj.posy+50)
			obj.drawing_coords[7] = (obj.posx, obj.posy+50, obj.posx-5, obj.posy+50-5)
			obj.drawing_coords[8] = (obj.posx, obj.posy+50, obj.posx+5, obj.posy+50-5)
		return 9

	if position==9:
		if draw:
			obj.drawing_coords[0] = (obj.posx, obj.posy, obj.posx+obj.widthx, obj.posy)
			obj.drawing_coords[1] = (obj.posx, obj.posy, obj.posx, obj.posy-obj.heighty)
			for i in range(0,obj.partsx):
				obj.drawing_coords[i+2] = (obj.posx+obj.widthx*i/obj.partsx, obj.posy-5, obj.posx+obj.widthx*i/obj.partsx, obj.posy+5)
			for i in range(0,obj.partsy):
				obj.drawing_coords[i+obj.partsx+2] = (obj.posx-5, obj.posy-obj.heighty*i/obj.partsy, obj.posx+5, obj.posy-obj.heighty*i/obj.partsy)
		return obj.partsx+obj.partsy+2

	if position==10:
		if draw:
			for i in range(0, len(obj.xvals)-1):
				obj.drawing_coords[i] = (obj.xvals[obj.items-i-1], obj.yvals[obj.items-i-1], obj.xvals[obj.items-i-2], obj.yvals[obj.items-i-2])

		return len(obj.xvals)

	if position==11:
		if draw:
			obj.drawing_coords[0] = (obj.posx, obj.posy)

		return 1

def initvars(obj, name, makecenter):

	obj.name = name
	obj.itr_move = 0
	obj.itr_rotate = 0
	obj.itr_scale = 0
	obj.no_of_itr_move = 0
	obj.no_of_itr_rotate = 0
	obj.no_of_itr_scale = 0
	obj.anglex = 0
	obj.angle = 0
	obj.deltax = 0
	obj.deltay = 0
	obj.deltax1 = 0
	obj.deltay1 = 0
	obj.scalex_int = 0
	obj.scaley_int = 0
	obj.scalex = 0
	obj.scaley = 0
	obj.appear = 0
	obj.scale_around_y = 0
	obj.scale_around_x = 0
	if obj.name == 'text':
		obj.anglem = 0
		obj.anglex = 0
		obj.angle = 0
		obj.deltax = 0
		obj.deltay = 0
		obj.deltaxm = 0
		obj.deltaym = 0
		obj.deltax1 = 0
		obj.deltay1 = 0
		obj.scale_int = 0
		obj.scalev = 0
	obj.items = drawings(obj, dict_of_elements[name], 0)
	obj.drawing = [0]*obj.items
	obj.drawing_coords = [0]*obj.items
	obj.drawing_coords_tmp = [0]*obj.items
	obj.drawing_coords_tmp2 = [0]*obj.items
	obj.drawing_coords_tmp3 = [0]*obj.items
	obj.xmin_tmp = [0]*obj.items
	obj.xmax_tmp = [0]*obj.items
	obj.ymin_tmp = [0]*obj.items
	obj.ymax_tmp = [0]*obj.items
	drawings(obj, dict_of_elements[name], 1)
	if makecenter:
		reCenter(obj)
	draw(obj, dict_of_elements[name])
	obj.center = getCenter(obj)
	return obj


def scrolate_kernal(obj):
	if obj.name == 'text':
		c.itemconfig(obj.drawing, font=("Arial", round(obj.scale_int)))
		c.itemconfig(obj.drawing, angle=obj.anglex)
		c.coords(obj.drawing, obj.deltax1, obj.deltay1)
		return

	if obj.name == 'plottrace':
		offset = 1

	for i in range(0,obj.items-offset):
		if obj.appear:
			obj.evl = (math_appear(obj.drawing_coords[i][0],obj.drawing_coords[i][1],obj.scale_around_x,obj.scale_around_y,obj.scalex_int,obj.scaley_int)[0],
					math_appear(obj.drawing_coords[i][0],obj.drawing_coords[i][1],obj.scale_around_x,obj.scale_around_y,obj.scalex_int,obj.scaley_int)[1],
					math_appear(obj.drawing_coords[i][2],obj.drawing_coords[i][3],obj.scale_around_x,obj.scale_around_y,obj.scalex_int,obj.scaley_int)[0],
					math_appear(obj.drawing_coords[i][2],obj.drawing_coords[i][3],obj.scale_around_x,obj.scale_around_y,obj.scalex_int,obj.scaley_int)[1])
		else:
			obj.evl = (math_scale(obj.drawing_coords[i][0],obj.drawing_coords[i][1],obj.scale_around_x,obj.scale_around_y,obj.scalex_int,obj.scaley_int)[0],
					math_scale(obj.drawing_coords[i][0],obj.drawing_coords[i][1],obj.scale_around_x,obj.scale_around_y,obj.scalex_int,obj.scaley_int)[1],
					math_scale(obj.drawing_coords[i][2],obj.drawing_coords[i][3],obj.scale_around_x,obj.scale_around_y,obj.scalex_int,obj.scaley_int)[0],
					math_scale(obj.drawing_coords[i][2],obj.drawing_coords[i][3],obj.scale_around_x,obj.scale_around_y,obj.scalex_int,obj.scaley_int)[1])

		c.coords(obj.drawing[i], obj.evl)

		obj.drawing_coords_tmp[i] = obj.evl


		obj.evl = (math_rotate(obj.drawing_coords_tmp[i][0], obj.drawing_coords_tmp[i][1], obj.rotate_around_x, obj.rotate_around_y, 1*obj.anglex/(1))[0],
					math_rotate(obj.drawing_coords_tmp[i][0], obj.drawing_coords_tmp[i][1], obj.rotate_around_x, obj.rotate_around_y, 1*obj.anglex/(1))[1],
					math_rotate(obj.drawing_coords_tmp[i][2], obj.drawing_coords_tmp[i][3], obj.rotate_around_x, obj.rotate_around_y, 1*obj.anglex/(1))[0],
					math_rotate(obj.drawing_coords_tmp[i][2], obj.drawing_coords_tmp[i][3], obj.rotate_around_x, obj.rotate_around_y, 1*obj.anglex/(1))[1])

		c.coords(obj.drawing[i], obj.evl)

		obj.drawing_coords_tmp2[i] = obj.evl

		obj.evl = (math_translate(obj.drawing_coords_tmp2[i][0], obj.drawing_coords_tmp2[i][1], obj.deltax1, obj.deltay1)[0],
					math_translate(obj.drawing_coords_tmp2[i][0], obj.drawing_coords_tmp2[i][1], obj.deltax1, obj.deltay1)[1],
					math_translate(obj.drawing_coords_tmp2[i][2], obj.drawing_coords_tmp2[i][3], obj.deltax1, obj.deltay1)[0],
					math_translate(obj.drawing_coords_tmp2[i][2], obj.drawing_coords_tmp2[i][3], obj.deltax1, obj.deltay1)[1])

		c.coords(obj.drawing[i], obj.evl)

		obj.drawing_coords_tmp3[i] = (obj.evl)

def scrolate_init(obj, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time):
	
	obj.itr_move = 0
	obj.deltax1 = 0
	obj.deltay1 = 0
	obj.animation_time = animation_time
	obj.deltax = newx - getCenter(obj)[0]
	obj.deltay = newy - getCenter(obj)[1]
	obj.no_of_itr_move = obj.animation_time * framerate/1000
	obj.no_of_itr_scale = obj.animation_time * framerate/1000
	obj.no_of_itr_rotate = obj.animation_time * framerate/1000

	obj.angle = rotateby
	obj.anglex = 0
	obj.itr_rotate = 0
	obj.itr_scale = 0
	obj.scalex = scalex
	obj.scaley = scaley

	obj.scale_around_x = around_sx
	obj.scale_around_y = around_sy

	obj.rotate_around_x = around_rx
	obj.rotate_around_y = around_ry

def scrolate_init_text(obj, newposx, newposy, angle, scale, animation_time):
	obj.animation_time = animation_time
	obj.anglem = angle
	obj.anglex = 0
	obj.itr_rotate = 0
	obj.itr_move = 0
	obj.deltax1 = 0
	obj.deltay1 = 0
	obj.deltaxm = newposx
	obj.deltaym = newposy
	obj.itr_scale = 0
	obj.scale_int = 0
	obj.scalev = scale
	obj.no_of_itr_rotate = obj.animation_time*framerate/1000
	obj.no_of_itr_scale = obj.animation_time*framerate/1000
	obj.no_of_itr_move = obj.animation_time*framerate/1000

def scrolate_trigger(obj):

	if obj.name == 'text':
		obj.itr_rotate = obj.itr_rotate + 1

		obj.anglex = obj.angle + obj.anglem * ease(obj.itr_rotate/obj.no_of_itr_rotate)

		obj.itr_move = obj.itr_move + 1

		obj.deltax1 = obj.posx + (obj.deltaxm - obj.posx) * ease(obj.itr_move/obj.no_of_itr_move)
		obj.deltay1 = obj.posy + (obj.deltaym - obj.posy) * ease(obj.itr_move/obj.no_of_itr_move)

		obj.itr_scale = obj.itr_scale + 1

		obj.scale_int = obj.fontsize + (obj.scalev - obj.fontsize) * ease(obj.itr_scale/ obj.no_of_itr_scale)

		if obj.itr_scale > obj.no_of_itr_scale:
			obj.fontsize = obj.scale_int
			obj.angle = obj.angle + obj.anglex
			obj.posx = obj.deltax1
			obj.posy = obj.deltay1
			return

	if obj.name != 'text':
		obj.itr_move = obj.itr_move + 1

		obj.deltax1 = obj.deltax * ease(obj.itr_move/obj.no_of_itr_move)
		obj.deltay1 = obj.deltay * ease(obj.itr_move/obj.no_of_itr_move)

		obj.itr_rotate = obj.itr_rotate + 1

		obj.anglex = obj.angle*ease(obj.itr_rotate/obj.no_of_itr_rotate)

		obj.itr_scale = obj.itr_scale + 1

		obj.scalex_int = obj.scalex * ease(obj.itr_scale/ obj.no_of_itr_scale)
		obj.scaley_int = obj.scaley * ease(obj.itr_scale/ obj.no_of_itr_scale)

		if obj.itr_scale > obj.no_of_itr_scale:
			obj.drawing_coords = obj.drawing_coords_tmp3[:]
			obj.appear = 0
			return

	scrolate_kernal(obj)

	master.after(int(1000./framerate), obj.scrolate_trigger)


class RES:
	def init(self, posx, posy, color='#00FF00', width=2):
		self.posx = posx
		self.posy = posy
		self.color = color
		self.width = width
		self = initvars(self, 'res', 1)
		self.appear = 1
		self.scrolate(getCenter(self)[0],getCenter(self)[1], 0, 1.3, 1.3, getCenter(self)[0],getCenter(self)[1],getCenter(self)[0],getCenter(self)[1],500)

	def terminals(self):
		return [[self.drawing_coords[0][0], self.drawing_coords[0][1]], [self.drawing_coords[11][2], self.drawing_coords[11][3]]]

	def scrolate(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time):
		scrolate_init(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time)
		self.scrolate_trigger()

	def scrolate_trigger(self):
		scrolate_trigger(self)

class NMOS:
	def init(self, posx, posy, flip=1, color='#00FF00', width=2):
		self.posx = posx
		self.posy = posy
		self.color = color
		self.width = width
		self.flip = flip
		self = initvars(self, 'nmos', 1)
		self.appear = 1
		self.scrolate(getCenter(self)[0],getCenter(self)[1], 0, 1.3, 1.3, getCenter(self)[0],getCenter(self)[1],getCenter(self)[0],getCenter(self)[1],500)

	def terminals(self):
		return [[self.drawing_coords[0][0], self.drawing_coords[0][1]], [self.drawing_coords[6][2], self.drawing_coords[6][3]], [self.drawing_coords[4][2], self.drawing_coords[4][3]]]

	def scrolate(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time):
		scrolate_init(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time)
		self.scrolate_trigger()

	def scrolate_trigger(self):
		scrolate_trigger(self)

class PMOS:
	def init(self, posx, posy, flip=1, color='#00FF00', width=2):
		self.posx = posx
		self.posy = posy
		self.color = color
		self.width = width
		self.flip = flip
		self = initvars(self, 'pmos', 1)
		self.appear = 1
		self.scrolate(getCenter(self)[0],getCenter(self)[1], 0, 1.3, 1.3, getCenter(self)[0],getCenter(self)[1],getCenter(self)[0],getCenter(self)[1],500)

	def terminals(self):
		return [[self.drawing_coords[0][0], self.drawing_coords[0][1]], [self.drawing_coords[6][2], self.drawing_coords[6][3]], [self.drawing_coords[4][2], self.drawing_coords[4][3]]]

	def scrolate(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time):
		scrolate_init(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time)
		self.scrolate_trigger()

	def scrolate_trigger(self):
		scrolate_trigger(self)

class WIRE:
	def init(self, startx, starty, endx, endy, width=2, color='#FFFFFF'):
		self.startx = startx
		self.starty = starty
		self.endx = endx
		self.endy = endy
		self.width = width
		self.color = color
		self = initvars(self, 'wire', 0)
		self.appear = 1
		self.scrolate(getCenter(self)[0],getCenter(self)[1], 0, 1, 1, getCenter(self)[0],getCenter(self)[1],getCenter(self)[0],getCenter(self)[1],500)
		return self

	def terminals(self):
		return [[self.drawing_coords[0][0], self.drawing_coords[0][1]], [self.drawing_coords[0][2], self.drawing_coords[0][3]]]

	def scrolate(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time):
		scrolate_init(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time)
		self.scrolate_trigger()

	def scrolate_trigger(self):
		scrolate_trigger(self)

class IDC:
	def init(self, posx, posy, width=2, color='#FFFFFF'):
		self.posx = posx
		self.posy = posy
		self.width = width
		self.color = color
		self = initvars(self, 'idc', 0)
		self.appear = 1
		self.scrolate(getCenter(self)[0],getCenter(self)[1], 0, 1.3, 1.3, getCenter(self)[0],getCenter(self)[1],getCenter(self)[0],getCenter(self)[1],500)

	def terminals(self):
		return [[self.drawing_coords[0][0], self.drawing_coords[0][1]], [self.drawing_coords[1][2], self.drawing_coords[1][3]]]

	def scrolate(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time):
		scrolate_init(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time)
		self.scrolate_trigger()

	def scrolate_trigger(self):
		scrolate_trigger(self)

class VDC:
	def init(self, posx, posy, width=2, color='#FFFFFF'):
		self.posx = posx
		self.posy = posy
		self.width = width
		self.color = color
		self = initvars(self, 'vdc', 0)
		self.appear = 1
		self.scrolate(getCenter(self)[0],getCenter(self)[1], 0, 1.3, 1.3, getCenter(self)[0],getCenter(self)[1],getCenter(self)[0],getCenter(self)[1],500)

	def terminals(self):
		return [[self.drawing_coords[0][0], self.drawing_coords[0][1]], [self.drawing_coords[1][2], self.drawing_coords[1][3]]]

	def scrolate(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time):
		scrolate_init(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time)
		self.scrolate_trigger()

	def scrolate_trigger(self):
		scrolate_trigger(self)

class VCVS:
	def init(self, posx, posy, width=2, color='#FFFFFF'):
		self.posx = posx
		self.posy = posy
		self.width = width
		self.color = color
		self = initvars(self, 'vcvs', 0)
		self.appear = 1
		self.scrolate(getCenter(self)[0],getCenter(self)[1], 0, 1.3, 1.3, getCenter(self)[0],getCenter(self)[1],getCenter(self)[0],getCenter(self)[1],500)

	def terminals(self):
		return [[self.drawing_coords[0][0], self.drawing_coords[0][1]], [self.drawing_coords[1][2], self.drawing_coords[1][3]]]

	def scrolate(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time):
		scrolate_init(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time)
		self.scrolate_trigger()

	def scrolate_trigger(self):
		scrolate_trigger(self)

class VCCS:
	def init(self, posx, posy, width=2, color='#FFFFFF'):
		self.posx = posx
		self.posy = posy
		self.width = width
		self.color = color
		self = initvars(self, 'vccs', 0)
		self.appear = 1
		self.scrolate(getCenter(self)[0],getCenter(self)[1], 0, 1.3, 1.3, getCenter(self)[0],getCenter(self)[1],getCenter(self)[0],getCenter(self)[1],500)

	def terminals(self):
		return [[self.drawing_coords[0][0], self.drawing_coords[0][1]], [self.drawing_coords[1][2], self.drawing_coords[1][3]]]

	def scrolate(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time):
		scrolate_init(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time)
		self.scrolate_trigger()

	def scrolate_trigger(self):
		scrolate_trigger(self)

class PLOT:
	def init(self, posx, posy, widthx, heighty, partsx, partsy, width=2, color='#FFFFFF'):
		self.posx = posx
		self.posy = posy
		self.widthx = widthx
		self.heighty = heighty
		self.partsx = partsx
		self.partsy = partsy
		self.width = width
		self.color = color
		self = initvars(self, 'plot', 1)
		self.appear = 1
		self.scrolate(getCenter(self)[0],getCenter(self)[1], 0, 1, 1, getCenter(self)[0],getCenter(self)[1],getCenter(self)[0],getCenter(self)[1],500)

	def scrolate(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time):
		scrolate_init(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time)
		self.scrolate_trigger()

	def scrolate_trigger(self):
		scrolate_trigger(self)

class PLOT_TRACE:
	def init(self, points, nTimes=50, width=2, color='#FFFFFF'):
		self.xvals, self.yvals = bezier_curve(points, nTimes)
		self.width = width
		self.color = color
		self = initvars(self, 'plottrace', 0)
		self.appear = 1
		self.scrolate(getCenter(self)[0],getCenter(self)[1], 0, 1, 1, getCenter(self)[0],getCenter(self)[1],getCenter(self)[0],getCenter(self)[1],500)

	def scrolate(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time):
		scrolate_init(self, newx, newy, rotateby, scalex, scaley, around_sx, around_sy, around_rx, around_ry, animation_time)
		self.scrolate_trigger()

	def scrolate_trigger(self):
		scrolate_trigger(self)

class TEXT:
	def init(self, posx, posy, text, fontsize=30, angle=0, color='#FFFFFF'):
		self.posx = posx
		self.posy = posy
		self.text = text
		self.fontsize = fontsize
		self.angle = angle
		self.color = color
		self = initvars(self, 'text', 0)
		self.appear = 0

	def scrolate(self, newposx, newposy, angle, targetfont, animation_time):
		scrolate_init_text(self, newposx, newposy, angle, targetfont, animation_time)
		self.scrolate_trigger()

	def scrolate_trigger(self):
		scrolate_trigger(self)

class Orchestra:

	def __init__(self):
		self.seq = 0

	def animate_kernal(self):
		if self.seq == 0:
			pl.init(100,100,'Hello World!')
		if self.seq == 1:
			pl.scrolate(200,200,-90,60,500)
		return

	def animate_trigger(self):
		self.animate_kernal()
		self.seq = self.seq + 1
		master.after(2000, self.animate_trigger)

	def animate(self):
	 	self.animate_trigger()

pl = TEXT()
o = Orchestra()
o.animate()

master.mainloop()
