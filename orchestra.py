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
framerate = 70

def cos(x):
    return math.cos(x)
def sin(x): 
    return math.sin(x)
def tan(x): 
    return math.tan(x)

c = Canvas(master, bg='#000000', bd=0, height=cheight, width=cwidth, highlightthickness=0, relief='ridge')
c.pack()

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

    polynomial_array = np.array([ bernstein_poly(i, nPoints-1, t) for i in range(0, nPoints)   ])

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


class DOT:
	def init(self, posx, posy):
		self.posx = posx
		self.posy = posy
		self.items = 1
		self.dot_p = [0]*self.items
		self.dot_p[0] = c.create_oval(posx-2, posy-2, posx+2, posy+2, outline='#FFFFFF',fill='#FFFFFF', width=2)

class BOX:
	def init(self, posx, posy, width, height):
		self.posx = posx;
		self.posy = posy;
		self.items = 1
		self.box_p = [0]*self.items
		self.box_p_coords = [0]*self.items
		self.box_p_coords_tmp = [0]*self.items
		self.box_p_coords_tmp2 = [0]*self.items
		self.box_p_coords_tmp3 = [0]*self.items
		self.xmin_tmp = [0]*self.items
		self.xmax_tmp = [0]*self.items
		self.ymin_tmp = [0]*self.items
		self.ymax_tmp = [0]*self.items

		self.itr_move = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.no_of_itr_move = 0
		self.no_of_itr_rotate = 0
		self.no_of_itr_scale = 0
		self.anglex = 0
		self.angle = 0
		self.deltax = 0
		self.deltay = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.scalex_int = 0
		self.scaley_int = 0
		self.deltay1 = 0
		self.scalex = 0
		self.scaley = 0

		self.box_p_coords[0] = (posx-width/2, posy-height/2, posx+width/2, posy-height/2, posx+width/2, posy+height/2, posx-width/2, posy+height/2)

		self.reCenter()

		self.box_p[0] = c.create_polygon(self.box_p_coords[0], fill='#FFFFFF', width=2)

		self.appear = 1
		self.scale(1,1,500)


	def reCenter(self):

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]

	def getCenter(self):	
		return [self.posx, self.posy]

	def rotate_kernal(self):

		for i in range(0, self.items):

			self.evl = (math_rotate(self.box_p_coords[i][0], self.box_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
						math_rotate(self.box_p_coords[i][0], self.box_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
						math_rotate(self.box_p_coords[i][2], self.box_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
						math_rotate(self.box_p_coords[i][2], self.box_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
						math_rotate(self.box_p_coords[i][4], self.box_p_coords[i][5], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
						math_rotate(self.box_p_coords[i][4], self.box_p_coords[i][5], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
						math_rotate(self.box_p_coords[i][6], self.box_p_coords[i][7], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
						math_rotate(self.box_p_coords[i][6], self.box_p_coords[i][7], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.box_p[i], self.evl)

			self.box_p_coords_tmp[i] = (self.evl)

	def rotate_trigger(self):
		
		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		if self.itr_rotate > self.no_of_itr_rotate:
			self.box_p_coords = self.box_p_coords_tmp[:]
			return

		self.rotate_kernal()
		master.after(int(1000/framerate), self.rotate_trigger)

	def rotate(self, angle, animation_time):
	 	self.animation_time = animation_time
	 	self.angle = angle
	 	self.anglex = 0
	 	self.itr_rotate = 0
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.no_of_itr_rotate = self.animation_time*framerate/1000
	 	self.rotate_trigger()

	def move_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_translate(self.box_p_coords[i][0], self.box_p_coords[i][1], self.deltax1, self.deltay1)[0],
						math_translate(self.box_p_coords[i][0], self.box_p_coords[i][1], self.deltax1, self.deltay1)[1],
						math_translate(self.box_p_coords[i][2], self.box_p_coords[i][3], self.deltax1, self.deltay1)[0],
						math_translate(self.box_p_coords[i][2], self.box_p_coords[i][3], self.deltax1, self.deltay1)[1],
						math_translate(self.box_p_coords[i][4], self.box_p_coords[i][5], self.deltax1, self.deltay1)[0],
						math_translate(self.box_p_coords[i][4], self.box_p_coords[i][5], self.deltax1, self.deltay1)[1],
						math_translate(self.box_p_coords[i][6], self.box_p_coords[i][7], self.deltax1, self.deltay1)[0],
						math_translate(self.box_p_coords[i][6], self.box_p_coords[i][7], self.deltax1, self.deltay1)[1])

			c.coords(self.box_p[i], self.evl)

			self.box_p_coords_tmp[i] = (self.evl)

	def move_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		if self.itr_move > self.no_of_itr_move:
			self.box_p_coords = self.box_p_coords_tmp[:]
			# for i in range(0, self.items):
			# 	self.box_p_coords_tmp[i] = (0,0,0,0)
			return
		self.move_kernal()
		master.after(int(1000./framerate), self.move_trigger)

	def translate(self, newposx, newposy, animation_time):
	 	self.itr_move = 0
	 	self.deltax1 = 0
	 	self.deltay1 = 0
	 	self.no_of_itr_move = self.animation_time * framerate/1000
	 	self.animation_time = animation_time
	 	self.deltax = newposx - self.getCenter()[0]
	 	self.deltay = newposy - self.getCenter()[1]
	 	self.move_trigger()


	def scrolate_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_scale(self.box_p_coords[i][0],self.box_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.box_p_coords[i][0],self.box_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
						math_scale(self.box_p_coords[i][2],self.box_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.box_p_coords[i][2],self.box_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
						math_scale(self.box_p_coords[i][4],self.box_p_coords[i][5],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.box_p_coords[i][4],self.box_p_coords[i][5],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
						math_scale(self.box_p_coords[i][6],self.box_p_coords[i][7],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.box_p_coords[i][6],self.box_p_coords[i][7],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

			c.coords(self.box_p[i], self.evl)

			self.box_p_coords_tmp[i] = self.evl

			self.evl = (math_rotate(self.box_p_coords_tmp[i][0], self.box_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
						math_rotate(self.box_p_coords_tmp[i][0], self.box_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
						math_rotate(self.box_p_coords_tmp[i][2], self.box_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
						math_rotate(self.box_p_coords_tmp[i][2], self.box_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
						math_rotate(self.box_p_coords_tmp[i][4], self.box_p_coords_tmp[i][5], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
						math_rotate(self.box_p_coords_tmp[i][4], self.box_p_coords_tmp[i][5], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
						math_rotate(self.box_p_coords_tmp[i][6], self.box_p_coords_tmp[i][7], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
						math_rotate(self.box_p_coords_tmp[i][6], self.box_p_coords_tmp[i][7], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.box_p[i], self.evl)

			self.box_p_coords_tmp2[i] = self.evl

			self.evl = (math_translate(self.box_p_coords_tmp2[i][0], self.box_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[0],
						math_translate(self.box_p_coords_tmp2[i][0], self.box_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[1],
						math_translate(self.box_p_coords_tmp2[i][2], self.box_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[0],
						math_translate(self.box_p_coords_tmp2[i][2], self.box_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[1],
						math_translate(self.box_p_coords_tmp2[i][4], self.box_p_coords_tmp2[i][5], self.deltax1, self.deltay1)[0],
						math_translate(self.box_p_coords_tmp2[i][4], self.box_p_coords_tmp2[i][5], self.deltax1, self.deltay1)[1],
						math_translate(self.box_p_coords_tmp2[i][6], self.box_p_coords_tmp2[i][7], self.deltax1, self.deltay1)[0],
						math_translate(self.box_p_coords_tmp2[i][6], self.box_p_coords_tmp2[i][7], self.deltax1, self.deltay1)[1])

			c.coords(self.box_p[i], self.evl)

			self.box_p_coords_tmp3[i] = (self.evl)

			
	def scrolate_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		# if self.itr_rotate > self.no_of_itr_rotate:
		# 	self.box_p_coords = self.box_p_coords_tmp2[:]
		# 	return

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.box_p_coords = self.box_p_coords_tmp3[:]
			# for i in range(0, self.items):
			# 	self.box_p_coords_tmp[i] = (0,0,0,0)
			return	
		
		self.scrolate_kernal()
		master.after(int(1000./framerate), self.scrolate_trigger)

	def scrolate(self, newposx, newposy, angle, scalex, scaley, animation_time):
		self.itr_move = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.animation_time = animation_time
		self.deltax = newposx - self.getCenter()[0]
		self.deltay = newposx - self.getCenter()[0]

		self.angle = angle
		self.anglex = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.no_of_itr_move = self.animation_time * framerate/1000
		self.no_of_itr_rotate = self.no_of_itr_move
		self.no_of_itr_scale = self.no_of_itr_rotate

		self.scaley = scaley
		self.scalex = scalex

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]
		self.scrolate_trigger()

	def scale_kernal(self):

		if self.appear:
			for i in range(0, self.items):

				self.evl =  (math_appear(self.box_p_coords[i][0],self.box_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.box_p_coords[i][0],self.box_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_appear(self.box_p_coords[i][2],self.box_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.box_p_coords[i][2],self.box_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_appear(self.box_p_coords[i][4],self.box_p_coords[i][5],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.box_p_coords[i][4],self.box_p_coords[i][5],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_appear(self.box_p_coords[i][6],self.box_p_coords[i][7],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.box_p_coords[i][6],self.box_p_coords[i][7],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.box_p[i],self.evl)

				self.box_p_coords_tmp[i] = self.evl
		else:
			for i in range(0, self.items):

				self.evl =  (math_scale(self.box_p_coords[i][0],self.box_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.box_p_coords[i][0],self.box_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_scale(self.box_p_coords[i][2],self.box_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.box_p_coords[i][2],self.box_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_scale(self.box_p_coords[i][4],self.box_p_coords[i][5],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.box_p_coords[i][4],self.box_p_coords[i][5],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_scale(self.box_p_coords[i][6],self.box_p_coords[i][7],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.box_p_coords[i][6],self.box_p_coords[i][7],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.box_p[i],self.evl)

				self.box_p_coords_tmp[i] = self.evl

	def scale_trigger(self):

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.box_p_coords = self.box_p_coords_tmp[:]
			self.appear = 0
			return

		self.scale_kernal()
		master.after(int(1000./framerate), self.scale_trigger)

	def scale(self, scalex, scaley, animation_time):
	 	self.animation_time = animation_time
	 	self.itr_scale = 0
	 	self.scaley = scaley
	 	self.scalex = scalex
	 	self.no_of_itr_scale = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.scale_trigger()

	
class RES:

	def init(self, posx, posy, color='#FFFFFF'):

		self.posx = posx
		self.posy = posy
		self.color = color

		self.items = 12
		self.res_p = [0]*self.items
		self.res_p_coords = [0]*self.items
		self.res_p_coords_tmp = [0]*self.items
		self.res_p_coords_tmp2 = [0]*self.items
		self.res_p_coords_tmp3 = [0]*self.items
		self.xmin_tmp = [0]*self.items
		self.xmax_tmp = [0]*self.items
		self.ymin_tmp = [0]*self.items
		self.ymax_tmp = [0]*self.items

		self.itr_move = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.no_of_itr_move = 0
		self.no_of_itr_rotate = 0
		self.no_of_itr_scale = 0
		self.anglex = 0
		self.angle = 0
		self.deltax = 0
		self.deltay = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.scalex_int = 0
		self.scaley_int = 0
		self.deltay1 = 0
		self.scalex = 0
		self.scaley = 0

		self.res_p_coords[0] = (posx, posy, posx+15, posy)
		self.res_p_coords[1] = (posx+15, posy, posx+15+4, posy+9)
		self.res_p_coords[2] = (posx+15+4, posy+9, posx+15+4+4, posy-9)
		self.res_p_coords[3] = (posx+15+4+4, posy-9, posx+15+4+4+4, posy+9)
		self.res_p_coords[4] = (posx+15+4+4+4, posy+9, posx+15+4+4+4+4, posy-9)
		self.res_p_coords[5] = (posx+15+4+4+4+4, posy-9, posx+15+4+4+4+4+4, posy+9)
		self.res_p_coords[6] = (posx+15+4+4+4+4+4, posy+9, posx+15+4+4+4+4+4+4, posy-9)
		self.res_p_coords[7] = (posx+15+4+4+4+4+4+4, posy-9, posx+15+4+4+4+4+4+4+4, posy+9)
		self.res_p_coords[8] = (posx+15+4+4+4+4+4+4+4, posy+9, posx+15+4+4+4+4+4+4+4+4, posy-9)
		self.res_p_coords[9] = (posx+15+4+4+4+4+4+4+4+4, posy-9, posx+15+4+4+4+4+4+4+4+4+4, posy+9)
		self.res_p_coords[10] = (posx+15+4+4+4+4+4+4+4+4+4, posy+9, posx+15+4+4+4+4+4+4+4+4+4+4, posy)
		self.res_p_coords[11] = (posx+15+4+4+4+4+4+4+4+4+4+4, posy, posx+15+4+4+4+4+4+4+4+4+4+4+15, posy)

		self.reCenter()

		self.res_p[0] = c.create_line(self.res_p_coords[0], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.res_p[1] = c.create_line(self.res_p_coords[1], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.res_p[2] = c.create_line(self.res_p_coords[2], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.res_p[3] = c.create_line(self.res_p_coords[3], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.res_p[4] = c.create_line(self.res_p_coords[4], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.res_p[5] = c.create_line(self.res_p_coords[5], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.res_p[6] = c.create_line(self.res_p_coords[6], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.res_p[7] = c.create_line(self.res_p_coords[7], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.res_p[8] = c.create_line(self.res_p_coords[8], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.res_p[9] = c.create_line(self.res_p_coords[9], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.res_p[10] = c.create_line(self.res_p_coords[10], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.res_p[11] = c.create_line(self.res_p_coords[11], arrow='none', fill=self.color, width=2, capstyle=ROUND)

		self.appear = 1
		self.scale(1,1,500)


	def getTerminals(self):
		return [[self.res_p_coords[0][0], self.res_p_coords[0][1]], [self.res_p_coords[11][2], self.res_p_coords[11][3]]]

	def reCenter(self):
		self.rotate_around_x_old = self.getCenter()[0]
		self.rotate_around_y_old = self.getCenter()[1]

		for i in range(0, self.items):
			self.res_p_coords[i] = (self.res_p_coords[i][0] + self.posx - self.rotate_around_x_old, self.res_p_coords[i][1] + self.posy - self.rotate_around_y_old, self.res_p_coords[i][2] + self.posx - self.rotate_around_x_old, self.res_p_coords[i][3] + self.posy -self.rotate_around_y_old)

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]

	def getCenter(self):
		for i in range(0, self.items):
			self.xmin_tmp[i] = min(self.res_p_coords[i][0], self.res_p_coords[i][2])
			self.xmax_tmp[i] = max(self.res_p_coords[i][0], self.res_p_coords[i][2])
			self.ymin_tmp[i] = min(self.res_p_coords[i][1], self.res_p_coords[i][3])
			self.ymax_tmp[i] = max(self.res_p_coords[i][1], self.res_p_coords[i][3])
		return [(min(self.xmin_tmp) + max(self.xmax_tmp))/2, (min(self.ymin_tmp) + max(self.ymax_tmp))/2]

	def rotate_kernal(self):

		for i in range(0, self.items):

			self.evl =  (math_rotate(self.res_p_coords[i][0], self.res_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
						math_rotate(self.res_p_coords[i][0], self.res_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
						math_rotate(self.res_p_coords[i][2], self.res_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
						math_rotate(self.res_p_coords[i][2], self.res_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])
			c.coords(self.res_p[i],self.evl)

			self.res_p_coords_tmp[i] = self.evl

	def rotate_trigger(self):

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		if self.itr_rotate > self.no_of_itr_rotate:
			self.res_p_coords = self.res_p_coords_tmp[:]
			return

		self.rotate_kernal()
		master.after(int(1000./framerate), self.rotate_trigger)

	def rotate(self, angle, animation_time):
	 	self.animation_time = animation_time
	 	self.angle = angle
	 	self.anglex = 0
	 	self.itr_rotate = 0
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.no_of_itr_rotate = self.animation_time*framerate/1000
	 	self.rotate_trigger()

	def move_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_translate(self.res_p_coords[i][0], self.res_p_coords[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.res_p_coords[i][0], self.res_p_coords[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.res_p_coords[i][2], self.res_p_coords[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.res_p_coords[i][2], self.res_p_coords[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.res_p[i], self.evl)

			self.res_p_coords_tmp[i] = (self.evl)

	def move_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		if self.itr_move > self.no_of_itr_move:
			self.res_p_coords = self.res_p_coords_tmp[:]
			# for i in range(0, self.items):
			# 	self.res_p_coords_tmp[i] = (0,0,0,0)
			return
		self.move_kernal()
		master.after(int(1000./framerate), self.move_trigger)

	def translate(self, newposx, newposy, animation_time):
	 	self.itr_move = 0
	 	self.deltax1 = 0
	 	self.deltay1 = 0
	 	self.animation_time = animation_time
	 	self.no_of_itr_move = self.animation_time * framerate/1000
	 	self.deltax = newposx - self.getCenter()[0]
	 	self.deltay = newposy - self.getCenter()[1]
	 	self.move_trigger()


	def scrolate_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_scale(self.res_p_coords[i][0],self.res_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.res_p_coords[i][0],self.res_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
						math_scale(self.res_p_coords[i][2],self.res_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.res_p_coords[i][2],self.res_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

			c.coords(self.res_p[i], self.evl)

			self.res_p_coords_tmp[i] = self.evl

			self.evl = (math_rotate(self.res_p_coords_tmp[i][0], self.res_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.res_p_coords_tmp[i][0], self.res_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.res_p_coords_tmp[i][2], self.res_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.res_p_coords_tmp[i][2], self.res_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.res_p[i], self.evl)

			self.res_p_coords_tmp2[i] = self.evl

			self.evl = (math_translate(self.res_p_coords_tmp2[i][0], self.res_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.res_p_coords_tmp2[i][0], self.res_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.res_p_coords_tmp2[i][2], self.res_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.res_p_coords_tmp2[i][2], self.res_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.res_p[i], self.evl)

			self.res_p_coords_tmp3[i] = (self.evl)

			
	def scrolate_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.res_p_coords = self.res_p_coords_tmp3[:]
			return	
		
		self.scrolate_kernal()
		master.after(int(1000./framerate), self.scrolate_trigger)

	def scrolate(self, newposx, newposy, angle, scalex, scaley, animation_time):
		self.itr_move = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.animation_time = animation_time
		self.deltax = newposx - self.getCenter()[0]
		self.deltay = newposx - self.getCenter()[0]
		self.no_of_itr_rotate = self.animation_time*framerate/1000
		self.no_of_itr_scale = self.animation_time*framerate/1000
		self.no_of_itr_move = self.animation_time*framerate/1000


		self.angle = angle
		self.anglex = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.scaley = scaley
		self.scalex = scalex
		self.no_of_itr_move = self.animation_time * framerate/1000

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]
		self.scrolate_trigger()

	def scale_kernal(self):

		if self.appear:

			for i in range(0, self.items):

				self.evl =  (math_appear(self.res_p_coords[i][0],self.res_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.res_p_coords[i][0],self.res_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_appear(self.res_p_coords[i][2],self.res_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.res_p_coords[i][2],self.res_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.res_p[i],self.evl)

				self.res_p_coords_tmp[i] = self.evl
		else:
			for i in range(0, self.items):

				self.evl =  (math_scale(self.res_p_coords[i][0],self.res_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.res_p_coords[i][0],self.res_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_scale(self.res_p_coords[i][2],self.res_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.res_p_coords[i][2],self.res_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.res_p[i],self.evl)

				self.res_p_coords_tmp[i] = self.evl

	def scale_trigger(self):

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.res_p_coords = self.res_p_coords_tmp[:]
			self.appear = 0
			return

		self.scale_kernal()
		master.after(int(1000./framerate), self.scale_trigger)

	def scale(self, scalex, scaley, animation_time):
	 	self.animation_time = animation_time
	 	self.itr_scale = 0
	 	self.scaley = scaley
	 	self.scalex = scalex
	 	self.no_of_itr_scale = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.scale_trigger()

class NMOS:

	def init(self, posx, posy, flip, color='#FFFFFF'):
		
		self.posx = posx
		self.posy = posy

		self.color = color

		if flip==0:
			self.flip = 1
		if flip==1:
			self.flip = -1

		self.items = 9
		self.mosfet_p = [0]*self.items
		self.mosfet_p_coords = [0]*self.items
		self.mosfet_p_coords_tmp = [0]*self.items
		self.mosfet_p_coords_tmp2 = [0]*self.items
		self.mosfet_p_coords_tmp3 = [0]*self.items
		self.xmin_tmp = [0]*self.items
		self.xmax_tmp = [0]*self.items
		self.ymin_tmp = [0]*self.items
		self.ymax_tmp = [0]*self.items

		self.itr_move = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.no_of_itr_move = 0
		self.no_of_itr_rotate = 0
		self.no_of_itr_scale = 0
		self.anglex = 0
		self.angle = 0
		self.deltax = 0
		self.deltay = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.scalex_int = 0
		self.scaley_int = 0
		self.deltay1 = 0
		self.scalex = 0
		self.scaley = 0

		self.mosfet_p_coords[0] = (posx, posy, posx, posy+15)
		self.mosfet_p_coords[1] = (posx, posy+15, posx-26*self.flip, posy+15)
		self.mosfet_p_coords[2] = (posx-26*self.flip, posy+15, posx-26*self.flip, posy+15+26)
		self.mosfet_p_coords[3] = (posx-26*self.flip, posy+15+26, posx, posy+15+26)
		self.mosfet_p_coords[4] = (posx, posy+15+26, posx, posy+15+26+15)
		self.mosfet_p_coords[5] = (posx-30*self.flip, posy+15, posx-30*self.flip, posy+15+26)
		self.mosfet_p_coords[6] = (posx-30*self.flip, posy+15+26/2, posx-30*self.flip-15*self.flip, posy+15+26/2)
		self.mosfet_p_coords[7] = (posx-3*self.flip, posy+15+26, posx-10*self.flip, posy+15+26-7)
		self.mosfet_p_coords[8] = (posx-3*self.flip, posy+15+26, posx-10*self.flip, posy+15+26+7)

		self.reCenter()

		self.mosfet_p[0]= c.create_line(self.mosfet_p_coords[0], arrow='none', fill=self.color, width=2, capstyle=PROJECTING)
		self.mosfet_p[1]= c.create_line(self.mosfet_p_coords[1], arrow='none', fill=self.color, width=2, capstyle=PROJECTING)
		self.mosfet_p[2]= c.create_line(self.mosfet_p_coords[2], arrow='none', fill=self.color, width=2, capstyle=PROJECTING)
		self.mosfet_p[3]= c.create_line(self.mosfet_p_coords[3], arrow='none', fill=self.color, width=2, capstyle=PROJECTING)
		self.mosfet_p[4]= c.create_line(self.mosfet_p_coords[4], arrow='none', fill=self.color, width=2, capstyle=PROJECTING)
		self.mosfet_p[5]= c.create_line(self.mosfet_p_coords[5], arrow='none', fill=self.color, width=2, capstyle=PROJECTING)
		self.mosfet_p[6]= c.create_line(self.mosfet_p_coords[6], arrow='none', fill=self.color, width=2, capstyle=PROJECTING)
		self.mosfet_p[7] = c.create_line(self.mosfet_p_coords[7], arrow='none', fill=self.color, width=2, capstyle=PROJECTING)
		self.mosfet_p[8] = c.create_line(self.mosfet_p_coords[8], arrow='none', fill=self.color, width=2, capstyle=PROJECTING)

		self.appear = 1
		self.scale(1,1,500)

		# self.mosfet_p[9] = c.create_oval(posx-2*self.flip, posy-2, posx+2*self.flip, posy+2, outline='#FFFFFF',fill='#FFFFFF', width=2)
		# self.mosfet_p[10] = c.create_oval(posx-2*self.flip, posy+15+26+15-2, posx+2*self.flip, posy+15+26+15+2, outline='#FFFFFF',fill='#FFFFFF', width=2)
		# self.mosfet_p[11] = c.create_oval(posx-30*self.flip-15*self.flip-2*self.flip, posy+15+26/2-2, posx-30*self.flip-15*self.flip+2*self.flip, posy+15+26/2+2, outline='#FFFFFF',fill='#FFFFFF', width=2)

	def getTerminals(self):
		return [[self.mosfet_p_coords[0][0], self.mosfet_p_coords[0][1]], [self.mosfet_p_coords[6][2], self.mosfet_p_coords[6][3]], [self.mosfet_p_coords[4][2], self.mosfet_p_coords[4][3]]]

	def reCenter(self):
		self.rotate_around_x_old = self.getCenter()[0]
		self.rotate_around_y_old = self.getCenter()[1]

		for i in range(0, self.items):
			self.mosfet_p_coords[i] = (self.mosfet_p_coords[i][0] + self.posx - self.rotate_around_x_old, self.mosfet_p_coords[i][1] + self.posy - self.rotate_around_y_old, self.mosfet_p_coords[i][2] + self.posx - self.rotate_around_x_old, self.mosfet_p_coords[i][3] + self.posy -self.rotate_around_y_old)

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]

	def getCenter(self):
		for i in range(0, self.items):
			self.xmin_tmp[i] = min(self.mosfet_p_coords[i][0], self.mosfet_p_coords[i][2])
			self.xmax_tmp[i] = max(self.mosfet_p_coords[i][0], self.mosfet_p_coords[i][2])
			self.ymin_tmp[i] = min(self.mosfet_p_coords[i][1], self.mosfet_p_coords[i][3])
			self.ymax_tmp[i] = max(self.mosfet_p_coords[i][1], self.mosfet_p_coords[i][3])
		return [(min(self.xmin_tmp) + max(self.xmax_tmp))/2, (min(self.ymin_tmp) + max(self.ymax_tmp))/2]

	def rotate_kernal(self):

		for i in range(0, self.items):

			self.evl = (math_rotate(self.mosfet_p_coords[i][0], self.mosfet_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.mosfet_p_coords[i][0], self.mosfet_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.mosfet_p_coords[i][2], self.mosfet_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.mosfet_p_coords[i][2], self.mosfet_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.mosfet_p[i], self.evl)

			self.mosfet_p_coords_tmp[i] = (self.evl)

	def rotate_trigger(self):

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		if self.itr_rotate > self.no_of_itr_rotate:
			self.mosfet_p_coords = self.mosfet_p_coords_tmp[:]
			return

		self.rotate_kernal()
		master.after(int(1000./framerate), self.rotate_trigger)

	def rotate(self, angle, animation_time):
	 	self.animation_time = animation_time
	 	self.angle = angle
	 	self.anglex = 0
	 	self.itr_rotate = 0
	 	self.no_of_itr_rotate = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.rotate_trigger()

	def move_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_translate(self.mosfet_p_coords[i][0], self.mosfet_p_coords[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.mosfet_p_coords[i][0], self.mosfet_p_coords[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.mosfet_p_coords[i][2], self.mosfet_p_coords[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.mosfet_p_coords[i][2], self.mosfet_p_coords[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.mosfet_p[i], self.evl)

			self.mosfet_p_coords_tmp[i] = (self.evl)

	def move_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		if self.itr_move > self.no_of_itr_move:
			self.mosfet_p_coords = self.mosfet_p_coords_tmp[:]
			# for i in range(0, self.items):
			# 	self.mosfet_p_coords_tmp[i] = (0,0,0,0)
			return
		self.move_kernal()
		master.after(int(1000./framerate), self.move_trigger)

	def translate(self, newposx, newposy, animation_time):
	 	self.itr_move = 0
	 	self.deltax1 = 0
	 	self.deltay1 = 0
	 	self.no_of_itr_move = self.animation_time * framerate/1000
	 	self.animation_time = animation_time
	 	self.deltax = newposx - self.getCenter()[0]
	 	self.deltay = newposy - self.getCenter()[1]
	 	self.move_trigger()


	def scrolate_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_scale(self.mosfet_p_coords[i][0],self.mosfet_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.mosfet_p_coords[i][0],self.mosfet_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
						math_scale(self.mosfet_p_coords[i][2],self.mosfet_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.mosfet_p_coords[i][2],self.mosfet_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

			c.coords(self.mosfet_p[i], self.evl)

			self.mosfet_p_coords_tmp[i] = self.evl

			self.evl = (math_rotate(self.mosfet_p_coords_tmp[i][0], self.mosfet_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.mosfet_p_coords_tmp[i][0], self.mosfet_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.mosfet_p_coords_tmp[i][2], self.mosfet_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.mosfet_p_coords_tmp[i][2], self.mosfet_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.mosfet_p[i], self.evl)

			self.mosfet_p_coords_tmp2[i] = self.evl

			self.evl = (math_translate(self.mosfet_p_coords_tmp2[i][0], self.mosfet_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.mosfet_p_coords_tmp2[i][0], self.mosfet_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.mosfet_p_coords_tmp2[i][2], self.mosfet_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.mosfet_p_coords_tmp2[i][2], self.mosfet_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.mosfet_p[i], self.evl)

			self.mosfet_p_coords_tmp3[i] = (self.evl)

			
	def scrolate_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)


		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)


		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.mosfet_p_coords = self.mosfet_p_coords_tmp3[:]
			return	
		
		self.scrolate_kernal()
		master.after(int(1000./framerate), self.scrolate_trigger)

	def scrolate(self, newposx, newposy, angle, scalex, scaley, animation_time):
		self.itr_move = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.animation_time = animation_time
		self.deltax = newposx - self.getCenter()[0]
		self.deltay = newposx - self.getCenter()[0]
		self.no_of_itr_move = self.animation_time * framerate/1000
		self.no_of_itr_rotate = self.animation_time*framerate/1000
		self.no_of_itr_scale = self.animation_time*framerate/1000

		self.angle = angle
		self.anglex = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.scaley = scaley
		self.scalex = scalex

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]
		self.scrolate_trigger()

	def scale_kernal(self):

		if self.appear:

			for i in range(0, self.items):

				self.evl =  (math_appear(self.mosfet_p_coords[i][0],self.mosfet_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.mosfet_p_coords[i][0],self.mosfet_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_appear(self.mosfet_p_coords[i][2],self.mosfet_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.mosfet_p_coords[i][2],self.mosfet_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.mosfet_p[i],self.evl)

				self.mosfet_p_coords_tmp[i] = self.evl
		else:
			for i in range(0, self.items):

				self.evl =  (math_scale(self.mosfet_p_coords[i][0],self.mosfet_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.mosfet_p_coords[i][0],self.mosfet_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_scale(self.mosfet_p_coords[i][2],self.mosfet_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.mosfet_p_coords[i][2],self.mosfet_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.mosfet_p[i],self.evl)

				self.mosfet_p_coords_tmp[i] = self.evl

	def scale_trigger(self):

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.mosfet_p_coords = self.mosfet_p_coords_tmp[:]
			self.appear = 0
			return

		self.scale_kernal()
		master.after(int(1000./framerate), self.scale_trigger)

	def scale(self, scalex, scaley, animation_time):
	 	self.animation_time = animation_time
	 	self.itr_scale = 0
	 	self.scaley = scaley
	 	self.scalex = scalex
	 	self.no_of_itr_scale = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.scale_trigger()


class PMOS:

	def init(self, posx, posy, flip, color='#FFFFFF'):
		
		self.posx = posx
		self.posy = posy

		self.color = color

		if flip==0:
			self.flip = 1
		if flip==1:
			self.flip = -1

		self.items = 9
		self.mosfet_p = [0]*self.items
		self.mosfet_p_coords = [0]*self.items
		self.mosfet_p_coords_tmp = [0]*self.items
		self.mosfet_p_coords_tmp2 = [0]*self.items
		self.mosfet_p_coords_tmp3 = [0]*self.items
		self.xmin_tmp = [0]*self.items
		self.xmax_tmp = [0]*self.items
		self.ymin_tmp = [0]*self.items
		self.ymax_tmp = [0]*self.items

		self.itr_move = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.no_of_itr_move = 0
		self.no_of_itr_rotate = 0
		self.no_of_itr_scale = 0
		self.anglex = 0
		self.angle = 0
		self.deltax = 0
		self.deltay = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.scalex_int = 0
		self.scaley_int = 0
		self.deltay1 = 0
		self.scalex = 0
		self.scaley = 0

		self.mosfet_p_coords[0] = (posx, posy, posx, posy+15)
		self.mosfet_p_coords[1] = (posx, posy+15, posx-26*self.flip, posy+15)
		self.mosfet_p_coords[2] = (posx-26*self.flip, posy+15, posx-26*self.flip, posy+15+26)
		self.mosfet_p_coords[3] = (posx-26*self.flip, posy+15+26, posx, posy+15+26)
		self.mosfet_p_coords[4] = (posx, posy+15+26, posx, posy+15+26+15)
		self.mosfet_p_coords[5] = (posx-30*self.flip, posy+15, posx-30*self.flip, posy+15+26)
		self.mosfet_p_coords[6] = (posx-30*self.flip, posy+15+26/2, posx-30*self.flip-15*self.flip, posy+15+26/2)
		self.mosfet_p_coords[7] = (posx-10*self.flip, posy+15, posx-3*self.flip, posy+15-7)
		self.mosfet_p_coords[8] = (posx-10*self.flip, posy+15, posx-3*self.flip, posy+15+7)

		self.reCenter()

		self.mosfet_p[0]= c.create_line(self.mosfet_p_coords[0], arrow='none', fill=self.color, width=2, capstyle=PROJECTING)
		self.mosfet_p[1]= c.create_line(self.mosfet_p_coords[1], arrow='none', fill=self.color, width=2, capstyle=PROJECTING)
		self.mosfet_p[2]= c.create_line(self.mosfet_p_coords[2], arrow='none', fill=self.color, width=2, capstyle=PROJECTING)
		self.mosfet_p[3]= c.create_line(self.mosfet_p_coords[3], arrow='none', fill=self.color, width=2, capstyle=PROJECTING)
		self.mosfet_p[4]= c.create_line(self.mosfet_p_coords[4], arrow='none', fill=self.color, width=2, capstyle=PROJECTING)
		self.mosfet_p[5]= c.create_line(self.mosfet_p_coords[5], arrow='none', fill=self.color, width=2, capstyle=PROJECTING)
		self.mosfet_p[6]= c.create_line(self.mosfet_p_coords[6], arrow='none', fill=self.color, width=2, capstyle=PROJECTING)
		self.mosfet_p[7] = c.create_line(self.mosfet_p_coords[7], arrow='none', fill=self.color, width=2, capstyle=PROJECTING)
		self.mosfet_p[8] = c.create_line(self.mosfet_p_coords[8], arrow='none', fill=self.color, width=2, capstyle=PROJECTING)

		self.appear = 1
		self.scale(1,1,500)

		# self.mosfet_p[9] = c.create_oval(posx-2*self.flip, posy-2, posx+2*self.flip, posy+2, outline='#FFFFFF',fill='#FFFFFF', width=2)
		# self.mosfet_p[10] = c.create_oval(posx-2*self.flip, posy+15+26+15-2, posx+2*self.flip, posy+15+26+15+2, outline='#FFFFFF',fill='#FFFFFF', width=2)
		# self.mosfet_p[11] = c.create_oval(posx-30*self.flip-15*self.flip-2*self.flip, posy+15+26/2-2, posx-30*self.flip-15*self.flip+2*self.flip, posy+15+26/2+2, outline='#FFFFFF',fill='#FFFFFF', width=2)

	def getTerminals(self):
		return [[self.mosfet_p_coords[0][0], self.mosfet_p_coords[0][1]], [self.mosfet_p_coords[6][2], self.mosfet_p_coords[6][3]], [self.mosfet_p_coords[4][2], self.mosfet_p_coords[4][3]]]

	def reCenter(self):
		self.rotate_around_x_old = self.getCenter()[0]
		self.rotate_around_y_old = self.getCenter()[1]

		for i in range(0, self.items):
			self.mosfet_p_coords[i] = (self.mosfet_p_coords[i][0] + self.posx - self.rotate_around_x_old, self.mosfet_p_coords[i][1] + self.posy - self.rotate_around_y_old, self.mosfet_p_coords[i][2] + self.posx - self.rotate_around_x_old, self.mosfet_p_coords[i][3] + self.posy -self.rotate_around_y_old)

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]

	def getCenter(self):
		for i in range(0, self.items):
			self.xmin_tmp[i] = min(self.mosfet_p_coords[i][0], self.mosfet_p_coords[i][2])
			self.xmax_tmp[i] = max(self.mosfet_p_coords[i][0], self.mosfet_p_coords[i][2])
			self.ymin_tmp[i] = min(self.mosfet_p_coords[i][1], self.mosfet_p_coords[i][3])
			self.ymax_tmp[i] = max(self.mosfet_p_coords[i][1], self.mosfet_p_coords[i][3])
		return [(min(self.xmin_tmp) + max(self.xmax_tmp))/2, (min(self.ymin_tmp) + max(self.ymax_tmp))/2]

	def rotate_kernal(self):

		for i in range(0, self.items):

			self.evl = (math_rotate(self.mosfet_p_coords[i][0], self.mosfet_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.mosfet_p_coords[i][0], self.mosfet_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.mosfet_p_coords[i][2], self.mosfet_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.mosfet_p_coords[i][2], self.mosfet_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.mosfet_p[i], self.evl)

			self.mosfet_p_coords_tmp[i] = (self.evl)

	def rotate_trigger(self):

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		if self.itr_rotate > self.no_of_itr_rotate:
			self.mosfet_p_coords = self.mosfet_p_coords_tmp[:]
			return

		self.rotate_kernal()
		master.after(int(1000./framerate), self.rotate_trigger)

	def rotate(self, angle, animation_time):
	 	self.animation_time = animation_time
	 	self.angle = angle
	 	self.anglex = 0
	 	self.itr_rotate = 0
	 	self.no_of_itr_rotate = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.rotate_trigger()

	def move_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_translate(self.mosfet_p_coords[i][0], self.mosfet_p_coords[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.mosfet_p_coords[i][0], self.mosfet_p_coords[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.mosfet_p_coords[i][2], self.mosfet_p_coords[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.mosfet_p_coords[i][2], self.mosfet_p_coords[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.mosfet_p[i], self.evl)

			self.mosfet_p_coords_tmp[i] = (self.evl)

	def move_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		if self.itr_move > self.no_of_itr_move:
			self.mosfet_p_coords = self.mosfet_p_coords_tmp[:]
			return
		self.move_kernal()
		master.after(int(1000./framerate), self.move_trigger)

	def translate(self, newposx, newposy, animation_time):
	 	self.itr_move = 0
	 	self.deltax1 = 0
	 	self.deltay1 = 0
	 	self.no_of_itr_move = self.animation_time * framerate/1000
	 	self.animation_time = animation_time
	 	self.deltax = newposx - self.getCenter()[0]
	 	self.deltay = newposy - self.getCenter()[1]
	 	self.move_trigger()


	def scrolate_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_scale(self.mosfet_p_coords[i][0],self.mosfet_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.mosfet_p_coords[i][0],self.mosfet_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
						math_scale(self.mosfet_p_coords[i][2],self.mosfet_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.mosfet_p_coords[i][2],self.mosfet_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

			c.coords(self.mosfet_p[i], self.evl)

			self.mosfet_p_coords_tmp[i] = self.evl

			self.evl = (math_rotate(self.mosfet_p_coords_tmp[i][0], self.mosfet_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.mosfet_p_coords_tmp[i][0], self.mosfet_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.mosfet_p_coords_tmp[i][2], self.mosfet_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.mosfet_p_coords_tmp[i][2], self.mosfet_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.mosfet_p[i], self.evl)

			self.mosfet_p_coords_tmp2[i] = self.evl

			self.evl = (math_translate(self.mosfet_p_coords_tmp2[i][0], self.mosfet_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.mosfet_p_coords_tmp2[i][0], self.mosfet_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.mosfet_p_coords_tmp2[i][2], self.mosfet_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.mosfet_p_coords_tmp2[i][2], self.mosfet_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.mosfet_p[i], self.evl)

			self.mosfet_p_coords_tmp3[i] = (self.evl)

			
	def scrolate_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)


		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		# if self.itr_rotate > self.no_of_itr_rotate:
		# 	self.mosfet_p_coords = self.mosfet_p_coords_tmp2[:]
		# 	return

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.mosfet_p_coords = self.mosfet_p_coords_tmp3[:]
			# for i in range(0, self.items):
			# 	self.mosfet_p_coords_tmp[i] = (0,0,0,0)
			return	
		
		self.scrolate_kernal()
		master.after(int(1000./framerate), self.scrolate_trigger)

	def scrolate(self, newposx, newposy, angle, scalex, scaley, animation_time):
		self.itr_move = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.animation_time = animation_time
		self.deltax = newposx - self.getCenter()[0]
		self.deltay = newposx - self.getCenter()[0]
		self.no_of_itr_move = self.animation_time * framerate/1000
		self.no_of_itr_scale = self.animation_time * framerate/1000
		self.no_of_itr_rotate = self.animation_time * framerate/1000

		self.angle = angle
		self.anglex = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.scaley = scaley
		self.scalex = scalex

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]
		self.scrolate_trigger()

	def scale_kernal(self):

		if self.appear:
			for i in range(0, self.items):

				self.evl =  (math_appear(self.mosfet_p_coords[i][0],self.mosfet_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.mosfet_p_coords[i][0],self.mosfet_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_appear(self.mosfet_p_coords[i][2],self.mosfet_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.mosfet_p_coords[i][2],self.mosfet_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.mosfet_p[i],self.evl)

				self.mosfet_p_coords_tmp[i] = self.evl
		else:
			for i in range(0, self.items):

				self.evl =  (math_scale(self.mosfet_p_coords[i][0],self.mosfet_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.mosfet_p_coords[i][0],self.mosfet_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_scale(self.mosfet_p_coords[i][2],self.mosfet_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.mosfet_p_coords[i][2],self.mosfet_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.mosfet_p[i],self.evl)

				self.mosfet_p_coords_tmp[i] = self.evl

	def scale_trigger(self):

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.mosfet_p_coords = self.mosfet_p_coords_tmp[:]
			self.appear = 0
			# for i in range(0, self.items):
			# 	self.mosfet_p_coords_tmp[i] = (0,0,0,0)
			return

		self.scale_kernal()
		master.after(int(1000./framerate), self.scale_trigger)

	def scale(self, scalex, scaley, animation_time):
	 	self.animation_time = animation_time
	 	self.itr_scale = 0
	 	self.scaley = scaley
	 	self.scalex = scalex
	 	self.no_of_itr_scale = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.scale_trigger()


class WIRE:
	def init(self, startx, starty, endx, endy, width=2, color='#FFFFFF'):
		self.startx = startx
		self.starty = starty
		self.endx = endx
		self.endy = endy
		self.color = color

		self.posx = startx
		self.posy = starty


		self.items = 1
		self.wire_p = [0]*self.items
		self.wire_p_coords = [0]*self.items
		self.wire_p_coords_tmp = [0]*self.items
		self.wire_p_coords_tmp2 = [0]*self.items
		self.wire_p_coords_tmp3 = [0]*self.items
		self.xmin_tmp = [0]*self.items
		self.xmax_tmp = [0]*self.items
		self.ymin_tmp = [0]*self.items
		self.ymax_tmp = [0]*self.items

		self.itr_move = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.no_of_itr_move = 0
		self.no_of_itr_rotate = 0
		self.no_of_itr_scale = 0
		self.anglex = 0
		self.angle = 0
		self.deltax = 0
		self.deltay = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.scalex_int = 0
		self.scaley_int = 0
		self.deltay1 = 0
		self.scalex = 0
		self.scaley = 0

		self.wire_p_coords[0] = (self.startx, self.starty, self.endx, self.endy)

		self.wire_p[0] = c.create_line(self.wire_p_coords[0], arrow='none', fill=self.color, width=width, joinstyle='bevel')

		self.appear = 1
		self.scale(1,1,500)

	def reCenter(self):
		self.rotate_around_x_old = self.getCenter()[0]
		self.rotate_around_y_old = self.getCenter()[1]

		for i in range(0, self.items):
			self.wire_p_coords[i] = (self.wire_p_coords[i][0] + self.posx - self.rotate_around_x_old, self.wire_p_coords[i][1] + self.posy - self.rotate_around_y_old, self.wire_p_coords[i][2] + self.posx - self.rotate_around_x_old, self.wire_p_coords[i][3] + self.posy -self.rotate_around_y_old)

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]

	def getCenter(self):
		for i in range(0, self.items):
			self.xmin_tmp[i] = min(self.wire_p_coords[i][0], self.wire_p_coords[i][2])
			self.xmax_tmp[i] = max(self.wire_p_coords[i][0], self.wire_p_coords[i][2])
			self.ymin_tmp[i] = min(self.wire_p_coords[i][1], self.wire_p_coords[i][3])
			self.ymax_tmp[i] = max(self.wire_p_coords[i][1], self.wire_p_coords[i][3])
		return [(min(self.xmin_tmp) + max(self.xmax_tmp))/2, (min(self.ymin_tmp) + max(self.ymax_tmp))/2]

	def rotate_kernal(self):

		for i in range(0, self.items):

			self.evl = (math_rotate(self.wire_p_coords[i][0], self.wire_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.wire_p_coords[i][0], self.wire_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.wire_p_coords[i][2], self.wire_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.wire_p_coords[i][2], self.wire_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.wire_p[i], self.evl)

			self.wire_p_coords_tmp[i] = (self.evl)

	def rotate_trigger(self):

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		if self.itr_rotate > self.no_of_itr_rotate:
			self.wire_p_coords = self.wire_p_coords_tmp[:]
			return

		self.rotate_kernal()
		master.after(int(1000./framerate), self.rotate_trigger)

	def rotate(self, angle, animation_time):
	 	self.animation_time = animation_time
	 	self.angle = angle
	 	self.anglex = 0
	 	self.itr_rotate = 0
	 	self.no_of_itr_rotate = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.rotate_trigger()

	def move_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_translate(self.wire_p_coords[i][0], self.wire_p_coords[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.wire_p_coords[i][0], self.wire_p_coords[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.wire_p_coords[i][2], self.wire_p_coords[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.wire_p_coords[i][2], self.wire_p_coords[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.wire_p[i], self.evl)

			self.wire_p_coords_tmp[i] = (self.evl)

	def move_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		if self.itr_move > self.no_of_itr_move:
			self.wire_p_coords = self.wire_p_coords_tmp[:]
			return
		self.move_kernal()
		master.after(int(1000./framerate), self.move_trigger)

	def translate(self, newposx, newposy, animation_time):
	 	self.itr_move = 0
	 	self.deltax1 = 0
	 	self.deltay1 = 0
	 	self.no_of_itr_move = self.animation_time * framerate/1000
	 	self.animation_time = animation_time
	 	self.deltax = newposx - self.getCenter()[0]
	 	self.deltay = newposy - self.getCenter()[1]
	 	self.move_trigger()


	def scrolate_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_scale(self.wire_p_coords[i][0],self.wire_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.wire_p_coords[i][0],self.wire_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
						math_scale(self.wire_p_coords[i][2],self.wire_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.wire_p_coords[i][2],self.wire_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

			c.coords(self.wire_p[i], self.evl)

			self.wire_p_coords_tmp[i] = self.evl

			self.evl = (math_rotate(self.wire_p_coords_tmp[i][0], self.wire_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.wire_p_coords_tmp[i][0], self.wire_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.wire_p_coords_tmp[i][2], self.wire_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.wire_p_coords_tmp[i][2], self.wire_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.wire_p[i], self.evl)

			self.wire_p_coords_tmp2[i] = self.evl

			self.evl = (math_translate(self.wire_p_coords_tmp2[i][0], self.wire_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.wire_p_coords_tmp2[i][0], self.wire_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.wire_p_coords_tmp2[i][2], self.wire_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.wire_p_coords_tmp2[i][2], self.wire_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.wire_p[i], self.evl)

			self.wire_p_coords_tmp3[i] = (self.evl)

			
	def scrolate_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.wire_p_coords = self.wire_p_coords_tmp3[:]
			return	
		
		self.scrolate_kernal()
		master.after(int(1000./framerate), self.scrolate_trigger)

	def scrolate(self, newposx, newposy, angle, scalex, scaley, animation_time):
		self.itr_move = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.animation_time = animation_time
		self.deltax = newposx - self.getCenter()[0]
		self.deltay = newposx - self.getCenter()[0]
		self.no_of_itr_move = self.animation_time * framerate/1000
		self.no_of_itr_scale = self.animation_time * framerate/1000
		self.no_of_itr_rotate = self.animation_time * framerate/1000

		self.angle = angle
		self.anglex = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.scaley = scaley
		self.scalex = scalex

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]
		self.scrolate_trigger()

	def scale_kernal(self):

		if self.appear:
			for i in range(0, self.items):

				self.evl =  (math_appear(self.wire_p_coords[i][0],self.wire_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.wire_p_coords[i][0],self.wire_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_appear(self.wire_p_coords[i][2],self.wire_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.wire_p_coords[i][2],self.wire_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.wire_p[i],self.evl)

				self.wire_p_coords_tmp[i] = self.evl
		else:
			for i in range(0, self.items):

				self.evl = (math_scale(self.wire_p_coords[i][0],self.wire_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.wire_p_coords[i][0],self.wire_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_scale(self.wire_p_coords[i][2],self.wire_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.wire_p_coords[i][2],self.wire_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.wire_p[i],self.evl)

				self.wire_p_coords_tmp[i] = self.evl			

	def scale_trigger(self):

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.wire_p_coords = self.wire_p_coords_tmp[:]
			return

		self.scale_kernal()
		master.after(int(1000./framerate), self.scale_trigger)

	def scale(self, scalex, scaley, animation_time):
	 	self.animation_time = animation_time
	 	self.itr_scale = 0
	 	self.scaley = scaley
	 	self.scalex = scalex
	 	self.no_of_itr_scale = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.scale_trigger()
		

class IDC:
	def init(self, posx, posy):
		self.posx = posx;
		self.posy = posy;
		self.items = 7
		self.idc_p = [0]*self.items
		self.idc_p_coords = [0]*self.items
		self.idc_p_coords_tmp = [0]*self.items
		self.idc_p_coords_tmp2 = [0]*self.items
		self.idc_p_coords_tmp3 = [0]*self.items
		self.xmin_tmp = [0]*self.items
		self.xmax_tmp = [0]*self.items
		self.ymin_tmp = [0]*self.items
		self.ymax_tmp = [0]*self.items

		self.itr_move = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.no_of_itr_move = 0
		self.no_of_itr_rotate = 0
		self.no_of_itr_scale = 0
		self.anglex = 0
		self.angle = 0
		self.deltax = 0
		self.deltay = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.scalex_int = 0
		self.scaley_int = 0
		self.deltay1 = 0
		self.scalex = 0
		self.scaley = 0

		self.idc_p_coords[0] = (posx, posy, posx, posy+20)
		self.idc_p_coords[1] = (posx, posy+60, posx, posy+80)
		self.idc_p_coords[2] = (posx-20, posy+20, posx+20, posy+60)
		self.idc_p_coords[3] = (posx-19, posy+21, posx+19, posy+59)
		self.idc_p_coords[4] = (posx, posy+20+10, posx, posy+20+30)
		self.idc_p_coords[5] = (posx, posy+20+30, posx-8, posy+20+30-8)
		self.idc_p_coords[6] = (posx, posy+20+30, posx+8, posy+20+30-8)

		self.reCenter()

		self.idc_p[0] = c.create_line(self.idc_p_coords[0], arrow='none', fill='#FFFFFF', width=2)
		self.idc_p[1] = c.create_line(self.idc_p_coords[1], arrow='none', fill='#FFFFFF', width=2)
		self.idc_p[2] = c.create_oval(self.idc_p_coords[2], outline='#FFFFFF',fill='#FFFFFF', width=2)
		self.idc_p[3] = c.create_oval(self.idc_p_coords[3], outline='#FFFFFF',fill='#000000', width=2)
		self.idc_p[4] = c.create_line(self.idc_p_coords[4], arrow='none', fill='#FFFFFF', width=2)
		self.idc_p[5] = c.create_line(self.idc_p_coords[5], arrow='none', fill='#FFFFFF', width=2)
		self.idc_p[6] = c.create_line(self.idc_p_coords[6], arrow='none', fill='#FFFFFF', width=2)

		self.appear = 1
		self.scale(1,1,500)


	def getTerminals(self):
		return [[self.idc_p_coords[0][0], self.idc_p_coords[0][1]], [self.idc_p_coords[1][2], self.idc_p_coords[1][3]]]

	def reCenter(self):
		self.rotate_around_x_old = self.getCenter()[0]
		self.rotate_around_y_old = self.getCenter()[1]

		for i in range(0, self.items):
			self.idc_p_coords[i] = (self.idc_p_coords[i][0] + self.posx - self.rotate_around_x_old, self.idc_p_coords[i][1] + self.posy - self.rotate_around_y_old, self.idc_p_coords[i][2] + self.posx - self.rotate_around_x_old, self.idc_p_coords[i][3] + self.posy -self.rotate_around_y_old)

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]

	def getCenter(self):
		for i in range(0, self.items):
			self.xmin_tmp[i] = min(self.idc_p_coords[i][0], self.idc_p_coords[i][2])
			self.xmax_tmp[i] = max(self.idc_p_coords[i][0], self.idc_p_coords[i][2])
			self.ymin_tmp[i] = min(self.idc_p_coords[i][1], self.idc_p_coords[i][3])
			self.ymax_tmp[i] = max(self.idc_p_coords[i][1], self.idc_p_coords[i][3])
		return [(min(self.xmin_tmp) + max(self.xmax_tmp))/2, (min(self.ymin_tmp) + max(self.ymax_tmp))/2]


	def rotate_kernal(self):

		for i in range(0, self.items):

			self.evl = (math_rotate(self.idc_p_coords[i][0], self.idc_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.idc_p_coords[i][0], self.idc_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.idc_p_coords[i][2], self.idc_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.idc_p_coords[i][2], self.idc_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.idc_p[i], self.evl)

			self.idc_p_coords_tmp[i] = (self.evl)

	def rotate_trigger(self):

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		if self.itr_rotate > self.no_of_itr_rotate:
			self.idc_p_coords = self.idc_p_coords_tmp[:]
			return

		self.rotate_kernal()
		master.after(int(1000./framerate), self.rotate_trigger)

	def rotate(self, angle, animation_time):
	 	self.animation_time = animation_time
	 	self.angle = angle
	 	self.anglex = 0
	 	self.itr_rotate = 0
	 	self.no_of_itr_rotate = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.rotate_trigger()

	def move_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_translate(self.idc_p_coords[i][0], self.idc_p_coords[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.idc_p_coords[i][0], self.idc_p_coords[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.idc_p_coords[i][2], self.idc_p_coords[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.idc_p_coords[i][2], self.idc_p_coords[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.idc_p[i], self.evl)

			self.idc_p_coords_tmp[i] = (self.evl)

	def move_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		if self.itr_move > self.no_of_itr_move:
			self.idc_p_coords = self.idc_p_coords_tmp[:]
			return
		self.move_kernal()
		master.after(int(1000./framerate), self.move_trigger)

	def translate(self, newposx, newposy, animation_time):
	 	self.itr_move = 0
	 	self.deltax1 = 0
	 	self.deltay1 = 0
	 	self.animation_time = animation_time
	 	self.no_of_itr_move = self.animation_time * framerate/1000
	 	self.deltax = newposx - self.getCenter()[0]
	 	self.deltay = newposy - self.getCenter()[1]
	 	self.move_trigger()


	def scrolate_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_scale(self.idc_p_coords[i][0],self.idc_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.idc_p_coords[i][0],self.idc_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
						math_scale(self.idc_p_coords[i][2],self.idc_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.idc_p_coords[i][2],self.idc_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

			c.coords(self.idc_p[i], self.evl)

			self.idc_p_coords_tmp[i] = self.evl

			self.evl = (math_rotate(self.idc_p_coords_tmp[i][0], self.idc_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.idc_p_coords_tmp[i][0], self.idc_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.idc_p_coords_tmp[i][2], self.idc_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.idc_p_coords_tmp[i][2], self.idc_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.idc_p[i], self.evl)

			self.idc_p_coords_tmp2[i] = self.evl

			self.evl = (math_translate(self.idc_p_coords_tmp2[i][0], self.idc_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.idc_p_coords_tmp2[i][0], self.idc_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.idc_p_coords_tmp2[i][2], self.idc_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.idc_p_coords_tmp2[i][2], self.idc_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.idc_p[i], self.evl)

			self.idc_p_coords_tmp3[i] = (self.evl)

			
	def scrolate_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		if self.itr_rotate == 0:
			self.no_of_itr_rotate = self.animation_time*framerate/1000

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.idc_p_coords = self.idc_p_coords_tmp3[:]
			return	
		
		self.scrolate_kernal()
		master.after(int(1000./framerate), self.scrolate_trigger)

	def scrolate(self, newposx, newposy, angle, scalex, scaley, animation_time):
		self.itr_move = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.animation_time = animation_time
		self.deltax = newposx - self.getCenter()[0]
		self.deltay = newposx - self.getCenter()[0]
		self.no_of_itr_move = self.animation_time * framerate/1000
		self.no_of_itr_scale = self.animation_time * framerate/1000
		self.no_of_itr_rotate = self.animation_time * framerate/1000

		self.angle = angle
		self.anglex = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.scaley = scaley
		self.scalex = scalex

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]
		self.scrolate_trigger()

	def scale_kernal(self):

		if self.appear:
			for i in range(0, self.items):

				self.evl =  (math_appear(self.idc_p_coords[i][0],self.idc_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.idc_p_coords[i][0],self.idc_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_appear(self.idc_p_coords[i][2],self.idc_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.idc_p_coords[i][2],self.idc_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.idc_p[i],self.evl)

				self.idc_p_coords_tmp[i] = self.evl
		else:
			for i in range(0, self.items):

				self.evl =  (math_scale(self.idc_p_coords[i][0],self.idc_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.idc_p_coords[i][0],self.idc_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_scale(self.idc_p_coords[i][2],self.idc_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.idc_p_coords[i][2],self.idc_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.idc_p[i],self.evl)

				self.idc_p_coords_tmp[i] = self.evl

	def scale_trigger(self):

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.idc_p_coords = self.idc_p_coords_tmp[:]
			return

		self.scale_kernal()
		master.after(int(1000./framerate), self.scale_trigger)

	def scale(self, scalex, scaley, animation_time):
	 	self.animation_time = animation_time
	 	self.itr_scale = 0
	 	self.scaley = scaley
	 	self.scalex = scalex
	 	self.no_of_itr_scale = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.scale_trigger()
		

class VDC:
	def init(self, posx, posy):
		self.posx = posx;
		self.posy = posy;
		self.items = 7
		self.vdc_p = [0]*self.items
		self.vdc_p_coords = [0]*self.items
		self.vdc_p_coords_tmp = [0]*self.items
		self.vdc_p_coords_tmp2 = [0]*self.items
		self.vdc_p_coords_tmp3 = [0]*self.items
		self.xmin_tmp = [0]*self.items
		self.xmax_tmp = [0]*self.items
		self.ymin_tmp = [0]*self.items
		self.ymax_tmp = [0]*self.items

		self.itr_move = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.no_of_itr_move = 0
		self.no_of_itr_rotate = 0
		self.no_of_itr_scale = 0
		self.anglex = 0
		self.angle = 0
		self.deltax = 0
		self.deltay = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.scalex_int = 0
		self.scaley_int = 0
		self.deltay1 = 0
		self.scalex = 0
		self.scaley = 0

		self.vdc_p_coords[0] = (posx, posy, posx, posy+20)
		self.vdc_p_coords[1] = (posx, posy+60, posx, posy+80)
		self.vdc_p_coords[2] = (posx-20, posy+20, posx+20, posy+60)
		self.vdc_p_coords[3] = (posx-19, posy+21, posx+19, posy+59)
		self.vdc_p_coords[4] = (posx, posy+20+10, posx, posy+20+18)
		self.vdc_p_coords[5] = (posx-4, posy+20+10+4, posx+4, posy+20+18-4)
		self.vdc_p_coords[6] = (posx-4, posy+20+30, posx+4, posy+20+30)

		self.reCenter()

		self.vdc_p[0] = c.create_line(self.vdc_p_coords[0], arrow='none', fill='#FFFFFF', width=2)
		self.vdc_p[1] = c.create_line(self.vdc_p_coords[1], arrow='none', fill='#FFFFFF', width=2)
		self.vdc_p[2] = c.create_oval(self.vdc_p_coords[2], outline='#FFFFFF',fill='#FFFFFF', width=2)
		self.vdc_p[3] = c.create_oval(self.vdc_p_coords[3], outline='#FFFFFF',fill='#000000', width=2)
		self.vdc_p[4] = c.create_line(self.vdc_p_coords[4], arrow='none', fill='#FFFFFF', width=2)
		self.vdc_p[5] = c.create_line(self.vdc_p_coords[5], arrow='none', fill='#FFFFFF', width=2)
		self.vdc_p[6] = c.create_line(self.vdc_p_coords[6], arrow='none', fill='#FFFFFF', width=2)

		self.appear = 1
		self.scale(1,1,500)


	def getTerminals(self):
		return [[self.vdc_p_coords[0][0], self.vdc_p_coords[0][1]], [self.vdc_p_coords[1][2], self.vdc_p_coords[1][3]]]

	def reCenter(self):
		self.rotate_around_x_old = self.getCenter()[0]
		self.rotate_around_y_old = self.getCenter()[1]

		for i in range(0, self.items):
			self.vdc_p_coords[i] = (self.vdc_p_coords[i][0] + self.posx - self.rotate_around_x_old, self.vdc_p_coords[i][1] + self.posy - self.rotate_around_y_old, self.vdc_p_coords[i][2] + self.posx - self.rotate_around_x_old, self.vdc_p_coords[i][3] + self.posy -self.rotate_around_y_old)

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]

	def getCenter(self):
		for i in range(0, self.items):
			self.xmin_tmp[i] = min(self.vdc_p_coords[i][0], self.vdc_p_coords[i][2])
			self.xmax_tmp[i] = max(self.vdc_p_coords[i][0], self.vdc_p_coords[i][2])
			self.ymin_tmp[i] = min(self.vdc_p_coords[i][1], self.vdc_p_coords[i][3])
			self.ymax_tmp[i] = max(self.vdc_p_coords[i][1], self.vdc_p_coords[i][3])
		return [(min(self.xmin_tmp) + max(self.xmax_tmp))/2, (min(self.ymin_tmp) + max(self.ymax_tmp))/2]


	def rotate_kernal(self):

		for i in range(0, self.items):

			self.evl = (math_rotate(self.vdc_p_coords[i][0], self.vdc_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.vdc_p_coords[i][0], self.vdc_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.vdc_p_coords[i][2], self.vdc_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.vdc_p_coords[i][2], self.vdc_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.vdc_p[i], self.evl)

			self.vdc_p_coords_tmp[i] = (self.evl)

	def rotate_trigger(self):

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		if self.itr_rotate > self.no_of_itr_rotate:
			self.vdc_p_coords = self.vdc_p_coords_tmp[:]
			return

		self.rotate_kernal()
		master.after(int(1000./framerate), self.rotate_trigger)

	def rotate(self, angle, animation_time):
	 	self.animation_time = animation_time
	 	self.angle = angle
	 	self.anglex = 0
	 	self.itr_rotate = 0
	 	self.no_of_itr_rotate = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.rotate_trigger()

	def move_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_translate(self.vdc_p_coords[i][0], self.vdc_p_coords[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.vdc_p_coords[i][0], self.vdc_p_coords[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.vdc_p_coords[i][2], self.vdc_p_coords[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.vdc_p_coords[i][2], self.vdc_p_coords[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.vdc_p[i], self.evl)

			self.vdc_p_coords_tmp[i] = (self.evl)

	def move_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		if self.itr_move > self.no_of_itr_move:
			self.vdc_p_coords = self.vdc_p_coords_tmp[:]
			return
		self.move_kernal()
		master.after(int(1000./framerate), self.move_trigger)

	def translate(self, newposx, newposy, animation_time):
	 	self.itr_move = 0
	 	self.deltax1 = 0
	 	self.deltay1 = 0
	 	self.no_of_itr_move = self.animation_time * framerate/1000
	 	self.animation_time = animation_time
	 	self.deltax = newposx - self.getCenter()[0]
	 	self.deltay = newposy - self.getCenter()[1]
	 	self.move_trigger()


	def scrolate_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_scale(self.vdc_p_coords[i][0],self.vdc_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.vdc_p_coords[i][0],self.vdc_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
						math_scale(self.vdc_p_coords[i][2],self.vdc_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.vdc_p_coords[i][2],self.vdc_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

			c.coords(self.vdc_p[i], self.evl)

			self.vdc_p_coords_tmp[i] = self.evl

			self.evl = (math_rotate(self.vdc_p_coords_tmp[i][0], self.vdc_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.vdc_p_coords_tmp[i][0], self.vdc_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.vdc_p_coords_tmp[i][2], self.vdc_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.vdc_p_coords_tmp[i][2], self.vdc_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.vdc_p[i], self.evl)

			self.vdc_p_coords_tmp2[i] = self.evl

			self.evl = (math_translate(self.vdc_p_coords_tmp2[i][0], self.vdc_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.vdc_p_coords_tmp2[i][0], self.vdc_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.vdc_p_coords_tmp2[i][2], self.vdc_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.vdc_p_coords_tmp2[i][2], self.vdc_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.vdc_p[i], self.evl)

			self.vdc_p_coords_tmp3[i] = (self.evl)

			
	def scrolate_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)


		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)


		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.vdc_p_coords = self.vdc_p_coords_tmp3[:]
			return	
		
		self.scrolate_kernal()
		master.after(int(1000./framerate), self.scrolate_trigger)

	def scrolate(self, newposx, newposy, angle, scalex, scaley, animation_time):
		self.itr_move = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.animation_time = animation_time
		self.deltax = newposx - self.getCenter()[0]
		self.deltay = newposx - self.getCenter()[0]
		self.no_of_itr_move = self.animation_time * framerate/1000
		self.no_of_itr_rotate = self.animation_time * framerate/1000
		self.no_of_itr_scale = self.animation_time * framerate/1000

		self.angle = angle
		self.anglex = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.scaley = scaley
		self.scalex = scalex

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]
		self.scrolate_trigger()

	def scale_kernal(self):

		if self.appear:
			for i in range(0, self.items):

				self.evl =  (math_scale(self.vdc_p_coords[i][0],self.vdc_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.vdc_p_coords[i][0],self.vdc_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_scale(self.vdc_p_coords[i][2],self.vdc_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.vdc_p_coords[i][2],self.vdc_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.vdc_p[i],self.evl)

				self.vdc_p_coords_tmp[i] = self.evl
		else:
			for i in range(0, self.items):

				self.evl =  (math_scale(self.vdc_p_coords[i][0],self.vdc_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.vdc_p_coords[i][0],self.vdc_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_scale(self.vdc_p_coords[i][2],self.vdc_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.vdc_p_coords[i][2],self.vdc_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.vdc_p[i],self.evl)

				self.vdc_p_coords_tmp[i] = self.evl

	def scale_trigger(self):

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.vdc_p_coords = self.vdc_p_coords_tmp[:]
			self.appear = 0
			return

		self.scale_kernal()
		master.after(int(1000./framerate), self.scale_trigger)

	def scale(self, scalex, scaley, animation_time):
	 	self.animation_time = animation_time
	 	self.itr_scale = 0
	 	self.scaley = scaley
	 	self.scalex = scalex
	 	self.no_of_itr_scale = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.scale_trigger()
		
class VCVS:
	def init(self, posx, posy):
		self.posx = posx;
		self.posy = posy;
		self.items = 9
		self.vcvs_p = [0]*self.items
		self.vcvs_p_coords = [0]*self.items
		self.vcvs_p_coords_tmp = [0]*self.items
		self.vcvs_p_coords_tmp2 = [0]*self.items
		self.vcvs_p_coords_tmp3 = [0]*self.items
		self.xmin_tmp = [0]*self.items
		self.xmax_tmp = [0]*self.items
		self.ymin_tmp = [0]*self.items
		self.ymax_tmp = [0]*self.items

		self.itr_move = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.no_of_itr_move = 0
		self.no_of_itr_rotate = 0
		self.no_of_itr_scale = 0
		self.anglex = 0
		self.angle = 0
		self.deltax = 0
		self.deltay = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.scalex_int = 0
		self.scaley_int = 0
		self.deltay1 = 0
		self.scalex = 0
		self.scaley = 0

		self.vcvs_p_coords[0] = (posx, posy, posx, posy+20)
		self.vcvs_p_coords[1] = (posx, posy+60, posx, posy+80)
		self.vcvs_p_coords[2] = (posx, posy+20, posx-20, posy+40)
		self.vcvs_p_coords[3] = (posx, posy+20, posx+20, posy+40)
		self.vcvs_p_coords[4] = (posx-20, posy+40, posx, posy+60)
		self.vcvs_p_coords[5] = (posx+20, posy+40, posx, posy+60)
		self.vcvs_p_coords[6] = (posx, posy+28, posx, posy+38)
		self.vcvs_p_coords[7] = (posx-5, posy+28+5, posx+5, posy+38-5)
		self.vcvs_p_coords[8] = (posx-5, posy+30+5+15, posx+5, posy+40-5+15)

		self.reCenter()

		self.vcvs_p[0] = c.create_line(self.vcvs_p_coords[0], arrow='none', fill='#FFFFFF', width=2)
		self.vcvs_p[1] = c.create_line(self.vcvs_p_coords[1], arrow='none', fill='#FFFFFF', width=2)
		self.vcvs_p[2] = c.create_line(self.vcvs_p_coords[2], arrow='none', fill='#FFFFFF', width=2)
		self.vcvs_p[3] = c.create_line(self.vcvs_p_coords[3], arrow='none', fill='#FFFFFF', width=2)
		self.vcvs_p[4] = c.create_line(self.vcvs_p_coords[4], arrow='none', fill='#FFFFFF', width=2)
		self.vcvs_p[5] = c.create_line(self.vcvs_p_coords[5], arrow='none', fill='#FFFFFF', width=2)
		self.vcvs_p[6] = c.create_line(self.vcvs_p_coords[6], arrow='none', fill='#FFFFFF', width=2)
		self.vcvs_p[7] = c.create_line(self.vcvs_p_coords[7], arrow='none', fill='#FFFFFF', width=2)
		self.vcvs_p[8] = c.create_line(self.vcvs_p_coords[8], arrow='none', fill='#FFFFFF', width=2)

		self.appear = 1
		self.scale(1,1,500)


	def getTerminals(self):
		return [[self.vcvs_p_coords[0][0], self.vcvs_p_coords[0][1]], [self.vcvs_p_coords[1][2], self.vcvs_p_coords[1][3]]]

	def reCenter(self):
		self.rotate_around_x_old = self.getCenter()[0]
		self.rotate_around_y_old = self.getCenter()[1]

		for i in range(0, self.items):
			self.vcvs_p_coords[i] = (self.vcvs_p_coords[i][0] + self.posx - self.rotate_around_x_old, self.vcvs_p_coords[i][1] + self.posy - self.rotate_around_y_old, self.vcvs_p_coords[i][2] + self.posx - self.rotate_around_x_old, self.vcvs_p_coords[i][3] + self.posy -self.rotate_around_y_old)

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]

	def getCenter(self):
		for i in range(0, self.items):
			self.xmin_tmp[i] = min(self.vcvs_p_coords[i][0], self.vcvs_p_coords[i][2])
			self.xmax_tmp[i] = max(self.vcvs_p_coords[i][0], self.vcvs_p_coords[i][2])
			self.ymin_tmp[i] = min(self.vcvs_p_coords[i][1], self.vcvs_p_coords[i][3])
			self.ymax_tmp[i] = max(self.vcvs_p_coords[i][1], self.vcvs_p_coords[i][3])
		return [(min(self.xmin_tmp) + max(self.xmax_tmp))/2, (min(self.ymin_tmp) + max(self.ymax_tmp))/2]


	def rotate_kernal(self):

		for i in range(0, self.items):

			self.evl = (math_rotate(self.vcvs_p_coords[i][0], self.vcvs_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.vcvs_p_coords[i][0], self.vcvs_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.vcvs_p_coords[i][2], self.vcvs_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.vcvs_p_coords[i][2], self.vcvs_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.vcvs_p[i], self.evl)

			self.vcvs_p_coords_tmp[i] = (self.evl)

	def rotate_trigger(self):

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		if self.itr_rotate > self.no_of_itr_rotate:
			self.vcvs_p_coords = self.vcvs_p_coords_tmp[:]
			return

		self.rotate_kernal()
		master.after(int(1000./framerate), self.rotate_trigger)

	def rotate(self, angle, animation_time):
	 	self.animation_time = animation_time
	 	self.angle = angle
	 	self.anglex = 0
	 	self.itr_rotate = 0
	 	self.no_of_itr_rotate = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.rotate_trigger()

	def move_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_translate(self.vcvs_p_coords[i][0], self.vcvs_p_coords[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.vcvs_p_coords[i][0], self.vcvs_p_coords[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.vcvs_p_coords[i][2], self.vcvs_p_coords[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.vcvs_p_coords[i][2], self.vcvs_p_coords[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.vcvs_p[i], self.evl)

			self.vcvs_p_coords_tmp[i] = (self.evl)

	def move_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		if self.itr_move > self.no_of_itr_move:
			self.vcvs_p_coords = self.vcvs_p_coords_tmp[:]
			return
		self.move_kernal()
		master.after(int(1000./framerate), self.move_trigger)

	def translate(self, newposx, newposy, animation_time):
	 	self.itr_move = 0
	 	self.deltax1 = 0
	 	self.deltay1 = 0

	 	self.animation_time = animation_time
	 	self.deltax = newposx - self.getCenter()[0]
	 	self.deltay = newposy - self.getCenter()[1]
	 	self.move_trigger()


	def scrolate_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_scale(self.vcvs_p_coords[i][0],self.vcvs_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.vcvs_p_coords[i][0],self.vcvs_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
						math_scale(self.vcvs_p_coords[i][2],self.vcvs_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.vcvs_p_coords[i][2],self.vcvs_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

			c.coords(self.vcvs_p[i], self.evl)

			self.vcvs_p_coords_tmp[i] = self.evl

			self.evl = (math_rotate(self.vcvs_p_coords_tmp[i][0], self.vcvs_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.vcvs_p_coords_tmp[i][0], self.vcvs_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.vcvs_p_coords_tmp[i][2], self.vcvs_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.vcvs_p_coords_tmp[i][2], self.vcvs_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.vcvs_p[i], self.evl)

			self.vcvs_p_coords_tmp2[i] = self.evl

			self.evl = (math_translate(self.vcvs_p_coords_tmp2[i][0], self.vcvs_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.vcvs_p_coords_tmp2[i][0], self.vcvs_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.vcvs_p_coords_tmp2[i][2], self.vcvs_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.vcvs_p_coords_tmp2[i][2], self.vcvs_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.vcvs_p[i], self.evl)

			self.vcvs_p_coords_tmp3[i] = (self.evl)

			
	def scrolate_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.vcvs_p_coords = self.vcvs_p_coords_tmp3[:]
			return	
		
		self.scrolate_kernal()
		master.after(int(1000./framerate), self.scrolate_trigger)

	def scrolate(self, newposx, newposy, angle, scalex, scaley, animation_time):
		self.itr_move = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.animation_time = animation_time
		self.deltax = newposx - self.getCenter()[0]
		self.deltay = newposx - self.getCenter()[0]
		self.no_of_itr_move = self.animation_time * framerate/1000
		self.no_of_itr_scale = self.animation_time * framerate/1000
		self.no_of_itr_rotate = self.animation_time * framerate/1000

		self.angle = angle
		self.anglex = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.scaley = scaley
		self.scalex = scalex

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]
		self.scrolate_trigger()

	def scale_kernal(self):

		if self.appear:
			for i in range(0, self.items):

				self.evl =  (math_appear(self.vcvs_p_coords[i][0],self.vcvs_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.vcvs_p_coords[i][0],self.vcvs_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_appear(self.vcvs_p_coords[i][2],self.vcvs_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.vcvs_p_coords[i][2],self.vcvs_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.vcvs_p[i],self.evl)

				self.vcvs_p_coords_tmp[i] = self.evl
		else:
			for i in range(0, self.items):

				self.evl =  (math_scale(self.vcvs_p_coords[i][0],self.vcvs_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.vcvs_p_coords[i][0],self.vcvs_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_scale(self.vcvs_p_coords[i][2],self.vcvs_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.vcvs_p_coords[i][2],self.vcvs_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.vcvs_p[i],self.evl)

				self.vcvs_p_coords_tmp[i] = self.evl

	def scale_trigger(self):

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.vcvs_p_coords = self.vcvs_p_coords_tmp[:]
			self.appear = 0
			return

		self.scale_kernal()
		master.after(int(1000./framerate), self.scale_trigger)

	def scale(self, scalex, scaley, animation_time):
	 	self.animation_time = animation_time
	 	self.itr_scale = 0
	 	self.scaley = scaley
	 	self.scalex = scalex
	 	self.no_of_itr_scale = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.scale_trigger()
		
class VCCS:
	def init(self, posx, posy):
		self.posx = posx;
		self.posy = posy;
		self.items = 9
		self.vccs_p = [0]*self.items
		self.vccs_p_coords = [0]*self.items
		self.vccs_p_coords_tmp = [0]*self.items
		self.vccs_p_coords_tmp2 = [0]*self.items
		self.vccs_p_coords_tmp3 = [0]*self.items
		self.xmin_tmp = [0]*self.items
		self.xmax_tmp = [0]*self.items
		self.ymin_tmp = [0]*self.items
		self.ymax_tmp = [0]*self.items

		self.itr_move = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.no_of_itr_move = 0
		self.no_of_itr_rotate = 0
		self.no_of_itr_scale = 0
		self.anglex = 0
		self.angle = 0
		self.deltax = 0
		self.deltay = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.scalex_int = 0
		self.scaley_int = 0
		self.deltay1 = 0
		self.scalex = 0
		self.scaley = 0

		self.vccs_p_coords[0] = (posx, posy, posx, posy+20)
		self.vccs_p_coords[1] = (posx, posy+60, posx, posy+80)
		self.vccs_p_coords[2] = (posx, posy+20, posx-20, posy+40)
		self.vccs_p_coords[3] = (posx, posy+20, posx+20, posy+40)
		self.vccs_p_coords[4] = (posx-20, posy+40, posx, posy+60)
		self.vccs_p_coords[5] = (posx+20, posy+40, posx, posy+60)
		self.vccs_p_coords[6] = (posx, posy+28, posx, posy+50)
		self.vccs_p_coords[7] = (posx, posy+50, posx-5, posy+50-5)
		self.vccs_p_coords[8] = (posx, posy+50, posx+5, posy+50-5)

		self.reCenter()

		self.vccs_p[0] = c.create_line(self.vccs_p_coords[0], arrow='none', fill='#FFFFFF', width=2)
		self.vccs_p[1] = c.create_line(self.vccs_p_coords[1], arrow='none', fill='#FFFFFF', width=2)
		self.vccs_p[2] = c.create_line(self.vccs_p_coords[2], arrow='none', fill='#FFFFFF', width=2)
		self.vccs_p[3] = c.create_line(self.vccs_p_coords[3], arrow='none', fill='#FFFFFF', width=2)
		self.vccs_p[4] = c.create_line(self.vccs_p_coords[4], arrow='none', fill='#FFFFFF', width=2)
		self.vccs_p[5] = c.create_line(self.vccs_p_coords[5], arrow='none', fill='#FFFFFF', width=2)
		self.vccs_p[6] = c.create_line(self.vccs_p_coords[6], arrow='none', fill='#FFFFFF', width=2)
		self.vccs_p[7] = c.create_line(self.vccs_p_coords[7], arrow='none', fill='#FFFFFF', width=2)
		self.vccs_p[8] = c.create_line(self.vccs_p_coords[8], arrow='none', fill='#FFFFFF', width=2)

		self.appear = 1
		self.scale(1,1,500)


	def getTerminals(self):
		return [[self.vccs_p_coords[0][0], self.vccs_p_coords[0][1]], [self.vccs_p_coords[1][2], self.vccs_p_coords[1][3]]]

	def reCenter(self):
		self.rotate_around_x_old = self.getCenter()[0]
		self.rotate_around_y_old = self.getCenter()[1]

		for i in range(0, self.items):
			self.vccs_p_coords[i] = (self.vccs_p_coords[i][0] + self.posx - self.rotate_around_x_old, self.vccs_p_coords[i][1] + self.posy - self.rotate_around_y_old, self.vccs_p_coords[i][2] + self.posx - self.rotate_around_x_old, self.vccs_p_coords[i][3] + self.posy -self.rotate_around_y_old)

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]

	def getCenter(self):
		for i in range(0, self.items):
			self.xmin_tmp[i] = min(self.vccs_p_coords[i][0], self.vccs_p_coords[i][2])
			self.xmax_tmp[i] = max(self.vccs_p_coords[i][0], self.vccs_p_coords[i][2])
			self.ymin_tmp[i] = min(self.vccs_p_coords[i][1], self.vccs_p_coords[i][3])
			self.ymax_tmp[i] = max(self.vccs_p_coords[i][1], self.vccs_p_coords[i][3])
		return [(min(self.xmin_tmp) + max(self.xmax_tmp))/2, (min(self.ymin_tmp) + max(self.ymax_tmp))/2]


	def rotate_kernal(self):

		for i in range(0, self.items):

			self.evl = (math_rotate(self.vccs_p_coords[i][0], self.vccs_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.vccs_p_coords[i][0], self.vccs_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.vccs_p_coords[i][2], self.vccs_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.vccs_p_coords[i][2], self.vccs_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.vccs_p[i], self.evl)

			self.vccs_p_coords_tmp[i] = (self.evl)

	def rotate_trigger(self):

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		if self.itr_rotate > self.no_of_itr_rotate:
			self.vccs_p_coords = self.vccs_p_coords_tmp[:]
			return

		self.rotate_kernal()
		master.after(int(1000./framerate), self.rotate_trigger)

	def rotate(self, angle, animation_time):
	 	self.animation_time = animation_time
	 	self.angle = angle
	 	self.anglex = 0
	 	self.itr_rotate = 0
	 	self.no_of_itr_rotate = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.rotate_trigger()

	def move_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_translate(self.vccs_p_coords[i][0], self.vccs_p_coords[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.vccs_p_coords[i][0], self.vccs_p_coords[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.vccs_p_coords[i][2], self.vccs_p_coords[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.vccs_p_coords[i][2], self.vccs_p_coords[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.vccs_p[i], self.evl)

			self.vccs_p_coords_tmp[i] = (self.evl)

	def move_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		if self.itr_move > self.no_of_itr_move:
			self.vccs_p_coords = self.vccs_p_coords_tmp[:]
			return

		self.move_kernal()
		master.after(int(1000./framerate), self.move_trigger)

	def translate(self, newposx, newposy, animation_time):
	 	self.itr_move = 0
	 	self.deltax1 = 0
	 	self.deltay1 = 0
	 	self.no_of_itr_move = self.animation_time * framerate/1000
	 	self.animation_time = animation_time
	 	self.deltax = newposx - self.getCenter()[0]
	 	self.deltay = newposy - self.getCenter()[1]
	 	self.move_trigger()


	def scrolate_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_scale(self.vccs_p_coords[i][0],self.vccs_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.vccs_p_coords[i][0],self.vccs_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
						math_scale(self.vccs_p_coords[i][2],self.vccs_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.vccs_p_coords[i][2],self.vccs_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

			c.coords(self.vccs_p[i], self.evl)

			self.vccs_p_coords_tmp[i] = self.evl

			self.evl = (math_rotate(self.vccs_p_coords_tmp[i][0], self.vccs_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.vccs_p_coords_tmp[i][0], self.vccs_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.vccs_p_coords_tmp[i][2], self.vccs_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.vccs_p_coords_tmp[i][2], self.vccs_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.vccs_p[i], self.evl)

			self.vccs_p_coords_tmp2[i] = self.evl

			self.evl = (math_translate(self.vccs_p_coords_tmp2[i][0], self.vccs_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.vccs_p_coords_tmp2[i][0], self.vccs_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.vccs_p_coords_tmp2[i][2], self.vccs_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.vccs_p_coords_tmp2[i][2], self.vccs_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.vccs_p[i], self.evl)

			self.vccs_p_coords_tmp3[i] = (self.evl)

			
	def scrolate_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		if self.itr_rotate == 0:
			self.no_of_itr_rotate = self.animation_time*framerate/1000

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.vccs_p_coords = self.vccs_p_coords_tmp3[:]
			return	
		
		self.scrolate_kernal()
		master.after(int(1000./framerate), self.scrolate_trigger)

	def scrolate(self, newposx, newposy, angle, scalex, scaley, animation_time):
		self.itr_move = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.animation_time = animation_time
		self.deltax = newposx - self.getCenter()[0]
		self.deltay = newposx - self.getCenter()[0]
		self.no_of_itr_move = self.animation_time * framerate/1000
		self.no_of_itr_scale = self.animation_time * framerate/1000
		self.no_of_itr_rotate = self.animation_time * framerate/1000

		self.angle = angle
		self.anglex = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.scaley = scaley
		self.scalex = scalex

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]
		self.scrolate_trigger()

	def scale_kernal(self):

		if self.appear:
			for i in range(0, self.items):

				self.evl =  (math_appear(self.vccs_p_coords[i][0],self.vccs_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.vccs_p_coords[i][0],self.vccs_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_appear(self.vccs_p_coords[i][2],self.vccs_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.vccs_p_coords[i][2],self.vccs_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.vccs_p[i],self.evl)

				self.vccs_p_coords_tmp[i] = self.evl
		else:
			for i in range(0, self.items):

				self.evl =  (math_scale(self.vccs_p_coords[i][0],self.vccs_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.vccs_p_coords[i][0],self.vccs_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_scale(self.vccs_p_coords[i][2],self.vccs_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.vccs_p_coords[i][2],self.vccs_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.vccs_p[i],self.evl)

				self.vccs_p_coords_tmp[i] = self.evl

	def scale_trigger(self):

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.vccs_p_coords = self.vccs_p_coords_tmp[:]
			self.appear = 0
			return

		self.scale_kernal()
		master.after(int(1000./framerate), self.scale_trigger)

	def scale(self, scalex, scaley, animation_time):
	 	self.animation_time = animation_time
	 	self.itr_scale = 0
	 	self.scaley = scaley
	 	self.scalex = scalex
	 	self.no_of_itr_scale = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.scale_trigger()


class PLOT:
	def init(self, posx, posy, size=cheight/2):
		self.posx = posx;
		self.posy = posy;
		self.size = size

		self.items = 20
		self.plot_p = [0]*self.items
		self.plot_p_coords = [0]*self.items
		self.plot_p_coords_tmp = [0]*self.items
		self.plot_p_coords_tmp2 = [0]*self.items
		self.plot_p_coords_tmp3 = [0]*self.items
		self.xmin_tmp = [0]*self.items
		self.xmax_tmp = [0]*self.items
		self.ymin_tmp = [0]*self.items
		self.ymax_tmp = [0]*self.items

		self.itr_move = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.no_of_itr_move = 0
		self.no_of_itr_rotate = 0
		self.no_of_itr_scale = 0
		self.anglex = 0
		self.angle = 0
		self.deltax = 0
		self.deltay = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.scalex_int = 0
		self.scaley_int = 0
		self.deltay1 = 0
		self.scalex = 0
		self.scaley = 0

		self.plot_p_coords[0] = (posx, posy, posx+self.size, posy)
		self.plot_p_coords[1] = (posx, posy, posx, posy-self.size)

		for i in range(0, 9):
			self.plot_p_coords[i+2] = (posx+(self.size/10)*(i+1), posy-4, posx+(self.size/10)*(i+1), posy+4)
		for i in range(0, 9):
			self.plot_p_coords[9+i+2] = (posx-4, posy-(self.size/10)*(i+1), posx+4, posy-(self.size/10)*(i+1))

		for i in range(0, self.items):
			if (i == 0) or (i == 1):
				self.plot_p[i] = c.create_line(self.plot_p_coords[i], arrow='last', fill='#FFFFFF', width=1)
			else:	
				self.plot_p[i] = c.create_line(self.plot_p_coords[i], arrow='none', fill='#FFFFFF', width=1)

		self.appear = 1
		self.scale(1,1,500)


	def getTerminals(self):
		return [[self.plot_p_coords[0][0], self.plot_p_coords[0][1]], [self.plot_p_coords[1][2], self.plot_p_coords[1][3]]]

	def getCenter(self):
		return [self.plot_p_coords[0][0], self.plot_p_coords[0][1]]

	def rotate_kernal(self):

		for i in range(0, self.items):

			self.evl = (math_rotate(self.plot_p_coords[i][0], self.plot_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.plot_p_coords[i][0], self.plot_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.plot_p_coords[i][2], self.plot_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.plot_p_coords[i][2], self.plot_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.plot_p[i], self.evl)

			self.plot_p_coords_tmp[i] = (self.evl)

	def rotate_trigger(self):

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		if self.itr_rotate > self.no_of_itr_rotate:
			self.plot_p_coords = self.plot_p_coords_tmp[:]
			return

		self.rotate_kernal()
		master.after(int(1000./framerate), self.rotate_trigger)

	def rotate(self, angle, animation_time):
	 	self.animation_time = animation_time
	 	self.angle = angle
	 	self.anglex = 0
	 	self.itr_rotate = 0
	 	self.no_of_itr_rotate = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.rotate_trigger()

	def move_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_translate(self.plot_p_coords[i][0], self.plot_p_coords[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.plot_p_coords[i][0], self.plot_p_coords[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.plot_p_coords[i][2], self.plot_p_coords[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.plot_p_coords[i][2], self.plot_p_coords[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.plot_p[i], self.evl)

			self.plot_p_coords_tmp[i] = (self.evl)

	def move_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		if self.itr_move > self.no_of_itr_move:
			self.plot_p_coords = self.plot_p_coords_tmp[:]
			return

		self.move_kernal()
		master.after(int(1000./framerate), self.move_trigger)

	def translate(self, newposx, newposy, animation_time):
	 	self.itr_move = 0
	 	self.deltax1 = 0
	 	self.deltay1 = 0
	 	self.animation_time = animation_time
	 	self.no_of_itr_move = self.animation_time * framerate/1000
	 	self.deltax = newposx - self.getCenter()[0]
	 	self.deltay = newposy - self.getCenter()[1]
	 	self.move_trigger()

	def scrolate_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_scale(self.plot_p_coords[i][0],self.plot_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.plot_p_coords[i][0],self.plot_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
						math_scale(self.plot_p_coords[i][2],self.plot_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.plot_p_coords[i][2],self.plot_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

			c.coords(self.plot_p[i], self.evl)

			self.plot_p_coords_tmp[i] = self.evl

			self.evl = (math_rotate(self.plot_p_coords_tmp[i][0], self.plot_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.plot_p_coords_tmp[i][0], self.plot_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.plot_p_coords_tmp[i][2], self.plot_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.plot_p_coords_tmp[i][2], self.plot_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.plot_p[i], self.evl)

			self.plot_p_coords_tmp2[i] = self.evl

			self.evl = (math_translate(self.plot_p_coords_tmp2[i][0], self.plot_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.plot_p_coords_tmp2[i][0], self.plot_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.plot_p_coords_tmp2[i][2], self.plot_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.plot_p_coords_tmp2[i][2], self.plot_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.plot_p[i], self.evl)

			self.plot_p_coords_tmp3[i] = self.evl

			
	def scrolate_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.plot_p_coords = self.plot_p_coords_tmp3[:]
			return	
		
		self.scrolate_kernal()
		master.after(int(1000./framerate), self.scrolate_trigger)

	def scrolate(self, newposx, newposy, angle, scalex, scaley, animation_time):
		self.itr_move = 0
		self.animation_time = animation_time
		self.deltax = newposx - self.getCenter()[0]
		self.deltay = newposy - self.getCenter()[1]
		self.no_of_itr_move = self.animation_time * framerate/1000
		self.no_of_itr_scale = self.animation_time * framerate/1000
		self.no_of_itr_rotate = self.animation_time * framerate/1000

		self.angle = angle
		self.anglex = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.scaley = scaley
		self.scalex = scalex

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]
		self.scrolate_trigger()

	def scale_kernal(self):

		if self.appear:
			for i in range(0, self.items):

				self.evl =  (math_appear(self.plot_p_coords[i][0],self.plot_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.plot_p_coords[i][0],self.plot_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_appear(self.plot_p_coords[i][2],self.plot_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.plot_p_coords[i][2],self.plot_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.plot_p[i],self.evl)

				self.plot_p_coords_tmp[i] = self.evl
		else:
			for i in range(0, self.items):

				self.evl =  (math_scale(self.plot_p_coords[i][0],self.plot_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.plot_p_coords[i][0],self.plot_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_scale(self.plot_p_coords[i][2],self.plot_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.plot_p_coords[i][2],self.plot_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.plot_p[i],self.evl)

				self.plot_p_coords_tmp[i] = self.evl

	def scale_trigger(self):

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.plot_p_coords = self.plot_p_coords_tmp[:]
			self.appear = 0
			return

		self.scale_kernal()
		master.after(int(1000./framerate), self.scale_trigger)

	def scale(self, scalex, scaley, animation_time):
	 	self.animation_time = animation_time
	 	self.itr_scale = 0
	 	self.scaley = scaley
	 	self.scalex = scalex
	 	self.no_of_itr_scale = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.scale_trigger()


class TRACE:
	def init(self, points, nTimes=50, color='#FFFFFF'):
		xvals, yvals = bezier_curve(points, nTimes)

		self.items = len(xvals)
		self.trace_p = [0]*self.items
		self.trace_p_coords = [0]*self.items
		self.trace_p_coords_tmp = [0]*self.items
		self.trace_p_coords_tmp2 = [0]*self.items
		self.trace_p_coords_tmp3 = [0]*self.items
		self.xmin_tmp = [0]*self.items
		self.xmax_tmp = [0]*self.items
		self.ymin_tmp = [0]*self.items
		self.ymax_tmp = [0]*self.items

		for i in range(0, len(xvals)-1):
			self.trace_p_coords[i] = (xvals[self.items-i-1], yvals[self.items-i-1], xvals[self.items-i-2], yvals[self.items-i-2])

		for i in range(0,len(xvals)-1):
			self.trace_p[i] = c.create_line(self.trace_p_coords[i], arrow='none', fill=color, width=2)

		self.itr_move = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.no_of_itr_move = 0
		self.no_of_itr_rotate = 0
		self.no_of_itr_scale = 0
		self.anglex = 0
		self.angle = 0
		self.deltax = 0
		self.deltay = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.scalex_int = 0
		self.scaley_int = 0
		self.deltay1 = 0
		self.scalex = 0
		self.scaley = 0

		self.appear = 1
		self.scale(1,1,500)

	def getCenter(self):
		return [self.trace_p_coords[0][0], self.trace_p_coords[0][1]]

	def rotate_kernal(self):

		for i in range(0, self.items-1):

			self.evl = (math_rotate(self.trace_p_coords[i][0], self.trace_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.trace_p_coords[i][0], self.trace_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.trace_p_coords[i][2], self.trace_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.trace_p_coords[i][2], self.trace_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.trace_p[i], self.evl)

			self.trace_p_coords_tmp[i] = (self.evl)

	def rotate_trigger(self):

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		if self.itr_rotate > self.no_of_itr_rotate:
			self.trace_p_coords = self.trace_p_coords_tmp[:]
			return

		self.rotate_kernal()
		master.after(int(1000./framerate), self.rotate_trigger)

	def rotate(self, angle, animation_time):
	 	self.animation_time = animation_time
	 	self.angle = angle
	 	self.anglex = 0
	 	self.itr_rotate = 0
	 	self.no_of_itr_rotate = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.rotate_trigger()

	def move_kernal(self):
		for i in range(0,self.items-1):

			self.evl = (math_translate(self.trace_p_coords[i][0], self.trace_p_coords[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.trace_p_coords[i][0], self.trace_p_coords[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.trace_p_coords[i][2], self.trace_p_coords[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.trace_p_coords[i][2], self.trace_p_coords[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.trace_p[i], self.evl)

			self.trace_p_coords_tmp[i] = (self.evl)

	def move_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		if self.itr_move > self.no_of_itr_move:
			self.trace_p_coords = self.trace_p_coords_tmp[:]
			return
		self.move_kernal()
		master.after(int(1000./framerate), self.move_trigger)

	def translate(self, newposx, newposy, animation_time):
	 	self.itr_move = 0
	 	self.deltax1 = 0
	 	self.deltay1 = 0
	 	self.no_of_itr_move = self.animation_time * framerate/1000
	 	self.animation_time = animation_time
	 	self.deltax = newposx - self.getCenter()[0]
	 	self.deltay = newposy - self.getCenter()[1]
	 	self.move_trigger()

	def scrolate_kernal(self):
		for i in range(0,self.items-1):

			self.evl = (math_scale(self.trace_p_coords[i][0],self.trace_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.trace_p_coords[i][0],self.trace_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
						math_scale(self.trace_p_coords[i][2],self.trace_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.trace_p_coords[i][2],self.trace_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

			c.coords(self.trace_p[i], self.evl)

			self.trace_p_coords_tmp[i] = self.evl

			self.evl = (math_rotate(self.trace_p_coords_tmp[i][0], self.trace_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.trace_p_coords_tmp[i][0], self.trace_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.trace_p_coords_tmp[i][2], self.trace_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.trace_p_coords_tmp[i][2], self.trace_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.trace_p[i], self.evl)

			self.trace_p_coords_tmp2[i] = self.evl

			self.evl = (math_translate(self.trace_p_coords_tmp2[i][0], self.trace_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.trace_p_coords_tmp2[i][0], self.trace_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.trace_p_coords_tmp2[i][2], self.trace_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.trace_p_coords_tmp2[i][2], self.trace_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.trace_p[i], self.evl)

			self.trace_p_coords_tmp3[i] = self.evl

			
	def scrolate_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.trace_p_coords = self.trace_p_coords_tmp3[:]
			return	
		
		self.scrolate_kernal()
		master.after(int(1000./framerate), self.scrolate_trigger)

	def scrolate(self, newposx, newposy, angle, scalex, scaley, animation_time):
		self.itr_move = 0
		self.animation_time = animation_time
		self.deltax = newposx - self.getCenter()[0]
		self.deltay = newposy - self.getCenter()[1]
		self.no_of_itr_move = self.animation_time * framerate/1000
		self.no_of_itr_scale = self.animation_time * framerate/1000
		self.no_of_itr_rotate = self.animation_time * framerate/1000

		self.angle = angle
		self.anglex = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.scaley = scaley
		self.scalex = scalex

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]
		self.scrolate_trigger()


	def scale_kernal(self):

		if self.appear:
			for i in range(0, self.items-1):

				self.evl =  (math_appear(self.trace_p_coords[i][0],self.trace_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.trace_p_coords[i][0],self.trace_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_appear(self.trace_p_coords[i][2],self.trace_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.trace_p_coords[i][2],self.trace_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.trace_p[i],self.evl)

				self.trace_p_coords_tmp[i] = self.evl
		else:
			for i in range(0, self.items-1):

				self.evl =  (math_scale(self.trace_p_coords[i][0],self.trace_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.trace_p_coords[i][0],self.trace_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_scale(self.trace_p_coords[i][2],self.trace_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.trace_p_coords[i][2],self.trace_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.trace_p[i],self.evl)

				self.trace_p_coords_tmp[i] = self.evl

	def scale_trigger(self):

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.trace_p_coords = self.trace_p_coords_tmp[:]
			self.appear = 0
			return

		self.scale_kernal()
		master.after(int(1000./framerate), self.scale_trigger)

	def scale(self, scalex, scaley, animation_time):
	 	self.animation_time = animation_time
	 	self.itr_scale = 0
	 	self.scaley = scaley
	 	self.scalex = scalex
	 	self.no_of_itr_scale = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.scale_trigger()


class WAVE:
	def init(self, posx, posy, eqn_in_x, stepsx=cheight/2, step_size=2, color='#FFFFFF'):

		self.posx = posx
		self.posy = posy

		xvals, yvals = eqn_interpolator(eqn_in_x, stepsx - stepsx/6, step_size)

		self.items = len(xvals)
		self.wave_p = [0]*self.items
		self.wave_p_coords = [0]*self.items
		self.wave_p_coords_tmp = [0]*self.items
		self.wave_p_coords_tmp2 = [0]*self.items
		self.xmin_tmp = [0]*self.items
		self.xmax_tmp = [0]*self.items
		self.ymin_tmp = [0]*self.items
		self.ymax_tmp = [0]*self.items

		for i in range(0, len(xvals)-1):
			self.wave_p_coords[i] = (self.posx + xvals[i], self.posy - yvals[i], self.posx + xvals[i+1], self.posy - yvals[i+1])

		for i in range(0,len(xvals)-1):
			self.wave_p[i] = c.create_line(self.wave_p_coords[i], arrow='none', fill=color, width=2)

		self.itr_move = 0
		self.itr_rotate = 0
		self.no_of_itr_move = 0
		self.no_of_itr_rotate = 0
		self.anglex = 0
		self.angle = 0
		self.deltax = 0
		self.deltay = 0
		self.deltax1 = 0
		self.deltay1 = 0

		self.appear = 1
		self.scale(1,1,500)

	def getCenter(self):
		return [self.wave_p_coords[0][0], self.wave_p_coords[0][1]]

	def rotate_kernal(self):

		for i in range(0, self.items-1):

			self.evl = (math_rotate(self.wave_p_coords[i][0], self.wave_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.wave_p_coords[i][0], self.wave_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.wave_p_coords[i][2], self.wave_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.wave_p_coords[i][2], self.wave_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.wave_p[i], self.evl)

			self.wave_p_coords_tmp[i] = (self.evl)

	def rotate_trigger(self):

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		if self.itr_rotate > self.no_of_itr_rotate:
			self.wave_p_coords = self.wave_p_coords_tmp[:]
			return

		self.rotate_kernal()
		master.after(int(1000./framerate), self.rotate_trigger)

	def rotate(self, angle, animation_time):
	 	self.animation_time = animation_time
	 	self.angle = angle
	 	self.anglex = 0
	 	self.itr_rotate = 0
	 	self.no_of_itr_rotate = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.rotate_trigger()

	def move_kernal(self):
		for i in range(0,self.items-1):

			self.evl = (math_translate(self.wave_p_coords[i][0], self.wave_p_coords[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.wave_p_coords[i][0], self.wave_p_coords[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.wave_p_coords[i][2], self.wave_p_coords[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.wave_p_coords[i][2], self.wave_p_coords[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.wave_p[i], self.evl)

			self.wave_p_coords_tmp[i] = self.evl

	def move_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		if self.itr_move > self.no_of_itr_move:
			self.wave_p_coords = self.wave_p_coords_tmp[:]
			return
		self.move_kernal()
		master.after(int(1000./framerate), self.move_trigger)

	def translate(self, newposx, newposy, animation_time):
	 	self.itr_move = 0
	 	self.deltax1 = 0
	 	self.deltay1 = 0
	 	self.animation_time = animation_time
	 	self.no_of_itr_move = self.animation_time * framerate/1000
	 	self.deltax = newposx - self.getCenter()[0]
	 	self.deltay = newposy - self.getCenter()[1]
	 	self.move_trigger()

	def scale_kernal(self):

		if self.appear:
			for i in range(0, self.items-1):

				self.evl =  (math_appear(self.wave_p_coords[i][0],self.wave_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.wave_p_coords[i][0],self.wave_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_appear(self.wave_p_coords[i][2],self.wave_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.wave_p_coords[i][2],self.wave_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.wave_p[i],self.evl)

				self.wave_p_coords_tmp[i] = self.evl
		else:
			for i in range(0, self.items-1):

				self.evl =  (math_scale(self.wave_p_coords[i][0],self.wave_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.wave_p_coords[i][0],self.wave_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_scale(self.wave_p_coords[i][2],self.wave_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.wave_p_coords[i][2],self.wave_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.wave_p[i],self.evl)

				self.wave_p_coords_tmp[i] = self.evl

	def scale_trigger(self):

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.wave_p_coords = self.wave_p_coords_tmp[:]
			self.appear = 0
			return

		self.scale_kernal()
		master.after(int(1000./framerate), self.scale_trigger)

	def scale(self, scalex, scaley, animation_time):
	 	self.animation_time = animation_time
	 	self.itr_scale = 0
	 	self.scaley = scaley
	 	self.scalex = scalex
	 	self.no_of_itr_scale = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.scale_trigger()


class WAVE2:
	def init(self, posx, posy, eqn_in_x, stepsx=cheight/2, step_size=1, color='#FFFFFF'):

		self.posx = posx
		self.posy = posy

		xvals, yvals = eqn_interpolator(eqn_in_x, 0.9*stepsx, step_size)

		self.items = len(xvals)
		self.wave2_p = [0]*self.items
		self.wave2_p_coords = [0]*self.items
		self.wave2_p_coords_tmp = [0]*self.items
		self.wave2_p_coords_tmp2 = [0]*self.items
		self.xmin_tmp = [0]*self.items
		self.xmax_tmp = [0]*self.items
		self.ymin_tmp = [0]*self.items
		self.ymax_tmp = [0]*self.items

		for i in range(0, len(xvals)-1):
			self.wave2_p_coords[i] = (self.posx + xvals[i], self.posy - yvals[i], self.posx + xvals[i+1], self.posy - yvals[i+1])

		for i in range(0,len(xvals)-1):
			self.wave2_p[i] = c.create_line(self.wave2_p_coords[i], arrow='none', fill=color, width=2)

		self.itr_move = 0
		self.itr_rotate = 0
		self.no_of_itr_move = 0
		self.no_of_itr_rotate = 0
		self.anglex = 0
		self.angle = 0
		self.deltax = 0
		self.deltay = 0
		self.deltax1 = 0
		self.deltay1 = 0

		self.appear = 1
		self.scale(1,1,500)

	def reCenter(self):
		self.rotate_around_x_old = self.getCenter()[0]
		self.rotate_around_y_old = self.getCenter()[1]

		for i in range(0, self.items-1):
			self.wave2_p_coords[i] = (self.wave2_p_coords[i][0] + self.posx - self.rotate_around_x_old, self.wave2_p_coords[i][1] + self.posy - self.rotate_around_y_old, self.wave2_p_coords[i][2] + self.posx - self.rotate_around_x_old, self.wave2_p_coords[i][3] + self.posy -self.rotate_around_y_old)

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]

	def getCenter(self):
		for i in range(0, self.items-1):
			self.xmin_tmp[i] = min(self.wave2_p_coords[i][0], self.wave2_p_coords[i][2])
			self.xmax_tmp[i] = max(self.wave2_p_coords[i][0], self.wave2_p_coords[i][2])
			self.ymin_tmp[i] = min(self.wave2_p_coords[i][1], self.wave2_p_coords[i][3])
			self.ymax_tmp[i] = max(self.wave2_p_coords[i][1], self.wave2_p_coords[i][3])
		self.xmin_tmp[self.items-1] = self.xmin_tmp[self.items-2]
		self.xmax_tmp[self.items-1] = self.xmax_tmp[self.items-2]
		self.ymin_tmp[self.items-1] = self.ymin_tmp[self.items-2]
		self.ymax_tmp[self.items-1] = self.ymax_tmp[self.items-2]
		return [(min(self.xmin_tmp) + max(self.xmax_tmp))/2, (min(self.ymin_tmp) + max(self.ymax_tmp))/2]

	def rotate_kernal(self):

		for i in range(0, self.items-1):

			self.evl = ()
			c.coords(self.wave2_p[i], self.evl)

			self.wave2_p_coords_tmp[i] = (self.evl)

	def rotate_trigger(self):

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		if self.itr_rotate > self.no_of_itr_rotate:
			self.wave2_p_coords = self.wave2_p_coords_tmp[:]
			return

		self.rotate_kernal()
		master.after(int(1000./framerate), self.rotate_trigger)

	def rotate(self, angle, animation_time):
	 	self.animation_time = animation_time
	 	self.angle = angle
	 	self.anglex = 0
	 	self.itr_rotate = 0
	 	self.no_of_itr_rotate = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.rotate_trigger()

	def move_kernal(self):
		for i in range(0,self.items-1):

			self.evl = (math_translate(self.wave2_p_coords[i][0], self.wave2_p_coords[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.wave2_p_coords[i][0], self.wave2_p_coords[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.wave2_p_coords[i][2], self.wave2_p_coords[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.wave2_p_coords[i][2], self.wave2_p_coords[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.wave2_p[i], self.evl)

			self.wave2_p_coords_tmp[i] = (self.evl)

	def move_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		if self.itr_move > self.no_of_itr_move:
			self.wave2_p_coords = self.wave2_p_coords_tmp[:]
			return
		self.move_kernal()
		master.after(int(1000./framerate), self.move_trigger)

	def translate(self, newposx, newposy, animation_time):
	 	self.itr_move = 0
	 	self.deltax1 = 0
	 	self.deltay1 = 0
	 	self.no_of_itr_move = self.animation_time * framerate/1000
	 	self.animation_time = animation_time
	 	self.deltax = newposx - self.getCenter()[0]
	 	self.deltay = newposy - self.getCenter()[1]
	 	self.move_trigger()

	def scale_kernal(self):

		if self.appear:
			for i in range(0, self.items-1):

				self.evl =  (math_appear(self.wave2_p_coords[i][0],self.wave2_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.wave2_p_coords[i][0],self.wave2_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_appear(self.wave2_p_coords[i][2],self.wave2_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.wave2_p_coords[i][2],self.wave2_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.wave2_p[i],self.evl)

				self.wave2_p_coords_tmp[i] = self.evl
		else:
			for i in range(0, self.items-1):

				self.evl =  (math_scale(self.wave2_p_coords[i][0],self.wave2_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.wave2_p_coords[i][0],self.wave2_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_scale(self.wave2_p_coords[i][2],self.wave2_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.wave2_p_coords[i][2],self.wave2_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.wave2_p[i],self.evl)

				self.wave2_p_coords_tmp[i] = self.evl

	def scale_trigger(self):

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.wave2_p_coords = self.wave2_p_coords_tmp[:]
			self.appear = 0
			return

		self.scale_kernal()
		master.after(int(1000./framerate), self.scale_trigger)

	def scale(self, scalex, scaley, animation_time):
	 	self.animation_time = animation_time
	 	self.itr_scale = 0
	 	self.scaley = scaley
	 	self.scalex = scalex
	 	self.no_of_itr_scale = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.scale_trigger()


class SUPPLY:

	def init(self, posx, posy, color='#FFFFFF'):

		self.posx = posx
		self.posy = posy
		self.color = color

		self.items = 8
		self.supply_p = [0]*self.items
		self.supply_p_coords = [0]*self.items
		self.supply_p_coords_tmp = [0]*self.items
		self.supply_p_coords_tmp2 = [0]*self.items
		self.supply_p_coords_tmp3 = [0]*self.items
		self.xmin_tmp = [0]*self.items
		self.xmax_tmp = [0]*self.items
		self.ymin_tmp = [0]*self.items
		self.ymax_tmp = [0]*self.items

		self.itr_move = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.no_of_itr_move = 0
		self.no_of_itr_rotate = 0
		self.no_of_itr_scale = 0
		self.anglex = 0
		self.angle = 0
		self.deltax = 0
		self.deltay = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.scalex_int = 0
		self.scaley_int = 0
		self.deltay1 = 0
		self.scalex = 0
		self.scaley = 0

		self.supply_p_coords[0] = (posx-30, posy, posx+20, posy)
		self.supply_p_coords[1] = (posx, posy, posx, posy+20)
		self.supply_p_coords[2] = (posx-30, posy, posx-20, posy-15)
		self.supply_p_coords[3] = (posx-20, posy, posx-10, posy-15)
		self.supply_p_coords[4] = (posx-10, posy, posx, posy-15)
		self.supply_p_coords[5] = (posx, posy, posx+10, posy-15)
		self.supply_p_coords[6] = (posx+10, posy, posx+20, posy-15)
		self.supply_p_coords[7] = (posx+20, posy, posx+30, posy-15)

		self.reCenter()

		self.supply_p[0] = c.create_line(self.supply_p_coords[0], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.supply_p[1] = c.create_line(self.supply_p_coords[1], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.supply_p[2] = c.create_line(self.supply_p_coords[2], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.supply_p[3] = c.create_line(self.supply_p_coords[3], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.supply_p[4] = c.create_line(self.supply_p_coords[4], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.supply_p[5] = c.create_line(self.supply_p_coords[5], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.supply_p[6] = c.create_line(self.supply_p_coords[6], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.supply_p[7] = c.create_line(self.supply_p_coords[7], arrow='none', fill=self.color, width=2, capstyle=ROUND)

		self.appear = 1
		self.scale(1,1,500)


	def getTerminals(self):
		return [[self.supply_p_coords[0][0], self.supply_p_coords[0][1]], [self.supply_p_coords[11][2], self.supply_p_coords[11][3]]]

	def reCenter(self):
		self.rotate_around_x_old = self.getCenter()[0]
		self.rotate_around_y_old = self.getCenter()[1]

		for i in range(0, self.items):
			self.supply_p_coords[i] = (self.supply_p_coords[i][0] + self.posx - self.rotate_around_x_old, self.supply_p_coords[i][1] + self.posy - self.rotate_around_y_old, self.supply_p_coords[i][2] + self.posx - self.rotate_around_x_old, self.supply_p_coords[i][3] + self.posy -self.rotate_around_y_old)

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]

	def getCenter(self):
		for i in range(0, self.items):
			self.xmin_tmp[i] = min(self.supply_p_coords[i][0], self.supply_p_coords[i][2])
			self.xmax_tmp[i] = max(self.supply_p_coords[i][0], self.supply_p_coords[i][2])
			self.ymin_tmp[i] = min(self.supply_p_coords[i][1], self.supply_p_coords[i][3])
			self.ymax_tmp[i] = max(self.supply_p_coords[i][1], self.supply_p_coords[i][3])
		return [(min(self.xmin_tmp) + max(self.xmax_tmp))/2, (min(self.ymin_tmp) + max(self.ymax_tmp))/2]

	def rotate_kernal(self):

		for i in range(0, self.items):

			self.evl =  (math_rotate(self.supply_p_coords[i][0], self.supply_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.supply_p_coords[i][0], self.supply_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.supply_p_coords[i][2], self.supply_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.supply_p_coords[i][2], self.supply_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])
			c.coords(self.supply_p[i],self.evl)

			self.supply_p_coords_tmp[i] = self.evl

	def rotate_trigger(self):

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		if self.itr_rotate > self.no_of_itr_rotate:
			self.supply_p_coords = self.supply_p_coords_tmp[:]
			return

		self.rotate_kernal()
		master.after(int(1000./framerate), self.rotate_trigger)

	def rotate(self, angle, animation_time):
	 	self.animation_time = animation_time
	 	self.angle = angle
	 	self.anglex = 0
	 	self.itr_rotate = 0
	 	self.no_of_itr_rotate = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.rotate_trigger()

	def move_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_translate(self.supply_p_coords[i][0], self.supply_p_coords[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.supply_p_coords[i][0], self.supply_p_coords[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.supply_p_coords[i][2], self.supply_p_coords[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.supply_p_coords[i][2], self.supply_p_coords[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.supply_p[i], self.evl)

			self.supply_p_coords_tmp[i] = (self.evl)

	def move_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		if self.itr_move > self.no_of_itr_move:
			self.supply_p_coords = self.supply_p_coords_tmp[:]
			return
		self.move_kernal()
		master.after(int(1000./framerate), self.move_trigger)

	def translate(self, newposx, newposy, animation_time):
	 	self.itr_move = 0
	 	self.deltax1 = 0
	 	self.deltay1 = 0
	 	self.no_of_itr_move = self.animation_time * framerate/1000
	 	self.animation_time = animation_time
	 	self.deltax = newposx - self.getCenter()[0]
	 	self.deltay = newposy - self.getCenter()[1]
	 	self.move_trigger()


	def scrolate_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_scale(self.supply_p_coords[i][0],self.supply_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.supply_p_coords[i][0],self.supply_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
						math_scale(self.supply_p_coords[i][2],self.supply_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.supply_p_coords[i][2],self.supply_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

			c.coords(self.supply_p[i], self.evl)

			self.supply_p_coords_tmp[i] = self.evl

			self.evl = (math_rotate(self.supply_p_coords_tmp[i][0], self.supply_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.supply_p_coords_tmp[i][0], self.supply_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.supply_p_coords_tmp[i][2], self.supply_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.supply_p_coords_tmp[i][2], self.supply_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.supply_p[i], self.evl)

			self.supply_p_coords_tmp2[i] = self.evl

			self.evl = (math_translate(self.supply_p_coords_tmp2[i][0], self.supply_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.supply_p_coords_tmp2[i][0], self.supply_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.supply_p_coords_tmp2[i][2], self.supply_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.supply_p_coords_tmp2[i][2], self.supply_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.supply_p[i], self.evl)

			self.supply_p_coords_tmp3[i] = (self.evl)

			
	def scrolate_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)


		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.supply_p_coords = self.supply_p_coords_tmp3[:]
			return	
		
		self.scrolate_kernal()
		master.after(int(1000./framerate), self.scrolate_trigger)

	def scrolate(self, newposx, newposy, angle, scalex, scaley, animation_time):
		self.itr_move = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.no_of_itr_move = self.animation_time * framerate/1000
		self.no_of_itr_scale = self.animation_time * framerate/1000
		self.no_of_itr_rotate = self.animation_time * framerate/1000
		self.animation_time = animation_time
		self.deltax = newposx - self.getCenter()[0]
		self.deltay = newposx - self.getCenter()[0]

		self.angle = angle
		self.anglex = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.scaley = scaley
		self.scalex = scalex

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]
		self.scrolate_trigger()

	def scale_kernal(self):

		if self.appear:

			for i in range(0, self.items):

				self.evl =  (math_appear(self.supply_p_coords[i][0],self.supply_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.supply_p_coords[i][0],self.supply_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_appear(self.supply_p_coords[i][2],self.supply_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.supply_p_coords[i][2],self.supply_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.supply_p[i],self.evl)

				self.supply_p_coords_tmp[i] = self.evl
		else:
			for i in range(0, self.items):

				self.evl =  (math_scale(self.supply_p_coords[i][0],self.supply_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.supply_p_coords[i][0],self.supply_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_scale(self.supply_p_coords[i][2],self.supply_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.supply_p_coords[i][2],self.supply_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.supply_p[i],self.evl)

				self.supply_p_coords_tmp[i] = self.evl

	def scale_trigger(self):

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.supply_p_coords = self.supply_p_coords_tmp[:]
			self.appear = 0
			return

		self.scale_kernal()
		master.after(int(1000./framerate), self.scale_trigger)

	def scale(self, scalex, scaley, animation_time):
	 	self.animation_time = animation_time
	 	self.itr_scale = 0
	 	self.scaley = scaley
	 	self.scalex = scalex
	 	self.no_of_itr_scale = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.scale_trigger()

class GROUND:

	def init(self, posx, posy, color='#FFFFFF'):

		self.posx = posx
		self.posy = posy
		self.color = color

		self.items = 4
		self.ground_p = [0]*self.items
		self.ground_p_coords = [0]*self.items
		self.ground_p_coords_tmp = [0]*self.items
		self.ground_p_coords_tmp2 = [0]*self.items
		self.ground_p_coords_tmp3 = [0]*self.items
		self.xmin_tmp = [0]*self.items
		self.xmax_tmp = [0]*self.items
		self.ymin_tmp = [0]*self.items
		self.ymax_tmp = [0]*self.items

		self.itr_move = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.no_of_itr_move = 0
		self.no_of_itr_rotate = 0
		self.no_of_itr_scale = 0
		self.anglex = 0
		self.angle = 0
		self.deltax = 0
		self.deltay = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.scalex_int = 0
		self.scaley_int = 0
		self.deltay1 = 0
		self.scalex = 0
		self.scaley = 0

		self.ground_p_coords[0] = (posx, posy, posx, posy-20)
		self.ground_p_coords[1] = (posx-10, posy, posx+10, posy)
		self.ground_p_coords[2] = (posx-10, posy, posx, posy+20)
		self.ground_p_coords[3] = (posx+10, posy, posx, posy+20)

		self.reCenter()

		self.ground_p[0] = c.create_line(self.ground_p_coords[0], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.ground_p[1] = c.create_line(self.ground_p_coords[1], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.ground_p[2] = c.create_line(self.ground_p_coords[2], arrow='none', fill=self.color, width=2, capstyle=ROUND)
		self.ground_p[3] = c.create_line(self.ground_p_coords[3], arrow='none', fill=self.color, width=2, capstyle=ROUND)

		self.appear = 1
		self.scale(1,1,500)


	def getTerminals(self):
		return [[self.ground_p_coords[0][0], self.ground_p_coords[0][1]], [self.ground_p_coords[11][2], self.ground_p_coords[11][3]]]

	def reCenter(self):
		self.rotate_around_x_old = self.getCenter()[0]
		self.rotate_around_y_old = self.getCenter()[1]

		for i in range(0, self.items):
			self.ground_p_coords[i] = (self.ground_p_coords[i][0] + self.posx - self.rotate_around_x_old, self.ground_p_coords[i][1] + self.posy - self.rotate_around_y_old, self.ground_p_coords[i][2] + self.posx - self.rotate_around_x_old, self.ground_p_coords[i][3] + self.posy -self.rotate_around_y_old)

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]

	def getCenter(self):
		for i in range(0, self.items):
			self.xmin_tmp[i] = min(self.ground_p_coords[i][0], self.ground_p_coords[i][2])
			self.xmax_tmp[i] = max(self.ground_p_coords[i][0], self.ground_p_coords[i][2])
			self.ymin_tmp[i] = min(self.ground_p_coords[i][1], self.ground_p_coords[i][3])
			self.ymax_tmp[i] = max(self.ground_p_coords[i][1], self.ground_p_coords[i][3])
		return [(min(self.xmin_tmp) + max(self.xmax_tmp))/2, (min(self.ymin_tmp) + max(self.ymax_tmp))/2]

	def rotate_kernal(self):

		for i in range(0, self.items):

			self.evl =  (math_rotate(self.ground_p_coords[i][0], self.ground_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.ground_p_coords[i][0], self.ground_p_coords[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.ground_p_coords[i][2], self.ground_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.ground_p_coords[i][2], self.ground_p_coords[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])
			c.coords(self.ground_p[i],self.evl)

			self.ground_p_coords_tmp[i] = self.evl

	def rotate_trigger(self):

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		if self.itr_rotate > self.no_of_itr_rotate:
			self.ground_p_coords = self.ground_p_coords_tmp[:]
			return

		self.rotate_kernal()
		master.after(int(1000./framerate), self.rotate_trigger)

	def rotate(self, angle, animation_time):
	 	self.animation_time = animation_time
	 	self.angle = angle
	 	self.anglex = 0
	 	self.itr_rotate = 0
	 	self.no_of_itr_rotate = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.rotate_trigger()

	def move_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_translate(self.ground_p_coords[i][0], self.ground_p_coords[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.ground_p_coords[i][0], self.ground_p_coords[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.ground_p_coords[i][2], self.ground_p_coords[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.ground_p_coords[i][2], self.ground_p_coords[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.ground_p[i], self.evl)

			self.ground_p_coords_tmp[i] = (self.evl)

	def move_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		if self.itr_move > self.no_of_itr_move:
			self.ground_p_coords = self.ground_p_coords_tmp[:]
			return
		self.move_kernal()
		master.after(int(1000./framerate), self.move_trigger)

	def translate(self, newposx, newposy, animation_time):
	 	self.itr_move = 0
	 	self.deltax1 = 0
	 	self.deltay1 = 0
	 	self.no_of_itr_move = self.animation_time * framerate/1000
	 	self.animation_time = animation_time
	 	self.deltax = newposx - self.getCenter()[0]
	 	self.deltay = newposy - self.getCenter()[1]
	 	self.move_trigger()


	def scrolate_kernal(self):
		for i in range(0,self.items):

			self.evl = (math_scale(self.ground_p_coords[i][0],self.ground_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.ground_p_coords[i][0],self.ground_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
						math_scale(self.ground_p_coords[i][2],self.ground_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
						math_scale(self.ground_p_coords[i][2],self.ground_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

			c.coords(self.ground_p[i], self.evl)

			self.ground_p_coords_tmp[i] = self.evl

			self.evl = (math_rotate(self.ground_p_coords_tmp[i][0], self.ground_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.ground_p_coords_tmp[i][0], self.ground_p_coords_tmp[i][1], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1],
									math_rotate(self.ground_p_coords_tmp[i][2], self.ground_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[0],
									math_rotate(self.ground_p_coords_tmp[i][2], self.ground_p_coords_tmp[i][3], self.rotate_around_x, self.rotate_around_y, self.anglex/(1))[1])

			c.coords(self.ground_p[i], self.evl)

			self.ground_p_coords_tmp2[i] = self.evl

			self.evl = (math_translate(self.ground_p_coords_tmp2[i][0], self.ground_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[0],
									math_translate(self.ground_p_coords_tmp2[i][0], self.ground_p_coords_tmp2[i][1], self.deltax1, self.deltay1)[1],
									math_translate(self.ground_p_coords_tmp2[i][2], self.ground_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[0],
									math_translate(self.ground_p_coords_tmp2[i][2], self.ground_p_coords_tmp2[i][3], self.deltax1, self.deltay1)[1])

			c.coords(self.ground_p[i], self.evl)

			self.ground_p_coords_tmp3[i] = (self.evl)

			
	def scrolate_trigger(self):

		self.itr_move = self.itr_move + 1

		self.deltax1 = self.deltax * ease(self.itr_move/self.no_of_itr_move)
		self.deltay1 = self.deltay * ease(self.itr_move/self.no_of_itr_move)

		self.itr_rotate = self.itr_rotate + 1

		self.anglex = self.angle*ease(self.itr_rotate/self.no_of_itr_rotate)

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.ground_p_coords = self.ground_p_coords_tmp3[:]
			return	
		
		self.scrolate_kernal()
		master.after(int(1000./framerate), self.scrolate_trigger)

	def scrolate(self, newposx, newposy, angle, scalex, scaley, animation_time):
		self.itr_move = 0
		self.deltax1 = 0
		self.deltay1 = 0
		self.no_of_itr_move = self.animation_time * framerate/1000
		self.no_of_itr_scale = self.animation_time * framerate/1000
		self.no_of_itr_rotate = self.animation_time * framerate/1000
		self.animation_time = animation_time
		self.deltax = newposx - self.getCenter()[0]
		self.deltay = newposx - self.getCenter()[0]

		self.angle = angle
		self.anglex = 0
		self.itr_rotate = 0
		self.itr_scale = 0
		self.scaley = scaley
		self.scalex = scalex

		self.rotate_around_x = self.getCenter()[0]
		self.rotate_around_y = self.getCenter()[1]
		self.scrolate_trigger()

	def scale_kernal(self):

		if self.appear:

			for i in range(0, self.items):

				self.evl =  (math_appear(self.ground_p_coords[i][0],self.ground_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.ground_p_coords[i][0],self.ground_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_appear(self.ground_p_coords[i][2],self.ground_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_appear(self.ground_p_coords[i][2],self.ground_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.ground_p[i],self.evl)

				self.ground_p_coords_tmp[i] = self.evl
		else:
			for i in range(0, self.items):

				self.evl =  (math_scale(self.ground_p_coords[i][0],self.ground_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.ground_p_coords[i][0],self.ground_p_coords[i][1],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1],
							math_scale(self.ground_p_coords[i][2],self.ground_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[0],
							math_scale(self.ground_p_coords[i][2],self.ground_p_coords[i][3],self.rotate_around_x,self.rotate_around_y,self.scalex_int,self.scaley_int)[1])

				c.coords(self.ground_p[i],self.evl)

				self.ground_p_coords_tmp[i] = self.evl

	def scale_trigger(self):

		self.itr_scale = self.itr_scale + 1

		self.scalex_int = self.scalex * ease(self.itr_scale/ self.no_of_itr_scale)
		self.scaley_int = self.scaley * ease(self.itr_scale/ self.no_of_itr_scale)

		if self.itr_scale > self.no_of_itr_scale:
			self.ground_p_coords = self.ground_p_coords_tmp[:]
			self.appear = 0
			return

		self.scale_kernal()
		master.after(int(1000./framerate), self.scale_trigger)

	def scale(self, scalex, scaley, animation_time):
	 	self.animation_time = animation_time
	 	self.itr_scale = 0
	 	self.scaley = scaley
	 	self.scalex = scalex
	 	self.no_of_itr_scale = self.animation_time*framerate/1000
	 	self.rotate_around_x = self.getCenter()[0]
	 	self.rotate_around_y = self.getCenter()[1]
	 	self.scale_trigger()




class Orchestra:

	def __init__(self):
		self.seq = 0

	def animate_kernal(self):
		if self.seq == 0:
			re.init(100,100,100,100)
		if self.seq == 1:
			re.scrolate(200,200,90,1,1,1000)
		return

	def animate_trigger(self):
		self.animate_kernal()
		self.seq = self.seq + 1
		master.after(2000, self.animate_trigger)

	def animate(self):	 	
	 	self.animate_trigger()

re = BOX()
o = Orchestra()
o.animate()

master.mainloop()
